#' Create a Per-Column Fitting Worker
#'
#' Closure factory that returns a function applying \code{\link{fitlm_multi}} to
#' each column of a data chunk. Isolates candidate, zero-filtering, and plotting
#' parameters from the outer calling environment so the returned closure can be
#' safely exported to \pkg{future} workers without capturing large grid objects.
#'
#' @param candidates A list of distribution names to fit.
#' @param ignore_zeros Logical; exclude zeros if \code{TRUE}.
#' @param zero_threshold Numeric threshold for zero filtering.
#' @param diagnostic_plots Logical; produce diagnostic plots if \code{TRUE}.
#' @param order Named list of L-moment orders per candidate.
#'
#' @return A function that takes a data chunk and returns a list of
#'   \code{fitlm_multi} results, one per column.
#' @noRd
.make_col_worker <- function(candidates, ignore_zeros, zero_threshold, diagnostic_plots, order) {
  force(candidates); force(ignore_zeros); force(zero_threshold); force(diagnostic_plots); force(order)
  function(chunk) lapply(seq_len(ncol(chunk)), function(j)
    fitlm_multi(chunk[, j, drop = FALSE], candidates = candidates,
                ignore_zeros = ignore_zeros, zero_threshold = zero_threshold,
                diagnostic_plots = diagnostic_plots, order = order))
}

#' Create a Bigmemory-Backed Fitting Worker
#'
#' Closure factory that returns a function which attaches a bigmemory descriptor
#' inside each \pkg{future} worker, extracts the assigned columns from the
#' shared memory-mapped matrix, and passes them to a column-wise fitting function.
#' The worker releases the memory map after copying its slice so the master
#' process can safely unlink the backing file.
#'
#' @param run_chunk A function taking a data chunk and returning fitting results.
#' @param desc A bigmemory descriptor object.
#'
#' @return A function that takes a column-index vector and returns fitting
#'   results for those columns.
#' @noRd
.make_mmap_worker <- function(run_chunk, desc) {
  force(run_chunk); force(desc)
  function(cols) {
    bmw <- bigmemory::attach.big.matrix(desc)
    chunk <- bmw[, cols, drop = FALSE]   # copy this worker's columns out
    rm(bmw); gc()                        # release the mmap so master can unlink
    run_chunk(chunk)
  }
}

#' Fit Distributions to Multiple Time Series
#'
#' @description Fits a set of candidate probability distributions to every column
#' of an xts object using the L-moment method. The function wraps
#' \code{\link{fitlm_multi}} for each series and collects the estimated parameters,
#' goodness-of-fit metrics, and diagnostic plots (Q-Q, P-P) into a list.
#' When \code{parallel = TRUE}, the work is distributed across available cores
#' via the \pkg{future} framework, with either file-backed shared memory (bigmemory
#' memory-mapped matrices for single-machine parallelism) or serialised column chunks
#' (for multi-node clusters). File-backed shared memory can be enabled by \code{shared_memory = TRUE}. 
#' This is reccomended for single machine usage and especially windows for efficiency and reduced RAM consumption. 
#' Panel layouts of the per-series Q-Q and P-P plots are
#' assembled via \pkg{patchwork} in pages of \code{nrow} by \code{ncol} sub-plots.
#' This function is designed for large ensembles of grid cells or stations where
#' fitting individually would be impractical.
#'
#' @param ts An xts object containing the time series data (multiple columns).
#' @param candidates A list of distribution names to fit (e.g. \code{list('exp','gamma3','gengamma')}).
#' @param nrow Number of rows per panel page in the diagnostic-plot grid.
#' @param ncol Number of columns per panel page in the diagnostic-plot grid.
#' @param ignore_zeros Logical. If \code{TRUE}, zeros are excluded before fitting. Default \code{FALSE}.
#' @param zero_threshold Numeric. Values below this threshold are treated as zero. Default 0.01.
#' @param parallel Logical. If \code{TRUE}, use parallel processing via \pkg{future}. Default \code{FALSE}.
#' @param ncores Number of worker processes when \code{parallel = TRUE} and no
#'   pre-existing \code{future::plan()} is set. Default 2.
#' @param shared_memory Logical. When \code{parallel = TRUE}, share the data grid
#'   with workers via a file-backed bigmatrix (memory-mapped, single machine only).
#'   Set \code{FALSE} for multi-node \code{plan(cluster)} setups. Default \code{TRUE}.
#' @param diagnostic_plots Logical. If \code{TRUE}, produce Q-Q and P-P diagnostic
#'   plots. Default \code{TRUE}.
#' @param order Optional named list mapping candidate distribution names to the
#'   vector of L-moment orders used by their optimiser, e.g.
#'   \code{list(gengamma = 1:5, expweibull = 1:3)}. Passed through to
#'   \code{\link{fitlm_multi}}. Default \code{NULL}.
#'
#' @return A list with elements \code{params} (per-column parameter tables and
#'   GoF summaries), and, if \code{diagnostic_plots = TRUE}, \code{diagnostic_plots},
#'   \code{QQ_plots}, \code{PP_plots}, \code{QQ_panels}, and \code{PP_panels}.
#'
#' @examples
#' # Synthetic daily data: 3 stations, 10 years
#' set.seed(123)
#' n <- 3650
#' ts <- xts::xts(cbind(station1 = rgamma(n, shape = 2, scale = 5),
#'                 station2 = rgamma(n, shape = 3, scale = 7),
#'                 station3 = rgamma(n, shape = 2.5, scale = 6)),
#'           order.by = seq.Date(as.Date("2000-01-01"), by = "day", length.out = n))
#'
#' fits <- fitlm_nxts(ts, candidates = list('exp','gamma3'),
#'                    nrow = 2, ncol = 2, ignore_zeros = FALSE)
#' fits$params[[1]]
#'
#' @export
fitlm_nxts <- function(ts, candidates, nrow = 5, ncol = 4, ignore_zeros = FALSE,
                       zero_threshold = 0.01, parallel = FALSE, ncores = 2,
                       shared_memory = TRUE, diagnostic_plots = TRUE, order = NULL){
  variables <- colnames(ts)

  run_chunk <- .make_col_worker(candidates, ignore_zeros, zero_threshold, diagnostic_plots, order)

  if (!parallel) {
    multi_fits <- run_chunk(ts)
  } else {
    # Respect a user-set HPC plan (cluster/batchtools/...); otherwise default per-OS
    # sized by `ncores`, and restore the previous plan on exit.
    if (inherits(future::plan(), "sequential")) {
      # Resolve the strategy to a value first: future::plan() evaluates its strategy
      # argument non-standardly, so passing the inline `if (...) ... else ...` expression
      # makes it try to tweak the `if` primitive itself (ls(envir=) error).
      strategy <- if (.Platform$OS.type == "windows") future::multisession else future::multicore
      oplan <- future::plan(strategy, workers = min(ncores, ncol(ts)))
      on.exit(future::plan(oplan), add = TRUE)
    }

    # Degree of parallelism follows the ACTIVE plan, not `ncores`, so a user who set
    # plan(cluster, workers = 64) uses all 64. One chunk per worker.
    nw <- future::nbrOfWorkers(); if (!is.finite(nw)) nw <- min(ncores, ncol(ts))
    idx <- parallel::splitIndices(ncol(ts), min(nw, ncol(ts)))

    # shared memory is a LOCAL mmap; fall back to chunks if the plan spans more than
    # the local machine's cores (its temp file is invisible to remote nodes).
    if (shared_memory && nw > parallel::detectCores()) {
      warning("shared_memory = TRUE is single-machine; plan has ", nw,
              " workers (> local cores) -> using serialized chunks instead.")
      shared_memory <- FALSE
    }

    if (shared_memory) {
      bk <- basename(tempfile("grid_")); dsc <- paste0(bk, ".desc")
      bm <- bigmemory::as.big.matrix(coredata(ts), type = "double",
              backingpath = tempdir(), backingfile = bk, descriptorfile = dsc)
      # rm + gc releases the mmap handle so the backing file can be deleted
      # (a still-mapped file is undeletable on Windows).
      on.exit({ rm(bm); gc(); unlink(file.path(tempdir(), c(bk, dsc))) }, add = TRUE)
      desc <- bigmemory::describe(bm)
      worker <- .make_mmap_worker(run_chunk, desc)
      fits <- future.apply::future_lapply(idx, worker, future.seed = TRUE,
                future.packages = c("anyFit", "bigmemory"))
    } else {
      # Chunks can exceed future's default 500 MB global limit on large grids.
      oopt <- options(future.globals.maxSize = Inf); on.exit(options(oopt), add = TRUE)
      chunks <- lapply(idx, function(cols) ts[, cols, drop = FALSE])
      fits <- future.apply::future_lapply(chunks, run_chunk, future.seed = TRUE,
                future.packages = "anyFit")
    }
    multi_fits <- unlist(fits, recursive = FALSE)   # back to per-column order
  }


  # params <- list()
  # for (i in 1:ncol(ts)){
  #   params <- c(params,list(list(params = multi_fits[[i]]$parameter_list,GoF_summary = multi_fits[[i]]$GoF_summary)))
  # }

  params <- lapply(1:ncol(ts), FUN = function(i) {list(params = multi_fits[[i]]$parameter_list,GoF_summary = multi_fits[[i]]$GoF_summary)})
  names(params) <- variables

  if (diagnostic_plots){
    diagnostic_plots <- list()
    QQ_plots <- list()
    PP_plots <- list()
    for (i in 1:ncol(ts)){
      diagnostic_plots <- c(diagnostic_plots,list(multi_fits[[i]]$diagnostics))
      QQ_plots <- c(QQ_plots,l =  list(multi_fits[[i]]$QQplot + ggtitle(variables[i])))
      PP_plots <- c(PP_plots,l =  list(multi_fits[[i]]$PPplot + ggtitle(variables[i])))
    }

    names(diagnostic_plots) <- variables
    names(QQ_plots) <- variables
    names(PP_plots) <- variables

    n_plots <- nrow*ncol
    panels <- pracma::ceil(ncol(ts)/n_plots)
    QQ_panels <- list()
    PP_panels <- list()
    for (i in 1:panels){
      if ((i*(n_plots)) < ncol(ts)){
        temp_QQ <- QQ_plots[(1 + (i-1)*n_plots):(i*(n_plots))]
        temp_PP <- PP_plots[(1 + (i-1)*n_plots):(i*(n_plots))]
      }else{
        temp_QQ <- QQ_plots[(1 + (i-1)*n_plots):ncol(ts)]
        temp_PP <- PP_plots[(1 + (i-1)*n_plots):ncol(ts)]
      }
      QQ_panels <- c(QQ_panels, list(patchwork::wrap_plots(plotlist = temp_QQ, nrow = nrow, ncol = ncol)+
                                       plot_layout(guides = "collect") & theme(legend.position = 'bottom')))

      PP_panels <- c(PP_panels, list(patchwork::wrap_plots(plotlist = temp_PP, nrow = nrow, ncol = ncol)+
                                       plot_layout(guides = "collect") & theme(legend.position = 'bottom')))
    }

    list_out <- list(params = params, diagnostic_plots = diagnostic_plots, QQ_plots = QQ_plots,
                     PP_plots = PP_plots, QQ_panels = QQ_panels, PP_panels = PP_panels)
  }else{
    list_out <- list(params = params)
  }

  return(list_out)
}
