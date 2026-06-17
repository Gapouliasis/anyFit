# Worker factories live at NAMESPACE top level (not inside fitlm_nxts) so the
# closures they return have the package namespace as their parent environment and
# therefore do NOT capture the caller's grid `ts`. A factory defined inside
# fitlm_nxts would chain back to that frame and drag the whole grid into every
# future export (gigabytes) -- the exact bug this design exists to avoid.
# They precede the roxygen block below on purpose: a roxygen block attaches to the
# NEXT object, so keeping these above it (with no roxygen tags) leaves fitlm_nxts'
# documentation correctly bound to fitlm_nxts.
.make_col_worker <- function(candidates, ignore_zeros, zero_threshold, diagnostic_plots) {
  force(candidates); force(ignore_zeros); force(zero_threshold); force(diagnostic_plots)
  function(chunk) lapply(seq_len(ncol(chunk)), function(j)
    fitlm_multi(chunk[, j, drop = FALSE], candidates = candidates,
                ignore_zeros = ignore_zeros, zero_threshold = zero_threshold,
                diagnostic_plots = diagnostic_plots))
}

.make_mmap_worker <- function(run_chunk, desc) {
  force(run_chunk); force(desc)
  function(cols) {
    bmw <- bigmemory::attach.big.matrix(desc)
    chunk <- bmw[, cols, drop = FALSE]   # copy this worker's columns out
    rm(bmw); gc()                        # release the mmap so master can unlink
    run_chunk(chunk)
  }
}

#' @title fitlm_nxts
#'
#' @description This function fits a list of candidate distributions using the L-moments method
#' to an arbitrary number of timeseries in xts format. Additionally to the list of the fitted parameters,
#'  goodness-of-fit metric, PP and QQ plots.
#'
#' @param ts A xts object containing the time series data.
#' @param candidates A list of distribution to fit.
#' @param nrow Number of rows for plotting.
#' @param ncol Number of columns for plotting.
#' @param ignore_zeros A logical value, if TRUE zeros will be ignored. Default is FALSE.
#' @param zero_threshold The threshold below which values are considered zero. Default is 0.01.
#' @param parallel Logical, whether to use parallel processing.
#' @param ncores Number of cores to use in the case of parallel computations
#' @param shared_memory Logical, when parallel, share the grid with workers via a
#'   filebacked big.matrix (mmap, single machine) instead of serializing column
#'   chunks. Set FALSE for multi-node \code{plan(cluster)} setups. Default TRUE.
#' @param diagnostic_plots A logical value, controls the output of diagnostic plots
#'
#' @return A list with the estimated parameters, diagnostic plots, QQ plots and PP plots.
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#' candidates <- list('exp','expweibull', 'gamma3')
#'
#' fits_all <- fitlm_nxts(data, candidates, nrow = 5, ncol = 4, ignore_zeros = TRUE)
#' fits_all$QQ_panels[[1]]
#'
#' @export
fitlm_nxts <- function(ts, candidates, nrow = 5, ncol = 4, ignore_zeros = FALSE,
                       zero_threshold = 0.01, parallel = FALSE, ncores = 2,
                       shared_memory = TRUE, diagnostic_plots = TRUE){
  variables <- colnames(ts)

  run_chunk <- .make_col_worker(candidates, ignore_zeros, zero_threshold, diagnostic_plots)

  if (!parallel) {
    multi_fits <- run_chunk(ts)
  } else {
    # Respect a user-set HPC plan (cluster/batchtools/...); otherwise default per-OS
    # sized by `ncores`, and restore the previous plan on exit.
    if (inherits(future::plan(), "sequential")) {
      oplan <- future::plan(
        if (.Platform$OS.type == "windows") future::multisession else future::multicore,
        workers = min(ncores, ncol(ts)))
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

