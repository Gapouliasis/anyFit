#' @title delim2xts
#'
#' @description Reads a delimited text file and returns an xts object with
#'   auto-detected timestamps. The date column is parsed via
#'   \code{\link[lubridate]{parse_date_time}} using a comprehensive set of
#'   order formats (ymd, dmy, mdy, and variants with hours, minutes, and
#'   seconds). When \code{strict_step = TRUE}, the function detects the
#'   regular time step from the first difference of parsed dates, constructs
#'   a complete index spanning the full date range, and fills missing
#'   positions with \code{NA} entries. Leap days (29 February) may be
#'   excluded from leap years via \code{exc_leaps}. A sentinel
#'   no-data value may be specified via \code{no_value} and replaced with
#'   \code{NA}. The resulting xts can optionally be saved to a
#'   tab-delimited text file.
#'
#' @param file_path Character. Path to the delimited data file.
#' @param time_zone Character. Time zone for date parsing (e.g. \code{"UTC"}).
#' @param date_index Integer. Column position of the date/time variable.
#'   Default \code{1}.
#' @param strict_step Logical. If \code{TRUE} (default), gaps are filled
#'   with \code{NA} at the detected regular interval.
#' @param delim Character. Field delimiter. Default \code{"\t"}.
#' @param skip_rows Integer. Number of header rows to skip. Default \code{0}.
#' @param col_names Logical. If \code{TRUE} (default), the first non-skipped
#'   row contains column names.
#' @param no_value Value to be treated as missing (\code{NA}). Default
#'   \code{NA}.
#' @param exc_leaps Logical. If \code{TRUE} (default), 29 February is
#'   excluded from leap years.
#' @param save_Xts Logical. If \code{TRUE}, the xts is saved as a
#'   tab-delimited text file. Default \code{FALSE}.
#' @param filename Character. Output filename (without extension) used when
#'   \code{save_Xts = TRUE}. Default \code{NA}.
#'
#' @return An \code{\link[xts]{xts}} object with parsed timestamps as the
#'   index and the remaining columns as the data matrix.
#'
#' @examples
#' # Create a temporary delimited file
#' tf <- tempfile(fileext = ".csv")
#' writeLines(c("date,value", "2020-01-01,10.5", "2020-01-02,12.3",
#'              "2020-01-03,11.0", "2020-01-04,9.8"), tf)
#' x <- delim2xts(tf, time_zone = "UTC", delim = ",", strict_step = TRUE)
#' head(x)
#'
#' @importFrom lubridate parse_date_time
#' @importFrom xts xts
#' @importFrom zoo write.zoo
#' @importFrom stats time
#' @importFrom utils read.table
#'
#' @export

#Function that reads data with dates from a delimited file and produces an xts object. With arguments for leap years
# and saving the xts
delim2xts <- function(file_path,time_zone,date_index = 1,strict_step = TRUE,
                      delim = '\t', skip_rows = 0, col_names = TRUE,
                       no_value = NA, exc_leaps = TRUE,save_Xts = FALSE, filename = NA){

  RawData <- read.table(file_path, skip = skip_rows, header = col_names,
                        sep = delim,  strip.white = FALSE)

  if (length(which(RawData == no_value))>0){
    Rawdata[RawData == no_value] <- NA
    I <- which(!RawData[,date_index] == no_value)
    Data <- RawData[I,]
  }else{
    Data <- RawData
  }

  funs = c("ymd", "ydm", "mdy", "myd", "dmy", "dym", "ymd H", "dmy H", "mdy H",
           "ydm H", "ymd HM", "dmy HM", "mdy HM", "ydm HM", "ymd HMS", "dmy HMS",
           "mdy HMS", "ydm HMS")

  dates = parse_date_time(Data[, date_index], orders = funs, tz = time_zone, train = FALSE)

  if (strict_step){
    starttime = min(dates)
    endtime = max(dates)
    # starttime=as.POSIXct(Data[1, date_index, drop = TRUE],format=date_format, tz = time_zone)
    # endtime=as.POSIXct(Data[nrow(Data), date_index, drop = TRUE],format=date_format, tz = time_zone)

    steps = diff(dates)
    if (min(steps)<max(steps)){
      stop('Time step is not strict')
    }else{
      time_step = steps[1]
      print(time_step)
    }

    datesNA=seq(from= starttime, to=endtime, by=time_step)

    ncols <- ncol(Data) - length(date_index)
    xtsNA=xts(x = matrix(NA,ncol=ncols,nrow=length(datesNA)), order.by = datesNA)

    # dates<-Data[,date_index, drop = TRUE]
    # dates<-as.POSIXct(dates,format=date_format, tz = time_zone)

    xtsInpData=xts(Data[,-date_index], order.by = dates)

    xtsNA[time(xtsInpData),] <- xtsInpData

    # exclude 29th of February from leap years
    if (exc_leaps == TRUE){
      myDates <- !(format(datesNA,'%m') == '02' & format(datesNA, '%d') == '29')
      leapdayspos<-which(myDates==FALSE)
      if (length(leapdayspos)>0){
        xtsNA<-xtsNA[-leapdayspos]
      }
    }

    colnames(xtsNA) <- colnames(Data)[-date_index]
    Dataxts<-xtsNA
  }else{
    Dataxts=xts(Data[,-date_index], order.by = dates)
    colnames(Dataxts) <- colnames(Data)[-date_index]
  }


  if (save_Xts == TRUE) {
    write.zoo(Dataxts,file=paste0(filename,'.txt'),sep='\t')
  }
  return(Dataxts)
}
