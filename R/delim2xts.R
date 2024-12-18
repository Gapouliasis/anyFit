#' @title delim2xts
#'
#' @description A function to convert data from delimited files to xts objects.
#' The function takes a file path, time zone, date index, time step, delimiter,
#' skip rows, column names, no value, and exc_leaps (exclude 29th of February from leap years) as arguments.
#' If the time step of the timeseries is strict then the function will try to identify potential gaps in the data and indicate them with NANs.
#' Otherwise the function will return the raw dataset.
#' The user is given the option to save the xts dataset in a text file.
#'
#' @param file_path Path to the file with the data.
#' @param time_zone Time zone of the data.
#' @param date_index Position of the date column in the data file.
#' @param strict_step Logical that indicates if the timeseries are recorded in strict timestep.
#' @param delim Delimiter used in the data file.
#' @param skip_rows Number of rows to be skipped in the data file.
#' @param col_names Logical that indicates if the data file has column names.
#' @param no_value Value used as NA.
#' @param exc_leaps Logical that indicates if leap years shall be excluded.
#' @param save_Xts Logical that indicates if the xts shall be saved as a txt file.
#' @param filename Name of the file where the xts will be saved.
#'
#' @return An xts object.
#'
#' @examples
#'file_path <- system.file("extdata", "KNMI_Daily.csv", package = "anyFit")
#'time_zone <- "UTC"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", strict_step = TRUE)
#'
#' @import xts
#' @import lubridate
#'
#' @export

#Function that reads data with dates from a delimited file and produces an xts object. With arguments for leap years
# and saving the xts
delim2xts <- function(file_path,time_zone,date_index = 1,strict_step = TRUE,
                      delim = '\t', skip_rows = 0, col_names = TRUE,
                       no_value = NA, exc_leaps = TRUE,save_Xts = FALSE, filename = NA){

  RawData <- readr::read_delim(file_path, skip = skip_rows, col_names = col_names,
                          delim = delim, escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

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

  dates = parse_date_time(Data[date_index, drop = TRUE], orders = funs, tz = time_zone)

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

