#' @title delim2xts
#'
#' @description A function to convert data from delimited files to xts objects.The function takes a file path, time zone, date index, time step, delimiter, skip rows, column names, no value, and exc_leaps (exclude 29th of February from leap years) as arguments. The use is given the option to save the xts dataset in a text file.
#'
#' @param file_path Path to the file with the data.
#' @param time_zone Time zone of the data.
#' @param date_index Position of the date column in the data file.
#' @param time_step Time step of the timeseries.
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
#'time_step <- "1 day"
#'
#'data <- delim2xts(file_path = file_path,
#'                  time_zone = "UTC", delim = " ", time_step = time_step)
#'
#' @import xts
#' @import lubridate
#'
#' @export

#Function that reads data with dates from a delimited file and produces an xts object. With arguments for leap years
# and saving the xts
delim2xts <- function(file_path,time_zone,date_index = 1,time_step = '1 day',delim = '\t', skip_rows = 0, col_names = TRUE,
                       no_value = NA, exc_leaps = TRUE,save_Xts = FALSE, filename = NA){
  # library(xts)
  # library(lubridate)
  # library(rdwd)
  # library(readr)

  RawData <- readr::read_delim(file_path, skip = skip_rows, col_names = col_names,
                          delim = delim, escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

  if (length(which(RawData == no_value))>0){
    Rawdata[RawData == no_value] <- NA
    I <- which(!RawData[,date_index] == no_value)
    Data <- RawData[I,]
  }else{
    Data <- RawData
  }


  # create empty xts full series to put variables
  funs = c(ymd, ydm, mdy, myd, dmy, dym,
           ymd_h, dmy_h, mdy_h, ydm_h,
           ymd_hm, dmy_hm, mdy_hm, ydm_hm,
           ymd_hms, dmy_hms, mdy_hms, ydm_hms)
  for (tfun in funs){
    param_list = list(data = Data[, date_index, drop = TRUE])
    param_list$tz = 'UTC'
    dates = tryCatch({do.call(tfun,param_list)},
                     warning = function(w) {})
    if (!is.null(dates)){
      break
    }
  }

  starttime = dates[1]
  endtime = dates[nrow(Data)]
  # starttime=as.POSIXct(Data[1, date_index, drop = TRUE],format=date_format, tz = time_zone)
  # endtime=as.POSIXct(Data[nrow(Data), date_index, drop = TRUE],format=date_format, tz = time_zone)

  datesNA=seq(from= starttime, to=endtime, by=time_step)

  ncols <- ncol(Data) - length(date_index)
  xtsNA=xts(x = matrix(NA,ncol=ncols,nrow=length(datesNA)), order.by = datesNA)

  # dates<-Data[,date_index, drop = TRUE]
  # dates<-as.POSIXct(dates,format=date_format, tz = time_zone)

  xtsInpData=xts(Data[,-1], order.by = dates)

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

  if (save_Xts == TRUE) {
    write.zoo(Dataxts,file=paste0(filename,'.txt'),sep='\t')
  }
  return(Dataxts)
}

