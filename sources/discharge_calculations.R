#' function to calculate discharge for events

# Initialization ----------------------------------------------------------
# load packages
library(tidyverse)
library(lubridate)
library(zoo)

# load data (only for testing)
#ev_samples <- read_csv("data/Event_samples.csv")
#ev_dates <- read_csv("ev_dates.csv")
#wb_data <- read_csv("data/WB_data.csv") # data from waterboard (Ha)
#cr6 <- read_csv("data/CR6.csv", col_types = "TddddddddddddddD") # from cr6 (Hb)
#crest_nap <- 79.331 # m + NAP as measured by waterboard

# Functions ----------------------------------------------------------------
#'Flow correction formula for discharge from a 2 feet Parshall flume
Qcor <- function(Ha,w,K) {
  Haf <- Ha * 3.28084 #m to feet
  wf <- w * 3.28084
  C <- ((Haf/(((1.8/K)^1.8)-2.45))^(4.57-3.14*K)+0.093*K)*(wf^0.815) #correction for submergence
  Qcor <- C / 35.3147 #cubic feet to m3
  return(Qcor)
}

#'Free flow formula for discharge from a 2 feet Parshall flume
Qparshall <- function(Ha,w) {
  Haf <- Ha * 3.28084 #m to feet
  wf <- w * 3.28084
  Q <- (4*wf*Haf^(1.522*wf^0.026)) / 35.3147
  return(Q)
}

# calculate discharge variable for events with K
discharge_vars <- function(event, crest_nap) {
  x <- event %>%
    mutate(Ha = Wh - crest_nap,
           K = (wat_level-0.015) / Ha,
           K = if_else(K > 0.90, 0.90, K),  #K should not be above 0.9, because the formula will result in negative values
           K_int = na.approx(K, rule = 2),  #interpolate the minutes between the WB_data
           Ha_int = (wat_level - 0.015) / K_int, #predict Ha for minutes between WB_data
           Qcor = Qcor(Ha, 0.609, K),
           Qold = Qparshall(Ha, 0.609),
           H_raw = if_else(wat_level < -0.014, 0, wat_level + 0.015),
           Q_raw = Qparshall(H_raw, 0.609),
           Q = if_else(K < 0.5, Qold, Qold-Qcor),
           Q = if_else(Q < 0, lag(Q), Q),
           Q_int = na.approx(Q, rule = 2))
  return(x)
}

# Calculate discharge -----------------------------------------------------

#' timestamps require special attention! All data is saved in UTC, so is now also loaded in UTC.
#' in the raw data corrections are done for DST errors, but there are also small syncing errors between
#' wb_data and cr6. The wb_data will be used as leading timeseries.

#' Calculate discharge from waterheight, using original formula by Parshall (1928).
#' all data is collected in metric units, the formula works with feet, so reshape all.
#' calculate corrected Q for all events in 2019 and 2020
#' time synchronization is done based on max water height for both cr6 and wb_data time series. not sure if this is 
#' the best approach. Other option is to synchronize based on the rising limb...

Q_calc <- function(dates, dat1, dat2, crest_nap) {
dates1 <- dates$date
name <- str_remove_all(dates$date, "-")
df_list <- vector("list", length = length(name))
df_tcor <- vector("list", length = length(name))
for (i in seq_along(dates1)) {
  ev_start <- dates$start_t[i]
  ev_end   <- dates$end_t[i]
  ev_wb    <-  dat1 %>%
    filter(timestamp > ev_start & timestamp < ev_end)
  ev_cr6   <- dat2 %>%
    filter(timestamp > ev_start & timestamp < ev_end)
  #calculate time difference based on max water height.
  wb_max <- ev_wb %>%
    filter(Wh == max(Wh))
  cr6_max <- ev_cr6 %>%
    mutate(h = if_else(is.na(wat_level), level_mm / 1000, wat_level)) %>%
    filter(h == max(h))
  t_cor <- wb_max$timestamp[1] - cr6_max$timestamp[1]
  df_tcor[[i]] <- t_cor
  ev_cr6 <- ev_cr6 %>%
    mutate(timestamp = timestamp + t_cor)
  df_list[[i]] <- ev_cr6 %>%
    left_join(dat1, by = "timestamp")
}

#calculate discharge for events with Hb data
dates_no_k <- dates %>%
  filter(date < "2019-08-01")
s <- length(dates_no_k$date) # events before this date do not have Hb
f <- length(dates$date)
for(i in (s+1):f) {
  #calculate discharge
  df_list[[i]] <- discharge_vars(df_list[[i]], crest_nap)
}

#' events before 2019-07 do not have Hb (wat_level) measurements. So the correction for Q cannot directly be calculated.
#' To make an estimate of free flow conditions, the submergence characteristics after 08-2019 are studied.
ev_data <- bind_rows(df_list) %>%
  filter(date(timestamp) > "2019-08-01") %>%
  filter(K > 0)
m <- lm(K ~ Ha + I(Ha^2) + I(Ha^3) + I(Ha^4), data = ev_data)

# calculate Qcor for early events
for(i in 1:s) {
  #calculate discharge
  df_list[[i]] <- df_list[[i]] %>%
    mutate(Ha = Wh - crest_nap)  #filter(!is.na(Ha))
  newdat <- data.frame(df_list[[i]]$Ha)
  names(newdat) <- "Ha"
  K_pred <-  predict(m, newdat)
  df_list[[i]] <- bind_cols(df_list[[i]], enframe(K_pred, name = NULL, value = "K_pred")) %>%
    mutate(Qold = Qparshall(Ha, 0.609),
           level_m = level_mm / 1000,
           Q_raw = Qparshall(level_m, 0.609)) %>%
    mutate(Qcor = Qcor(Ha, 0.609, K_pred),
           Q = if_else(K_pred < 0.5, Qold, Qold-Qcor),
           Q_int = na.approx(Q, rule = 2))
}
Q_tcor <- vector("list", length = 2)
Q_tcor[[1]] <- df_list
Q_tcor[[2]] <- df_tcor
return(Q_tcor)
}

#test
#df_dat <- Q_calc(ev_dates, wb_data, cr6, crest_nap)
