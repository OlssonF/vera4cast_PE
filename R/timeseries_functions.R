#' Resample timeseries
#' 
#' @param ts a vector of data 
#' @param length.out number - how long should the resampled timeseries be
#' @param n number - how many times should the sampling occur
#' @param doy vector of days associates to the ts (optional)
#'
#' @returns a dataframe 

resample <- function(ts, length.out = 100, n = 100, doy = NA) {
  
  length_ts <- length(ts)
  
  resample_list <- NULL
  
  start_index <- sample(1:(length_ts - length.out), size = n, replace = F)
  
  
  for (i in 1:n) {
  
  
    ts_subset <- data.frame(observation = ts[start_index[i]:(start_index[i]+length.out-1)],
                            index = seq(1, length.out, by=1)) |> 
      mutate(n = i)
    
    if (is.na(doy[1])) {
      resample_list <- bind_rows(ts_subset, resample_list)
    } else {
      ts_subset <- ts_subset |> 
        mutate(doy = doy[start_index[i]])
      resample_list <- bind_rows(ts_subset, resample_list)
    }
    
    # print(i)
  }
  return(resample_list)
}

#' Downsample timeseries
#' 
#' @param ts a dataframe of raw timeseries data (single site-variable-depth combination)
#' Needs a datetime column that is regularly sampled at the timestep smaller than the out_freq
#' @param out_freq string - what should the temporal frequency of the output be options are 'hourly', 'daily', 'weekly'
#' @param method string - how should the downsampling occur - 'aggregate' or 'sample'
#' @param target_out for 'sample' what is the target time to sample 
#' MM:SS for hourly - string, 
#' #' HH:MM:SS for daily - string, 
#' time and day of week for weekly - vector with the time and the day e.g. c("12:00:00", 2) - noon on Tuesday
#' @param max_distance how far from the target is okay? numeric vector length 1-2
#' hourly in minutes
#' for daily in hours
#' for weekly this is a vector that relates to the time and then day of the week e.g. c(1,1) - 1 hour and 1 day max distance 
#' Default is any distance - Inf
#' @param min_sample for aggregate what is the minimum number of samples required to calculate the average. numeric. Default = 1
#' @returns a dataframe

downsample <- function(ts, 
                       out_freq = 'daily', 
                       target_out = '12:00:00',
                       max_distance = Inf,
                       min_sample = 1,
                       method = 'sample') {
  
  ts <- ts |> 
    tsibble::as_tsibble(key = any_of(c("site_id", "depth_m", "variable")),
                        index = "datetime") |>
    tsibble::fill_gaps() 
  
  if (method == 'aggregate') {
    if (out_freq == 'hourly') {
      target_out <- '00:00'
      max_distance <- Inf
    }
    
    if (out_freq == 'daily') {
      target_out <- '12:00:00'
      max_distance <- Inf
      
    }
    
    if (out_freq == 'weekly') {
      target_out <- c('12:00:00', '1')
      max_distance <- c(Inf, Inf)
    }
  }
  
  
  # Do the downsampling
  if (str_detect(out_freq, pattern = 'hour')) {
    ts_downsample <- 
      ts |> 
      mutate(target = as_datetime(paste0(as_date(datetime), ' ',hour(datetime),':', target_out)), 
             hour = as_datetime(paste0(as_date(datetime), ' ', hour(datetime), ':00:00'))) |> 
      tsibble::index_by(hour) |> 
      tsibble::group_by_key() |> 
      # na.omit() |> 
      dplyr::summarise(observation_aggregate = mean(observation, na.rm = T),
                       observation_sample = nth(observation, which.min(abs(as.numeric(datetime - target)))), # takes the observation closest
                       n = n(),
                       distance = min(abs(as.numeric(datetime - target)))) |> 
      dplyr::mutate(observation_aggregate = ifelse(method == 'aggregate' & n > 3,
                                                   observation_aggregate, NA),
                    observation_sample = ifelse(method == 'sample' & distance <= max_distance,
                                                observation_sample, NA)) |> # removes if there are less than 3 obs in any one hour
      as_tibble()  |> 
      pivot_longer(cols = observation_aggregate:observation_sample,
                   values_to='observation', names_prefix = 'observation_') |> 
      filter(name == method) |> 
      rename(method = name,
             datetime = hour) |> 
      select(any_of(c("datetime", "observation", "site_id", "depth_m", "variable", "method")))
  }
  
  
  if (str_detect(out_freq, pattern = c('daily')) | 
      str_detect(out_freq, pattern = c('weekly'))) {
    
    target_time <- ifelse(out_freq == 'weekly', target_out[1], target_out) 
    max_distance_t <- ifelse(out_freq == 'weekly', max_distance[1], max_distance) 
    
    ts_downsample <- 
      ts |> 
      mutate(target = as_datetime(paste0(as_date(datetime), ' ', target_time)), 
             date = as_date(datetime),
             distance = abs(as.numeric(datetime - target))) |> 
      tsibble::index_by(date) |> 
      tsibble::group_by_key() |> 
      # na.omit() |> 
      dplyr::summarise(observation_aggregate = mean(observation, na.rm = T),
                       observation_sample = nth(observation, which.min(distance)), # takes the observation closest
                       n = n(),
                       distance = min(distance)) |> 
      dplyr::mutate(observation_aggregate = ifelse(method == 'aggregate' & n > min_sample,
                                                   observation_aggregate, NA),
                    observation_sample = ifelse(method == 'sample' & distance <= max_distance_t,
                                                observation_sample, NA)) |> 
      as_tibble()  |> 
      pivot_longer(cols = observation_aggregate:observation_sample,
                   values_to='observation', names_prefix = 'observation_') |> 
      filter(name == method) |> 
      rename(method = name) |> 
      select(any_of(c("date", "observation", "site_id", "depth_m", "variable", "method")))
  }
  
  
  if (str_detect(out_freq, pattern = 'week')) {
    target_day <- target_out[2] 
    max_distance_d <- max_distance[2]
    
    ts_downsample <- ts_downsample |> 
      tsibble::as_tsibble(key = any_of(c("site_id", "depth_m", "variable")),
                          index = "date") |>
      tsibble::fill_gaps()  |> 
      mutate(dow = lubridate::wday(date, week_start = 1)) |> 
      tsibble::group_by_key() |>
      tsibble::index_by(year_week = ~ yearweek(.)) |>
      dplyr::summarise(observation_aggregate = mean(observation, na.rm = T),
                       observation_sample = nth(observation, which.min(abs(as.numeric(dow - as.numeric(target_day))))), # takes the observation closest
                       n = n(),
                       distance = min(abs(as.numeric(dow - as.numeric(target_day), 'hours')))) |>
      dplyr::filter((method == 'aggregate' & n > min_sample) | (method == 'sample' & distance <= max_distance_d)) |> # removes if there are less than 3 obs in any one hour
      as_tibble()  |>
      pivot_longer(cols = observation_aggregate:observation_sample,
                   values_to='observation', names_prefix = 'observation_') |>
      filter(name == method) |>
      rename(method = name) |>
      select(any_of(c("year_week", "observation", "site_id", "depth_m", "variable", "method")))

  }
  
  return(ts_downsample)
}

#' Shuffle a timeseries
#' 
#' @param ts a vector of data 
#' @param times how many shuffled realisations should be generated
#' @returns a dataframe 

shuffle <- function(ts, times = 1) {
  
  resample_list <- NULL
  
  message('Generating ', times , ' shuffled realisation(s)')
  for (i in 1:times) {
    
    ts_shuffle <- data.frame(observation = sample(ts)) |> 
      mutate(n = i)
    
    resample_list <- bind_rows(ts_shuffle, resample_list)
    
  }
  return(resample_list)
}


#' Get season from date
#' 
#' @param ts a vector of data 
#' @param times how many shuffled realisations should be generated
#' @returns a dataframe
season <- function(dates) {
  WS <- as.Date("2012-12-15", format = "%Y-%m-%d") # Winter Solstice
  SE <- as.Date("2012-3-15",  format = "%Y-%m-%d") # Spring Equinox
  SS <- as.Date("2012-6-15",  format = "%Y-%m-%d") # Summer Solstice
  FE <- as.Date("2012-9-15",  format = "%Y-%m-%d") # Fall Equinox
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(dates, format="2012-%m-%d"))
  
  ifelse (d >= WS | d < SE, "Winter",
          ifelse (d >= SE & d < SS, "Spring",
                  ifelse (d >= SS & d < FE, "Summer", "Autumn")))
}
