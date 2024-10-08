library(tidyverse)
library(zoo)
library(RCurl)

source('R/get_targets.R')
source('R/PE_functions.R')
source('R/timeseries_functions.R')


# Get targets ---------------------------------------
# Get observational target data

fcre_EDI <- "https://pasta.lternet.edu/package/data/eml/edi/271/8/fbb8c7a0230f4587f1c6e11417fe9dce"
# bvre_EDI <- "https://pasta.lternet.edu/package/data/eml/edi/725/4/9adadd2a7c2319e54227ab31a161ea12"
bvre_EDI <- "https://pasta-s.lternet.edu/package/data/eml/edi/157/31/9adadd2a7c2319e54227ab31a161ea12" 
# fcre_met <-  "https://pasta.lternet.edu/package/data/eml/edi/389/8/d4c74bbb3b86ea293e5c52136347fbb0" 

fcre_depths <- c(1.6, 9)
bvre_depths <- c(1.5, 13)

targets <- get_targets(infiles = c(fcre_EDI, bvre_EDI),
                       is_met = c(F,F),
                       interpolate = T, maxgap = 12) |> 
  mutate(depth_m = ifelse(is.na(depth_m), 'met',
                          ifelse(depth_m < 5, 'surface', ifelse(depth_m > 5, 'bottom', NA))))

targets_PT1H <- downsample(ts = targets, 
                           out_freq = 'hourly', 
                           method = 'sample', 
                           target_out = '00:00',
                           max_distance = 20) # 20 minutes either side 

targets_P1D <- downsample(ts = targets, 
                          out_freq = 'daily', 
                          method = 'sample', 
                          target_out = '12:00:00',
                          max_distance = 2) # 2 hours either side



targets_P1W <- downsample(ts = targets, 
                          out_freq = 'weekly', 
                          method = 'sample', 
                          target_out = c('12:00:00', 2), # 2nd DOW ie Tuesday
                          max_distance = c(2, 1)) # 2 hours from 12 and 1 day either side

# interpolation?
targets_P1D_interp <- targets_P1D |> 
  na.omit() |> 
  tsibble::as_tsibble(index = date, key = c(site_id, variable, depth_m)) |> 
  tsibble::fill_gaps() |> 
  as_tibble() |> 
  group_by(variable, site_id, depth_m) |> 
  arrange(date) |> 
  mutate(observation = zoo::na.approx(observation, na.rm = T, maxgap = 3))


# Get temperature profiles
temp_profiles <- 
  map2(c('none', 'none'), c(fcre_EDI, bvre_EDI), get_temp_profiles) |>
  list_rbind() |> 
  filter(year(datetime) %in% c(2019, 2020, 2021, 2022, 2023))


strat_dates <- calc_strat_dates(density_diff = 0.1, temp_profiles = temp_profiles) |> 
  mutate(start = as_date(start),
         end = as_date(end))

# Plot the observations
targets_PT1H |> 
  ggplot(aes(x=datetime, y=observation, colour = depth_m)) +
  geom_line() +
  facet_wrap(site_id~variable, scales = 'free_y', nrow = 3)

targets_P1D |> 
  ggplot(aes(x=date, y=observation, colour = depth_m)) +
  geom_line() +
  facet_grid(variable~site_id, scales = 'free')

targets_P1W |> 
  filter(variable == 'Temp_C') |> 
  ggplot(aes(x=as_date(year_week), y=observation)) +
  geom_point() +
  facet_wrap(site_id~depth_m, scales = 'free', nrow = 3)


# Resample timeseries
resample_length <- 50 # how long are the sections?
resample_n <- 500 # how many times should the timeseries be sampled
targets_P1D_resample <- targets_P1D |>
  filter(year(date) >= 2021) # only the years where there are the same data
  na.omit() |> 
  tsibble::as_tsibble(index = date, key = c(variable, site_id, depth_m)) |> 
  arrange(variable, depth_m, site_id, date) |> 
  reframe(.by = c(variable, site_id, depth_m),
          resample(ts = observation, 
                   length.out = resample_length, 
                   n = resample_n, 
                   doy = date))

# Generate shuffles

targets_P1D_shuffled <- targets_P1D |>
  na.omit() |> 
  tsibble::as_tsibble(index = date, key = c(variable, site_id, depth_m)) |> 
  arrange(variable, depth_m, site_id, date) |>
  reframe(.by = c(variable, site_id, depth_m),
          shuffle(ts = observation, times = 500)) 
