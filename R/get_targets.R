#' get daily targets
#' 
#' @param infiles vector of EDI data files
#' @param interpolate should interpolation be carried out
#' @return a targets dataframe in VERA format

get_targets <- function(infiles, interpolate = T, maxgap = 12) {

  message('processing wq file(s)')
  
  standard_names <- data.frame(variable_new = c('Tw_C','SpCond_uScm', 'Chla_ugL', 'fDOM_QSU'),
                               variable = c('EXOTemp_C_1', 'EXOSpCond_uScm_1', 'EXOChla_ugL_1', 'EXOfDOM_QSU_1'))
  targets_wq <- NULL
  
  # Load data
  message('Reading WQ data from EDI...')
  for (i in 1:length(infiles)) {
    df <- read_csv(infiles[i], show_col_types = F, progress = F) |> 
      filter(Site == 50) |> 
      # mutate(site_id = ifelse(Reservoir == 'BVR', 'bvre', ifelse(Reservoir == 'FCR', 'fcre', Reservoir))) |> 
      rename(datetime = DateTime,
             site_id = Reservoir) |> 
      filter(site_id %in% c('FCR', 'BVR')) |> 
      select(-Site)
    
    df_flags <- df |> 
      select(any_of(c('datetime', 'site_id')) | contains('Flag') & contains(standard_names$variable)) |> 
      pivot_longer(cols = contains('Flag'),
                   names_to = 'variable', 
                   values_to = 'flag_value', 
                   names_prefix = 'Flag_')
    
    df_observations <- df |> 
      select(any_of(c('datetime', 'site_id')) | contains(standard_names$variable) & !contains('Flag')) |> 
      pivot_longer(cols = -any_of(c('site_id', 'datetime')),
                   names_to = 'variable', 
                   values_to = 'observation', 
                   names_prefix = 'Flag_')
    
    
    # Filter any flagged data
    message('Filtering flags...')
    df_long <- inner_join(df_observations, df_flags, 
                          # by = c('site_id', 'datetime', 'depth_m', 'variable'),
                          relationship = 'many-to-many') |> 
      na.omit() |> 
      filter(!flag_value %in% c(9, 7, 2, 1, 5, 3)) |> 
      select(-contains('flag')) |> 
      mutate(depth_m = as.numeric(str_split_i(variable, "_", 3)), 
             depth_m = ifelse(str_detect(variable, 'EXO') & site_id == 'FCR', 1.6, 
                              ifelse(str_detect(variable, 'EXO') & site_id == 'BVR', 1.5, depth_m)),
             variable = str_split_i(variable, "\\.", 1)) |> 
      full_join(standard_names) |> 
      mutate(variable = variable_new) |> 
      select(-variable_new)
    
    # get DO
    df_DO <- df |>
      select(datetime, site_id, any_of(c('RDO_mgL_9_adjusted', 'RDO_mgL_13', 'EXODO_mgL_1', 'EXODO_mgL_1.5'))) |>
      pivot_longer(cols = contains('DO'), names_to = 'variable', values_to = 'observation') |>
      mutate(depth_m = as.numeric(str_split_i(variable, "_", 3)), 
             depth_m = ifelse(str_detect(variable, 'EXO') & site_id == 'FCR', 1.6, 
                              ifelse(str_detect(variable, 'EXO') & site_id == 'BVR', 1.5, depth_m)),
             variable = 'DO_mgL') 
    
    #get bottom temp
    df_bottomT <- df |>
      select(datetime, site_id, starts_with(c('ThermistorTemp', 'Flag_ThermistorTemp'))) |>
      pivot_longer(cols = contains('Thermistor'), names_to = 'variable', values_to = 'observation') |>
      mutate(depth_m = str_split_i(variable, "_", -1),
             depth_m = as.numeric(ifelse(depth_m == 'surface', 0, depth_m))) |> 
      filter(depth_m == max(depth_m)) |> 
      pivot_wider(names_from = variable, values_from = observation) |> 
      filter(if_any(starts_with("Flag"),  ~.x == 0)) |> 
      select(-starts_with('Flag')) |> 
      rename(observation = starts_with('Thermistor')) |> 
      mutate(variable = 'Tw_C')
    
    # combine
    df_all <- df_long |>
      dplyr::bind_rows(df_DO) |> 
      dplyr::bind_rows(df_bottomT) |>
      na.omit() |> 
      dplyr::select(datetime, site_id, depth_m, observation, variable) |>
      dplyr::mutate(observation = ifelse(!is.finite(observation),NA,observation)) |> 
      tsibble::as_tsibble(key = any_of(c("site_id", "depth_m", "variable")),
                          index = "datetime") |>
      tsibble::fill_gaps() |> 
      as_tibble()
    
    
    targets_wq <- bind_rows(targets_wq, df_all) 
  }
 
  targets <- targets_wq
  
  if (interpolate == T) {
    test <- targets |> 
      group_by(site_id, depth_m, variable) |> 
      arrange(datetime) |> 
      mutate(observation = zoo::na.approx(observation, na.rm = T,
                                          maxgap = maxgap, rule = 1))
  }
  
  return(targets)
}

# ================================================================#


# =================== Temperature profiles ==================

get_temp_profiles <- function(current_file = 'none', historic_file){
  source('R/find_depths.R')
 
  if (current_file != 'none') {
    message('reading ', current_file)
    current_df <- readr::read_csv(current_file, show_col_types = F) |>
      dplyr::filter(Site == 50) |>
      dplyr::select(Reservoir, DateTime,
                    dplyr::starts_with('ThermistorTemp'))
    
    if (current_df$Reservoir[1] == 'BVR') {
      bvr_depths <- find_depths(data_file = current_file,
                                depth_offset = "https://raw.githubusercontent.com/FLARE-forecast/BVRE-data/bvre-platform-data-qaqc/BVR_Depth_offsets.csv",
                                output <- NULL,
                                date_offset <- "2021-04-05",
                                offset_column1<- "Offset_before_05APR21",
                                offset_column2 <- "Offset_after_05APR21") |>
        dplyr::filter(variable == 'ThermistorTemp') |>
        dplyr::select(Reservoir, DateTime, variable, depth_bin, Position)
      
      current_df_1 <- current_df  |>
        tidyr::pivot_longer(cols = starts_with('ThermistorTemp'),
                            names_to = c('variable','Position'),
                            names_sep = '_C_',
                            values_to = 'observation') |>
        dplyr::mutate(date = lubridate::as_date(DateTime),
                      Position = as.numeric(Position)) |>
        na.omit() |>
        dplyr::left_join(bvr_depths,
                         by = c('Position', 'DateTime', 'Reservoir', 'variable')) |>
        dplyr::group_by(date, Reservoir, depth_bin) |>
        dplyr::summarise(observation = mean(observation, na.rm = T),
                         n = dplyr::n(),
                         .groups = 'drop') |>
        dplyr::mutate(observation = ifelse(n < 144/3, NA, observation), # 144 = 24(hrs) * 6(10 minute intervals/hr)
                      Reservoir = 'BVR') |>
        
        dplyr::rename(site_id = Reservoir,
                      datetime = date,
                      depth = depth_bin) |>
        dplyr::select(-n) |> 
        dplyr::mutate(depth = as.character(depth))
    }
    
    # read in differently for FCR
    if (current_df$Reservoir[1] == 'FCR') {
      current_df_1 <- current_df |>
        tidyr::pivot_longer(cols = starts_with('ThermistorTemp'),
                            names_to = 'depth',
                            names_prefix = 'ThermistorTemp_C_',
                            values_to = 'observation') |>
        dplyr::mutate(#Reservoir = ifelse(Reservoir == 'FCR',
                        #                 'fcre',
                         #                ifelse(Reservoir == 'BVR',
                          #                      'bvre', NA)),
                      date = lubridate::as_date(DateTime)) |>
        na.omit() |>
        dplyr::group_by(date, Reservoir, depth) |>
        dplyr::summarise(observation = mean(observation, na.rm = T),
                         n = dplyr::n(),
                         .groups = 'drop') |>
        dplyr::mutate(observation = ifelse(n < 144/2, NA, observation),
                      depth = as.character(depth)) |> # 144 = 24(hrs) * 6(10 minute intervals/hr)
        dplyr::rename(site_id = Reservoir,
                      datetime = date) |>
        dplyr::select(-n)
    }
    message('Current file ready')
  } else {
    current_df_1 <- NULL
    message('No current file')
  }
  
  # read in historical data file
  # EDI
  # infile <- tempfile()
  # try(download.file(historic_file, infile, method="curl"))
  # if (is.na(file.size(infile))) download.file(historic_file,infile,method="auto")
  
  historic_df <- readr::read_csv(historic_file, show_col_types = FALSE) |>
    dplyr::filter(Site == 50) |>
    dplyr::select(Reservoir, DateTime,
                  dplyr::starts_with('ThermistorTemp'))
  
  # Extract depths for BVR
  if (historic_df$Reservoir[1] == 'BVR') {
    bvr_depths <- find_depths(data_file = historic_file,
                              depth_offset = "https://raw.githubusercontent.com/FLARE-forecast/BVRE-data/bvre-platform-data-qaqc/BVR_Depth_offsets.csv",
                              output <- NULL,
                              date_offset <- "2021-04-05",
                              offset_column1<- "Offset_before_05APR21",
                              offset_column2 <- "Offset_after_05APR21") |>
      dplyr::filter(variable == 'ThermistorTemp') |>
      dplyr::select(Reservoir, DateTime, variable, depth_bin, Position)
    
    historic_df_1 <- historic_df |>
      tidyr::pivot_longer(cols = starts_with('ThermistorTemp'),
                          names_to = c('variable','Position'),
                          names_sep = '_C_',
                          values_to = 'observation') |>
      dplyr::mutate(date = lubridate::as_date(DateTime),
                    Position = as.numeric(Position)) |>
      na.omit() |>
      dplyr::left_join(bvr_depths,
                       by = c('Position', 'DateTime', 'Reservoir', 'variable')) |>
      dplyr::group_by(date, Reservoir, depth_bin) |>
      dplyr::summarise(observation = mean(observation, na.rm = T),
                       n = dplyr::n(),
                       .groups = 'drop') |>
      dplyr::mutate(observation = ifelse(n < 144/3, NA, observation), # 144 = 24(hrs) * 6(10 minute intervals/hr)
                    Reservoir = 'BVR') |>
      dplyr::rename(site_id = Reservoir,
                    datetime = date,
                    depth = depth_bin) |>
      dplyr::select(-n) |> 
      dplyr::mutate(depth = as.character(depth))
  }
  
  if (historic_df$Reservoir[1] == 'FCR') {
    historic_df_1 <- historic_df |>
      tidyr::pivot_longer(cols = starts_with('ThermistorTemp'),
                          names_to = 'depth',
                          names_prefix = 'ThermistorTemp_C_',
                          values_to = 'observation') |>
      dplyr::mutate(#Reservoir = ifelse(Reservoir == 'FCR',
                     #                  'fcre',
                      #                 ifelse(Reservoir == 'BVR',
                       #                       'bvre', NA)),
                    date = lubridate::as_date(DateTime)) |>
      dplyr::group_by(date, Reservoir, depth)  |>
      dplyr::summarise(observation = mean(observation, na.rm = T),
                       n = dplyr::n(),
                       .groups = 'drop') |>
      dplyr::mutate(observation = ifelse(n < 6/2, NA, observation)) |> # 6 = 6(10 minute intervals/hr)
      dplyr::rename(site_id = Reservoir,
                    datetime = date)|>
      dplyr::select(-n) |> 
      dplyr::mutate(depth = as.character(depth))
  }
  
  message('EDI file ready')
  
  ## manipulate the data files to match each other
  
  
  ## bind the two files using row.bind()
  final_df <- dplyr::bind_rows(historic_df_1, current_df_1) |>
    dplyr::mutate(variable = 'Temp_C_mean',
                  depth = as.numeric(ifelse(depth == "surface", 0.1, depth))) |>
    rename(depth_m = depth)
  
  final_df <- final_df |>
    mutate(observation = ifelse(is.nan(observation), NA, observation)) |>
    drop_na(depth_m)
  ## Match data to flare targets file
  # Use pivot_longer to create a long-format table
  # for time specific - use midnight UTC values for daily
  # for hourly
  
  ## return dataframe formatted to match FLARE targets
  return(final_df)
}


calc_strat_dates <- function(density_diff = 0.1,
                             temp_profiles) {
  
  ## extract the depths that will be used to calculate the density difference (surface, bottom)
  depths_use <- temp_profiles |>
    na.omit() |> 
    dplyr::group_by(datetime, site_id) |>
    dplyr::summarise(top = min(as.numeric(depth_m, na.rm = T)),
                     bottom = max(as.numeric(depth_m, na.rm = T)),.groups = 'drop') |>
    tidyr::pivot_longer(cols = top:bottom, 
                        names_to = 'location',
                        values_to = 'depth_m')
  
  sites <- distinct(depths_use, site_id) |> pull()
  
  strat_dates <- NULL
  
  for (site in sites) {
    temp_profile_site <- filter(temp_profiles, site_id == site)
    # need a full timeseries
    all_dates <- data.frame(datetime = seq.Date(min(temp_profile_site$datetime), 
                                                max(temp_profile_site$datetime),
                                                'day'))
    density_obs <-
      filter(depths_use, site_id == site) |> 
      inner_join(na.omit(temp_profile_site), by = join_by(datetime, site_id, depth_m)) |> 
      mutate(density = rLakeAnalyzer::water.density(observation)) |> 
      select(datetime, site_id, density, observation, location) |> 
      pivot_wider(values_from = c(density, observation), names_from = location, id_cols = c(datetime, site_id)) |> 
      full_join(all_dates, by = 'datetime') |> 
      mutate(dens_diff = density_bottom - density_top,
             strat = ifelse(abs(dens_diff > 0.1) & observation_top > observation_bottom, 1, 0),
             strat = zoo::na.approx(strat, na.rm = F))
    
    
    # extract the dates of the stratified periods
    #using a loop function to go through each year and do the rle function
    
    strat <- data.frame(year = unique(year(density_obs$datetime)), 
                        length = NA,
                        start = NA,
                        end = NA)
    
    for (i in 1:nrow(strat)) {
      year_use <- strat$year[i]
      
      temp.dens <- density_obs %>%
        filter(year(datetime) == year_use)
      
      if (nrow(temp.dens) >= 300) {
        #run length encoding according to the strat var
        temp.rle <- rle(temp.dens$strat)
        
        #what is the max length for which the value is "norm"
        strat$length[i] <- max(temp.rle$lengths[temp.rle$values==1], 
                               na.rm = T)
        
        #stratification dates
        rle.strat <- data.frame(strat = temp.rle$values, 
                                lengths = temp.rle$lengths)
        
        # Get the end of ech run
        rle.strat$end <- cumsum(rle.strat$lengths)
        # Get the start of each run
        rle.strat$start <- rle.strat$end - rle.strat$lengths + 1
        
        # Sort rows by whehter it is stratified or not
        rle.strat <- rle.strat[order(rle.strat$strat), ]
        
        start.row <- rle.strat$start[which(rle.strat$length == max(rle.strat$lengths)
                                           & rle.strat$strat == 1)] 
        #gets the row with the start date
        #of the run which has the max length and is 1
        
        end.row <- rle.strat$end[which(rle.strat$length == max(rle.strat$lengths)
                                       & rle.strat$strat == 1)] 
        #gets the row with the end date
        #of the run which has the max length and is TRuE
        
        strat$start[which(strat$year == year_use)] <- as.character(temp.dens$datetime[start.row])
        strat$end[which(strat$year == year_use)] <- as.character(temp.dens$datetime[end.row])
        
        strat$site_id <- site
      }
     
    } 
    strat_dates <- bind_rows(strat, strat_dates)
    message(site)
  }
  
  return(na.omit(strat_dates))
}


get_ice_continuous <- function(infile) {
  
  historic_df <- readr::read_csv(infile, show_col_types = F) |> 
    filter(Site == 50)
  
  # apply to each reservoir
  all_ice <- NULL
  for (site_id in unique(historic_df$Reservoir)) {
    
    # dates
    period <- historic_df |>
      dplyr::filter(Reservoir == site_id) |> 
      dplyr::reframe(.by = Reservoir, 
                     first = min(Date),
                     last = max(Date))
    
    # get all the days to fill in with 0
    all_dates <- expand.grid(Date = seq.Date(period$first,
                                             period$last, by = 'day'),
                             Reservoir = site_id)
    
    
    ice_subset <- historic_df |>
      dplyr::filter(Reservoir == site_id) |> 
      right_join(all_dates, by = join_by(Reservoir, Date)) |> 
      arrange(Date) |> 
      mutate(IceOn = zoo::na.locf(IceOn),
             IceOff = zoo::na.locf(IceOff)) |> 
      select(Reservoir, Date, IceOn, IceOff)
    
    all_ice <- bind_rows(all_ice, ice_subset)
  }
  
  return(all_ice)
}
