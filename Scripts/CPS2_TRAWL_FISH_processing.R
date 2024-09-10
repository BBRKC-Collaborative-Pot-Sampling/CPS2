# PURPOSE -------------------------------------------------------------------------------------------------------------------
  # 1) To automate processing trawl data specimen tables and catch summaries from Collaborative Pot 
  #    Sampling II (CPS2) 2024 for groundfish species
  
  # Authors: Shannon Hennessey, Emily Ryznar, NOAA-AFSC


# INSTALL PACKAGES ----------------------------------------------------------------------------------------------------------
    #install.packages("tidyverse")


# LOAD PACKAGES -------------------------------------------------------------------------------------------------------------
    library(tidyverse)


# LOAD FISH DATA -----------------------------------------------------------------------------------------------------------------

  # Set trawl data filepath
      path <- "Y:/KOD_Survey/CPS2/Data/Trawl Data/Groundfish/" 

  # Load summary catch tables to get full species list
      catch <- list.files(paste0(path, "Catch - FTP/")) %>%
               purrr::map_df(~read.csv(paste0(path, "Catch - FTP/", .x))) %>%
               dplyr::select(SPECIES_CODE, SPECIES_NAME) %>%
               distinct()

  # Load raw data for processing below
      raw_sample <- list.files(paste0(path, "Raw - FTP/"), pattern = "_SAMPLE_0") %>% # RECORDS of SAMPLE INFO
                    purrr::map_df(~read.csv(paste0(path, "Raw - FTP/", .x))) #E.G. SEX, SPECIES

      raw_sample_values <- list.files(paste0(path, "Raw - FTP/"), pattern = "_CATCH_VALUE") %>% # RECORDS OF WEIGHT TOSSED
                           purrr::map_df(~read.csv(paste0(path, "Raw - FTP/", .x))) 

      raw_haul <- list.files(paste0(path, "Raw - FTP/"), pattern = "_HAUL_0") %>%
                  purrr::map_df(~ read.csv(paste0(path, "Raw - FTP/", .x))) %>%
                  select(-c(WEIGHT, WEIGHT_UNITS))
      
      
  # Read in haul data
      hauls <- read.csv("Y:/KOD_Survey/CPS2/Data/Trawl Data/VA_HAULS.csv") %>% # Vesteraalen
               dplyr::filter(SAMPLED == 1)


# PROCESS FISH DATA ----------------------------------------------------------------------------------------------------------------
  # Join raw_sample_values and raw_sample to get total weight per species per haul
      samples <- left_join(raw_sample_values %>% select(CATCH_SAMPLE_VALUE_ID, CATCH_SAMPLE_ID, SAMPLE_VALUE_WEIGHT, SAMPLE_VALUE_WEIGHT_UNITS,
                                                          SAMPLE_VALUE_COUNT, KEEP_FLAG, RECORDING_DEVICE, RECORDER), 
                             raw_sample %>% select(HAUL_ID, CATCH_SAMPLE_ID, SPECIES_CODE, RECORDING_DEVICE, RECORDER), 
                             by = c("CATCH_SAMPLE_ID", "RECORDING_DEVICE", "RECORDER")) %>%
                   left_join(., raw_haul %>% select(HAUL_ID, CRUISE_ID, HAUL, STATION, RECORDING_DEVICE)) %>%
                   # remove crab from samples
                   dplyr::filter(!SPECIES_CODE %in% c(68560, 69322, 69400)) %>%
                   # join with catch to get species names
                   left_join(., catch) %>%  
                   # filter out weight corrections
                   filter(if_all(c(HAUL_ID, SPECIES_CODE), complete.cases)) %>%
                   # calculate weight totals for each species
                   dplyr::group_by(HAUL, HAUL_ID, CATCH_SAMPLE_ID, SPECIES_CODE, SPECIES_NAME, RECORDING_DEVICE) %>%
                   dplyr::summarise(WEIGHT = sum(SAMPLE_VALUE_WEIGHT))
      
                   # *Might have to re-incorporate this bit of code to deal with counts??
                   # distinct() %>%
                   # dplyr::left_join(samples %>% select(-c(SAMPLE_VALUE_WEIGHT, CATCH_SAMPLE_VALUE_ID))) %>%
                   # distinct() %>%
                   # dplyr::select(HAUL, HAUL_ID, CATCH_SAMPLE_ID, WEIGHT, SPECIES_CODE, SPECIES_NAME, RECORDING_DEVICE) %>%
                    

  # Calculate sampling factor from samples, join back with samples file to 
  # get specimen information, join with raw_haul file to get haul and station #s, join with hauls
  # file to get lat/lon, set/haul date and time for each haul (with positive catch)
      specimen_table <- samples %>%
                        # dplyr::group_by(HAUL, HAUL_ID, CATCH_SAMPLE_ID, RECORDING_DEVICE) %>%
                        dplyr::reframe(FLATFISH = ifelse(SPECIES_CODE %in% c(10210, 10220, 10261, 10285, 10110, 10130,
                                                                             10269, 10211, 10112, 10270, 10200), 1, 0)) %>%
                        right_join(., samples) %>%
                        dplyr::group_by(HAUL, HAUL_ID, RECORDING_DEVICE) %>%
                        dplyr::mutate(IDENT = sum(WEIGHT[FLATFISH == 1]), across(),
                                      UNIDENT = sum(WEIGHT[SPECIES_NAME == "flatfish unid."]), across(),
                                      SAMPLING_FACTOR = (IDENT + UNIDENT)/IDENT) %>%
                        # remove unidentified flatfish bin
                        dplyr::filter(!SPECIES_CODE == 10001) %>% 
                        # set sampling factor to 1 if not a flatfish
                        dplyr::mutate(SAMPLING_FACTOR = ifelse(!SPECIES_CODE %in% c(10210, 10220, 10261, 10285, 10110, 10130,
                                                                                   10269, 10211, 10112, 10270, 10200), 1, SAMPLING_FACTOR)) %>%
                        dplyr::left_join(., raw_haul, by = c("HAUL", "HAUL_ID")) %>%
                        dplyr::left_join(., hauls %>% select(STATION, LAT_DD, LON_DD)) %>%
                        dplyr::mutate(VESSEL = "Vesteraalen", 
                                      CRUISE = 202301) %>% 
                        dplyr::select(CRUISE, VESSEL, HAUL, STATION, LAT_DD, LON_DD, #DATE_HAUL, TIME_HAUL, SOAK_TIME, DEPTH_F,
                                      SPECIES_CODE, SPECIES_NAME, WEIGHT, SAMPLING_FACTOR) %>%
                        dplyr::mutate(CALC_WEIGHT = WEIGHT*SAMPLING_FACTOR)
      
    # Save .csv
        write.csv(specimen_table, "./Data/CPS2_2024_trawl_fishcatch.csv", row.names = FALSE)
      
        
    # Reframe and aggregate calculated weights by main bycatch species (to match POT bycatch dataframe)
        fish_sum <- specimen_table %>%
                    dplyr::mutate(SPP_LABS = dplyr::case_when(SPECIES_NAME == "Pacific cod" ~ "PacificCod",
                                                              SPECIES_NAME == "Pacific halibut" ~ "Halibut",
                                                              SPECIES_NAME == "great sculpin" ~ "GreatSculpin", 
                                                              SPECIES_NAME == "yellowfin sole" ~ "YellowfinSole", 
                                                              SPECIES_NAME == "walleye pollock" ~ "Pollock", 
                                                              SPECIES_NAME == "starry flounder" ~ "StarryFlounder", 
                                                              SPECIES_NAME == "northern rock sole" ~ "RockSole",
                                                              TRUE ~ "Other")) %>%
                    dplyr::group_by(VESSEL, STATION, LAT_DD, LON_DD, SPP_LABS) %>%
                    dplyr::summarise(WEIGHT = sum(CALC_WEIGHT))
        
        bycatch_spp <- unique(fish_sum$SPP_LABS)
        
        fish_sum %>%
          dplyr::right_join(expand_grid(SPP_LABS = bycatch_spp,
                                        hauls)) %>%
          replace_na(list(WEIGHT = 0, VESSEL = "Vesteraalen")) %>%
          dplyr::select(VESSEL, STATION, LAT_DD, LON_DD, SPP_LABS, WEIGHT) %>%
          pivot_wider(., id_cols = c(VESSEL, STATION, LAT_DD, LON_DD,),
                      names_from = "SPP_LABS", values_from = "WEIGHT") %>%
          select(VESSEL, STATION, LAT_DD, LON_DD, PacificCod, Halibut, 
                 GreatSculpin, YellowfinSole, Pollock, StarryFlounder, RockSole, Other) %>%
          write.csv("./Data/CPS2_2024_trawl_fish_bycatch.csv", row.names = FALSE)
      
  