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
      path <- "Y:/KOD_Survey/CPS2/Trawl Data/Groundfish/" 

  # Load summary catch tables to get full species list
      # basically just using 'catch' files to assign species code because there aren't
      # any haul identifiers in each file -- unless you go haul by haul with all files
      # this doesn't work for any batch processing...
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
      hauls <- read.csv("Y:/KOD_Survey/CPS2/Trawl Data/VA_HAULS.csv") %>% # Vesteraalen
               dplyr::filter(SAMPLED == 1)


# PROCESS FISH DATA ----------------------------------------------------------------------------------------------------------------
  # Join raw_sample_values and raw_sample to get total weight per species per haul
  # Expand unidentified/subsampled flatfish weights to ID'd flatfish
  # Estimate counts based on basket count
  # Expand counts and weights by sampling factor
  # Join with haul info, format final output
      specimen_table <- left_join(raw_sample_values %>% select(CATCH_SAMPLE_VALUE_ID, CATCH_SAMPLE_ID, SAMPLE_VALUE_WEIGHT, SAMPLE_VALUE_WEIGHT_UNITS,
                                                               SAMPLE_VALUE_COUNT, KEEP_FLAG, RECORDING_DEVICE, RECORDER), 
                                  raw_sample %>% select(HAUL_ID, CATCH_SAMPLE_ID, SPECIES_CODE, RECORDING_DEVICE, RECORDER), 
                                  by = c("CATCH_SAMPLE_ID", "RECORDING_DEVICE", "RECORDER")) %>%
                        left_join(., raw_haul %>% select(HAUL_ID, CRUISE_ID, HAUL, STATION, RECORDING_DEVICE)) %>%
                        # remove crab from samples
                        dplyr::filter(!SPECIES_CODE %in% c(68560, 69322, 69400)) %>%
                        # join with catch to get species names
                        left_join(., catch) %>%  
                        # rename "Pacific cod (adult)" --> "Pacific cod"
                        mutate(SPECIES_NAME = ifelse(SPECIES_NAME == "Pacific cod (adult)", "Pacific cod", SPECIES_NAME),
                               SPECIES_CODE = ifelse(SPECIES_CODE == 21722, 21720, SPECIES_CODE)) %>%
                        # filter out weight corrections
                        # not sure how to handle these b/c they have no species identifier
                        filter(if_all(c(HAUL_ID, SPECIES_CODE), complete.cases)) %>%
                        # parse unidentified flatfish weights to respective species
                        # arrowtooth, Kam, flathead, rex, yellowfin, longhead, starry, rock, butter, plaice 
                        group_by(STATION) %>%
                        mutate(FLATFISH = ifelse(SPECIES_CODE %in% c(10210, 10220, 10261, 10285, 10110, 10130,
                                                                     10211, 10112, 10270, 10200), 1, 0),
                               IDENT = sum(SAMPLE_VALUE_WEIGHT[FLATFISH == 1]), 
                               UNIDENT = sum(SAMPLE_VALUE_WEIGHT[SPECIES_NAME == "flatfish unid."]),
                               SAMPLING_FACTOR = (IDENT + UNIDENT)/IDENT, 
                               # set sampling factor to 1 if not a flatfish
                               SAMPLING_FACTOR = ifelse(!SPECIES_CODE %in% c(10210, 10220, 10261, 10285, 10110, 10130,
                                                                             10211, 10112, 10270, 10200), 1, SAMPLING_FACTOR)) %>%
                        # remove unidentified flatfish bin
                        dplyr::filter(!SPECIES_CODE == 10001) %>% 
                        select(-IDENT, -UNIDENT, -FLATFISH) %>%
                        # expand basket count to rest of sample
                        mutate(IND_WEIGHT_CALC = SAMPLE_VALUE_WEIGHT/SAMPLE_VALUE_COUNT) %>%
                        group_by(HAUL, HAUL_ID, CATCH_SAMPLE_ID, RECORDING_DEVICE, STATION, 
                                 SPECIES_CODE, SPECIES_NAME, SAMPLING_FACTOR) %>%
                        nest() %>%
                        mutate(data = map(data, function(df){
                          # pull individual fish weight
                          ind_weight <- unique(na.omit(df$IND_WEIGHT_CALC))
                          # if more than one basket counted, take mean of the calculated individual weight
                          if(length(ind_weight) > 1) ind_weight <- mean(ind_weight)
                          # calculate estimated counts from weight
                          counts <- df %>%
                                    {if(length(ind_weight) == 1) mutate(., COUNT = SAMPLE_VALUE_WEIGHT/ind_weight,
                                                                           IND_WEIGHT = ind_weight) 
                                      else(mutate(., COUNT = NA, IND_WEIGHT = NA))}
                          left_join(df, counts)
                        })) %>%
                        unnest(cols = c(data)) %>%
                        # calculate weight and count totals for each species
                        summarise(IND_WEIGHT = mean(IND_WEIGHT),
                                  WEIGHT = sum(SAMPLE_VALUE_WEIGHT), 
                                  COUNT = sum(COUNT)) %>%
                        left_join(., hauls %>% select(STATION, LAT_DD, LON_DD)) %>%
                        mutate(VESSEL = "Vesteraalen", 
                               CRUISE = 202301) %>%
                        mutate(CALC_WEIGHT = WEIGHT*SAMPLING_FACTOR, 
                               CALC_COUNT = COUNT*SAMPLING_FACTOR) %>%
                        ungroup() %>%
                        select(CRUISE, VESSEL, HAUL, STATION, LAT_DD, LON_DD, #DATE_HAUL, TIME_HAUL, SOAK_TIME, DEPTH_F,
                               SPECIES_CODE, SPECIES_NAME, IND_WEIGHT, WEIGHT, COUNT,
                               SAMPLING_FACTOR, CALC_WEIGHT, CALC_COUNT)

        
    # Save .csv
        write.csv(specimen_table, "./Data/CPS2_2024_trawl_fishcatch.csv", row.names = FALSE)
      
        
    # Reframe and aggregate calculated weights by main bycatch species (to match POT bycatch dataframe)
        fish_sum <- specimen_table %>%
                    # remove non-fish taxa
                    filter(!SPECIES_NAME %in% c("Octopus sp.", "octopus unid.", "Hanasaki king crab, Male, All Sizes")) %>%
                    dplyr::mutate(SPP_LABS = dplyr::case_when(SPECIES_NAME == "Pacific cod" ~ "PacificCod",
                                                              SPECIES_NAME == "Pacific halibut" ~ "Halibut",
                                                              SPECIES_NAME == "great sculpin" ~ "GreatSculpin", 
                                                              SPECIES_NAME == "yellowfin sole" ~ "YellowfinSole", 
                                                              SPECIES_NAME == "walleye pollock" ~ "Pollock", 
                                                              SPECIES_NAME == "starry flounder" ~ "StarryFlounder", 
                                                              SPECIES_NAME == "northern rock sole" ~ "RockSole",
                                                              TRUE ~ "Other")) %>%
                    dplyr::group_by(VESSEL, STATION, LAT_DD, LON_DD, SPP_LABS) %>%
                    dplyr::summarise(WEIGHT = sum(CALC_WEIGHT),
                                     COUNT = sum(CALC_COUNT))
        
        write.csv(fish_sum, "./Data/CPS2_2024_trawl_fish_station_weightcount.csv", row.names = FALSE)
        
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
      
  