# PURPOSE -------------------------------------------------------------------------------------------------------------------
  # 1) To automate processing trawl data specimen tables and catch summaries from Collaborative Pot 
  #    Sampling II (CPS2) 2024 for BBRKC
  # 2) To calculate cpue by haul and maturity/sex category for BBRKC
  # 3) To run error checks on processed specimen and catch summaries
  
  # Authors: Shannon Hennessey, Emily Ryznar, NOAA-AFSC


# INSTALL PACKAGES ----------------------------------------------------------------------------------------------------------
  # install.packages("tidyverse")


# LOAD PACKAGES -------------------------------------------------------------------------------------------------------------
  library(tidyverse)


# LOAD CRAB DATA -----------------------------------------------------------------------------------------------------------------

  # Set trawl data filepath
      path <- "Y:/KOD_Survey/CPS2/Trawl Data/Crab/"  

  # Load summary catch and specimen tables
      # **delete catch and specimen files 0015, 0019, 0021, 0022, 0028, 0030, 0031, 0040, 0053, 0064, 0068, 0088, 0128 because no crab catch
      catch <- list.files(paste0(path, "Catch - FTP/")) %>% 
               purrr::map_df(~read.csv(paste0(path, "Catch - FTP/", .x)))
  
      specimen <- list.files(paste0(path, "Specimen - FTP/")) %>%
                  purrr::map_df(~read.csv(paste0(path, "Specimen - FTP/", .x))) 
  
  # Load raw data for processing below   
      raw_sample <- list.files(paste0(path, "Raw Data - FTP/"), pattern = "_SAMPLE_0") %>% # RECORDS of SAMPLE INFO
                    purrr::map_df(~read.csv(paste0(path, "Raw Data - FTP/", .x))) #E.G. SEX, SPECIES
  
      raw_sample_values <- list.files(paste0(path, "Raw Data - FTP/"), pattern = "_SAMPLE_VALUES") %>% #RECORDS OF # TOSSED
                           purrr::map_df(~read.csv(paste0(path, "Raw Data - FTP/", .x))) 
  
      raw_specimen <- list.files(paste0(path, "Raw Data - FTP/"), pattern = "_SPECIMEN_0") %>% 
                      purrr::map_df(~read.csv(paste0(path, "Raw Data - FTP/", .x))) 
  
      raw_specimen_bio <- list.files(paste0(path, "Raw Data - FTP/"), pattern = "_SPECIMEN_BIOMETRICS") %>% 
                          purrr::map_df(~read.csv(paste0(path, "Raw Data - FTP/", .x))) %>%
                          right_join(., catch %>% select(HAUL, HAUL_ID, RECORDING_DEVICE), 
                                     relationship = "many-to-many") %>%
                          distinct()
  
  # Read in haul data
      hauls <- read.csv("Y:/KOD_Survey/CPS2/Trawl Data/VA_HAULS.csv") %>% # Vesteraalen
               dplyr::filter(SAMPLED == 1)
    
      
# PROCESS CRAB DATA ----------------------------------------------------------------------------------------------------------------
    
      samples <- left_join(raw_sample_values %>% select(HAUL_ID, CATCH_SAMPLE_VALUE_ID, CATCH_SAMPLE_ID, WEIGHT, KEEP, 
                                                        RECORDING_DEVICE, RECORDER), 
                           raw_sample %>% select(HAUL_ID, CATCH_SAMPLE_ID, MIX, SAMPLED_ALL, PARENT_SAMPLE_ID, SPECIES_CODE, SPECIES_NAME,
                                                 SEX, RECORDING_DEVICE, RECORDER), 
                           by = c("HAUL_ID", "CATCH_SAMPLE_ID", "RECORDING_DEVICE", "RECORDER")) %>%
                 left_join(., distinct(catch %>% select(HAUL, HAUL_ID, RECORDING_DEVICE))) 
    
  # if 2 baskets of same spp/sex in a CATCH_SAMPLE_ID, sum weights and replace for all obs....
      samples2 <- samples %>%  
                  dplyr::group_by(HAUL, HAUL_ID, CATCH_SAMPLE_ID, SEX, SPECIES_CODE, RECORDING_DEVICE) %>%
                  dplyr::summarise(WEIGHT = sum(WEIGHT)) %>%
                  # distinct() %>%
                  dplyr::left_join(samples %>% select(-c(WEIGHT, CATCH_SAMPLE_VALUE_ID))) %>%
                  distinct() %>%
                  dplyr::select(HAUL, HAUL_ID, CATCH_SAMPLE_ID, MIX, SAMPLED_ALL, PARENT_SAMPLE_ID,
                                SPECIES_CODE, SPECIES_NAME, SEX, WEIGHT, KEEP, RECORDING_DEVICE) #
    
  # Expand specimen biometric table, join to raw_specimen table to get catch sample ID, join with samples file to get 
  # number tossed
      specimen_sum <- raw_specimen_bio %>%
                      dplyr::select(HAUL, HAUL_ID, SPECIMEN_ID, BIOMETRIC_NAME, VALUE, RECORDING_DEVICE) %>%
                      pivot_wider(., id_cols = c(HAUL, HAUL_ID, SPECIMEN_ID, RECORDING_DEVICE), 
                                  names_from = "BIOMETRIC_NAME", values_from = "VALUE") %>%
                      dplyr::rename(SHELL_CONDITION = CRAB_SHELL_CONDITION, EGG_COLOR = CRAB_EGG_COLOR,
                                    EGG_CONDITION = CRAB_EGG_CONDITION, CLUTCH_SIZE = CRAB_EGG_CLUTCH_SIZE,
                                    LENGTH = CARAPACE_LENGTH, WIDTH = CARAPACE_WIDTH) %>%
                      left_join(., raw_specimen) %>%
                      dplyr::select(-c(WEIGHT)) %>% # remove weight
                      right_join(samples2, ., by = c("HAUL", "HAUL_ID", "SEX", "RECORDING_DEVICE", "SPECIES_CODE", "CATCH_SAMPLE_ID"), 
                                 relationship = "many-to-many") %>%
                      distinct()

  # Calculate sampling factor from specimen summary table, join back with specimen_sum file to 
  # get specimen information, join with catch file to get vessel and station #s, join with hauls
  # file to get lat/lon, set/haul date and time for each haul (with positive catch)
      specimen_table <- specimen_sum %>%
                        dplyr::group_by(HAUL, HAUL_ID, CATCH_SAMPLE_ID, SEX, RECORDING_DEVICE) %>%
                        # # dplyr::reframe(KEPT = n(),
                        # #                TOSSED = TOSSED,
                        # #                SAMPLING_FACTOR = (KEPT + TOSSED)/KEPT) %>%
                        # distinct() %>%
                        # dplyr::right_join(specimen_sum, multiple = "all") %>%
                        dplyr::right_join(catch %>% select(-c(WEIGHT)), ., by = c("HAUL", "HAUL_ID", "SPECIES_CODE", "RECORDING_DEVICE")) %>%
                        dplyr::left_join(specimen) %>%
                        dplyr::right_join(hauls, ., by = c("STATION"), relationship = "many-to-many") %>%
                        dplyr::mutate(VESSEL = ifelse(VESSEL == 176, "Vesteraalen", ""),
                                      SAMPLING_FACTOR = 1) %>% ## NEED TO DOUBLE CHECK SF (make sure no >1 ever)..this doesn't join right
                        dplyr::select(CRUISE, VESSEL, HAUL, STATION, LAT_DD, LON_DD, #DATE_HAUL, TIME_HAUL, SOAK_TIME, DEPTH_F,
                                      SPECIES_CODE, SPECIES_NAME, SEX, LENGTH, WIDTH, SAMPLING_FACTOR, SHELL_CONDITION, EGG_COLOR, EGG_CONDITION, 
                                      CLUTCH_SIZE, WEIGHT, DISEASE_CODE, DISEASE_DORSAL, DISEASE_VENTRAL, DISEASE_LEGS,  
                                      CHELA_HEIGHT, MERUS_LENGTH, COMMENTS, NOTES)
    

  # Process specimen table with all haul data, save
      specimen_table %>% dplyr::filter(SPECIES_CODE == 69322) %>% # filter just RKC
        write.csv("./Data/CPS2_2024_Processed_Trawl_Specimen_Data.csv", row.names = FALSE)
    
  # Update catch summary table with new crab #s from sampling factor
      catch_summary <- specimen_table %>%
                       dplyr::group_by(CRUISE, VESSEL, HAUL, STATION, SPECIES_CODE) %>%
                       dplyr::reframe(NUMBER_CRAB = sum(SAMPLING_FACTOR)) %>%
                       dplyr::right_join(catch %>% dplyr::rename(N_ENTRIES = NUMBER_CRAB) %>%
                                                   dplyr::mutate(VESSEL = ifelse(VESSEL == 176, "Vesteraalen", ""))) %>%
                       dplyr::select(CRUISE, VESSEL, HAUL, STATION, SPECIES_CODE, NUMBER_CRAB, N_ENTRIES) %>%
                       na.omit() # bad/gear testing potlifts will have NA for SPN and # crab
      
  # Print lines where N_CRAB =/= N_ENTRIES
      catch_summary %>% filter(NUMBER_CRAB != N_ENTRIES)
      # R6 - 1 less tanner than should be?
                   

# CALCULATE BBRKC CPUE -------------------------------------------------------------------------------------------------------     
    
  # Make non-overlapping maturity/sex and legal/sublegal categories, bind together
      maturity <- specimen_table %>%
                  dplyr::filter(SPECIES_CODE == 69322) %>% # filter just RKC
                  dplyr::mutate(MAT_SEX = dplyr::case_when((SPECIES_CODE == 69322 & SEX == 1 & LENGTH >= 120) ~ "Mature male",
                                                           (SPECIES_CODE == 69322 & SEX == 1 & LENGTH < 120) ~ "Immature male",
                                                           (SPECIES_CODE == 69322 & SEX == 2 & CLUTCH_SIZE >= 1) ~ "Mature female",
                                                           (SPECIES_CODE == 69322 & SEX == 2 & CLUTCH_SIZE == 0) ~ "Immature female"))
    
      legal <- specimen_table %>%
               dplyr::filter(SPECIES_CODE == 69322) %>% # filter just RKC
               dplyr::mutate(MAT_SEX = dplyr::case_when((SPECIES_CODE == 69322 & SEX == 1 & LENGTH >= 135) ~ "Legal male",
                                                        (SPECIES_CODE == 69322 & SEX == 1 & LENGTH < 135) ~ "Sublegal male")) %>%
               dplyr::filter(is.na(MAT_SEX) == "FALSE")
           
      mat_spec <- rbind(maturity, legal) #bind
    
  
  # Summarize counts by sex/maturity category per station
      mat_spec_sum <- mat_spec %>%
                      dplyr::group_by(VESSEL, HAUL, STATION, LAT_DD, LON_DD, MAT_SEX) %>%
                      dplyr::reframe(COUNT = sum(SAMPLING_FACTOR))
    
      positive_haul_cpue <- mat_spec_sum
  

  # Expand potlifts file to all mat-sex categories and potlifts, join to positive catch file to get zeros 
      mat_sex_combos <- c("Mature male", "Immature male", "Mature female", "Immature female", "Legal male", "Sublegal male")
    
      trawl_cpue <- positive_haul_cpue %>%
                    dplyr::right_join(expand_grid(MAT_SEX = mat_sex_combos,
                                                  hauls)) %>%
                    replace_na(list(COUNT = 0, VESSEL = "Vesteraalen")) %>%
                    dplyr::select(VESSEL, HAUL, STATION, LAT_DD, LON_DD,
                                  MAT_SEX, COUNT) #, CATCH_PER_HOUR) 

  # Save csv
      write.csv(trawl_cpue, "./Data/CPS2_2024_trawlcatch.csv", row.names = FALSE)
    
    
# CALCULATE CRAB BYCATCH ----------------------------------------------------------------------------------------------------
    
  # Summarize counts by spp/sex per station
      crab_sum <- specimen_table %>%
                  dplyr::group_by(VESSEL, HAUL, STATION, LAT_DD, LON_DD, SEX, SPECIES_NAME, SPECIES_CODE) %>%
                  dplyr::reframe(COUNT = sum(SAMPLING_FACTOR)) %>%
                  dplyr::filter(!SPECIES_CODE %in% c(69322)) %>%
                  dplyr::mutate(SPP_SEX = dplyr::case_when((SPECIES_CODE == 69400 & (SEX == 1 | SEX == 2)) ~ "HairCrab",
                                                           (SPECIES_CODE == 68560 & SEX == 1) ~ "MaleTanner",
                                                           (SPECIES_CODE == 68560 & SEX == 2) ~ "FemaleTanner")) %>%
                  dplyr::group_by(VESSEL, HAUL, STATION, LAT_DD, LON_DD, SPP_SEX) %>%
                  dplyr::summarise(COUNT = sum(COUNT))
                
 
      spp_sex_combos <- c("MaleTanner", "FemaleTanner", "MaleSnow", "FemaleSnow", "MaleHybrid", "FemaleHybrid", "HairCrab")
    
      crab_sum %>%
        dplyr::right_join(expand_grid(SPP_SEX = spp_sex_combos,
                                      hauls)) %>%
        replace_na(list(COUNT = 0, VESSEL = "Vesteraalen")) %>%
        dplyr::select(VESSEL, STATION,LAT_DD, LON_DD, SPP_SEX, COUNT) %>%
        pivot_wider(., id_cols = c(VESSEL, STATION, LAT_DD, LON_DD,),
                    names_from = "SPP_SEX", values_from = "COUNT") %>%
        # sum across male and female crabs to get species-level counts
        dplyr::mutate(Tanner = MaleTanner + FemaleTanner,
                      Snow = MaleSnow + FemaleSnow,
                      Hybrid = MaleHybrid + FemaleHybrid) %>%
        write.csv("./Data/CPS2_2024_trawl_crab_bycatch.csv", row.names = FALSE)
    
      
# ERROR CHECKING -----------------------------------------------------------------------      
      
  # Load error checking function
      source("./Scripts/error_check.R")
      
  # Run error checks
      error_chk(method = "TRAWL", 
                specimen_table = specimen_table[which(specimen_table$SPECIES_CODE == 69322),], 
                catch_summary = catch_summary[which(catch_summary$SPECIES_CODE == 69322),], 
                potlifts = NULL, 
                haul = hauls, 
                cpue = trawl_cpue)
      
    
