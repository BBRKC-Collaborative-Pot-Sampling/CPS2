
# LOAD DATA -----------------------------------------------------------------------------------------------------------------

  # Load mapping layers
      source("./Scripts/map_setup.R")
  
# Read in POT and TRAWL data -----
      pot_cpue <- read.csv("./Data/CPS2_2024_potcatch.csv") %>%
                  dplyr::rename(HAUL = SPN,
                                STATION = POT_ID) %>%
                  filter(!nchar(STATION) > 3) # filter out CAM, COFFIN, and BAIT POT_IDs
        
      trawl_cpue <- read.csv("./Data/CPS2_2024_trawlcatch.csv")

  # Combine pot and trawl counts into one dataframe
      cpue <- trawl_cpue %>% 
              rbind(., pot_cpue %>% dplyr::select(-c(BUOY, DATE_SET, TIME_SET, DATE_HAUL, TIME_HAUL, SOAK_TIME, CATCH_PER_HOUR)))
      
    # Combine sex/mat categories for total BBRKC
      tot_cpue <- cpue %>%
                  dplyr::filter(MAT_SEX %in% c("Mature male", "Immature male", "Mature female", "Immature female")) %>%
                  dplyr::group_by(VESSEL, STATION, LAT_DD, LON_DD) %>%
                  dplyr::summarise(COUNT = sum(COUNT))
      
  # Transform catch data to correct crs
      cpue_mapdat <- cpue %>%
                     sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
                     sf::st_transform(crs = map.crs)
                  
      tot_cpue_mapdat <- tot_cpue %>%
                         sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs)
      

# Read in bycatch data -----
    # Read in POT and TRAWL data
      pot_bycatch <- read.csv("./Data/CPS2_2024_POT_bycatch.csv") %>%
                     dplyr::rename(HAUL = SPN,
                                   STATION = POT_ID) %>%
                     dplyr::mutate(RockSole = 0) %>%
                     filter(!nchar(STATION) > 3) %>% # filter out CAM, COFFIN, and BAIT POT_IDs
                     dplyr::select(VESSEL, STATION, LAT_DD, LON_DD, MaleTanner, FemaleTanner, Tanner,
                                   MaleSnow, FemaleSnow, Snow, MaleHybrid, FemaleHybrid, Hybrid, HairCrab,
                                   PacificCod, Halibut, GreatSculpin, YellowfinSole, Pollock, StarryFlounder,
                                   RockSole, Other)
      
      trawl_fish_bycatch <- read.csv("./Data/CPS2_2024_TRAWL_FISH_bycatch.csv") 
      trawl_crab_bycatch <- read.csv("./Data/CPS2_2024_TRAWL_CRAB_bycatch.csv") 
      
      
  # Combine crab and fish trawl bycatch into one dataframe
      trawl_bycatch <- trawl_fish_bycatch %>% 
                       left_join(., trawl_crab_bycatch) %>%
                       dplyr::select(VESSEL, STATION, LAT_DD, LON_DD, MaleTanner, FemaleTanner, Tanner, 
                                     MaleSnow, FemaleSnow, Snow, MaleHybrid, FemaleHybrid, Hybrid, HairCrab, 
                                     PacificCod, Halibut, GreatSculpin, YellowfinSole, Pollock, StarryFlounder, 
                                     RockSole, Other)
       
  # Combine pot and trawl bycatch into one dataframe
      bycatch <- trawl_bycatch %>% 
                 rbind(., pot_bycatch)
      
  # Transform bycatch data to correct crs
      bycatch <- bycatch %>%
                 sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
                 sf::st_transform(crs = map.crs)
  
      
# Read in specimen data ------
  # Read in POT and TRAWL data
      pot_specimen <- read.csv("./Data/CPS2_2024_Processed_Pot_Specimen_Data.csv") %>%
                      dplyr::rename(HAUL = SPN,
                                    STATION = POT_ID) %>%
                      filter(!nchar(STATION) > 3) # filter out CAM, COFFIN, and BAIT POT_IDs
      
      trawl_specimen <- read.csv("./Data/CPS2_2024_Processed_Trawl_Specimen_Data.csv")
      
  # Combine pot and trawl specimen data into one dataframe
      specimen <- trawl_specimen %>% dplyr::select(-c(SPECIES_NAME)) %>%
                  rbind(., pot_specimen %>% dplyr::select(-c(BUOY, DATE_HAUL, TIME_HAUL, SOAK_TIME, DEPTH_F))) %>% 
                  dplyr::mutate(MAT_SEX_SC = dplyr::case_when((SPECIES_CODE == 69322 & SEX == 1 & LENGTH >= 120 & SHELL_CONDITION == 1) ~ "Mature male SC1",
                                                              (SPECIES_CODE == 69322 & SEX == 1 & LENGTH >= 120 & SHELL_CONDITION == 2) ~ "Mature male SC2",
                                                              (SPECIES_CODE == 69322 & SEX == 1 & LENGTH >= 120 & SHELL_CONDITION == 3) ~ "Mature male SC3",
                                                              (SPECIES_CODE == 69322 & SEX == 2 & CLUTCH_SIZE >= 1 & SHELL_CONDITION == 1) ~ "Mature female SC1",
                                                              (SPECIES_CODE == 69322 & SEX == 2 & CLUTCH_SIZE >= 1 & SHELL_CONDITION == 2) ~ "Mature female SC2",
                                                              (SPECIES_CODE == 69322 & SEX == 2 & CLUTCH_SIZE >= 1 & SHELL_CONDITION == 3) ~ "Mature female SC3")) %>%
                  dplyr::filter(!is.na(MAT_SEX_SC))
      
      # now...need to wrangle specific shell conditions and make counts...
      # Calculate COUNT and CATCH_PER_HOUR per pot, CHANGE VESSEL TO ACTUAL NAME
      positive_counts <- specimen %>%
                            dplyr::group_by(VESSEL, STATION, MAT_SEX_SC) %>%
                            dplyr::reframe(COUNT = sum(SAMPLING_FACTOR))
      
      mat_sex_sc_combos <- c("Mature male SC1", "Mature male SC2", "Mature male SC3", 
                             "Mature female SC1", "Mature female SC2", "Mature female SC3")
      
      sc_counts <- positive_counts %>%
                   dplyr::right_join(., expand_grid(MAT_SEX_SC = mat_sex_sc_combos,
                                                    tot_cpue %>% select(-c(COUNT)))) %>%
                   replace_na(list(COUNT = 0)) 
      
      
  # Transform specimen data to correct crs
      sc_counts <- sc_counts %>%
                    sf::st_as_sf(coords = c(x = "LON_DD", y = "LAT_DD"), crs = sf::st_crs(4326)) %>%
                    sf::st_transform(crs = map.crs)
      
      
# INTERPOLATE TEMPERATURE: Temperature map with bathy, CPS1 outline, and RKCSA -------
      temp <- read.csv("./Data/Temperature Data/CPS2_TempData_Combined.csv")
      # ggplot(temp) + geom_point(aes(x = longitude, y = latitude, color = temperature_C))
  
    # Set up interpolation raster
      sp_interp.raster <- raster::raster(st_as_sf(CPS1_bound), res = 1000)
      
      BB_rast <- st_read(survey_gdb, layer = "BristolBaySurveyStrata") %>%
                 st_transform(map.crs)
  
  # Transform data for interpolation
      in.crs <- "+proj=longlat +datum=NAD83"
      sp_interp.df <- unique(temp)
      sp::coordinates(sp_interp.df) <- c(x = "longitude", y = "latitude")
      sp::proj4string(sp_interp.df) <- sp::CRS(in.crs)
      sp_interp.df <- sp::spTransform(sp_interp.df, sp::CRS(map.crs))
  
  # Set up a new IDW for ordinary kriging 
      idw_vgm_fit <- gstat::gstat(formula = temperature_C ~ 1, 
                                  locations = sp_interp.df, 
                                  nmax = Inf)
  
  # Ordinary Kriging: VGM
      ste.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit), 
                                        vgm(c("Bes")))
      
      ste_fit <- gstat::gstat(formula = temperature_C ~ 1, 
                              locations = sp_interp.df, 
                              model = ste.vgfit, 
                              nmax = Inf)
      
      ste.predict <- predict(ste_fit, as(sp_interp.raster, "SpatialGrid"))
      
  # write unmasked surfaces to raster, stacked by year
      temp_rast <- ste.predict %>%
                   raster::raster(.) %>%
                   mask(CPS1_bound)
  
  # extract interpolated data from raster to data frame
      coords <- sp::coordinates(temp_rast)
  
      temp_df <- na.omit(data.frame(coords, temperature = temp_rast@data@values))
      
      temp_breaks <- c(-Inf, seq(-1,4,1), Inf)
      viridis_option <- "H" # viridis turbo palette
      n_temp_breaks <- length(temp_breaks)-1
  


# GENERAL MAPS --------------------------------------------------------------------------------------------------------
  # Temperature -----
  # Plot
      temp_map_continuous <- ggplot() +
                             geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature))+
                             scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                                                guide = guide_colorbar(title.position = "top"))+
                             geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
                             geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1)+
                             #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                             geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha = 0.9, linewidth = 1) +
                             geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha = 0.5, linewidth = 1) +
                             # geom_sf(data = temploggers, color = "black", size = 2)+
                             geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
                             scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
                             scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
                             labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Temperature")+
                             scale_color_manual(values = c("black", "red"), 
                                                labels = c("CPS survey boundary", "Red King Crab Savings Area"),
                                                name = "") +
                             coord_sf(xlim = plot.boundary$x,
                                      ylim = plot.boundary$y) +
                             #ggplot2::scale_fill_manual(name = "Temperature (?C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
                             #labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                             #"4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
                             geom_sf_text(sf::st_as_sf(data.frame(lab = c("50m", "100m"), 
                                                                  x = c(-161.5, -165), y = c(58.3, 56.1)),
                                                       coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                                            sf::st_transform(crs = map.crs),
                                          mapping = aes(label = lab))+
                             guides(color = guide_legend(nrow = 2)) +
                             theme_bw() +
                             theme(axis.title = element_blank(),
                                   axis.text = element_text(size = 10),
                                   legend.text = element_text(size = 10),
                                   legend.title = element_text(size = 10),
                                   legend.position = "bottom",
                                   legend.direction = "horizontal",
                                   plot.title = element_text(face = "bold", size = 15),
                                   plot.subtitle = element_text(size = 12))
      
      ggsave(plot = temp_map_continuous, "./Figures/tempcont.png", height = 7, width = 10, units = "in")


# ABUNDANCE MAPS --------------------------------------------------------------------------------------------------------
  # Specify palette
      pot_pal <- viridis::mako(10) 
      trawl_pal <- viridis::rocket(10)
      
  # Total BBRKC with temperature -------
    # Plot total counts
      total.map <- ggplot() +
                    # plot temperature
                    geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
                    scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                                       guide = guide_colorbar(title.position = "top", 
                                                              discrete = FALSE)) +
                    guides(alpha = "none") +
                    new_scale("fill") +
                    # new_scale("shape") +
                    # new_scale("size") +
        
                    # plot mapping layers
                    geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
                    # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                    geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
                    geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha = 0.9, linewidth = 1) +
                    geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
                    scale_color_manual(values = c("black", "red"),
                                       labels = c("CPS survey boundary", "Red King Crab Savings Area"),
                                       name = "") +
                    guides(color = guide_legend(nrow = 2)) +
                    geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +

                    # add POT points
                    geom_sf(data = tot_cpue_mapdat %>% filter(VESSEL %in% c("Arctic Lady", "Seabrooke")),
                            mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), 
                            alpha = 0.5, colour = "black") +
                    scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
                    scale_size_continuous(range = c(2, 10), limits = c(0, 175), breaks = seq(0, 175, by = 25), guide = "none") + 
                    scale_fill_gradientn(limits = c(0, 175), breaks = seq(0, 175, by = 25),
                                         colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
                    
                    # POT legend
                    guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                               override.aes = list(shape = c(4, rep(21, 7)),
                                                                   size = seq(2, 10, by = ((10-2)/(8-1))))),
                           size = guide_legend(),
                           color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
                    
                    # start a new scale
                    new_scale("fill") +
                    new_scale("shape") +
                    new_scale("size") +
                    
                    # add TRAWL points
                    geom_sf(data = tot_cpue_mapdat %>% filter(VESSEL %in% c("Vesteraalen")),
                            mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), 
                            alpha = 0.5, colour = "black") +
                    scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none") +
                    scale_size_continuous(range = c(2, 10), limits = c(0, 154),
                                          breaks = c(seq(0, 25, by = 5), 98, 154)) +
                    scale_fill_gradientn(breaks = c(seq(0, 25, by = 5), 98, 154),
                                         limits = c(0, 154), 
                                         colors = c("gray", rev(trawl_pal[5:length(trawl_pal)]))) +
                    scale_x_continuous(breaks = map_layers$lon.breaks) +
                    scale_y_continuous(breaks = map_layers$lat.breaks) +
                    # scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
                    # scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
                    labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Total BBRKC") + 
                    
                    # TRAWL legend
                    guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                               override.aes = list(shape = c(8, rep(24, 7)))),
                           fill = guide_legend(title = "Trawl Count", title.position = "top"),
                           color = guide_legend(nrow = 2)) +
                    
                    # crop spatial extent, add bathymetry, set theme
                    coord_sf(xlim = plot.boundary$x,
                             ylim = plot.boundary$y) +
                    geom_sf_text(sf::st_as_sf(data.frame(lab = c("50m", "100m"), 
                                                         x = c(-161.5, -165), y = c(58.3, 56.1)),
                                              coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                                   sf::st_transform(crs = map.crs),
                                 mapping = aes(label = lab)) +
                    theme_bw() +
                    theme(axis.title = element_blank(),
                          axis.text = element_text(size = 10),
                          legend.text = element_text(size = 10),
                          legend.title = element_text(size = 10),
                          legend.position = "right",
                          legend.direction = "horizontal",
                          plot.title = element_text(face = "bold", size = 15),
                          plot.subtitle = element_text(size = 12))
      
      ggsave(plot = total.map, "./Figures/BBRKC_TOTAL_temp.png", height = 7, width = 10, units = "in")
      
      
  # BBRKC mature male with temperature -------
      BBRKC.matmale.map <- ggplot() +
                          # plot temperature
                          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
                          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                                             guide = guide_colorbar(title.position = "top", 
                                                                    discrete = FALSE)) +
                          guides(alpha = "none") +
                          new_scale("fill") +

                          # plot mapping layers
                          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
                          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
                          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
                          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
                          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
                          scale_color_manual(values = c("black", "red"), 
                                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                                             name = "") +
                          guides(color = guide_legend(nrow = 2)) +
                          
                          # add POT points
                          geom_sf(data = filter(cpue_mapdat, MAT_SEX == "Mature male",
                                                VESSEL %in% c("Arctic Lady", "Seabrooke")),
                                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black") +
                          scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
                          scale_size_continuous(range = c(2, 10), limits = c(0, 125), guide = "none",
                                                breaks = seq(0, 125, by = 25))+ 
                          scale_fill_gradientn(breaks = seq(0, 125, by = 25),
                                               limits = c(0, 125), 
                                               colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
                          
                          # POT legend
                          guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                                     override.aes = list(shape = c(4, rep(21, 5)),
                                                                         size = seq(2, 10, by = ((10-2)/(6-1))))),
                                 size = guide_legend(),
                                 color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
                          
                          # start a new scale
                          new_scale_fill() +
                          new_scale("shape") +
                          new_scale("size") +
                          
                          # add TRAWL points
                          geom_sf(data = filter(cpue_mapdat, MAT_SEX == "Mature male",
                                                VESSEL %in% c("Vesteraalen")),
                                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
                          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
                          scale_size_continuous(range = c(2, 10), limits = c(0, 25), 
                                                breaks = c(seq(0, 20, by = 5), 25))+ 
                          scale_fill_gradientn(breaks = c(seq(0, 20, by = 5), 25),
                                               limits = c(0, 25), 
                                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
                          # TRAWL legend
                          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                                     override.aes = list(shape = c(8, rep(24,5)))),
                                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                                 color = guide_legend(nrow = 2)) +
                          
                          scale_x_continuous(breaks = map_layers$lon.breaks)+
                          scale_y_continuous(breaks = map_layers$lat.breaks)+
                          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature male (>= 120mm)")+
                          coord_sf(xlim = plot.boundary$x,
                                   ylim = plot.boundary$y) +
                          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                                         sf::st_transform(crs = map.crs),
                                       mapping = aes(label = lab))+
                          theme_bw() +
                          theme(axis.title = element_blank(),
                                axis.text = element_text(size = 10),
                                legend.text = element_text(size = 10),
                                legend.title = element_text(size = 10),
                                legend.position = "right",
                                legend.direction = "horizontal",
                                plot.title = element_text(face = "bold", size = 15),
                                plot.subtitle = element_text(size = 12))

        ggsave(plot = BBRKC.matmale.map, "./Figures/BBRKC_matmale_temp.png", height = 7, width = 10, units = "in")

        
  # BBRKC immature male with temperature -------
        BBRKC.immale.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # add POT points
          geom_sf(data = filter(cpue_mapdat, MAT_SEX == "Immature male",
                                VESSEL %in% c("Arctic Lady", "Seabrooke")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black") +
          scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 125), guide = "none",
                                breaks = seq(0, 125, by = 25))+ 
          scale_fill_gradientn(breaks = seq(0, 125, by = 25),
                               limits = c(0, 125), 
                               colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
          
          # POT legend
          guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                     override.aes = list(shape = c(4, rep(21, 5)),
                                                         size = seq(2, 10, by = ((10-2)/(6-1))))),
                 size = guide_legend(),
                 color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(cpue_mapdat, MAT_SEX == "Immature male",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 59), 
                                breaks = c(seq(0, 20, by = 5), 59))+ 
          scale_fill_gradientn(breaks = c(seq(0, 20, by = 5), 59),
                               limits = c(0, 59), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24,5)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Immature male (< 120mm)")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.immale.map, "./Figures/BBRKC_immale_temp.png", height = 7, width = 10, units = "in")
        
        
  # BBRKC mature female with temperature -------
        BBRKC.matfemale.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # add POT points
          geom_sf(data = filter(cpue_mapdat, MAT_SEX == "Mature female",
                                VESSEL %in% c("Arctic Lady", "Seabrooke")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black") +
          scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 125), guide = "none",
                                breaks = seq(0, 125, by = 25))+ 
          scale_fill_gradientn(breaks = seq(0, 125, by = 25),
                               limits = c(0, 125), 
                               colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
          
          # POT legend
          guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                     override.aes = list(shape = c(4, rep(21, 5)),
                                                         size = seq(2, 10, by = ((10-2)/(6-1))))),
                 size = guide_legend(),
                 color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(cpue_mapdat, MAT_SEX == "Mature female",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 147), 
                                breaks = c(seq(0, 20, by = 5), 147))+ 
          scale_fill_gradientn(breaks = c(seq(0, 20, by = 5), 147),
                               limits = c(0, 147), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24,5)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature female")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matfemale.map, "./Figures/BBRKC_matfem_temp.png", height = 7, width = 10, units = "in")
        
        
  # BBRKC immature female with temperature -------
        BBRKC.imfemale.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # add POT points
          geom_sf(data = filter(cpue_mapdat, MAT_SEX == "Immature female",
                                VESSEL %in% c("Arctic Lady", "Seabrooke")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black") +
          scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 125), guide = "none",
                                breaks = seq(0, 125, by = 25))+ 
          scale_fill_gradientn(breaks = seq(0, 125, by = 25),
                               limits = c(0, 125), 
                               colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
          
          # POT legend
          guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                     override.aes = list(shape = c(4, rep(21, 5)),
                                                         size = seq(2, 10, by = ((10-2)/(6-1))))),
                 size = guide_legend(),
                 color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(cpue_mapdat, MAT_SEX == "Immature female",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 25), 
                                breaks = c(seq(0, 20, by = 5), 25))+ 
          scale_fill_gradientn(breaks = c(seq(0, 20, by = 5), 25),
                               limits = c(0, 25), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24,5)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Immature female")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.imfemale.map, "./Figures/BBRKC_imfem_temp.png", height = 7, width = 10, units = "in")
        
        
  # Tanner male with temperature -------
        species = "MaleTanner"
        
        # Plot bycatch counts and weights (pot vs. trawl)
        tanner.male.map <- ggplot() +
                          # plot temperature
                          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
                          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                                             guide = guide_colorbar(title.position = "top", 
                                                                    discrete = FALSE)) +
                          guides(alpha = "none") +
                          new_scale("fill") +
                          
                          # plot mapping layers
                          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
                          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
                          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha = 0.9, linewidth = 1) +
                          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
                          scale_color_manual(values = c("black", "red"),
                                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                                             name = "") +
                          guides(color = guide_legend(nrow = 2)) +            
                          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
                          
                          # add POT points
                          geom_sf(data = filter(bycatch, VESSEL %in% c("Arctic Lady", "Seabrooke")),
                                  mapping = aes(size = get(species), fill = get(species), shape = get(species) == 0), 
                                  alpha = 0.5, colour = "black") +
                          scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
                          scale_size_continuous(range = c(2, 10), limits = c(0, 30), 
                                                breaks = seq(0, 30, by = 10), guide = "none") + 
                          scale_fill_gradientn(limits = c(0, 30), 
                                               breaks = seq(0, 30, by = 10),
                                               colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
                          
                          # POT legend
                          guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                                     override.aes = list(shape = c(4, rep(21, 3)),
                                                                         size = seq(2, 10, by = ((10-2)/(4-1))))),
                                 size = guide_legend(),
                                 color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
                          
                          # start a new scale
                          new_scale("fill") +
                          new_scale("shape") +
                          new_scale("size") +
                          
                          # add TRAWL points
                          geom_sf(data = bycatch %>% filter(VESSEL %in% c("Vesteraalen")),
                                  mapping = aes(size = get(species), fill = get(species), shape = get(species) == 0), 
                                  alpha = 0.5, colour = "black") +
                          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none") +
                          scale_size_continuous(range = c(2, 10), limits = c(0, 175),
                                                breaks = c(seq(0, 175, by = 25))) +
                          scale_fill_gradientn(breaks = c(seq(0, 175, by = 25)),
                                               limits = c(0, 175), 
                                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)]))) +
                          scale_x_continuous(breaks = map_layers$lon.breaks) +
                          scale_y_continuous(breaks = map_layers$lat.breaks) +
                          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Male Tanner crab") +
                          
                          # TRAWL legend
                          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                                     override.aes = list(shape = c(8, rep(24, 7)))),
                                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                                 color = guide_legend(nrow = 2)) +
                          
                          # crop spatial extent, add bathymetry, set theme
                          coord_sf(xlim = plot.boundary$x,
                                   ylim = plot.boundary$y) +
                          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                                         sf::st_transform(crs = map.crs),
                                       mapping = aes(label = lab))+
                          theme_bw() +
                          theme(axis.title = element_blank(),
                                axis.text = element_text(size = 10),
                                legend.text = element_text(size = 10),
                                legend.title = element_text(size = 10),
                                legend.position = "right",
                                legend.direction = "horizontal",
                                plot.title = element_text(face = "bold", size = 15),
                                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = tanner.male.map, paste0("./Figures/", species, "_temp.png"), height = 7, width = 10, units = "in")
        
        
  # Tanner female with temperature -------
        species = "FemaleTanner"
        
        # Plot bycatch counts and weights (pot vs. trawl)
        tanner.female.map <- ggplot() +
                            # plot temperature
                            geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
                            scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                                               guide = guide_colorbar(title.position = "top", 
                                                                      discrete = FALSE)) +
                            guides(alpha = "none") +
                            new_scale("fill") +
                            
                            # plot mapping layers
                            geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
                            # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                            geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
                            geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha = 0.9, linewidth = 1) +
                            geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
                            scale_color_manual(values = c("black", "red"),
                                               labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                                               name = "") +
                            guides(color = guide_legend(nrow = 2)) +            
                            geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
                            
                            # add POT points
                            geom_sf(data = filter(bycatch, VESSEL %in% c("Arctic Lady", "Seabrooke")),
                                    mapping = aes(size = get(species), fill = get(species), shape = get(species) == 0), 
                                    alpha = 0.5, colour = "black") +
                            scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
                            scale_size_continuous(range = c(2, 10), limits = c(0, 5), 
                                                  breaks = seq(0, 5, by = 1), guide = "none") + 
                            scale_fill_gradientn(limits = c(0, 5), 
                                                 breaks = seq(0, 5, by = 1),
                                                 colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
                            
                            # POT legend
                            guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                                       override.aes = list(shape = c(4, rep(21, 5)),
                                                                           size = seq(2, 10, by = ((10-2)/(6-1))))),
                                   size = guide_legend(),
                                   color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
                            
                            # start a new scale
                            new_scale("fill") +
                            new_scale("shape") +
                            new_scale("size") +
                            
                            # add TRAWL points
                            geom_sf(data = bycatch %>% filter(VESSEL %in% c("Vesteraalen")),
                                    mapping = aes(size = get(species), fill = get(species), shape = get(species) == 0), 
                                    alpha = 0.5, colour = "black") +
                            scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none") +
                            scale_size_continuous(range = c(2, 10), limits = c(0, 50),
                                                  breaks = c(seq(0, 50, by = 10))) +
                            scale_fill_gradientn(breaks = c(seq(0, 50, by = 10)),
                                                 limits = c(0, 50), 
                                                 colors = c("gray", rev(trawl_pal[5:length(trawl_pal)]))) +
                            scale_x_continuous(breaks = map_layers$lon.breaks) +
                            scale_y_continuous(breaks = map_layers$lat.breaks) +
                            labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Female Tanner crab") +
                            
                            # TRAWL legend
                            guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                                       override.aes = list(shape = c(8, rep(24, 5)))),
                                   fill = guide_legend(title = "Trawl Count", title.position = "top"),
                                   color = guide_legend(nrow = 2)) +
                            
                            # crop spatial extent, add bathymetry, set theme
                            coord_sf(xlim = plot.boundary$x,
                                     ylim = plot.boundary$y) +
                            geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                                                 x = c(-161.5, -165), y = c(58.3, 56.1)),
                                                      coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                                           sf::st_transform(crs = map.crs),
                                         mapping = aes(label = lab))+
                            theme_bw() +
                            theme(axis.title = element_blank(),
                                  axis.text = element_text(size = 10),
                                  legend.text = element_text(size = 10),
                                  legend.title = element_text(size = 10),
                                  legend.position = "right",
                                  legend.direction = "horizontal",
                                  plot.title = element_text(face = "bold", size = 15),
                                  plot.subtitle = element_text(size = 12))
                          
        ggsave(plot = tanner.female.map, paste0("./Figures/", species, "_temp.png"), height = 7, width = 10, units = "in")
        
        
  # Pacific Cod with temperature -------
        species = "PacificCod"
        
        # Plot bycatch counts and weights (pot vs. trawl)
        pcod.map <- ggplot() +
                    # plot temperature
                    geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
                    scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                                       guide = guide_colorbar(title.position = "top", 
                                                              discrete = FALSE)) +
                    guides(alpha = "none") +
                    new_scale("fill") +
                    
                    # plot mapping layers
                    geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
                    # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                    geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
                    geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha = 0.9, linewidth = 1) +
                    geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
                    scale_color_manual(values = c("black", "red"),
                                       labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                                       name = "") +
                    guides(color = guide_legend(nrow = 2)) +            
                    geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
                    
                    # add POT points
                    geom_sf(data = filter(bycatch, VESSEL %in% c("Arctic Lady", "Seabrooke")),
                            mapping = aes(size = get(species), fill = get(species), shape = get(species) == 0), 
                            alpha = 0.5, colour = "black") +
                    scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
                    scale_size_continuous(range = c(2, 10), limits = c(0, 20), 
                                          breaks = seq(0, 20, by = 5), guide = "none") + 
                    scale_fill_gradientn(limits = c(0, 20), 
                                         breaks = seq(0, 20, by = 5),
                                         colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
                    
                    # POT legend
                    guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                               override.aes = list(shape = c(4, rep(21, 4)),
                                                                   size = seq(2, 10, by = ((10-2)/(5-1))))),
                           size = guide_legend(),
                           color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
                    
                    # start a new scale
                    new_scale("fill") +
                    new_scale("shape") +
                    new_scale("size") +
                    
                    # add TRAWL points
                    geom_sf(data = bycatch %>% filter(VESSEL %in% c("Vesteraalen")),
                            mapping = aes(size = get(species), fill = get(species), shape = get(species) == 0), 
                            alpha = 0.5, colour = "black") +
                    scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none") +
                    scale_size_continuous(range = c(2, 10), limits = c(0, 100),
                                          breaks = c(seq(0, 100, by = 25))) +
                    scale_fill_gradientn(breaks = c(seq(0, 100, by = 25)),
                                         limits = c(0, 100), 
                                         colors = c("gray", rev(trawl_pal[5:length(trawl_pal)]))) +
                    scale_x_continuous(breaks = map_layers$lon.breaks) +
                    scale_y_continuous(breaks = map_layers$lat.breaks) +
                    labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Pacific cod") +
                    
                    # TRAWL legend
                    guides(size = guide_legend(title = "Trawl Weight (kg)", title.position = "top", nrow = 2, 
                                               override.aes = list(shape = c(8, rep(24, 4)))),
                           fill = guide_legend(title = "Trawl Weight (kg)", title.position = "top"),
                           color = guide_legend(nrow = 2)) +
                    
                    # crop spatial extent, add bathymetry, set theme
                    coord_sf(xlim = plot.boundary$x,
                             ylim = plot.boundary$y) +
                    geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                                         x = c(-161.5, -165), y = c(58.3, 56.1)),
                                              coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                                   sf::st_transform(crs = map.crs),
                                 mapping = aes(label = lab))+
                    theme_bw() +
                    theme(axis.title = element_blank(),
                          axis.text = element_text(size = 10),
                          legend.text = element_text(size = 10),
                          legend.title = element_text(size = 10),
                          legend.position = "right",
                          legend.direction = "horizontal",
                          plot.title = element_text(face = "bold", size = 15),
                          plot.subtitle = element_text(size = 12))
        
        ggsave(plot = pcod.map, paste0("./Figures/", species, "_map.png"), height = 7, width = 10, units = "in")

# SHELL CONDITION MAPS --------------------------------------------------------------------------------------------------
    
  # BBRKC mature male SC 2 with temperature -------
        BBRKC.matmale.sc2.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # add POT points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature male SC2",
                                VESSEL %in% c("Arctic Lady", "Seabrooke")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black") +
          scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 75), guide = "none",
                                breaks = seq(0, 75, by = 25))+ 
          scale_fill_gradientn(breaks = seq(0, 75, by = 25),
                               limits = c(0, 75), 
                               colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
          
          # POT legend
          guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                     override.aes = list(shape = c(4, rep(21, 3)),
                                                         size = seq(2, 10, by = ((10-2)/(4-1))))),
                 size = guide_legend(),
                 color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature male SC2",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 20), 
                                breaks = c(seq(0, 20, by = 5)))+ 
          scale_fill_gradientn(breaks = c(seq(0, 20, by = 5)),
                               limits = c(0, 20), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24,4)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature male (>= 120mm), shell condition 2")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matmale.sc2.map, "./Figures/BBRKC_matmale_SC2_temp.png", height = 7, width = 10, units = "in")
        
        
  # BBRKC mature male SC 3 with temperature -------
        BBRKC.matmale.sc3.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # add POT points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature male SC3",
                                VESSEL %in% c("Arctic Lady", "Seabrooke")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black") +
          scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 75), guide = "none",
                                breaks = seq(0, 75, by = 25))+ 
          scale_fill_gradientn(breaks = seq(0, 75, by = 25),
                               limits = c(0, 75), 
                               colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
          
          # POT legend
          guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                     override.aes = list(shape = c(4, rep(21, 3)),
                                                         size = seq(2, 10, by = ((10-2)/(4-1))))),
                 size = guide_legend(),
                 color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature male SC3",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 20), 
                                breaks = c(seq(0, 20, by = 5)))+ 
          scale_fill_gradientn(breaks = c(seq(0, 20, by = 5)),
                               limits = c(0, 20), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24,4)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature male (>= 120mm), shell condition 3")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matmale.sc3.map, "./Figures/BBRKC_matmale_SC3_temp.png", height = 7, width = 10, units = "in")
        
        
  # BBRKC mature female SC 2 with temperature -------
        BBRKC.matfemale.sc2.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # add POT points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature female SC2",
                                VESSEL %in% c("Arctic Lady", "Seabrooke")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black") +
          scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 75), guide = "none",
                                breaks = seq(0, 75, by = 25))+ 
          scale_fill_gradientn(breaks = seq(0, 75, by = 25),
                               limits = c(0, 75), 
                               colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
          
          # POT legend
          guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                     override.aes = list(shape = c(4, rep(21, 3)),
                                                         size = seq(2, 10, by = ((10-2)/(4-1))))),
                 size = guide_legend(),
                 color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature female SC2",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 143), 
                                breaks = c(seq(0, 15, by = 5), 143))+ 
          scale_fill_gradientn(breaks = c(seq(0, 15, by = 5), 143),
                               limits = c(0, 143), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24, 4)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature female, shell condition 2")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matfemale.sc2.map, "./Figures/BBRKC_matfem_SC2_temp.png", height = 7, width = 10, units = "in")
        
        
  # BBRKC mature female SC 3 with temperature ------- 
        BBRKC.matfemale.sc3.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # add POT points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature female SC3",
                                VESSEL %in% c("Arctic Lady", "Seabrooke")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black") +
          scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 75), guide = "none",
                                breaks = seq(0, 75, by = 25))+ 
          scale_fill_gradientn(breaks = seq(0, 75, by = 25),
                               limits = c(0, 75), 
                               colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
          
          # POT legend
          guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                     override.aes = list(shape = c(4, rep(21, 3)),
                                                         size = seq(2, 10, by = ((10-2)/(4-1))))),
                 size = guide_legend(),
                 color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature female SC3",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 20), 
                                breaks = c(seq(0, 20, by = 5)))+ 
          scale_fill_gradientn(breaks = c(seq(0, 20, by = 5)),
                               limits = c(0, 20), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24,4)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature female, shell condition 3")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matfemale.sc3.map, "./Figures/BBRKC_matfem_SC3_temp.png", height = 7, width = 10, units = "in")
        
        
        
  # BBRKC mature male SC 1 (trawl only) with temperature -------
        BBRKC.matmale.sc1.trawl.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature male SC1",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 20), 
                                breaks = c(seq(0, 20, by = 5)))+ 
          scale_fill_gradientn(breaks = c(seq(0, 20, by = 5)),
                               limits = c(0, 20), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24,4)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature male (>= 120mm), shell condition 1")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matmale.sc1.trawl.map, "./Figures/BBRKC_matmale_SC1_trawl_temp.png", height = 7, width = 10, units = "in")
        
  # BBRKC mature male SC 2 (trawl only) with temperature -------
        BBRKC.matmale.sc2.trawl.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature male SC2",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 20), 
                                breaks = c(seq(0, 20, by = 5)))+ 
          scale_fill_gradientn(breaks = c(seq(0, 20, by = 5)),
                               limits = c(0, 20), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24,4)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature male (>= 120mm), shell condition 2")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matmale.sc2.trawl.map, "./Figures/BBRKC_matmale_SC2_trawl_temp.png", height = 7, width = 10, units = "in")
        
        
  # BBRKC mature male SC 3 (trawl only) with temperature -------
        BBRKC.matmale.sc3.trawl.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature male SC3",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 20), 
                                breaks = c(seq(0, 20, by = 5)))+ 
          scale_fill_gradientn(breaks = c(seq(0, 20, by = 5)),
                               limits = c(0, 20), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24,4)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature male (>= 120mm), shell condition 3")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matmale.sc3.trawl.map, "./Figures/BBRKC_matmale_SC3_trawl_temp.png", height = 7, width = 10, units = "in")
        
        
  # BBRKC mature female SC 1 (trawl only) with temperature -------
        BBRKC.matfemale.sc1.trawl.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature female SC1",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 20), 
                                breaks = c(seq(0, 20, by = 5)))+ 
          scale_fill_gradientn(breaks = c(seq(0, 20, by = 5)),
                               limits = c(0, 20), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24,4)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature female, shell condition 1")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matfemale.sc1.trawl.map, "./Figures/BBRKC_matfem_SC1_trawl_temp.png", height = 7, width = 10, units = "in")
        
        
  # BBRKC mature female SC 2 (trawl only) with temperature -------
        BBRKC.matfemale.sc2.trawl.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature female SC2",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 143), 
                                breaks = c(seq(0, 15, by = 5), 143))+ 
          scale_fill_gradientn(breaks = c(seq(0, 15, by = 5), 143),
                               limits = c(0, 143), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24,4)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature female, shell condition 2")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matfemale.sc2.trawl.map, "./Figures/BBRKC_matfem_SC2_trawl_temp.png", height = 7, width = 10, units = "in")
        
        
  # BBRKC mature female SC 3 (trawl only) with temperature -------     
        BBRKC.matfemale.sc3.trawl.map <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE)) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # start a new scale
          new_scale_fill() +
          new_scale("shape") +
          new_scale("size") +
          
          # add TRAWL points
          geom_sf(data = filter(sc_counts, MAT_SEX_SC == "Mature female SC3",
                                VESSEL %in% c("Vesteraalen")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 20), 
                                breaks = c(seq(0, 20, by = 5)))+ 
          scale_fill_gradientn(breaks = c(seq(0, 20, by = 5)),
                               limits = c(0, 20), 
                               colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # TRAWL legend
          guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                                     override.aes = list(shape = c(8, rep(24,4)))),
                 fill = guide_legend(title = "Trawl Count", title.position = "top"),
                 color = guide_legend(nrow = 2)) +
          
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature female, shell condition 3")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matfemale.sc3.trawl.map, "./Figures/BBRKC_matfem_SC3_trawl_temp.png", height = 7, width = 10, units = "in")
        
      
# CPS1 / CPS2 COMPARISONS ----------------------------------------------------------------------------------------------
  # CPS1 BBRKC mature female with temperature -------
        CPS1_temp_df <- readRDS("./Data/CPS1_TEMP.rds")
        CPS1_cpue_mapdat <- readRDS("./Data/CPS1_COUNTS.rds")
        
        BBRKC.matfemale.cps1 <- ggplot() +
                                # plot temperature
                                geom_tile(data = CPS1_temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
                                scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                                                   guide = guide_colorbar(title.position = "top", 
                                                                          discrete = FALSE),
                                                   limits = c(min(c(min(temp_df$temperature), min(CPS1_temp_df$temperature))), 
                                                                  max(c(max(temp_df$temperature), max(CPS1_temp_df$temperature))))) +
                                guides(alpha = "none") +
                                new_scale("fill") +
                                
                                # plot mapping layers
                                geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
                                # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                                geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
                                geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
                                geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
                                geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
                                scale_color_manual(values = c("black", "red"), 
                                                   labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                                                   name = "") +
                                guides(color = guide_legend(nrow = 2)) +
                                
                                # add POT points
                                geom_sf(data = filter(CPS1_cpue_mapdat, MAT_SEX == "Mature female"),
                                        mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black") +
                                scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
                                scale_size_continuous(range = c(2, 10), limits = c(0, 125), guide = "none",
                                                      breaks = seq(0, 125, by = 25))+ 
                                scale_fill_gradientn(breaks = seq(0, 125, by = 25),
                                                     limits = c(0, 125), 
                                                     colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
                                
                                # POT legend
                                guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                                           override.aes = list(shape = c(4, rep(21, 5)),
                                                                               size = seq(2, 10, by = ((10-2)/(6-1))))),
                                       # size = guide_legend(),
                                       color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
                                
                                scale_x_continuous(breaks = map_layers$lon.breaks)+
                                scale_y_continuous(breaks = map_layers$lat.breaks)+
                                labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "Mature female")+
                                coord_sf(xlim = plot.boundary$x,
                                         ylim = plot.boundary$y) +
                                geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                                                     x = c(-161.5, -165), y = c(58.3, 56.1)),
                                                          coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                                               sf::st_transform(crs = map.crs),
                                             mapping = aes(label = lab))+
                                theme_bw() +
                                theme(axis.title = element_blank(),
                                      axis.text = element_text(size = 10),
                                      legend.text = element_text(size = 10),
                                      legend.title = element_text(size = 10),
                                      legend.position = "right",
                                      legend.direction = "horizontal",
                                      plot.title = element_text(face = "bold", size = 15),
                                      plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matfemale.cps1, "./Figures/BBRKC_CPS1matfem_temp.png", height = 7, width = 10, units = "in")
        
  
  # CPS2 BBRKC mature female with temperature -------      
        BBRKC.matfemale.cps2 <- ggplot() +
          # plot temperature
          geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
          scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                             guide = guide_colorbar(title.position = "top", 
                                                    discrete = FALSE),
                             limits = c(min(c(min(temp_df$temperature), min(CPS1_temp_df$temperature))), 
                                        max(c(max(temp_df$temperature), max(CPS1_temp_df$temperature))))) +
          guides(alpha = "none") +
          new_scale("fill") +
          
          # plot mapping layers
          geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
          # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
          geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
          geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
          scale_color_manual(values = c("black", "red"), 
                             labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                             name = "") +
          guides(color = guide_legend(nrow = 2)) +
          
          # add POT points
          geom_sf(data = filter(cpue_mapdat, MAT_SEX == "Mature female",
                                VESSEL %in% c("Arctic Lady", "Seabrooke")),
                  mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black") +
          scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
          scale_size_continuous(range = c(2, 10), limits = c(0, 125), guide = "none",
                                breaks = seq(0, 125, by = 25))+ 
          scale_fill_gradientn(breaks = seq(0, 125, by = 25),
                               limits = c(0, 125), 
                               colors = c("gray", rev(pot_pal[5:length(pot_pal)]))) +
          
          # POT legend
          guides(fill = guide_legend(title = "Pot Count",  title.position = "top",
                                     override.aes = list(shape = c(4, rep(21, 5)),
                                                         size = seq(2, 10, by = ((10-2)/(6-1))))),
                 # size = guide_legend(),
                 color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2)) +
          
          # # start a new scale
          # new_scale_fill() +
          # new_scale("shape") +
          # new_scale("size") +
          
          # # add TRAWL points
          # geom_sf(data = filter(cpue_mapdat, MAT_SEX == "Mature female",
          #                       VESSEL %in% c("Vesteraalen")),
          #         mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), alpha = 0.5, colour = "black")+
          # scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none")+
          # scale_size_continuous(range = c(2, 10), limits = c(0, 147), 
          #                       breaks = c(seq(0, 20, by = 5), 147))+ 
          # scale_fill_gradientn(breaks = c(seq(0, 20, by = 5), 147),
          #                      limits = c(0, 147), 
          #                      colors = c("gray", rev(trawl_pal[5:length(trawl_pal)])))+
          # # TRAWL legend
          # guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
          #                            override.aes = list(shape = c(8, rep(24,5)))),
          #        fill = guide_legend(title = "Trawl Count", title.position = "top"),
          #        color = guide_legend(nrow = 2)) +
          # 
          scale_x_continuous(breaks = map_layers$lon.breaks)+
          scale_y_continuous(breaks = map_layers$lat.breaks)+
          labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Mature female")+
          coord_sf(xlim = plot.boundary$x,
                   ylim = plot.boundary$y) +
          geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                               x = c(-161.5, -165), y = c(58.3, 56.1)),
                                    coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs),
                       mapping = aes(label = lab))+
          theme_bw() +
          theme(axis.title = element_blank(),
                axis.text = element_text(size = 10),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 10),
                legend.position = "right",
                legend.direction = "horizontal",
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 12))
        
        ggsave(plot = BBRKC.matfemale.cps2, "./Figures/BBRKC_CPS2matfem_temp.png", height = 7, width = 10, units = "in")
        
     
    # CPS1 rescaled temperature -------
        temp.scaled.cps1 <- ggplot() +
                        # plot temperature
                        geom_tile(data = CPS1_temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
                        scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                                           guide = guide_colorbar(title.position = "top", 
                                                                  discrete = FALSE),
                                           limits = c(min(c(min(temp_df$temperature), min(CPS1_temp_df$temperature))), 
                                                      max(c(max(temp_df$temperature), max(CPS1_temp_df$temperature))))) +
                        guides(alpha = "none") +
                        new_scale("fill") +
                        
                        # plot mapping layers
                        geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
                        # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                        geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
                        geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
                        geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
                        geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
                        scale_color_manual(values = c("black", "red"), 
                                           labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                                           name = "") +
                        guides(color = guide_legend(nrow = 2)) +
                        scale_x_continuous(breaks = map_layers$lon.breaks)+
                        scale_y_continuous(breaks = map_layers$lat.breaks)+
                        labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = "Temperature")+
                        coord_sf(xlim = plot.boundary$x,
                                 ylim = plot.boundary$y) +
                        geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                                             x = c(-161.5, -165), y = c(58.3, 56.1)),
                                                  coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                                       sf::st_transform(crs = map.crs),
                                     mapping = aes(label = lab))+
                        theme_bw() +
                        theme(axis.title = element_blank(),
                              axis.text = element_text(size = 10),
                              legend.text = element_text(size = 10),
                              legend.title = element_text(size = 10),
                              legend.position = "bottom",
                              legend.direction = "horizontal",
                              plot.title = element_text(face = "bold", size = 15),
                              plot.subtitle = element_text(size = 12))
                      
                      ggsave(plot = temp.scaled.cps1, "./Figures/CPS1_scaled_temp.png", height = 7, width = 10, units = "in")
        
        
        # CPS2 BBRKC mature female with temperature -------      
        temp.scaled.cps2 <- ggplot() +
                            # plot temperature
                            geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature, alpha = 0.8)) +
                            scale_fill_viridis(name = "Temperature (°C)", option = "plasma", 
                                               guide = guide_colorbar(title.position = "top", 
                                                                      discrete = FALSE),
                                               limits = c(min(c(min(temp_df$temperature), min(CPS1_temp_df$temperature))), 
                                                          max(c(max(temp_df$temperature), max(CPS1_temp_df$temperature))))) +
                            guides(alpha = "none") +
                            new_scale("fill") +
                            
                            # plot mapping layers
                            geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
                            # geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                            geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
                            geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha= 0.9, linewidth = 1) +
                            geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
                            geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
                            scale_color_manual(values = c("black", "red"), 
                                               labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                                               name = "") +
                            guides(color = guide_legend(nrow = 2)) +
                            scale_x_continuous(breaks = map_layers$lon.breaks)+
                            scale_y_continuous(breaks = map_layers$lat.breaks)+
                            labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Temperature")+
                            coord_sf(xlim = plot.boundary$x,
                                     ylim = plot.boundary$y) +
                            geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                                                 x = c(-161.5, -165), y = c(58.3, 56.1)),
                                                      coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                                           sf::st_transform(crs = map.crs),
                                         mapping = aes(label = lab))+
                            theme_bw() +
                            theme(axis.title = element_blank(),
                                  axis.text = element_text(size = 10),
                                  legend.text = element_text(size = 10),
                                  legend.title = element_text(size = 10),
                                  legend.position = "bottom",
                                  legend.direction = "horizontal",
                                  plot.title = element_text(face = "bold", size = 15),
                                  plot.subtitle = element_text(size = 12))
        
        ggsave(plot = temp.scaled.cps2, "./Figures/CPS2_scaled_temp.png", height = 7, width = 10, units = "in")
        
      
      