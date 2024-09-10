
# INSTALL PACKAGES ----------------------------------------------------------------------------------------------------------
# install.packages(c("tidyverse", "gsubfn", "terra", "rgdal", "colorRamps", "sf", "viridis", "grid", "shadowtext",
#                     "ggnewscale"))


# LOAD PACKAGES -------------------------------------------------------------------------------------------------------------
  library(tidyverse)


# LOAD DATA -----------------------------------------------------------------------------------------------------------------

  # Load mapping layers
    source("./Scripts/map_setup.R")

  # Read in POT and TRAWL data
    # I think I need to reframe this to make long rather than wide df...
      pot_crab_bycatch <- read.csv("./Data/CPS2_2024_pot_crab_bycatch.csv") 
      pot_fish_bycatch <- read.csv("./Data/CPS2_2024_pot_fish_bycatch.csv")

      # pot_bycatch <- left_join(pot_crab_bycatch, pot_fish_bycatch) %>%
      #                dplyr::rename(HAUL = SPN,
      #                              STATION = POT_ID) %>%
      #                dplyr::mutate(RockSole = 0) %>%
      #                dplyr::select(VESSEL, STATION, LAT_DD, LON_DD, MaleTanner, FemaleTanner, Tanner, 
      #                              MaleSnow, FemaleSnow, Snow, MaleHybrid, FemaleHybrid, Hybrid, HairCrab, 
      #                              PacificCod, Halibut, GreatSculpin, YellowfinSole, Pollock, StarryFlounder, 
      #                              RockSole, Other)
      pot_bycatch <- read.csv("./Data/CPS2_2024_POT_bycatch.csv") %>%
                     dplyr::rename(HAUL = SPN,
                                   STATION = BUOY) %>%
                     dplyr::mutate(RockSole = 0) %>%
                     filter(!nchar(STATION) > 3) %>% # filter out CAM, COFFIN, and BAIT POT_IDs
                     dplyr::select(VESSEL, STATION, LAT_DD, LON_DD, MaleTanner, FemaleTanner, Tanner,
                                   MaleSnow, FemaleSnow, Snow, MaleHybrid, FemaleHybrid, Hybrid, HairCrab,
                                   PacificCod, Halibut, GreatSculpin, YellowfinSole, Pollock, Starry.Flounder,
                                   RockSole, Other)

        
    trawl_fish_bycatch <- read.csv("./Data/CPS2_2024_TRAWL_FISH_bycatch.csv") 
    trawl_crab_bycatch <- read.csv("./Data/CPS2_2024_TRAWL_CRAB_bycatch.csv") 
    
    
# MAP POT AND TRAWL BYCATCH --------------------------------------------------------------------------------------------------------

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
  

  # Set up bycatch spp labels and crab/fish categories
    bycatch_labs <- data.frame(lab = c("Male Tanner crab", "Female Tanner crab", "Tanner crab", "Male snow crab",
                                       "Female snow crab", "Snow crab", "Male hybrid crab", "Female hybrid crab", 
                                       "Hybrid crab", "Hair crab", "Pacific cod", "Pacific halibut", "Great sculpin",
                                       "Yellowfin sole", "Pollock", "Starry flounder", "Rock sole"),
                               SPP_LABS = c("MaleTanner", "FemaleTanner", "Tanner", "MaleSnow", "FemaleSnow", "Snow", 
                                            "MaleHybrid", "FemaleHybrid", "Hybrid", "HairCrab", "PacificCod", "Halibut", 
                                            "GreatSculpin", "YellowfinSole", "Pollock", "StarryFlounder", "RockSole"))
    
    crab <- c("MaleTanner", "FemaleTanner", "Tanner", "MaleSnow", "FemaleSnow", "Snow", 
              "MaleHybrid", "FemaleHybrid", "Hybrid", "HairCrab")
    fish <- c("Hybrid", "HairCrab", "PacificCod", "Halibut", "GreatSculpin", 
              "YellowfinSole", "Pollock", "StarryFlounder", "RockSole")  
  

  
    
    # method = "POT"; data = bycatch; species = "Pollock"; species_lab = "Pollock"; out_name = "pollock"; breaks = 4
    # method = "POT"; data = bycatch; species = "PacificCod"; species_lab = "Pacific Cod"; out_name = "pcod"; breaks = 4
    # method = "POT"; data = bycatch; species = "YellowfinSole"; species_lab = "Yellowfin Sole"; out_name = "yfs"; breaks = 4
    # method = "POT"; data = bycatch; species = "RockSole"; species_lab = "Rock Sole"; out_name = "rs"; breaks = 4
    # 
    # bycatch_map(data = bycatch, species = "Pollock", species_lab = "Pollock")
    # bycatch_map(data = bycatch, species = "Tanner", species_lab = "Tanner crab")
    
    for(i in 1:nrow(bycatch_labs)){
      bycatch_map(data = bycatch, species = bycatch_labs[i, "SPP_LABS"], species_lab = bycatch_labs[i, "lab"])
    }
    

  bycatch_map <- function(data, species, species_lab){
   
    # Set up shape mapping
    shapes <- c(0, 15) # set shape mapping
    names(shapes) <- c("n=1", "n>1")
    
    # Specify palette
    pot_pal <- viridis::mako(10) 
    trawl_pal <- viridis::rocket(10)
    
    # Make function for rounding 
    roundUP <- function(x, m){x + m - x %% m}
    
    # calculate count limits for plotting
    pot_lim <- data %>% 
      dplyr::filter(VESSEL %in% c("Arctic Lady", "Seabrooke")) %>%
      select(all_of(species)) %>%
      st_drop_geometry() 
    
    pot_lim <- roundUP(max(pot_lim[,1]), 5)
    
    
    trawl_lim <- data %>% 
      dplyr::filter(VESSEL %in% c("Vesteraalen")) %>%
      select(all_of(species)) %>%
      st_drop_geometry() 
    
    trawl_lim <- roundUP(max(trawl_lim[,1]), 5)
    
    trawl_legend <- ifelse(species %in% crab, "Trawl Count", "Trawl Weight (kg)")
    
    
    
    # Plot bycatch counts and weights (pot vs. trawl)
    bycatch.maps <- ggplot() +
                    geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
                    geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
                    geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha = 0.9, linewidth = 1) +
                    geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
                    scale_color_manual(values = c("black", "red"),
                                       labels = c("EBS Summer Survey Boundary", "Red King Crab Savings Area"),
                                       name = "") +
                    guides(color = guide_legend(nrow = 2)) +            
                    geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
                    
                    
                    # geom_sf(data = data,
                    #         mapping = aes(size = get(species), fill = get(species), shape = get(species) == 0), 
                    #         alpha = 0.5, colour = "black") +
                    # scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
                    # scale_size_continuous(range = c(2, 10)) + 
                    # scale_fill_gradientn(colors = c("gray", rev(pal[5:length(pal)]))) +
                    
                    # add POT points
                    # ** automate the scaling to ID max count and do 5 breaks or something? rounded?
                  geom_sf(data = data %>% filter(VESSEL %in% c("Arctic Lady", "Seabrooke")),
                          mapping = aes(size = get(species), fill = get(species), shape = get(species) == 0), 
                          alpha = 0.5, colour = "black") +
                    scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
                    scale_size_continuous(range = c(2, 10), limits = c(0, pot_lim), breaks = seq(0, pot_lim, by = pot_lim/5), guide = "none") + 
                    scale_fill_gradientn(limits = c(0, pot_lim), breaks = seq(0, pot_lim, by = pot_lim/5),
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
                    # **need an ifelse to ID fish vs. crab to plot weight or count in legend...
                    geom_sf(data = data %>% filter(VESSEL %in% c("Vesteraalen")),
                            mapping = aes(size = get(species), fill = get(species), shape = get(species) == 0), 
                            alpha = 0.5, colour = "black") +
                    scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none") +
                    scale_size_continuous(range = c(2, 10), limits = c(0, trawl_lim),
                                          breaks = c(seq(0, trawl_lim, by = trawl_lim/5))) +
                    scale_fill_gradientn(breaks = c(seq(0, trawl_lim, by = trawl_lim/5)),
                                         limits = c(0, trawl_lim), 
                                         colors = c("gray", rev(trawl_pal[5:length(trawl_pal)]))) +
                    scale_x_continuous(breaks = map_layers$lon.breaks) +
                    scale_y_continuous(breaks = map_layers$lat.breaks) +
                    labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = species_lab) +
                    
                    # TRAWL legend
                    guides(size = guide_legend(title = trawl_legend, title.position = "top", nrow = 2, 
                                               override.aes = list(shape = c(8, rep(24, 5)))),
                           fill = guide_legend(title = trawl_legend, title.position = "top"),
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
    
    # Save figures
    ggsave(plot = bycatch.maps, paste0("./Figures/bycatch/", species, "_map.png"),
           height = 7, width = 10, units = "in")
  }
  
  # # Save plots
  #   plot_labs <- c("matmale", "immale", "matfem", "imfem", "legalmale", "sublegalmale")
  #   
  #   for(i in 1:4){
  #     ggsave(plot = BBRKC.maps[[i]], paste0("./Figures/BBRKC_", plot_labs[i], ".png"), 
  #            height = 7, width = 10, units = "in")
  #   }
  
