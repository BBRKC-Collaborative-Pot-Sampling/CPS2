# CPS2 Report Figures 
# Category:
# Script 1/X 
# Madi Heller-Shipley BSFRF 

# Package Requirements 
#devtools::install_github("AFSC-Shellfish-Assessment-Program/crabpack")
library(crabpack)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggthemes)
library(gridExtra)
library(here)
library(ggpattern)


#This script has dependencies from scripts provided by NOAA

## Source plot function scripts 
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/abund_ts_plots.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/male_1mm_sc_plots.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/female_1mm_shell_plots.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/female_1mm_egg_plots.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/female_1mm_clutch_plots.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/female_1mm_maturity_plots.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/weighted_crds_maps.R")


## Function to automate size-sex category labels for abundance timeseries plots
bioabund_labs <- function(species, stock){
  
  if(species %in% c("RKC", "BKC")){
    # if(stock %in% c("BBRKC", "PribRKC", "PribBKC", "StMattBKC", "NorthRKC", "NSRKC", "BKCNBS")){
    # Specify table and plot labels
    labs <- c(paste0("Immature male (<", mat_lookup$cutline[mat_lookup$stock == stock], " mm)"),
              paste0("Mature male (\u2265", mat_lookup$cutline[mat_lookup$stock == stock], " mm)"),
              paste0("Legal male (\u2265", mat_lookup$legal[mat_lookup$stock == stock], " mm)"),
              "Immature female",
              "Mature female")
    
    names(labs) <- c("Immature Male", "Mature Male", "Legal Male",  "Immature Female", "Mature Female") 
    
  } 
  
  if(species %in% c("SNOW", "TANNER", "HYBRID")){
    # if (stock %in% c("TannerE", "TannerW", "Snow", "Hybrid", "SnowNBS")){
    # Specify table and plot labels
    labs <- c(paste0("Small male (<", mat_lookup$cutline[mat_lookup$stock == stock], " mm)"),
              paste0("Large male (\u2265", mat_lookup$cutline[mat_lookup$stock == stock], " mm)"),
              paste0("Legal male (\u2265", mat_lookup$legal[mat_lookup$stock == stock], " mm)"),
              paste0("Industry preferred male (\u2265", mat_lookup$preferred[mat_lookup$stock == stock], " mm)"),
              "Immature female",
              "Mature female")
    
    names(labs) <- c("Small Male", "Large Male", "Legal Male", "Industry Preferred", "Immature Female", "Mature Female")
  } 
  
  if(species == "HAIR"){
    # if(stock %in% c("Hair", "HairNBS")){
    # Specify plot labels
    labs <- c(paste0("Sublegal male (<", mat_lookup$legal[mat_lookup$stock == stock], " mm)"),
              paste0("Legal male (\u2265", mat_lookup$legal[mat_lookup$stock == stock], " mm)"),
              "Female") 
    
    names(labs) <- c("Sublegal Male", "Legal Male", "Female") 
  }
  
  return(labs)
}

#specify AKFIN API
channel <- "API"

#======================
#Set-up Work Flow 
#======================

fig_dir <- paste0("./Figures_Report/")
out_dir <- paste0("./Outputs_Report/")

#Figure 1: BBRKC Abundance 
#==========

## NMFS 
RKC_EBS <- crabpack::get_specimen_data(species = "RKC",
                                       region = "EBS",
                                       years = c(1979:2024),
                                       channel = channel)

# Specify current year and recent (last 5) years
current_year <- 2024
recent_years <- c(2019, 2021:2024)


# Mature female abundance
matFem_BBRKC<- crabpack::calc_bioabund(crab_data = RKC_EBS,
                        species = "RKC",
                        region = "EBS",
                        district = "BB",
                        crab_category = "mature_female",
                        female_maturity = "morphometric") %>%
  mutate(STOCK = "BBRKC")

# Mature male abundance 
matMal_BBRKC<- crabpack::calc_bioabund(crab_data = RKC_EBS,
                                       species = "RKC",
                                       region = "EBS",
                                       district = "BB",
                                       crab_category = "mature_male"
                                       ) %>%
  mutate(STOCK = "BBRKC")


Fig1_a <- abund_ts_plot(matMal_BBRKC)
Fig1_b <- abund_ts_plot(matFem_BBRKC)

Fig1<-grid.arrange(Fig1_a, Fig1_b, nrow = 2, ncol = 1)
ggsave("C:/GitHub/CPS2/Figures_Report/Fig1.png", plot = Fig1, width = 12, height = 8, dpi = 300)


# Figure 2: Closure area map
#===============

# MAP OF TRAWL CLOSURE AREAS W/CPS2 SURVEY AREA 
# Load boundaries
st_read(survey_gdb,layer="Area516") ->  area516

st_read(survey_gdb,layer="BycatchZone1") ->  zone1

st_read(survey_gdb,layer="NBBTCA") ->  nbbtca
st_read(survey_gdb,layer="TogiakTrawlArea") ->  togtrawl
st_read("C:/GitHub/CPS2/Data/BBDistrict.gdb", layer = "BB_District") -> bb_dist

# Transform plot boundary
plot.boundary.untrans <- data.frame(y = c(54, 59.5), 
                                    x = c(-168, -158)) 

plot.boundary <-  plot.boundary.untrans %>%
  sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
  sf::st_transform(crs = map.crs) %>%
  sf::st_coordinates() %>%
  as.data.frame() %>%
  dplyr::rename(x = X, y = Y) # plot boundary projected

breaks.x <- map_layers$lon.breaks[(map_layers$lon.breaks >= plot.boundary.untrans$x[1] &  # set lon breaks
                                     map_layers$lon.breaks < plot.boundary.untrans$x[2]) == TRUE]

breaks.y <- map_layers$lat.breaks[(map_layers$lat.breaks > plot.boundary.untrans$y[1] & # set lat breaks
                                     map_layers$lat.breaks < plot.boundary.untrans$y[2]) == TRUE]
# Plot
Fig2<- ggplot() +
  geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color=alpha("grey70")) +
  geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(fill = "firebrick3"), color = "black", alpha= 0.5, linewidth = 1) +
  geom_sf(data = st_as_sf(RKCSA), aes(fill = "firebrick2"),  color = "black", alpha =0.5, linewidth = 0.5) +
  geom_sf(data = st_as_sf(bb_dist), fill = NA, aes(colour = "black"), linewidth = 1) +
  #geom_sf(data = st_as_sf(zone1$Shape), fill = "blue", alpha = 0.15, aes(color = "black"), linewidth = 0.5) +
  geom_sf_pattern(data = st_as_sf(zone1$Shape),
                  aes(pattern_type = "stripe", pattern_angle = 30), fill = NA, color = "black", pattern_alpha = 0.15)+
  geom_sf(data = st_as_sf(area516$Shape), fill = NA, aes(color = "darkblue"), linewidth = 1) +
  geom_sf(data = st_as_sf(togtrawl$Shape), aes(fill = "gold2"), alpha = 0.5) +
  geom_sf_pattern(data = st_as_sf(nbbtca$Shape),
                  aes(pattern_type = "stripe", pattern_angle = 120), fill = NA, color = "black", pattern_alpha = 0.15)+
  geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "olivedrab2"), linewidth = 1)+
  geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
  scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
  scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
  labs(title = "BBRKC Collaborative Pot Sampling", subtitle = "Survey extent and closure areas")+
  scale_color_manual(values = c("black", "darkblue", "olivedrab2"), 
                     labels = c("Bristol Bay management boundary", "Area 516", "CPS survey extent"),
                     name = "") +
  scale_fill_manual(values = c(alpha("firebrick3", 0.5), alpha("indianred1", 0.25), alpha("gold2", 0.5)),
                    labels = c("Red King Crab Savings Area", "Red King Crab Savings Subarea", "Nearshore Bristol Bay Trawl Area"),
                    name = "")+
  
  scale_pattern_type_manual(values = c("stripe", "stripe"),
                            labels = c("Bycatch Limitation Zone 1", "Northern Bristol Bay Trawl Closure Area"),
                            name = "")+
  scale_pattern_angle_continuous(range = c(30, 120),
                                 breaks = c(30, 120),
                                 labels = c("Bycatch Limitation Zone 1", "Nearshore Bristol Bay Trawl Closure Area"),
                                 name = "")+
  coord_sf(xlim = plot.boundary$x,
           ylim = plot.boundary$y) +
  geom_sf_text(sf::st_as_sf(data.frame(lab= c("50m", "100m"), 
                                       x = c(-165.8, -166.2), y = c(58.3, 56.5)),
                            coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                 sf::st_transform(crs = map.crs),
               mapping = aes(label = lab))+
  guides(color = guide_legend(nrow = 3), fill = guide_legend(nrow = 3), pattern_angle = guide_legend(nrow = 2), 
         pattern_type = "none") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 12),
        panel.grid.major = element_blank())

ggsave(plot = Fig2, "./Figures_Report/Fig2.png", height=7, width=10, units="in")


# Figure 3: CPS Planned v actual stations 
#===============
# Load planned pot IDS and coordinates csv, transform lat and lon to mapping crs

## YOU NEED TO ADD THE TRAWL STATIONS STILL!!!!

surv_plan <- read.csv("./Data/Pot IDs and Coordinates.csv") %>%
  mutate(Longitude = Longitude*-1)

surv_effort <-read.csv("./Data/CPS2_2024_Potlifts.csv")

surv_plan <- surv_plan %>%
  rename(POT_ID = Pot.ID, LAT_DD = Latitude, LON_DD = Longitude)

unsurveyed_pots <- dplyr::anti_join(surv_plan,surv_effort, by = "POT_ID") %>%
  mutate(VESSEL = "Not surveyed")

surv_all <- dplyr::bind_rows(surv_effort %>% dplyr::select(VESSEL, POT_ID, LAT_DD, LON_DD),
                             unsurveyed_pots %>% dplyr::select(VESSEL, POT_ID, LAT_DD, LON_DD))

surv_all <- surv_all %>%
  sf::st_as_sf(coords = c("LON_DD", "LAT_DD"), crs = 4326) %>%
  sf::st_transform(crs = map.crs) %>%
  distinct()

surv_all <- surv_all %>%
  mutate(LAT_DD = st_coordinates(.)[,2],  # Extract new lat
         LON_DD = st_coordinates(.)[,1])  # Extract new lon

summary(surv_all %>% filter(VESSEL == "Seabrooke") %>% select(LAT_DD, LON_DD))

# Transform plot boundary
plot.boundary.untrans <- data.frame(y = c(54.5, 58.5), 
                                    x = c(-164.8, -159))

plot.boundary <-  plot.boundary.untrans %>%
  sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
  sf::st_transform(crs = map.crs) %>%
  sf::st_coordinates() %>%
  as.data.frame() %>%
  dplyr::rename(x = X, y = Y) # plot boundary projected

breaks.x <- map_layers$lon.breaks[(map_layers$lon.breaks >= plot.boundary.untrans$x[1] &  # set lon breaks
                                     map_layers$lon.breaks < plot.boundary.untrans$x[2]) == TRUE]

breaks.y <- map_layers$lat.breaks[(map_layers$lat.breaks > plot.boundary.untrans$y[1] & # set lat breaks
                                     map_layers$lat.breaks < plot.boundary.untrans$y[2]) == TRUE]

# Make row and column labels
sta_labs <- data.frame(lab = c("A", "B", "C", "D", "E", "F", "G", "H", "I",
                               "J", "K", "10", "20", "30", "40", "50", "60", "70", "80"),
                       x = c(-164.2729, -163.8046, -163.3358, -162.8671, -162.3984, -161.9298, -161.4614,-160.9909, -160.5225,
                             -160.0544, -159.583, -164.6662, -164.6588, -164.6516, -164.6446, -164.6376, -164.6308, -164.6241,
                             -164.6175),
                       y = c(57.99736, 58.01606, 58.01643, 58.01504, 58.0119, 58.007, 58.00035, 58.02509, 58.01492,
                             58.003, 58.02244, 57.49888, 57.16714, 56.83529, 56.50332, 56.17121, 55.83898, 55.50659,
                             55.17405)) %>%
  sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
  sf::st_transform(crs = map.crs)


# Plot
Fig3a <- ggplot() +
  #geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "olivedrab2"), linewidth = 1)+
  geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
  geom_sf(data = surv_all,
          mapping = aes(fill = VESSEL), shape = 21, size = 3)+   
  scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
  scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
  labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Survey vessel effort")+
  coord_sf(xlim = plot.boundary$x,
           ylim = plot.boundary$y) +
  scale_fill_manual(values = c("cyan4", "white", "darkgoldenrod2"),
                    labels = c("Arctic Lady", "Not surveyed","Seabrooke"),
                    name = "")+
  geom_sf_text(data = sta_labs,
               mapping = aes(label = lab))+
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 12),
        panel.grid.major = element_blank())

ggsave(plot = Fig3a, "./Figures_Report/Fig3a.png", height=7, width=10, units="in")

## YOU NEED TO ADD THE TRAWL STATIONS STILL!!!!



# Fig 16: Krieged temps 
#============

# Fig 17: Ice Ice baby
#============

Ice_3.17.24 <- st_read("./Data/Spatial Layers/Daily Ice Extent/3.17.24/nic_autoc2024077n_pl_a.shp") %>%
  sf::st_transform(crs = map.crs)

Ice_3.26.24 <- st_read("./Data/Spatial Layers/Daily Ice Extent/3.26.24/nic_autoc2024086n_pl_a.shp") %>%
  sf::st_transform(crs = map.crs)

Ice_4.5.24 <- st_read("./Data/Spatial Layers/Daily Ice Extent/4.5.24/nic_autoc2024096n_pl_a.shp") %>%
  sf::st_transform(crs = map.crs)

# Transform plot boundary
plot.boundary.untrans <- data.frame(y = c(54, 59.5), 
                                    x = c(-168, -158)) 

plot.boundary <-  plot.boundary.untrans %>%
  sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
  sf::st_transform(crs = map.crs) %>%
  sf::st_coordinates() %>%
  as.data.frame() %>%
  dplyr::rename(x = X, y = Y) # plot boundary projected

breaks.x <- map_layers$lon.breaks[(map_layers$lon.breaks >= plot.boundary.untrans$x[1] &  # set lon breaks
                                     map_layers$lon.breaks < plot.boundary.untrans$x[2]) == TRUE]

breaks.y <- map_layers$lat.breaks[(map_layers$lat.breaks > plot.boundary.untrans$y[1] & # set lat breaks
                                     map_layers$lat.breaks < plot.boundary.untrans$y[2]) == TRUE]
# Plot
Fig17 <- ggplot() +
  geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
  
  # Outline Red King Crab Savings Area without fill
  geom_sf(data = st_as_sf(RKCSA_sub), fill = NA, color = "grey50", linewidth = 0.5, aes(linetype = "RKCSA Sub")) +
  geom_sf(data = st_as_sf(RKCSA), fill = NA, color = "grey50", linewidth = 0.5, aes(linetype = "RKCSA")) +
  
  geom_sf(data = st_as_sf(CPS1_bound), fill = NA, aes(color = "black"), linewidth = 1) +
  geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
  
  # Add ice extent layers with transparency and borders
  geom_sf(data = Ice_3.17.24, aes(fill = "Ice Extent 3/17/24"), color = alpha("lightblue", 0.5), size = 0.1) +
  geom_sf(data = Ice_3.26.24, aes(fill = "Ice Extent 3/26/24"), color = alpha("deepskyblue", 0.5), size = 0.1) +
  geom_sf(data = Ice_4.5.24, aes(fill = "Ice Extent 4/5/24"), color = alpha("darkblue", 0.5), size = 0.1) +
  
  scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W")) +
  scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N")) +
  
  labs(title = "2024 BBRKC Collaborative Pot Sampling", 
       subtitle = "Ice extent") +
  
  scale_color_manual(values = c("black"), 
                     labels = c("CPS2 Survey extent"),
                     name = "Survey Boundaries") +
  
  scale_fill_manual(values = c(alpha("lightblue", 0.3), alpha("deepskyblue", 0.3), alpha("darkblue", 0.3)),
                    labels = c("Ice Extent 3/17/24",
                               "Ice Extent 3/26/24",
                               "Ice Extent 4/5/24"),
                    name = "Ice Extent") +
  
  scale_linetype_manual(values = c("solid", "dashed"),
                        labels = c("Red King Crab Savings Area", 
                                   "Red King Crab Savings Subarea"),
                        name = "Management Areas") +
  
  coord_sf(xlim = plot.boundary$x, ylim = plot.boundary$y) +
  
  guides(color = guide_legend(order = 1, nrow = 3), 
         fill = guide_legend(order = 3, nrow = 3), 
         linetype = guide_legend(order = 2, nrow = 2)) +
  
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "vertical",
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 12),
        panel.grid.major = element_blank())

ggsave(plot = Fig17, "./Figures_Report/Fig17.png", height=7, width=10, units="in")


