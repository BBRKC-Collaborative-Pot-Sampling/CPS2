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
library(Cairo)
library(sp)
library(ggh4x)
library(sf)
library(gstat)
library(raster)
library(sf)
library(viridis)
library(shadowtext) # for geom_shadowtext
library(patchwork)
#This script has dependencies from scripts provided by NOAA

## Source plot function scripts 
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/lookup_tables.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/abund_ts_plots.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/male_1mm_sc_plots.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/female_1mm_shell_plots.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/female_1mm_egg_plots.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/female_1mm_clutch_plots.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/female_1mm_maturity_plots.R")
source("C:/GitHub/NOAA Scripts/scripts_for_Madi/weighted_crds_maps.R")
source("C:/GitHub/CPS2/Scripts/map_setup.R")
#source("C:/GitHub/CPS1/Scripts/TechReport_2023_processing.R")


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

sizeComp_dataPrep <- function(df, gear, species, sex, metric) {
  
  if(gear == "trawl") {
  df <- df %>%
    filter(SPECIES_NAME == species, SEX == sex)%>%
    mutate(SIZE_1MM = round(LENGTH))
  } else if( gear == "pot"){
    df <- df %>%
      filter(SEX == sex) %>%
      mutate(SIZE_1MM = round(LENGTH))
      }
  
  if (metric == "shell condition") {
    if ("SHELL_CONDITION" %in% names(df)) {
      df <- df %>%
        mutate(METRIC_TEXT = case_when(
          SHELL_CONDITION %in% 0:1 ~ "Soft Molting",
          SHELL_CONDITION == 2 ~ "New Hard",
          SHELL_CONDITION == 3 ~ "Old",
          SHELL_CONDITION %in% 4:5 ~ "Very Old",
          TRUE ~ "Unknown"
        ))
    }
  } else if (metric == "egg condition" && sex == 2) {
    if ("EGG_CONDITION" %in% names(df)) {
      df <- df %>%
        mutate(METRIC_TEXT = case_when(
          EGG_CONDITION == 0 ~ "No Eggs",
          EGG_CONDITION == 1 ~ "Uneyed Eggs",
          EGG_CONDITION == 2 ~ "Eyed Eggs",
          EGG_CONDITION == 3 ~ "Dead Eggs",
          EGG_CONDITION == 4 ~ "Empty Egg Cases",
          EGG_CONDITION == 5 ~ "Hatching",
          ((CLUTCH_SIZE == 999 & EGG_CONDITION == 0) | is.na(EGG_CONDITION)) ~ "Unknown",
          TRUE ~ "Unknown"
        ))
    }
  } else if (metric == "clutch fullness" && sex == 2) {
    if ("CLUTCH_SIZE" %in% names(df)) {
      df <- df %>%
        mutate(METRIC_TEXT = case_when(
          CLUTCH_SIZE == 0 ~ "Immature",
          CLUTCH_SIZE == 1 ~ "Mature Barren",
          CLUTCH_SIZE == 2 ~ "Trace",
          CLUTCH_SIZE == 3 ~ "Quarter Full",
          CLUTCH_SIZE == 4 ~ "Half Full",
          CLUTCH_SIZE == 5 ~ "Three Quarter Full",
          CLUTCH_SIZE == 6 ~ "Full",
          CLUTCH_SIZE == 999 | is.na(CLUTCH_SIZE) ~ "Unknown",
          TRUE ~ "Unknown"
        ))
    }
  } else if (metric == "maturity status" && sex == 2) {
    if ("CLUTCH_SIZE" %in% names(df)) {
      df <- df %>%
        mutate(METRIC_TEXT = case_when(
          CLUTCH_SIZE == 0 ~ "Immature",
          CLUTCH_SIZE %in% 1:6 ~ "Mature",
          is.na(CLUTCH_SIZE) | CLUTCH_SIZE == 999 ~ "Unknown",
          TRUE ~ "Unknown"
        ))
    }
  }
  return(df)
}


#specify AKFIN API
channel <- "API"

#======================
#Set-up Work Flow 
#======================

fig_dir <- paste0("./Figures_Report/")
out_dir <- paste0("./Outputs_Report/")

# Data =================
## NMFS 
RKC_EBS <- crabpack::get_specimen_data(species = "RKC",
                                       region = "EBS",
                                       years = c(1979:2024),
                                       channel = channel)
NMFS_Hauls<-RKC_EBS$haul

## BSFRF 

## Trawl
TrawlStations<-read.csv(here("data/CPS2_NephropsTows.csv"))
Vest_Crab<-read.csv(here("data/CPS2_2024_Processed_Trawl_Specimen_Data.csv"))

### Pots 
Pots_2024<-read.csv(here("data/CPS2_2024_Processed_Pot_Specimen_Data.csv"))
Pots_2023 <- read.csv("C:/GitHub/CPS1/Output/CPS1_2023_Processed_Specimen_Data.csv")


#Fig 1: BBRKC Abundance ===========

# Specify current year and recent (last 5) years
current_year <- 2024
recent_years <- c(2022:2024)


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
#ggsave("C:/GitHub/CPS2/Figures_Report/Fig1.png", plot = Fig1, width = 12, height = 8, dpi = 300)


# Fig 2: Closure area map ========

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
                    labels = c("Red King Crab Savings Area", "Red King Crab Savings Subarea", "Northern Bristol Bay Trawl Area"),
                    name = "")+
  
  scale_pattern_type_manual(values = c("black", "black"),
                            labels = c("Bycatch Limitation Zone 1", "Nearshore Bristol Bay Trawl Closure Area"),
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

#Issue with legend rendering when using ggsave, save from consol
#ggsave(plot = Fig2, "./Figures_Report/Fig2.2.png", height=7, width=10, units="in", dpi = 300)


# Fig 3: CPS Planned v actual stations ============

# Load planned station IDs and coordinates csv, transform lat and lon to mapping crs

#read in trawl station data
TrawlStations<-read.csv(here("data/CPS2_NephropsTows.csv"))

TrawlStat <- TrawlStations %>%
  sf::st_as_sf(coords = c("Lon", "Lat"), crs = 4326) %>%
  sf::st_transform(crs = map.crs) %>%
  distinct()

#Extract lat and long
TrawlStat <- TrawlStat %>%
  mutate(Lat = st_coordinates(.)[,2],  # Extract new lat
         Lon = st_coordinates(.)[,1])  # Extract new lon


# read in pot station data
surv_plan <- read.csv("./Data/Pot IDs and Coordinates.csv") %>%
  mutate(Longitude = Longitude*-1)

# read in completed pot stations
surv_effort <-read.csv("./Data/CPS2_2024_Potlifts.csv")

surv_plan <- surv_plan %>%
  rename(POT_ID = Pot.ID, LAT_DD = Latitude, LON_DD = Longitude)

#Create Not surveyed category for vessel
unsurveyed_pots <- dplyr::anti_join(surv_plan,surv_effort, by = "POT_ID") %>%
  mutate(VESSEL = "Not surveyed")

#Bind together
surv_all <- dplyr::bind_rows(surv_effort %>% dplyr::select(VESSEL, POT_ID, LAT_DD, LON_DD),
                             unsurveyed_pots %>% dplyr::select(VESSEL, POT_ID, LAT_DD, LON_DD))
# Make into a spatial layer 
surv_all <- surv_all %>%
  sf::st_as_sf(coords = c("LON_DD", "LAT_DD"), crs = 4326) %>%
  sf::st_transform(crs = map.crs) %>%
  distinct()

#Extract lat and long
surv_all <- surv_all %>%
  mutate(LAT_DD = st_coordinates(.)[,2],  # Extract new lat
         LON_DD = st_coordinates(.)[,1])  # Extract new lon

#summary(surv_all %>% filter(VESSEL == "Seabrooke") %>% select(LAT_DD, LON_DD))

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
  geom_sf(data = TrawlStat,
          mapping = aes(fill = Fished), shape = 24, size = 2)+
  scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
  scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
  labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Survey vessel effort")+
  coord_sf(xlim = plot.boundary$x,
           ylim = plot.boundary$y) +
  scale_fill_manual(values = c("cyan4","white", "white","darkgoldenrod2", "red4" ),
                    labels = c("Arctic Lady", "Not surveyed (trawl)", "Not surveyed (pots)", "Seabrooke", "Vesteraalen"),
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

# Fig 15: Soak Time Distribution ===========

Fig15<- ggplot()+
  geom_histogram(surv_effort, mapping = aes(x = SOAK_TIME), fill = "cyan4", bins = 35, color = "black")+
  theme_bw() +
  ylab("Frequency")+
  xlab("Soak time (hours)")+
  scale_x_continuous(breaks = seq(0,max(surv_effort$SOAK_TIME), 2)) +
  labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Soak time")+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 12),
        panel.grid.major = element_blank()) 

#ggsave(plot = Fig15, "./Figures_Report/Fig15.png", height=7, width=10, units="in", dpi = 300)

# Fig 16: Krieged temps =========

# Start with NOAA and years 

# Load function to generate NMFS temperature maps for Bristol Bay
#NOTE!! MODIFIED TO OVERLAP CPUE AS WELL

temp_map_ebs_nbs <- function(haul_ebs, years, cpue_data = NULL, mat_sex = NULL){
  
  haul_ebs %>%
    dplyr::filter(YEAR %in% years, HAUL_TYPE != 17) %>%
    dplyr::select("REGION", "YEAR", "STATION_ID", 
                  "MID_LATITUDE", "MID_LONGITUDE", "GEAR_TEMPERATURE") %>%
    dplyr::filter(!is.na(GEAR_TEMPERATURE), 
                  !is.na(MID_LATITUDE), 
                  !is.na(MID_LONGITUDE)) -> temp
  
  # load EBS-NBS survey extent for masking
  interpolation.crs <- map.crs
  
  # Make raster for interpolation
  cell.resolution = 1000
  in.crs = "+proj=longlat +datum=NAD83"
  extrap.box = c(xmn = -165, xmx = -157, ymn = 50, ymx = 59)
  
  # Transform plot boundary
  plot.boundary.untrans <- data.frame(y = c(50, 59), 
                                      x = c(-165, -157))
  
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
  
  
  n_dim <- floor(abs(plot.boundary$x[1] - plot.boundary$x[2]))/cell.resolution
  
  sp_interp.raster <- raster::raster(xmn = plot.boundary$x[1], 
                                     xmx = plot.boundary$x[2], 
                                     ymn = plot.boundary$y[1], 
                                     ymx = plot.boundary$y[2], 
                                     nrow = n_dim, 
                                     ncol = n_dim)
  
  raster::projection(sp_interp.raster) <- interpolation.crs
  
  # Transform data for interpolation ----
  sp_interp.df <- unique(temp)
  sp::coordinates(sp_interp.df) <- c(x = "MID_LONGITUDE", y = "MID_LATITUDE")
  sp::proj4string(sp_interp.df) <- sp::CRS(in.crs)
  sp_interp.df <- sp::spTransform(sp_interp.df, sp::CRS(interpolation.crs))
  
  # Set up a new IDW for ordinary kriging ----
  idw_vgm_fit <- gstat::gstat(formula = GEAR_TEMPERATURE ~ 1, 
                              locations = sp_interp.df, 
                              nmax = Inf)
  
  # Ordinary Kriging: Stein's Matern VGM----
  ste.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit), 
                                    vgm(c("Ste")))
  
  ste_fit <- gstat::gstat(formula = GEAR_TEMPERATURE ~ 1, 
                          locations = sp_interp.df, 
                          model = ste.vgfit, 
                          nmax = Inf)
  
  ste.predict <- predict(ste_fit, as(sp_interp.raster, "SpatialGrid"))
  
  # write unmasked surfaces to raster, stacked by year
  ste.predict %>%
    raster::raster(.) %>%
    mask(st_transform(map_layers$survey.area, map.crs)) ->  temp_rast
  
  print("Finished Kriging")
  
  # extract interpolated data from raster to data frame
  coords<-coordinates(temp_rast)
  
  temp_df_nmfs <-na.omit(data.frame(coords, temperature = temp_rast@data@values, year = years))
  
  
  temp_breaks <- c(-Inf, seq(-1,8,1), Inf)
  viridis_option <- "H" # viridis turbo palette
  n_temp_breaks <- length(temp_breaks)-1
  
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
  # Year labels
  sf::st_as_sf(data.frame(lab= paste("", years, "\nNMFS"), 
                          x = c(-160), y = c(55.2)),
               coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    cbind(years, st_coordinates(.)) -> year_lab
  
  # Start of Plotting
  
  if(is.null(cpue_data)){# Plot interpolated data 
    print("Plotting only discrete temperature")
    
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = temp_df_nmfs, 
                       aes(x = x, 
                           y = y,
                           fill = cut(temperature, 
                                      breaks = temp_breaks))) +
    geom_sf(data = st_transform(map_layers$survey.area, map.crs), fill = NA, linewidth = 0.5) +
    geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1)+
    geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
    scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
    scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y) +
    ggplot2::scale_fill_manual(name = "Temperature (°C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
                               labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                                          "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
    geom_shadowtext(year_lab,
                    mapping = aes(label = lab, x = X, y = Y), size = 8,  color = "black", bg.color = "white")+
    theme_bw() +
    theme(axis.title = element_blank(),
           axis.text = element_text(size = 15),
          # legend.text = element_text(size = 15),
          # legend.title = element_text(size = 15),
           legend.position = "none",
          # legend.direction = "horizontal",
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 12))  -> map
  
  }else{
  
    print("Plotting temperature and CPUE")
  
    # 1. Prep filtered specimen data
    nmfs_filtered <- cpue_data$specimen %>%
      filter(YEAR %in% years, REGION == "EBS", SPECIES == "RKC") %>%
      mutate(
        legal = SEX == 1 & SIZE_1MM >= 135,
        mature = (SEX == 1 & SIZE_1MM >= 120) | (SEX == 2 & SIZE_1MM >= 90),
        immature = (SEX == 1 & SIZE_1MM < 120) | (SEX == 2 & SIZE_1MM < 90)
      ) %>%
      mutate(
        MAT_SEX_LIST = case_when(
          legal ~ list(c("Mature Male", "Legal Male")),
          mature & SEX == 1 & SIZE_1MM < 135 ~ list("Mature Male"),
          immature & SEX == 1 ~ list("Immature Male"),
          mature & SEX == 2 ~ list("Mature Female"),
          immature & SEX == 2 ~ list("Immature Female"),
          TRUE ~ list(NA_character_)
        )
      ) %>%
      unnest(MAT_SEX_LIST) %>%
      rename(MAT_SEX = MAT_SEX_LIST) %>%
      filter(!is.na(MAT_SEX))
      
    
      # 2 Summarize crab counts per haul 
      crab_counts <- nmfs_filtered %>%
        group_by(HAULJOIN, MAT_SEX) %>%
        summarise(TOTAL = sum(SAMPLING_FACTOR), .groups = "drop")
      
      
      haul_info <- cpue_data$haul %>%
        dplyr::select(HAULJOIN, STATION_ID, AREA_SWEPT, MID_LATITUDE, MID_LONGITUDE)
      
      cpue_summary <- crab_counts %>%
        left_join(haul_info, by = "HAULJOIN") %>%
        mutate(CPUE = (TOTAL / AREA_SWEPT))
    
      cpue_station <- cpue_summary %>%
       group_by(STATION_ID, MAT_SEX) %>%
       summarise(CPUE = mean(CPUE), .groups = "drop") # IS IT MEAN OR SUM 
     
      all_hauls <- read.csv("C:/GitHub/CPS1/Data/haul_newtimeseries.csv")
      
      #stations <- na.omit(pull(sta, "BBRKC")) # or whatever your calc_factor is
      
      BB_Stations <- all_hauls %>%
        filter(SURVEY_YEAR %in% 2023, HAUL_TYPE != 17) %>%
        filter(STATIONID %in% stations) %>%
        group_by(STATIONID) %>%
        summarise(
          MID_LATITUDE = first(MID_LATITUDE),
          MID_LONGITUDE = first(MID_LONGITUDE),
          AREA_SWEPT = mean(DISTANCE_FISHED * NET_WIDTH / 10, na.rm = TRUE),  # or whatever your CPUE calc uses
          .groups = "drop"
        ) %>%
        rename(STATION_ID = STATIONID) 

    # 3. Create all station × sex combos
    mat_sex_combos <- unique(nmfs_filtered$MAT_SEX)

    full_station_combos <- expand_grid(
      STATION_ID = BB_Stations$STATION_ID,
      MAT_SEX = mat_sex_combos,
    ) %>%
      left_join(cpue_station, by = c("STATION_ID", "MAT_SEX")) %>%
      left_join(BB_Stations %>% dplyr::select(STATION_ID, MID_LATITUDE, MID_LONGITUDE), by = "STATION_ID") %>%
      mutate(CPUE = replace_na(CPUE, 0),
             zero_cpue = CPUE == 0)

    # 4. Filter CPUE to the relevant sex
    
    cpue_layer <- full_station_combos %>%
      filter(MAT_SEX == mat_sex) %>%
      mutate(CPUE = replace_na(CPUE, 0),
             zero_cpue = CPUE == 0) %>%
      st_as_sf(coords = c("MID_LONGITUDE", "MID_LATITUDE"), crs = 4326) %>%
      st_transform(crs = map.crs)
    
    # Number of size breaks
    if(!mat_sex %in% c("Mature Female", "All crab")){
      num = 500
    } else{
      num = 1500
    }  
 
    
    # Year labeling 
    
    sf::st_as_sf(data.frame(lab= paste("", years, "\nNMFS"), 
                            x = c(-160), y = c(55.2)),
                 coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
      sf::st_transform(crs = map.crs) %>%
      cbind(years, st_coordinates(.)) -> year_lab
    
    # Plot
    ggplot() +
      ggplot2::geom_tile(data = temp_df_nmfs, 
                         aes(x = x, 
                             y = y,
                             fill = cut(temperature, 
                                        breaks = temp_breaks))) +
      geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 0.5)+
      geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
      #CPUE 
      geom_sf(data = cpue_layer,
              aes(size=CPUE, shape = zero_cpue), 
              color = "black",
              fill = "black",
              alpha = 0.5,
              stroke = 0.3)+
      # Shape legend override
      scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none")+
      # Map boundaries
      geom_sf(data = st_transform(map_layers$survey.area, map.crs), fill = NA, linewidth = 0.5) +
      scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
      scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
      #labs(title = paste(mat_sex))+
      coord_sf(xlim = plot.boundary$x,
               ylim = plot.boundary$y) +
      # Fill legend
      ggplot2::scale_fill_manual(name = "Temperature (°C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
                                 labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                                            "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), 
                                 drop = FALSE) +
      geom_shadowtext(year_lab,
                      mapping = aes(label = lab, x = X, y = Y), size = 8,  color = "black", bg.color = "white")+
      scale_size_continuous(range = c(2, 10), limits = c(0, max(cpue_layer$CPUE)), 
                            breaks =seq(0, max(cpue_layer$CPUE), by = num))+
      guides(size = guide_legend(
        title.position = "top",
        ncol = 1, 
        order = 1, 
        override.aes = list(
          shape = c(4, rep(21, length(seq(0, max(cpue_layer$CPUE), by = num))-1)))),
          fill = "none")+
      geom_shadowtext(year_lab,
                      mapping = aes(label = lab, x = X, y = Y), size = 8,  color = "black", bg.color = "white")+
      theme_bw() +
      theme(axis.title = element_blank(),
            axis.text = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.position = "right",
            legend.direction = "horizontal") -> map
            #plot.title = "none",
            #plot.subtitle = "none")  
    
    
  }
    return(map)
}

krigeTemp_2022<-temp_map_ebs_nbs(NMFS_Hauls, 2022)
krigeTemp_2023<-temp_map_ebs_nbs(NMFS_Hauls, 2023)
krigeTemp_2024<-temp_map_ebs_nbs(NMFS_Hauls, 2024)

# CPS1 

# Data 
CPS1_TEMP <- readRDS("C:/GitHub/CPS2/Data/CPS1_TEMP.rds")
temploggers <- read.csv("C:/GitHub/CPS1/Data/2023_BBRKC_ALL_TEMPS_SAL/CPS1_2023BB_RKC_LoggerData_ADFG.csv")

temploggers <- sf::st_as_sf(temploggers, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  sf::st_transform(crs = map.crs)

# Generate temperature map for CPS1
# Year labels
sf::st_as_sf(data.frame(lab= paste("", 2023, "\nCPS1"), 
                        x = c(-160), y = c(55.2)),
             coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
  sf::st_transform(crs = map.crs) %>%
  cbind(2023, st_coordinates(.)) -> CPS1year_lab

# Temp breaks
temp_breaks <- c(-Inf, seq(-1,8,1), Inf)
viridis_option <- "H" # viridis turbo palette
n_temp_breaks <- length(temp_breaks)-1

# Plot
ggplot() +
  ggplot2::geom_tile(data = CPS1_TEMP, 
                     aes(x = x, 
                         y = y,
                         fill = cut(temperature, 
                                    breaks = temp_breaks))) +
  geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1)+
  geom_sf(data = CPS1_bathy, color=alpha("white")) +
  #geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
  geom_sf(data = temploggers, color = "black", size = 1)+
  geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
  scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W"))+
  scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N"))+
  coord_sf(xlim = plot.boundary$x,
           ylim = plot.boundary$y) +
  ggplot2::scale_fill_manual(name = "Temperature (°C)", values = viridis_pal(option = viridis_option)(n_temp_breaks),
                             labels = c(expression(""<=-1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                                        "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"), drop = FALSE) +
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  geom_shadowtext(CPS1year_lab,
                  mapping = aes(label = lab, x = X, y = Y), size = 8,  color = "black", bg.color = "white")+
  geom_shadowtext(data = (sf::st_as_sf(data.frame(lab= c("35m", "45m", "55m", "65m", "75m", "85m"), 
                                                  x = c(-160.6, -161.6, -161.6, -161.2, -163.5, -163.6), 
                                                  y = c(56.1, 56, 57.5, 57.2, 56.76, 56.38)),
                                       coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                            sf::st_transform(crs = map.crs) %>%
                            cbind(st_coordinates(.))),
                  mapping = aes(label = lab, x = X, y = Y), color = "black", bg.color = "white", size = 3.5)+
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 15),
        #legend.key.width = unit(12, "mm"),
        legend.position = "none",
        #legend.direction = "horizontal",
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 12)) -> CPS1_tempdiscrete

#CPS2 (WARNING, MADI ONLY KIND OF KNOWS WHAT SHE IS DOING HERE)

# Data 
CPS2_Temp<-read.csv(here("Data/Temperature Data/CPS2_TempData_Combined.csv"))

# Clean and convert to sf
CPS2_clean <- CPS2_Temp %>%
  filter(!is.na(latitude), !is.na(longitude), !is.na(temperature_C))

CPS2_sf <- st_as_sf(CPS2_clean, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = map.crs)

# Convert to Spatial for autoKrige
CPS2_sp <- as(CPS2_sf, "Spatial")

# Build raster grid from bounding box of data
bb <- st_bbox(CPS1_bound)
cell_size <- 1000  # 10km

# Create raster from bounding box
r <- raster(
  xmn = bb["xmin"], xmx = bb["xmax"],
  ymn = bb["ymin"], ymx = bb["ymax"],
  res = cell_size,
  crs = map.crs
)

# Convert to SpatialGrid
grid <- as(r, "SpatialGrid")

# Run kriging USING AUTOKRIG BECAUSE ISSUES WITH NOAA METHOD AND CONVERGENCE
auto <- autoKrige(temperature_C ~ 1,
                  input_data = CPS2_sp,
                  new_data = grid)

# Convert to df for plotting
kriged_raster <- raster(auto$krige_output)
kriged_df <- as.data.frame(kriged_raster, xy = TRUE)
names(kriged_df)[3] <- "temperature"

# # Init Plot
# ggplot(kriged_df, aes(x = x, y = y, fill = temperature)) +
#   geom_tile() +
#   scale_fill_viridis_c(name = "Temperature (°C)") +
#   coord_equal() +
#   theme_minimal()


# Set up temp breaks to mimic CPS1 style 

# Same temp breaks and style as CPS1
temp_breaks <- c(-Inf, seq(-1, 8, 1), Inf)
viridis_option <- "H"
n_temp_breaks <- length(temp_breaks) - 1

# CPS2 year label (adjust x/y as needed for label placement)
CPS2year_lab <- sf::st_as_sf(data.frame(lab = "2024\nCPS2",
                                        x = c(-160), y = c(55.2)),
                             coords = c("x", "y"), crs = 4326) %>%
  sf::st_transform(crs = map.crs) %>%
  cbind(st_coordinates(.))

# Trim boundaries 
CPS2_bound_trans <- st_transform(CPS1_bound, crs = crs(kriged_raster))

# 2. Convert to Spatial for masking
CPS2_bound_sp <- as(CPS2_bound_trans, "Spatial")

# 3. Mask the kriged raster using the CPS2 boundary
masked_raster <- raster::mask(kriged_raster, CPS2_bound_sp)

# 4. Convert the masked raster to dataframe for ggplot
masked_df <- as.data.frame(masked_raster, xy = TRUE)
names(masked_df)[3] <- "temperature"
masked_df <- na.omit(masked_df)

# 5. Temp df so all colors show up in legend
temp_levels <- c("<=-1", "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                 "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1")

masked_df$temperature_bin <- cut(masked_df$temperature, 
                                 breaks = temp_breaks, 
                                 labels = temp_levels,
                                 include.lowest = TRUE, 
                                 right = TRUE)
masked_df$temperature_bin <- factor(masked_df$temperature_bin, levels = temp_levels)

# Create dummy rows — one for each *missing* level
missing_levels <- setdiff(temp_levels, unique(masked_df$temperature_bin))
dummy_rows <- data.frame(
  x = NA, y = NA, temperature_bin = factor(missing_levels, levels = temp_levels)
)

# Combine
plot_df <- rbind(masked_df[, c("x", "y", "temperature_bin")], dummy_rows)
krig_cps2<-plot_df
# Plot
ggplot() +
  geom_tile(data = plot_df,
            aes(x = x, y = y, fill = temperature_bin)) +
            #aes(x = x, y = y, fill = cut(temperature, breaks = temp_breaks))) +
  geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1) +
  geom_sf(data = CPS1_bathy, color = alpha("white")) +
  geom_sf(data = CPS2_sf, color = "black", size = 1) +
  geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
  scale_x_continuous(breaks = c(-165, -160), labels = paste0(c(165, 160), "°W")) +
  scale_y_continuous(breaks = c(56, 58), labels = paste0(c(56, 58), "°N")) +
  coord_sf(xlim = plot.boundary$x, ylim = plot.boundary$y) +
  scale_fill_manual(name = "Temperature (°C)",
                    values = viridis::viridis(length(temp_levels), option = viridis_option),
                    labels = temp_levels,
                    drop = FALSE) +
                    #values = viridis::viridis(n_temp_breaks, option = viridis_option),
                    #labels = c(expression("" <= -1), "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                    #           "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1"),
                    #drop = FALSE) +
  guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 2)) +
  geom_shadowtext(data = CPS2year_lab,
                  aes(label = lab, x = X, y = Y),
                  size = 8, color = "black", bg.color = "white") +
  geom_shadowtext(data = sf::st_as_sf(data.frame(
    lab = c("35m", "45m", "55m", "65m", "75m", "85m"),
    x = c(-160.6, -161.6, -161.6, -161.2, -163.5, -163.6),
    y = c(56.1, 56, 57.5, 57.2, 56.76, 56.38)),
    coords = c("x", "y"), crs = 4326) %>%
      sf::st_transform(crs = map.crs) %>%
      cbind(st_coordinates(.)),
    aes(label = lab, x = X, y = Y),
    color = "black", bg.color = "white", size = 3.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 15),
        legend.key.width = unit(4, "mm"),
        legend.position = "right",
        legend.direction = "vertical",
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 12)) -> CPS2_tempdiscrete


# Finalize the figure 
# Combine the three NMFS plots into a row

top_row <- krigeTemp_2022 + krigeTemp_2023 + krigeTemp_2024 + 
  plot_layout(ncol = 3)

# Combine CPS plots with an empty spacer to help center them
bottom_row <- plot_spacer() + CPS1_tempdiscrete + CPS2_tempdiscrete + plot_spacer() +
  plot_layout(ncol = 4, widths = c(0.25, 1.5, 1.5, 0.25))

# Stack top and bottom
final_layout <- top_row / bottom_row + 
  plot_layout(heights = c(1, 1.4))  # Equal height rows; tweak if needed

#ggsave("Figures_Report/Fig16.png", plot = final_layout, width = 15, height = 10, dpi = 300)

# Fig 17: Ice Ice baby ==============

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

#ggsave(plot = Fig17, "./Figures_Report/Fig17.png", height=7, width=10, units="in")



# Figs 18-22: Male size-frequency distribution by shell condition ==========

# Step 1: parce the data 
# NMFS Data
recent_years <- c(2022:2024)

male_abundance_NMFS <- calc_bioabund(crab_data = RKC_EBS,species = "RKC", district = "BB", years = recent_years, sex = "male", shell_condition = "all_categories", bin_1mm = TRUE) %>%
  mutate(STOCK = "BBRKC") 

female_abundance_NMFS <- calc_bioabund(crab_data = RKC_EBS,species = "RKC", district = "BB", years = recent_years, sex = "female", female_maturity = "all_categories", shell_condition = "all_categories", egg_condition = "all_categories", clutch_size = "all_categories", bin_1mm = TRUE) %>%
  mutate(STOCK = "BBRKC") 

pot_23_male_sc  <- sizeComp_dataPrep(Pots_2023, gear = "pot", species = "red king crab", sex = 1, metric = "shell condition")
pot_24_male_sc  <- sizeComp_dataPrep(Pots_2024, gear = "pot", species = "red king crab", sex = 1, metric = "shell condition")
trawl_24_male_sc <- sizeComp_dataPrep(Vest_Crab, gear = "trawl", species = "red king crab", sex = 1, metric = "shell condition")

pot_23_fem_sc  <- sizeComp_dataPrep(Pots_2023, gear = "pot", species = "red king crab", sex = 2, metric = "shell condition")
pot_24_fem_sc  <- sizeComp_dataPrep(Pots_2024, gear = "pot", species = "red king crab", sex = 2, metric = "shell condition")
trawl_24_fem_sc <- sizeComp_dataPrep(Vest_Crab, gear = "trawl", species = "red king crab", sex = 2, metric = "shell condition")

pot_23_fem_clutch  <- sizeComp_dataPrep(Pots_2023, gear = "pot", species = "red king crab", sex = 2, metric = "clutch fullness")
pot_24_fem_clutch  <- sizeComp_dataPrep(Pots_2024, gear = "pot", species = "red king crab", sex = 2, metric = "clutch fullness")
trawl_24_fem_clutch <- sizeComp_dataPrep(Vest_Crab, gear = "trawl", species = "red king crab", sex = 2, metric = "clutch fullness")

pot_23_fem_eggs  <- sizeComp_dataPrep(Pots_2023, gear = "pot", species = "red king crab", sex = 2, metric = "egg condition")
pot_24_fem_eggs  <- sizeComp_dataPrep(Pots_2024, gear = "pot", species = "red king crab", sex = 2, metric = "egg condition")
trawl_24_fem_eggs <- sizeComp_dataPrep(Vest_Crab, gear = "trawl", species = "red king crab", sex = 2, metric = "egg condition")

pot_23_fem_mat  <- sizeComp_dataPrep(Pots_2023, gear = "pot", species = "red king crab", sex = 2, metric = "maturity status")
pot_24_fem_mat  <- sizeComp_dataPrep(Pots_2024, gear = "pot", species = "red king crab", sex = 2, metric = "maturity status")
trawl_24_fem_mat <- sizeComp_dataPrep(Vest_Crab, gear = "trawl", species = "red king crab", sex = 2, metric = "maturity status")


CompFigures <- function(NMFS_dat, Vest_dat, Pots_2023, Pots_2024, 
                        metric, sex_label, stock = "BBRKC") {
  
    # --- Define reusable color palettes by metric
  metric_colors <- list(
    "shell condition" = c("Very Old" = "#a6611a", 
                          "Old" = "#dfc27d", 
                          "New Hard" = "#80cdc1", 
                          "Soft Molting" = "#018571"),
    
    "clutch fullness" = c("Immature" = "#a6cee3", 
                          "Mature Barren" = "#1f78b4", 
                          "Trace" = "#b2df8a", 
                          "Quarter Full" = "#33a02c",
                          "Half Full" = "#fb9a99", 
                          "Three Quarter Full" = "#e31a1c", 
                          "Full" = "#fdbf6f", 
                          "Unknown" = "#999999"),
    
    "egg condition" = c("No Eggs" = "#8dd3c7", 
                        "Uneyed Eggs" = "#ffffb3", 
                        "Eyed Eggs" = "#bebada", 
                        "Dead Eggs" = "#fb8072",
                        "Empty Egg Cases" = "#80b1d3", 
                        "Hatching" = "#fdb462", 
                        "Unknown" = "#b3b3b3"),
    
    "maturity status" = c("Immature" = "#7570b3", 
                          "Mature" = "#d95f02", 
                          "Unknown" = "#b3b3b3")
  )
  
  # --- Harmonize NMFS METRIC_TEXT
  if (!"METRIC_TEXT" %in% names(NMFS_dat)) {
    if (metric == "shell condition" && "SHELL_TEXT" %in% names(NMFS_dat)) {
      NMFS_dat <- NMFS_dat %>% rename(METRIC_TEXT = SHELL_TEXT)
    } else if (metric == "egg condition" && "EGG_CONDITION_TEXT" %in% names(NMFS_dat)) {
      NMFS_dat <- NMFS_dat %>% rename(METRIC_TEXT = EGG_CONDITION_TEXT)
    } else if (metric == "clutch fullness" && "CLUTCH_TEXT" %in% names(NMFS_dat)) {
      NMFS_dat <- NMFS_dat %>% rename(METRIC_TEXT = CLUTCH_TEXT)
    } else if (metric == "maturity status" && "CLUTCH_TEXT" %in% names(NMFS_dat)) {
      NMFS_dat <- NMFS_dat %>% rename(METRIC_TEXT = CLUTCH_TEXT)
      #   NMFS_dat <- NMFS_dat %>%
    #     mutate(METRIC_TEXT = ifelse(CLUTCH_TEXT == "Immature", "Immature", "Mature"))
     }
  }
  
  #print(NMFS_dat)
  # --- Create standardized dataframes
  NMFS_plot <- NMFS_dat %>%
    group_by(YEAR, SIZE_1MM, METRIC_TEXT) %>%
    summarise(VALUE = sum(ABUNDANCE, na.rm = TRUE)/1e6, .groups = "drop") %>%
    #summarise(VALUE = sum(ABUNDANCE, na.rm = TRUE), .groups = "drop") %>%
    mutate(SOURCE = paste0("NMFS ", YEAR))

  
  
  #print(unique(NMFS_plot$METRIC_TEXT))
  
  Pots_2023_plot <- Pots_2023 %>%
    group_by(SIZE_1MM, METRIC_TEXT) %>%
    summarise(VALUE = n(), .groups = "drop") %>%
    mutate(SOURCE = "Pots 2023")
  
  Pots_2024_plot <- Pots_2024 %>%
    group_by(SIZE_1MM, METRIC_TEXT) %>%
    summarise(VALUE = n(), .groups = "drop") %>%
    mutate(SOURCE = "Pots 2024")
  
  Vest_plot <- Vest_dat %>%
    group_by(SIZE_1MM, METRIC_TEXT) %>%
    summarise(VALUE = n(), .groups = "drop") %>%
    mutate(SOURCE = "Trawl 2024")
  
  print(unique(Vest_plot$METRIC_TEXT))
  
  # --- Combine
  df_all <- bind_rows(NMFS_plot, Pots_2023_plot, Pots_2024_plot, Vest_plot)
  
  # --- Standardize labels
  if(metric == "shell condition"){
  df_all <- df_all %>%
    mutate(METRIC_TEXT = case_when(
      METRIC_TEXT %in% c("Soft Molting", "soft_molting") ~ "Soft Molting",
      METRIC_TEXT %in% c("New Hard", "new_hardshell") ~ "New Hard",
      METRIC_TEXT %in% c("Old", "oldshell") ~ "Old",
      METRIC_TEXT %in% c("Very Old", "very_oldshell") ~ "Very Old",
      TRUE ~ METRIC_TEXT
    ))} else if(metric == "clutch fullness"){
      df_all <-df_all %>%
        mutate(METRIC_TEXT = case_when(
          METRIC_TEXT %in% c("Full", "full") ~ "Full", #6 full
          METRIC_TEXT %in% ("three_quarter") ~ "Three Quarter Full", #5 3/4
          METRIC_TEXT %in% ("half") ~ "Half Full", #4 half
          METRIC_TEXT %in% ("quarter") ~ "Quarter Full", #3 quarter
          METRIC_TEXT %in% ("trace") ~ "Trace", #2 trace
          METRIC_TEXT %in% ("mature_barren") ~ "Mature Barren", #1 mature barren 
          METRIC_TEXT %in% ("immature") ~ "Immature",#0 immature
          METRIC_TEXT %in% c("unknown", "unknown" ,"999", NA) ~ "Unknown",
          TRUE ~ METRIC_TEXT
        ))} else if (metric == "egg condition") {
          df_all <- df_all %>%
            mutate(METRIC_TEXT = case_when(
              METRIC_TEXT %in% c("none", "No Eggs") ~ "No Eggs",
              METRIC_TEXT %in% c("uneyed", "Uneyed Eggs") ~ "Uneyed Eggs",
              METRIC_TEXT %in% c("eyed", "Eyed Eggs") ~ "Eyed Eggs",
              METRIC_TEXT %in% c("dead") ~ "Dead Eggs",
              METRIC_TEXT %in% c("hatching") ~ "Hatching",
              METRIC_TEXT %in% c("empty_cases", "empty cases") ~ "Empty Egg Cases",
              TRUE ~ METRIC_TEXT
            ))} else if (metric == "maturity status") {
              df_all <- df_all %>%
                mutate(METRIC_TEXT = case_when(
                  METRIC_TEXT %in% c("Immature", "immature", "0") ~ "Immature",
                  METRIC_TEXT %in% c("Mature Barren", "Trace", "Quarter Full", "Half Full", 
                                     "Three Quarter Full", "Full") ~ "Mature",
                  METRIC_TEXT %in% c("Unknown", "999", NA_character_) ~ "Unknown",
                  TRUE ~ "Mature"  # catch all for anything unlisted but not Immature or Unknown
                ))
            }
  
  # --- Factor ordering for facets and legend
  df_all$SOURCE <- factor(df_all$SOURCE, 
                          levels = c("NMFS 2022", "Pots 2023", "NMFS 2023", "Pots 2024", "NMFS 2024","Trawl 2024"))
  
  df_all$METRIC_TEXT <- factor(df_all$METRIC_TEXT, levels = names(metric_colors[[metric]]))
  
  # --- Plot
  p <- ggplot(df_all, aes(x = SIZE_1MM, y = VALUE, fill = METRIC_TEXT)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~SOURCE, ncol = 2, scales = "free_y") +
    scale_x_continuous(limits = c(0, 210), expand = c(0, 0)) +
    scale_fill_manual(values = metric_colors[[metric]], name = metric) +
    labs(
      x = paste0("Size (mm) – ", sex_label),
      y = "Abundance (millions) / Count",
      fill = metric,
      title = paste0(stock," ", sex_label, " Size Composition – ", metric)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = "black")
    )
  
  return(p)
}


Fig18<-CompFigures(
  NMFS_dat = male_abundance_NMFS,
  Vest_dat = trawl_24_male_sc,
  Pots_2023 = pot_23_male_sc,
  Pots_2024 = pot_24_male_sc,
  metric = "shell condition",
  sex_label = "Male",
  stock = "BBRKC"
)
#ggsave(plot = Fig18, "./Figures_Report/Fig18.png", height=7, width=10, units="in")

Fig19<-CompFigures(
  NMFS_dat = female_abundance_NMFS,
  Vest_dat = trawl_24_fem_mat,
  Pots_2023 = pot_23_fem_mat,
  Pots_2024 = pot_24_fem_mat,
  metric = "maturity status",
  sex_label = "Female",
  stock = "BBRKC"
)
#ggsave(plot = Fig19, "./Figures_Report/Fig19.png", height=7, width=10, units="in")

Fig20<-CompFigures(
  NMFS_dat = female_abundance_NMFS,
  Vest_dat = trawl_24_fem_sc,
  Pots_2023 = pot_23_fem_sc,
  Pots_2024 = pot_24_fem_sc,
  metric = "shell condition",
  sex_label = "Female",
  stock = "BBRKC"
)
#ggsave(plot = Fig20, "./Figures_Report/Fig20.png", height=7, width=10, units="in")

Fig21<-CompFigures(
  NMFS_dat = female_abundance_NMFS,
  Vest_dat = trawl_24_fem_eggs,
  Pots_2023 = pot_23_fem_eggs,
  Pots_2024 = pot_24_fem_eggs,
  metric = "egg condition",
  sex_label = "Female",
  stock = "BBRKC"
)
#ggsave(plot = Fig21, "./Figures_Report/Fig21.png", height=7, width=10, units="in")

Fig22<-CompFigures(
  NMFS_dat = female_abundance_NMFS,
  Vest_dat = trawl_24_fem_clutch,
  Pots_2023 = pot_23_fem_clutch,
  Pots_2024 = pot_24_fem_clutch,
  metric = "clutch fullness",
  sex_label = "Female",
  stock = "BBRKC"
)
#ggsave(plot = Fig22, "./Figures_Report/Fig22.png", height=7, width=10, units="in")


# Fig 23: Overall RKC distributions single panel (just CPS2, both gears). NO TEMP. ========
# Combine pot and trawl counts into one dataframe

# Read in POT and TRAWL data
pot_cpue <- read.csv("./Data/CPS2_2024_potcatch.csv") %>%
  dplyr::rename(HAUL = SPN,
                STATION = POT_ID) %>%
  filter(!nchar(STATION) > 3) # filter out CAM, COFFIN, and BAIT POT_IDs

# Combine pot and trawl counts into one dataframe

trawl_cpue <- read.csv("./Data/CPS2_2024_trawlcatch.csv")

cpue <- trawl_cpue %>% 
  rbind(., pot_cpue %>% dplyr::select(-c(BUOY, DATE_SET, TIME_SET, DATE_HAUL, TIME_HAUL, SOAK_TIME, CATCH_PER_HOUR)))

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

# max.date <- max(pot_cpue_mapdat$DATE_HAUL) # label for most recent pot haul date

# Set up shape mapping
shapes <- c(0,15) # set shape mapping
names(shapes) <- c("n=1", "n>1")

# Set up RKC labels
mat_labs <- c("Mature female", "Immature female", "Mature male (>= 120mm)", "Immature male (< 120mm)", "Legal male (>= 135mm)", "Sublegal male")
names(mat_labs) <- c("Mature female", "Immature female", "Mature male", "Immature male", "Legal male", "Sublegal male")

mat_labs <- data.frame(lab = c("Mature female", "Immature female", "Mature male (>= 120mm)", 
                               "Immature male (< 120mm)", "Legal male (>= 135mm)", "Sublegal male"),
                       MAT_SEX = c("Mature female", "Immature female", "Mature male", "Immature male","Legal male", "Sublegal male"))

mat_sex_combos <- c("Mature male", "Immature male", "Mature female", "Immature female", "Legal male", "Sublegal male")

# Specify palette
pot_pal <- viridis::mako(10) 
trawl_pal <- viridis::rocket(10)


# Plot total counts
total.map <-  ggplot() +
  # plot mapping layers
  geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
  geom_sf(data = st_as_sf(BB_strata), fill = NA, mapping = aes(color = "black"), linewidth = 1) +
  geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha = 0.9, linewidth = 1) +
  geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
  scale_color_manual(values = c("black", "red"),
                     labels = c("EBS Summer Survey Boundary", "Red King Crab Savings Area"),
                     name = "") +
  #guides(color = guide_legend(nrow = 2)) +
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
         #size = guide_legend(),
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
  scale_size_continuous(range = c(2, 10), 
                        limits = c(0, 154),
                        breaks = c(seq(0, 25, by = 5), 98, 154),
                        guide = "none") +
  scale_fill_gradientn(breaks = c(seq(0, 25, by = 5), 98, 154),
                       limits = c(0, 154), 
                       colors = c("gray", rev(trawl_pal[5:length(trawl_pal)]))) +
  scale_x_continuous(breaks = map_layers$lon.breaks) +
  scale_y_continuous(breaks = map_layers$lat.breaks) +
  labs(title = "2024 BBRKC Collaborative Pot Sampling", subtitle = "Total BBRKC") + #,
  # caption = "* preliminary information") +
  
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
        plot.subtitle = element_text(size = 12)) #, 
# plot.caption = element_text(hjust = 0, face = "italic"))

ggsave(plot = total.map, "./Figures_Report//Fig23.png", height = 7, width = 10, units = "in")

# Figs 24-29. RKC distributions placed on temperature surfaces, in 3-panel format =========
        #(CPS1, CPS2, 2024 summer trawl), by demographic: 24) legal males; 
                                                        # 25) mature-size males; 
                                                        # 26) immature-sized males; 
                                                        # 27) mature females; 
                                                        # 28) immature females.
# NMFS Data 

# This code is from TechReport_2023_processing.R
# I've put the necessary functions here and used a wrapper to create the above figures

# The original function for these figures is on line 433. 

# NOAA Plots First
# 2023
krigTemp_CPUE_2023_MM<-temp_map_ebs_nbs(NMFS_Hauls, 2023, cpue_data = RKC_EBS, mat_sex = "Mature Male")
krigTemp_CPUE_2023_LM<-temp_map_ebs_nbs(NMFS_Hauls, 2023, cpue_data = RKC_EBS, mat_sex = "Legal Male")
krigTemp_CPUE_2023_IM<-temp_map_ebs_nbs(NMFS_Hauls, 2023, cpue_data = RKC_EBS, mat_sex = "Immature Male")
krigTemp_CPUE_2023_MF<-temp_map_ebs_nbs(NMFS_Hauls, 2023, cpue_data = RKC_EBS, mat_sex = "Mature Female")
krigTemp_CPUE_2023_IF<-temp_map_ebs_nbs(NMFS_Hauls, 2023, cpue_data = RKC_EBS, mat_sex = "Immature Female")

# 2024
krigTemp_CPUE_2024_MM<-temp_map_ebs_nbs(NMFS_Hauls, 2024, cpue_data = RKC_EBS, mat_sex = "Mature Male")
krigTemp_CPUE_2024_LM<-temp_map_ebs_nbs(NMFS_Hauls, 2024, cpue_data = RKC_EBS, mat_sex = "Legal Male")
krigTemp_CPUE_2024_IM<-temp_map_ebs_nbs(NMFS_Hauls, 2024, cpue_data = RKC_EBS, mat_sex = "Immature Male")
krigTemp_CPUE_2024_MF<-temp_map_ebs_nbs(NMFS_Hauls, 2024, cpue_data = RKC_EBS, mat_sex = "Mature Female")
krigTemp_CPUE_2024_IF<-temp_map_ebs_nbs(NMFS_Hauls, 2024, cpue_data = RKC_EBS, mat_sex = "Immature Female")


# CPS Figures 

CPS1_COUNTS <- readRDS("C:/GitHub/CPS2/Data/CPS1_COUNTS.rds")

# CPS1 Function for plotting discrete temperatures with ONLY pot counts

cps1_temp_cpue_plots<-function(cpue_mapdat = CPS1_COUNTS, # CPUE sf object
                               plot_df = CPS1_TEMP, # Kriged temperature df
                               mat_sex, # e.g., "Mature male"
                               map_layers=map_layers, 
                               plot_boundary = plot.boundary,
                               raw_temp_sf = temploggers){
  
  # Filter Catch data to mat_sex designation
  mapdat<-cpue_mapdat %>%
    filter(MAT_SEX == mat_sex)
  
  # Set up temp breaks and labels
  temp_breaks <- c(-Inf, seq(-1, 8, 1), Inf)
  temp_levels <- c("<=-1", "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                   "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1")
  
  
  # Bin the temps to be discrete
  plot_df$temperature_bin <- cut(plot_df$temperature, breaks = temp_breaks, labels = temp_levels, include.lowest = TRUE)
  
  
  # Create year label
  sf::st_as_sf(data.frame(lab= paste("", 2023, "\nCPS1"), 
                          x = c(-160), y = c(55.2)),
               coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    cbind(years, st_coordinates(.)) -> year_lab
  
  
  # Plot total counts
  ggplot() +
    # Temp Tiles
    geom_tile(data = plot_df,
              aes(x = x, y = y, fill = temperature_bin), show.legend = FALSE) +
    scale_fill_manual(  # your discrete temp legend
      name = "Temperature (°C)",
      values = setNames(
        viridis::viridis(length(temp_levels), option = "H"),
        temp_levels),
      drop = FALSE) +
    new_scale("fill") +
    # Base map Layers
    geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
    geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1) +
    geom_sf(data = st_as_sf(RKCSA_sub), fill = NA, color = "red", alpha = 0.9, linewidth = 1) +
    geom_sf(data = st_as_sf(RKCSA), fill = NA, color = "red", alpha = 0.5, linewidth = 1) +
    # scale_color_manual(values = c("black", "red"),
    #                    labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
    #                    name = "") +
    #guides(color = guide_legend(nrow = 2)) +
    new_scale("color") +
    geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
    
    # add POT points
    geom_sf(data = mapdat %>% filter(VESSEL %in% c("Silver Spray", "Summer Bay")),
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
           #size = guide_legend(),
           color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2, order = 1)) +
    # start a new scale
    new_scale("fill") +
    #labs(title = "2023 BBRKC Collaborative Pot Sampling", subtitle = mat_sex) +
    # crop spatial extent, add bathymetry, set theme
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y) +
    geom_sf_text(sf::st_as_sf(data.frame(lab = c("50m", "100m"), 
                                         x = c(-161.5, -165), y = c(58.3, 56.1)),
                              coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                   sf::st_transform(crs = map.crs),
                 mapping = aes(label = lab)) +
    geom_shadowtext(year_lab,
                    mapping = aes(label = lab, x = X, y = Y), size = 8,  color = "black", bg.color = "white")+
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 10),
          # legend.text = element_text(size = 10),
          # legend.title = element_text(size = 10),
          legend.position = "none") -> map
          # legend.direction = "horizontal",
          #plot.title = element_text(face = "bold", size = 15),
          #plot.subtitle = element_text(size = 12)) -> map 
  
  return(map)
  
}

# CPS2 Function for plotting discrete temperatures with pot and trawl counts
# This function creates a figure with one legend()

cps2_temp_cpue_plots<-function(cpue_mapdat, # CPUE sf object
                               plot_df = plot_df, # Kriged temperature df
                               mat_sex, # e.g., "Mature male"
                               map_layers=map_layers, 
                               plot_boundary = plot.boundary,
                               raw_temp_sf = CPS2_sf){
  
 # Filter Catch data to mat_sex designation
   mapdat<-cpue_mapdat %>%
    filter(MAT_SEX == mat_sex)
   
   # Set up temp breaks and labels
   temp_breaks <- c(-Inf, seq(-1, 8, 1), Inf)
   temp_levels <- c("<=-1", "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                    "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1")
   
   # Create year label
   sf::st_as_sf(data.frame(lab= paste("", 2024, "\nCPS2"), 
                           x = c(-160), y = c(55.2)),
                coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
     sf::st_transform(crs = map.crs) %>%
     cbind(years, st_coordinates(.)) -> year_lab
   

  # Plot total counts
 ggplot() +
    # Temp Tiles
    geom_tile(data = plot_df,
             aes(x = x, y = y, fill = temperature_bin)) +
   scale_fill_manual(  # discrete temp legend
     name = "Temperature (°C)",
     values = setNames(
       viridis::viridis(length(temp_levels), option = "H"),
       temp_levels),
     drop = FALSE
   ) +
   guides(
     fill = guide_legend(
       title = "Temperature (°C)", title.position = "top", order =3)) +
   new_scale("fill") +
   # Base map Layers
    geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
    geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1) + 
    geom_sf(data = st_as_sf(CPS1_bound), fill = NA, mapping = aes(color = "black"), linewidth = 1) + #this is hacky, but removing it messes the map up.....
    geom_sf(data = st_as_sf(RKCSA_sub), mapping = aes(color = "red"), fill = NA, alpha = 0.9, linewidth = 1) +
    geom_sf(data = st_as_sf(RKCSA), fill = NA,  color = "red", alpha =0.5, linewidth = 1) +
    scale_color_manual(values = c("black", "red"),
                       labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                       name = "") +
    #guides(color = guide_legend(nrow = 2)) +
   new_scale("color") +
   geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
    
    # add POT points
    geom_sf(data = mapdat %>% filter(VESSEL %in% c("Arctic Lady", "Seabrooke")),
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
           #size = guide_legend(),
           color = guide_legend(title = "Pot Count", title.position = "top", nrow = 2, order = 1)) +
    # start a new scale
    new_scale("fill") +
    new_scale("shape") +
    new_scale("size") +
    
    # add TRAWL points
    geom_sf(data = mapdat %>% filter(VESSEL %in% c("Vesteraalen")),
            mapping = aes(size = COUNT, fill = COUNT, shape = COUNT == 0), 
            alpha = 0.5, colour = "black") +
    scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none") +
    scale_size_continuous(range = c(2, 10), 
                          limits = c(0, 154),
                          breaks = c(seq(0, 25, by = 5), 98, 154),
                          guide = "none") +
    scale_fill_gradientn(breaks = c(seq(0, 25, by = 5), 98, 154),
                         limits = c(0, 154), 
                         colors = c("gray", rev(trawl_pal[5:length(trawl_pal)]))) +
    scale_x_continuous(breaks = map_layers$lon.breaks) +
    scale_y_continuous(breaks = map_layers$lat.breaks) +
    #labs(title = mat_sex) + #,
    # caption = "* preliminary information") +
    
    # TRAWL legend
    guides(size = guide_legend(title = "Trawl Count", title.position = "top", nrow = 2, 
                               override.aes = list(shape = c(8, rep(24, 7)))),
           fill = guide_legend(title = "Trawl Count", title.position = "top"),
           color = guide_legend(nrow = 2,  order = 2)) +
   new_scale("fill") +
    # crop spatial extent, add bathymetry, set theme
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y) +
    geom_sf_text(sf::st_as_sf(data.frame(lab = c("50m", "100m"), 
                                         x = c(-161.5, -165), y = c(58.3, 56.1)),
                              coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                   sf::st_transform(crs = map.crs),
                 mapping = aes(label = lab)) +
   geom_shadowtext(data = st_drop_geometry(year_lab),
                   aes(x = X, y = Y, label = lab),
                   size = 8, color = "black", bg.color = "white", inherit.aes = FALSE)+
   theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.position = "right",
          legend.direction = "horizontal",
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 12)) -> map 
  
  return(map)
  
}


#Function Calls 
#CPS1
Temp_Count_CPS1_MM<-cps1_temp_cpue_plots(cpue_mapdat = CPS1_COUNTS, # CPUE sf object
                                         plot_df = CPS1_TEMP, # Kriged temperature df
                                         mat_sex = "Mature male", # e.g., "Mature male"
                                         map_layers=map_layers,
                                         plot_boundary = plot.boundary,
                                         raw_temp_sf = temploggers)

Temp_Count_CPS1_LM<-cps1_temp_cpue_plots(cpue_mapdat = CPS1_COUNTS, # CPUE sf object
                                         plot_df = CPS1_TEMP, # Kriged temperature df
                                         mat_sex = "Legal male", # e.g., "Mature male"
                                         map_layers=map_layers, 
                                         plot_boundary = plot.boundary,
                                         raw_temp_sf = temploggers)

Temp_Count_CPS1_IM<-cps1_temp_cpue_plots(cpue_mapdat = CPS1_COUNTS, # CPUE sf object
                                         plot_df = CPS1_TEMP, # Kriged temperature df
                                         mat_sex = "Immature male", # e.g., "Mature male"
                                         map_layers=map_layers,
                                         plot_boundary = plot.boundary,
                                         raw_temp_sf = temploggers)

Temp_Count_CPS1_MF<-cps1_temp_cpue_plots(cpue_mapdat = CPS1_COUNTS, # CPUE sf object
                                         plot_df = CPS1_TEMP, # Kriged temperature df
                                         mat_sex = "Mature female", # e.g., "Mature male"
                                         map_layers=map_layers,
                                         plot_boundary = plot.boundary,
                                         raw_temp_sf = temploggers)

Temp_Count_CPS1_IF<-cps1_temp_cpue_plots(cpue_mapdat = CPS1_COUNTS, # CPUE sf object
                                         plot_df = CPS1_TEMP, # Kriged temperature df
                                         mat_sex = "Immature female", # e.g., "Mature male"
                                         map_layers=map_layers,
                                         plot_boundary = plot.boundary,
                                         raw_temp_sf = temploggers)

#CPS2
Temp_Count_CPS2_MM<-cps2_temp_cpue_plots(
  cpue_mapdat=cpue_mapdat, # CPUE sf object
  plot_df = krig_cps2, # Kriged temperature df
  mat_sex = "Mature male", # e.g., "Mature male"
  map_layers=map_layers,
  plot_boundary = plot.boundary,
  raw_temp_sf = CPS2_sf)

Temp_Count_CPS2_LM<-cps2_temp_cpue_plots(
  cpue_mapdat=cpue_mapdat, # CPUE sf object
  plot_df = krig_cps2, # Kriged temperature df
  mat_sex = "Legal male", # e.g., "Mature male"
  map_layers=map_layers,
  plot_boundary = plot.boundary,
  raw_temp_sf = CPS2_sf)

Temp_Count_CPS2_IM<-cps2_temp_cpue_plots(
  cpue_mapdat=cpue_mapdat, # CPUE sf object
  plot_df = krig_cps2, # Kriged temperature df
  mat_sex = "Immature male", # e.g., "Mature male"
  map_layers=map_layers,
  plot_boundary = plot.boundary,
  raw_temp_sf = CPS2_sf)

Temp_Count_CPS2_MF<-cps2_temp_cpue_plots(
  cpue_mapdat=cpue_mapdat, # CPUE sf object
  plot_df = krig_cps2, # Kriged temperature df
  mat_sex = "Mature female", # e.g., "Mature male"
  map_layers=map_layers,
  plot_boundary = plot.boundary,
  raw_temp_sf = CPS2_sf)

Temp_Count_CPS2_IF<-cps2_temp_cpue_plots(
  cpue_mapdat=cpue_mapdat, # CPUE sf object
  plot_df = krig_cps2, # Kriged temperature df
  mat_sex = "Immature female", # e.g., "Mature male"
  map_layers=map_layers,
  plot_boundary = plot.boundary,
  raw_temp_sf = CPS2_sf)

make_4panel_plot <- function(p1, p2, p3, p4, title = NULL) {
  
  top_row <- p1 + p2 + p3 + 
    plot_layout(ncol = 3)
  
  # Combine CPS plots with an empty spacer to help center them
  bottom_row <- plot_spacer() + p4 + plot_spacer() +
    plot_layout(ncol = 3, widths = c(0.25, 2, 0.25))
  
  # Stack top and bottom
  final_layout <- top_row / bottom_row + 
    plot_layout(heights = c(1, 1.5)) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))
     
  return(final_layout)
}

MM<-make_4panel_plot(krigTemp_CPUE_2023_MM, krigTemp_CPUE_2024_MM, Temp_Count_CPS1_MM, Temp_Count_CPS2_MM, title = "Mature Males")
ggsave("Figures_Report/Fig25.png", plot = MM, width = 15, height = 10, dpi = 300)

LM<-make_4panel_plot(krigTemp_CPUE_2023_LM, krigTemp_CPUE_2024_LM, Temp_Count_CPS1_LM, Temp_Count_CPS2_LM, title = "Legal Males")
ggsave("Figures_Report/Fig24.png", plot = LM, width = 15, height = 10, dpi = 300)

IM<-make_4panel_plot(krigTemp_CPUE_2023_IM, krigTemp_CPUE_2024_IM, Temp_Count_CPS1_IM, Temp_Count_CPS2_IM, title = "Immature Males")
ggsave("Figures_Report/Fig26.png", plot = IM, width = 15, height = 10, dpi = 300)

MF<-make_4panel_plot(krigTemp_CPUE_2023_MF, krigTemp_CPUE_2024_MF, Temp_Count_CPS1_MF, Temp_Count_CPS2_MF, title = "Mature Females")
ggsave("Figures_Report/Fig27.png", plot = MF, width = 15, height = 10, dpi = 300)

IF<-make_4panel_plot(krigTemp_CPUE_2023_IF, krigTemp_CPUE_2024_IF, Temp_Count_CPS1_IF, Temp_Count_CPS2_IF, title = "Immature Females")
ggsave("Figures_Report/Fig28.png", plot = IF, width = 15, height = 10, dpi = 300)

# Fig 30.  (RKC geographic centers of distribution (no temperatures), by demographic ==============
          # 6 panels = legal males; sublegal males; mature-size males; immature-sized males; mature females; immature females.
          # Each panel contains 5 points = 2022 summer trawl, CPS1, 2023 summer trawl, CPS2, 2024 summer trawl.
          # I have computed CPS values and will endeavor to do the trawl-survey points, as well (noting that the latter will compute summer values only for the stations within the CPS footprint)



# Fig 31-35. Bycatch distributions on temperature surfaces, single-panle by species ==========
  #31) Tanners; 32) opies; 33) cod; 34) yellowfin sole; 35) rock sole.)))
# Fig. 36, perhaps? Length-frequency distributions for RKC captured in the two bait experiments====
    # 2-panel = hanging bait; bait bags.