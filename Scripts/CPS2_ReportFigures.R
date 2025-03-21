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


#Figure 1: BBRKC Abundance ===========

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


# Figure 2: Closure area map ========

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


# Figure 3: CPS Planned v actual stations ============

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

# Load function to generate NMFS temperature maps for Bristol Bay
temp_map_ebs_nbs <- function(haul_ebs, years){
  
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
  
  # Plot interpolated data
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
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.position = "bottom",
          legend.direction = "horizontal",
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 12))  -> temp_map
  return(temp_map)
  
}

krigeTemp_2022<-temp_map_ebs_nbs(NMFS_Hauls, 2022)
krigeTemp_2023<-temp_map_ebs_nbs(NMFS_Hauls, 2023)
krigeTemp_2024<-temp_map_ebs_nbs(NMFS_Hauls, 2024)

krigeTemp_22_23_24<-temp_map_ebs_nbs(NMFS_Hauls, 2022:2024)

Fig16<-grid.arrange(kreigTemp_2022, kreigTemp_2023,kreigTemp_2024, nrow = 1, ncol = 3)


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



#Figs 18-22: Male size-frequency distribution by shell condition ==========

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



male_shell_tabplot(male_abundance_NMFS)
