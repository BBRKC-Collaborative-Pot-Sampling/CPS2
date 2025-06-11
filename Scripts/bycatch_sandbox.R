
make_bycatch_map <- function(species,
                             species_lab,
                             pot_sf = bycatch_pots,
                             trawl_sf = bycatch_trawl,
                             temp_df = plot_df,
                             map_layers = map_layers,
                             plot_boundary = plot.boundary) {
  pot_pal   <- viridis(7, option = "C")
  trawl_pal <- viridis(7, option = "E")
  
  temp_breaks <- c(-Inf, seq(-1, 8, 1), Inf)
  temp_levels <- c("<=-1", "-0.9-0", "0.1-1", "1.1-2", "2.1-3", "3.1-4",
                   "4.1-5", "5.1-6", "6.1-7", "7.1-8", ">8.1")
  
  sf::st_as_sf(data.frame(lab= paste("", 2024, "\nCPS2"), 
                          x = c(-160), y = c(55.2)),
               coords = c("x", "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    cbind(years, st_coordinates(.)) -> year_lab
  
  pot_vals   <- pot_sf %>% dplyr::select(VESSEL, SPN, species, geometry)
  trawl_vals <- trawl_sf %>% dplyr::select(VESSEL, STATION, species, geometry)
  
  max_trawl <- max(trawl_vals[[species]], na.rm = TRUE)
  max_pot   <- max(pot_vals[[species]], na.rm = TRUE)
  
  # Breaks & aesthetics
  if (max_trawl == 0) {
    trawl_breaks <- 0
    trawl_shapes <- 8
    trawl_sizes  <- 1
    trawl_fills  <- "gray"
  } else {
    trawl_breaks <- round(seq(0, max_trawl, length.out = 6))
    if (length(trawl_breaks) < 6) trawl_breaks <- unique(c(trawl_breaks, max_trawl))
    trawl_shapes <- c(8, rep(24, length(trawl_breaks) - 1))
    vals <- seq(0, 1, length.out = length(trawl_breaks) - 1)
    trawl_sizes <- c(1.5, 3 + (vals^1.8) * 10)
    trawl_fills <- colorRampPalette(c("gray", rev(trawl_pal)))(length(trawl_breaks))
  }
  
  if (max_pot == 0) {
    pot_breaks <- 0
    pot_shapes <- 4
    pot_sizes  <- 1
    pot_fills  <- "gray"
  } else {
    pot_breaks <- round(seq(0, max_pot, length.out = 6))
    if (length(pot_breaks) < 6) pot_breaks <- unique(c(pot_breaks, max_pot))
    pot_shapes <- c(4, rep(21, length(pot_breaks) - 1))
    vals <- seq(0, 1, length.out = length(pot_breaks) - 1)
    pot_sizes <- c(1.5, 3 + (vals^1.8) * 10)
    pot_fills <- colorRampPalette(c("gray", rev(pot_pal)))(length(pot_breaks))
  }
  
  make_size_scale <- function(max_val, breaks) {
    if (max_val == 0) {
      scale_size_continuous(range = c(1.5, 1.5), limits = c(0, 1), guide = "legend",
                            breaks = NULL, labels = NULL)
    } else {
      scale_size_continuous(range = c(1.5, 10), limits = c(0, max_val),
                            breaks = breaks, guide = "none")
    }
  }
  
  # Build the plot
  map2 <- ggplot() +
    geom_tile(data = temp_df, aes(x = x, y = y, fill = temperature_bin)) +
    scale_fill_manual(name = "Temperature (°C)",
                      values = setNames(viridis::viridis(length(temp_levels), option = "H"), temp_levels),
                      drop = FALSE) +
    guides(fill = guide_legend(title = "Temperature (°C)", title.position = "top", order = 3)) +
    new_scale("fill") +
    geom_sf(data = st_transform(map_layers$bathymetry, map.crs), color = alpha("grey70")) +
    geom_sf(data = st_as_sf(CPS1_bound), fill = NA, color = "black", linewidth = 1) +
    geom_sf(data = st_as_sf(RKCSA_sub), fill = NA, color = "red", alpha = 0.9, linewidth = 1) +
    geom_sf(data = st_as_sf(RKCSA), fill = NA, color = "red", alpha = 0.5, linewidth = 1) +
    new_scale("color") +
    scale_color_manual(values = c("black", "red"),
                       labels = c("CPS Survey Boundary", "Red King Crab Savings Area"),
                       name = "") +
    new_scale("color") +
    geom_sf(data = st_transform(map_layers$akland, map.crs), fill = "grey80") +
    
    # Pot points
    geom_sf(data = pot_vals,
            mapping = aes(size = .data[[species]],
                          fill = .data[[species]],
                          shape = .data[[species]] == 0),
            alpha = 0.5, colour = "black") +
    scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 21), guide = "none") +
    make_size_scale(max_pot, pot_breaks) +
    scale_fill_gradientn(limits = c(0, max_pot), breaks = pot_breaks,
                         colors = pot_fills, guide = "none") +
    guides(fill = guide_legend(title = "Pot Count", title.position = "top",
                               override.aes = list(shape = pot_shapes,
                                                   size = pot_sizes,
                                                   fill = pot_fills))) +
    new_scale("fill") +
    new_scale("shape") +
    new_scale("size") +
    
    # Trawl points
    geom_sf(data = trawl_vals,
            mapping = aes(size = .data[[species]],
                          fill = .data[[species]],
                          shape = .data[[species]] == 0),
            alpha = 0.5, colour = "black") +
    scale_shape_manual(values = c('TRUE' = 8, 'FALSE' = 24), guide = "none") +
    make_size_scale(max_trawl, trawl_breaks) +
    scale_fill_gradientn(limits = c(0, max_trawl), breaks = trawl_breaks,
                         colors = trawl_fills, guide = "none") +
    new_scale("fill") +
    new_scale("shape") +
    new_scale("size") +
    
    scale_x_continuous(breaks = map_layers$lon.breaks) +
    scale_y_continuous(breaks = map_layers$lat.breaks) +
    coord_sf(xlim = plot_boundary$x, ylim = plot_boundary$y) +
    geom_sf_text(sf::st_as_sf(data.frame(lab = c("50m", "100m"),
                                         x = c(-161.5, -165), y = c(58.3, 56.1)),
                              coords = c("x", "y"), crs = sf::st_crs(4326)) %>%
                   sf::st_transform(crs = map.crs),
                 mapping = aes(label = lab)) +
    geom_shadowtext(data = st_drop_geometry(year_lab),
                    aes(x = X, y = Y, label = lab),
                    size = 8, color = "black", bg.color = "white", inherit.aes = FALSE) +
    labs(title = paste0("2024 CPS2 Bycatch: ", species_lab)) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.position = "right",
          legend.direction = "horizontal",
          plot.title = element_text(face = "bold", size = 15))
  
  # Cool, so the only way to get this 
  
  if (max_trawl > 0) {
    # 1) make a tiny legend‐only data.frame
    legend_df <- data.frame(
      x     = NA_real_,
      y     = NA_real_,
      count = trawl_breaks
    )
    
    map2 <- map2 +
      # 2) draw invisible points *mapped* to count for size+fill+shape
      geom_point(
        data        = legend_df,
        aes(
          x     = x,
          y     = y,
          size  = count,
          fill  = count,
          shape = factor(count)      # <-- map shape to count!
        ),
        alpha       = 0.5,              # invisible on the map
        show.legend = TRUE
      ) +
      # 3) tell ggplot which triangle (and circle) to use at each break
      scale_shape_manual(
        name   = "Trawl Count",
        breaks = trawl_breaks,
        values = trawl_shapes,        # your c(8,24,24,…)
        guide  = guide_legend(title.position = "top")
      ) +
      # 4) size legend—override its aes so that each key shows the right size & fill
      scale_size_continuous(
        name   = "Trawl Count",
        breaks = trawl_breaks,
        guide  = guide_legend(
          override.aes = list(
            shape = trawl_shapes,
            fill  = trawl_fills,
            size  = trawl_sizes
          )
        )
      ) +
      # 5) register a fill scale (but hide its legend)
      scale_fill_gradientn(
        colours = trawl_fills,
        breaks  = trawl_breaks,
        guide   = "none"
      )
  }
  
  
  return(map2)
}



Tanner_map <- make_bycatch_map(
  species        = "Tanner",
  species_lab    = "Tanner Crab",
  pot_sf         = bycatch_pots,
  trawl_sf       = bycatch_trawl,
  temp_df        = plot_df,
  map_layers     = map_layers,
  plot_boundary  = plot.boundary
)

Snow_map <- make_bycatch_map(
  species        = "Snow",
  species_lab    = "Snow Crab",
  pot_sf         = bycatch_pots,
  trawl_sf       = bycatch_trawl,
  temp_df        = plot_df,
  map_layers     = map_layers,
  plot_boundary  = plot.boundary
)

cod_map <- make_bycatch_map(
  species        = "PacificCod",
  species_lab    = "Pacific Cod",
  pot_sf         = bycatch_pots,
  trawl_sf       = bycatch_trawl,
  temp_df        = plot_df,
  map_layers     = map_layers,
  plot_boundary  = plot.boundary
)

yellow_map <- make_bycatch_map(
  species        = "YellowfinSole",
  species_lab    = "Yellowfin Sole",
  pot_sf         = bycatch_pots,
  trawl_sf       = bycatch_trawl,
  temp_df        = plot_df,
  map_layers     = map_layers,
  plot_boundary  = plot.boundary
)

rock_map <- make_bycatch_map(
  species        = "RockSole",
  species_lab    = "Rock Sole",
  pot_sf         = bycatch_pots,
  trawl_sf       = bycatch_trawl,
  temp_df        = plot_df,
  map_layers     = map_layers,
  plot_boundary  = plot.boundary
)




