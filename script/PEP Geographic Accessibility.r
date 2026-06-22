
# 1. LOAD PACKAGES
library(sf)
library(terra)
library(raster)
library(gdistance)
library(dplyr)
library(readr)
library(exactextractr)
library(ggplot2)
library(ggspatial)
library(scales)


# 2. USER INPUTS

level2_3_file         <- "~/Documents/RABIES_ANALYSIS/inputs/khfa_level_2_3_cleaned.csv"
level4_file           <- "~/Documents/RABIES_ANALYSIS/inputs_00/khfa_level_4_cleaned.gpkg"

population_file       <- "~/Documents/RABIES_ANALYSIS/inputs/KENYA POPULATION GRID3 DATA.tif"
urbanisation_file     <- "~/Documents/RABIES_ANALYSIS/inputs_00/DEGREE OF URBANISATION 2026 RECLASSIFIED.tif"

counties_file         <- "~/Documents/RABIES_ANALYSIS/inputs/counties.gpkg"
walk_friction_file    <- "~/Documents/RABIES_ANALYSIS/inputs/friction_walking_kenya.tif"
motor_friction_file   <- "~/Documents/RABIES_ANALYSIS/inputs/friction_motorized_kenya.tif"

output_dir <- "~/Documents/RABIES_ANALYSIS/outputs/paper 4 outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

points_id_col     <- "facility"
counties_name_col <- "county"
lon_col           <- "longitude"
lat_col           <- "latitude"
keph_col          <- "keph_level"

target_epsg <- "EPSG:32737"
output_epsg <- "EPSG:4326"

service_thresholds <- c(30, 60, 120)

# ===============================
# 3. HELPER FUNCTIONS
# ===============================

clean_map_theme <- theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text        = element_blank(),
  axis.ticks       = element_blank(),
  axis.title       = element_blank(),
  plot.title       = element_text(face = "bold", hjust = 0.5),
  legend.title     = element_text(face = "bold")
)

save_plot <- function(p, filename, width = 10, height = 8) {
  ggsave(
    filename = file.path(output_dir, filename),
    plot = p,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

write_wgs84_raster <- function(rast_obj, filename, method = "bilinear") {
  rast_wgs84 <- terra::project(rast_obj, output_epsg, method = method)
  terra::writeRaster(rast_wgs84, filename, overwrite = TRUE)
}

classify_time <- function(time_raster) {
  terra::ifel(
    time_raster <= 30, 1,
    terra::ifel(
      time_raster > 30 & time_raster <= 60, 2,
      terra::ifel(
        time_raster > 60 & time_raster <= 120, 3,
        terra::ifel(time_raster > 120, 4, NA)
      )
    )
  )
}

safe_exact_sum <- function(rast_obj, poly_obj) {
  if (is.null(poly_obj) || nrow(poly_obj) == 0) return(0)
  x <- exactextractr::exact_extract(rast_obj, poly_obj, "sum")
  x <- sum(as.numeric(x), na.rm = TRUE)
  ifelse(is.na(x), 0, x)
}

weighted_median <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  
  if (length(x) == 0) return(NA_real_)
  
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  
  cw <- cumsum(w) / sum(w)
  x[which(cw >= 0.5)[1]]
}

safe_weighted_median_time <- function(time_rast, weight_rast, poly_obj) {
  if (is.null(poly_obj) || nrow(poly_obj) == 0) return(NA_real_)
  
  exactextractr::exact_extract(
    time_rast,
    poly_obj,
    fun = function(values, coverage_fraction, weights) {
      weighted_median(values, weights * coverage_fraction)
    },
    weights = weight_rast
  )[[1]]
}

raster_to_plot_df <- function(rast_obj, value_name) {
  rast_wgs <- terra::project(rast_obj, output_epsg, method = "near")
  df <- as.data.frame(rast_wgs, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- value_name
  df
}

make_county_labels <- function(counties_layer) {
  suppressWarnings(st_point_on_surface(counties_layer))
}

# ===============================
# 4. OUTPUT FILES
# ===============================

walk_time_file        <- file.path(output_dir, "walk_time.tif")
motor_time_file       <- file.path(output_dir, "motor_time.tif")
level4_time_file      <- file.path(output_dir, "level4_motor_time.tif")

walk_class_file       <- file.path(output_dir, "walk_time_class.tif")
motor_class_file      <- file.path(output_dir, "motor_time_class.tif")

combined_county_csv   <- file.path(output_dir, "all_thresholds_county_population_by_access_zone.csv")
combined_national_csv <- file.path(output_dir, "all_thresholds_national_population_by_access_zone.csv")

# ===============================
# 5. LOAD DATA
# ===============================

level2_3 <- readr::read_csv(level2_3_file, show_col_types = FALSE)
level4_sf <- sf::st_read(level4_file, quiet = FALSE)

counties_sf <- sf::st_read(counties_file, quiet = FALSE)

pop_rast_original <- terra::rast(population_file)
urban_rast_original <- terra::rast(urbanisation_file)

walk_fric  <- terra::rast(walk_friction_file)
motor_fric <- terra::rast(motor_friction_file)

cat("Original population raster total:\n")
print(terra::global(pop_rast_original, "sum", na.rm = TRUE)[1, 1])

# ===============================
# 6. CHECK REQUIRED COLUMNS
# ===============================

required_point_cols <- c(points_id_col, lon_col, lat_col, keph_col)
missing_point_cols <- setdiff(required_point_cols, names(level2_3))

if (length(missing_point_cols) > 0) {
  stop("Missing columns in level2_3 CSV: ", paste(missing_point_cols, collapse = ", "))
}

if (!(counties_name_col %in% names(counties_sf))) {
  stop("County name column missing: ", counties_name_col)
}

# ===============================
# 7. CLEAN AND PROJECT DATA
# ===============================

level2_3 <- level2_3 %>%
  filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]]))

level2_3_sf <- sf::st_as_sf(
  level2_3,
  coords = c(lon_col, lat_col),
  crs = 4326,
  remove = FALSE
)

level2_3_sf <- st_transform(level2_3_sf, target_epsg)
level4_sf   <- st_transform(level4_sf, target_epsg)
counties_sf <- st_transform(counties_sf, target_epsg)

walk_fric  <- terra::project(walk_fric, target_epsg, method = "bilinear")
motor_fric <- terra::project(motor_fric, target_epsg, method = "bilinear")

motor_fric <- terra::resample(motor_fric, walk_fric, method = "bilinear")

counties_vect <- terra::vect(counties_sf)

walk_fric <- terra::crop(walk_fric, counties_vect)
walk_fric <- terra::mask(walk_fric, counties_vect)

motor_fric <- terra::crop(motor_fric, counties_vect)
motor_fric <- terra::mask(motor_fric, counties_vect)

inside_idx <- lengths(sf::st_within(level2_3_sf, counties_sf)) > 0
level2_3_sf <- level2_3_sf[inside_idx, ]

inside_l4_idx <- lengths(sf::st_within(level4_sf, counties_sf)) > 0
level4_sf <- level4_sf[inside_l4_idx, ]

if (nrow(level2_3_sf) == 0) {
  stop("No level2_3 facility points remain inside county boundaries.")
}

if (nrow(level4_sf) == 0) {
  stop("No level 4 facility points remain inside county boundaries.")
}

# ===============================
# 8. COMPUTE TRAVEL TIME SURFACES
# ===============================

walk_fric_r  <- raster::raster(walk_fric)
motor_fric_r <- raster::raster(motor_fric)

level2_3_sp <- methods::as(level2_3_sf, "Spatial")
level4_sp   <- methods::as(level4_sf, "Spatial")

walk_fric_r[walk_fric_r <= 0]   <- NA
motor_fric_r[motor_fric_r <= 0] <- NA

walk_conductance_r  <- 1 / walk_fric_r
motor_conductance_r <- 1 / motor_fric_r

walk_tr <- gdistance::transition(walk_conductance_r, mean, directions = 8)
walk_tr <- gdistance::geoCorrection(walk_tr, type = "c", multpl = FALSE, scl = FALSE)

motor_tr <- gdistance::transition(motor_conductance_r, mean, directions = 8)
motor_tr <- gdistance::geoCorrection(motor_tr, type = "c", multpl = FALSE, scl = FALSE)

cat("Computing walking travel time surface to level2_3 facilities...\n")
walk_time_r <- gdistance::accCost(walk_tr, level2_3_sp)

cat("Computing motorised travel time surface to level2_3 facilities...\n")
motor_time_r <- gdistance::accCost(motor_tr, level2_3_sp)

cat("Computing motorised travel time surface to nearest level 4 facility...\n")
level4_time_r <- gdistance::accCost(motor_tr, level4_sp)

walk_time   <- terra::rast(walk_time_r)
motor_time  <- terra::rast(motor_time_r)
level4_time <- terra::rast(level4_time_r)

names(walk_time)   <- "walk_time"
names(motor_time)  <- "motor_time"
names(level4_time) <- "level4_motor_time"

terra::crs(walk_time)   <- target_epsg
terra::crs(motor_time)  <- target_epsg
terra::crs(level4_time) <- target_epsg

walk_time   <- terra::mask(walk_time, counties_vect)
motor_time  <- terra::mask(motor_time, counties_vect)
level4_time <- terra::mask(level4_time, counties_vect)

write_wgs84_raster(walk_time, walk_time_file, method = "bilinear")
write_wgs84_raster(motor_time, motor_time_file, method = "bilinear")
write_wgs84_raster(level4_time, level4_time_file, method = "bilinear")

# ===============================
# 9. CLASSIFY WALKING AND MOTORISED TRAVEL TIME
# ===============================

walk_time_class  <- classify_time(walk_time)
motor_time_class <- classify_time(motor_time)

names(walk_time_class)  <- "walk_time_class"
names(motor_time_class) <- "motor_time_class"

write_wgs84_raster(walk_time_class, walk_class_file, method = "near")
write_wgs84_raster(motor_time_class, motor_class_file, method = "near")

# ===============================
# 10. PREPARE MAP DATA
# ===============================

counties_wgs <- st_transform(counties_sf, output_epsg)
county_labels_wgs <- make_county_labels(counties_wgs)

country_wgs <- st_union(counties_wgs) %>% st_as_sf()
level2_3_wgs <- st_transform(level2_3_sf, output_epsg)

walk_class_df <- raster_to_plot_df(walk_time_class, "class")
motor_class_df <- raster_to_plot_df(motor_time_class, "class")

class_labels <- c(
  "1" = "<=30 mins",
  "2" = "30-60 mins",
  "3" = "60-120 mins",
  "4" = ">120 mins"
)

travel_colors <- c(
  "1" = "#1a9850",
  "2" = "#fee08b",
  "3" = "#f46d43",
  "4" = "#a50026"
)

access_zone_colors <- c(
  "1" = "#1a9850",
  "2" = "#fdae61",
  "3" = "#d7191c"
)

# ===============================
# 11. MAP: LEVEL2_3 POINTS WITH COUNTIES
# ===============================

p_level2_3 <- ggplot() +
  geom_sf(data = counties_wgs, fill = NA, color = "grey30", linewidth = 0.3) +
  geom_sf_text(
    data = county_labels_wgs,
    aes(label = .data[[counties_name_col]]),
    size = 2.0,
    color = "black",
    check_overlap = TRUE
  ) +
  geom_sf(
    data = level2_3_wgs,
    aes(color = .data[[keph_col]]),
    size = 1,
    alpha = 0.85
  ) +
  annotation_north_arrow(location = "tr", which_north = "true") +
  labs(
    title = "Level 2 and 3 Facilities by KEPH Level",
    color = "KEPH level"
  ) +
  coord_sf(crs = st_crs(output_epsg)) +
  clean_map_theme

save_plot(p_level2_3, "level2_3_points_by_keph_level.png")

# ===============================
# 12. PLOT CLASSIFIED WALKING TRAVEL TIME MAP
# ===============================

p_walk_class <- ggplot() +
  geom_tile(
    data = walk_class_df,
    aes(x = x, y = y, fill = factor(class))
  ) +
  geom_sf(data = counties_wgs, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_sf_text(
    data = county_labels_wgs,
    aes(label = .data[[counties_name_col]]),
    size = 2.0,
    color = "black",
    check_overlap = TRUE
  ) +
  scale_fill_manual(
    values = travel_colors,
    labels = class_labels,
    name = "Walking time"
  ) +
  annotation_north_arrow(location = "tr", which_north = "true") +
  labs(title = "Classified Walking Travel Time to Nearest Level 2/3 Facility") +
  coord_sf(crs = st_crs(output_epsg)) +
  clean_map_theme

save_plot(p_walk_class, "classified_walking_travel_time.png")

# ===============================
# 13. PLOT CLASSIFIED MOTORISED TRAVEL TIME MAP
# ===============================

p_motor_class <- ggplot() +
  geom_tile(
    data = motor_class_df,
    aes(x = x, y = y, fill = factor(class))
  ) +
  geom_sf(data = counties_wgs, fill = NA, color = "grey30", linewidth = 0.2) +
  geom_sf_text(
    data = county_labels_wgs,
    aes(label = .data[[counties_name_col]]),
    size = 2.0,
    color = "black",
    check_overlap = TRUE
  ) +
  scale_fill_manual(
    values = travel_colors,
    labels = class_labels,
    name = "Motorised time"
  ) +
  annotation_north_arrow(location = "tr", which_north = "true") +
  labs(title = "Classified Motorised Travel Time to Nearest Level 2/3 Facility") +
  coord_sf(crs = st_crs(output_epsg)) +
  clean_map_theme

save_plot(p_motor_class, "classified_motorised_travel_time.png")

# ===============================
# 14. PREPARE POPULATION + URBANISATION EXTRACTION
# ===============================

pop_rast_r <- raster::raster(pop_rast_original)
names(pop_rast_r) <- "population"

urban_rast_original <- terra::project(
  urban_rast_original,
  pop_rast_original,
  method = "near"
)

rural_pop_rast <- terra::ifel(urban_rast_original == 1, pop_rast_original, NA)
urban_pop_rast <- terra::ifel(urban_rast_original == 2, pop_rast_original, NA)

rural_pop_rast_r <- raster::raster(rural_pop_rast)
urban_pop_rast_r <- raster::raster(urban_pop_rast)

# ===============================
# ALIGN LEVEL 4 TIME RASTER TO POPULATION GRID
# ===============================

level4_time_pop_crs <- terra::project(
  level4_time,
  pop_rast_original,
  method = "bilinear"
)


level4_time_pop_crs <- terra::resample(
  level4_time_pop_crs,
  pop_rast_original,
  method = "bilinear"
)

level4_time_pop_r <- raster::raster(level4_time_pop_crs)
pop_crs <- sf::st_crs(pop_rast_original)

if (is.na(pop_crs)) {
  stop("Population raster has no valid CRS.")
}

counties_pop_crs <- st_transform(counties_sf, pop_crs)
country_union_pop_crs <- st_union(counties_pop_crs) %>% st_as_sf()

total_pop_national <- safe_exact_sum(pop_rast_r, country_union_pop_crs)

cat("Total population within county boundaries:\n")
print(total_pop_national)

# ===============================
# 15. RUN ACCESS-ZONE SCENARIOS
# ===============================

all_national_results <- list()
all_county_results <- list()

for (threshold in service_thresholds) {
  
  cat("\nProcessing service threshold:", threshold, "minutes\n")
  
  threshold_label <- paste0(threshold, "min")
  
  threshold_dir <- file.path(output_dir, threshold_label)
  dir.create(threshold_dir, recursive = TRUE, showWarnings = FALSE)
  
  # -------------------------------
  # 15.1 Create access zones
  # -------------------------------
  
  access_zones <- terra::ifel(
    walk_time <= threshold,
    1,
    terra::ifel(
      walk_time > threshold & motor_time <= threshold,
      2,
      terra::ifel(
        walk_time > threshold & motor_time > threshold,
        3,
        NA
      )
    )
  )
  
  names(access_zones) <- "access_zone"
  
  access_zones_file <- file.path(threshold_dir, paste0(threshold_label, "_access_zones.tif"))
  write_wgs84_raster(access_zones, access_zones_file, method = "near")
  
  # -------------------------------
  # 15.2 Convert access zones to polygons
  # -------------------------------
  
  access_zone_poly <- terra::as.polygons(
    access_zones,
    dissolve = TRUE,
    na.rm = TRUE
  )
  
  access_zone_sf <- st_as_sf(access_zone_poly)
  st_crs(access_zone_sf) <- st_crs(target_epsg)
  
  access_zone_sf <- st_make_valid(access_zone_sf) %>%
    rename(zone = access_zone) %>%
    mutate(
      zone_label = case_when(
        zone == 1 ~ "Zone 1: Well-served",
        zone == 2 ~ "Zone 2: Motorised-dependent",
        zone == 3 ~ "Zone 3: Critically underserved",
        TRUE ~ NA_character_
      )
    )
  
  access_zone_pop_crs <- st_transform(access_zone_sf, pop_crs)
  access_zone_wgs <- st_transform(access_zone_sf, output_epsg)
  
  zone3_pop_crs <- access_zone_pop_crs %>%
    filter(zone == 3)
  
  zone3_wgs <- access_zone_wgs %>%
    filter(zone == 3)
  
  zone3_dissolved_wgs <- zone3_wgs %>%
    summarise(zone = 3, geometry = st_union(geometry)) %>%
    st_make_valid()
  
  # -------------------------------
  # 15.3 National population by zone
  # -------------------------------
  
  national_zone_list <- vector("list", nrow(access_zone_pop_crs))
  
  national_zone3_pop <- safe_exact_sum(pop_rast_r, zone3_pop_crs)
  national_zone3_rural_pop <- safe_exact_sum(rural_pop_rast_r, zone3_pop_crs)
  national_zone3_urban_pop <- safe_exact_sum(urban_pop_rast_r, zone3_pop_crs)
  
  national_zone3_median_time_level4 <- safe_weighted_median_time(
    level4_time_pop_r,
    pop_rast_r,
    zone3_pop_crs
  )
  
  for (z in seq_len(nrow(access_zone_pop_crs))) {
    
    zone_poly <- access_zone_pop_crs[z, ]
    pop_sum <- safe_exact_sum(pop_rast_r, zone_poly)
    
    national_zone_list[[z]] <- data.frame(
      threshold_min = threshold,
      zone = zone_poly$zone,
      zone_label = zone_poly$zone_label,
      population = pop_sum
    )
  }
  
  national_zone_df <- bind_rows(national_zone_list) %>%
    group_by(threshold_min, zone, zone_label) %>%
    summarise(population = sum(population, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      pct_total_population = (population / total_pop_national) * 100,
      
      zone3_rural_population = ifelse(zone == 3, national_zone3_rural_pop, NA_real_),
      zone3_urban_population = ifelse(zone == 3, national_zone3_urban_pop, NA_real_),
      
      pct_zone3_population_rural = ifelse(
        zone == 3 & national_zone3_pop > 0,
        (national_zone3_rural_pop / national_zone3_pop) * 100,
        NA_real_
      ),
      
      pct_zone3_population_urban = ifelse(
        zone == 3 & national_zone3_pop > 0,
        (national_zone3_urban_pop / national_zone3_pop) * 100,
        NA_real_
      ),
      
      median_motorised_time_to_nearest_level4_for_zone3_population = ifelse(
        zone == 3,
        national_zone3_median_time_level4,
        NA_real_
      )
    ) %>%
    arrange(zone)
  
  write_csv(
    national_zone_df,
    file.path(threshold_dir, paste0(threshold_label, "_national_population_by_access_zone.csv"))
  )
  
  all_national_results[[threshold_label]] <- national_zone_df
  
  # -------------------------------
  # 15.4 County population by zone
  # -------------------------------
  
  county_zone_list <- list()
  counter <- 1
  
  for (i in seq_len(nrow(counties_pop_crs))) {
    
    county_poly <- counties_pop_crs[i, ]
    county_name <- county_poly[[counties_name_col]][1]
    
    county_total_pop <- safe_exact_sum(pop_rast_r, county_poly)
    
    county_zone3_intersection <- suppressWarnings(st_intersection(county_poly, zone3_pop_crs))
    county_zone3_intersection <- st_make_valid(county_zone3_intersection)
    
    county_zone3_pop <- safe_exact_sum(pop_rast_r, county_zone3_intersection)
    county_zone3_rural_pop <- safe_exact_sum(rural_pop_rast_r, county_zone3_intersection)
    county_zone3_urban_pop <- safe_exact_sum(urban_pop_rast_r, county_zone3_intersection)
    
    county_zone3_median_time_level4 <- safe_weighted_median_time(
      level4_time_pop_r,
      pop_rast_r,
      county_zone3_intersection
    )
    
    for (z in c(1, 2, 3)) {
      
      zone_poly <- access_zone_pop_crs %>%
        filter(zone == z)
      
      zone_label <- case_when(
        z == 1 ~ "Zone 1: Well-served",
        z == 2 ~ "Zone 2: Motorised-dependent",
        z == 3 ~ "Zone 3: Critically underserved"
      )
      
      if (nrow(zone_poly) == 0) {
        pop_zone_county <- 0
      } else {
        intersection_poly <- suppressWarnings(st_intersection(county_poly, zone_poly))
        intersection_poly <- st_make_valid(intersection_poly)
        pop_zone_county <- safe_exact_sum(pop_rast_r, intersection_poly)
      }
      
      county_zone_list[[counter]] <- data.frame(
        threshold_min = threshold,
        county = county_name,
        zone = z,
        zone_label = zone_label,
        population = pop_zone_county,
        county_total_population = county_total_pop,
        pct_county_population = ifelse(
          county_total_pop > 0,
          (pop_zone_county / county_total_pop) * 100,
          NA
        ),
        
        zone3_rural_population = ifelse(z == 3, county_zone3_rural_pop, NA_real_),
        zone3_urban_population = ifelse(z == 3, county_zone3_urban_pop, NA_real_),
        
        pct_zone3_population_rural = ifelse(
          z == 3 & county_zone3_pop > 0,
          (county_zone3_rural_pop / county_zone3_pop) * 100,
          NA_real_
        ),
        
        pct_zone3_population_urban = ifelse(
          z == 3 & county_zone3_pop > 0,
          (county_zone3_urban_pop / county_zone3_pop) * 100,
          NA_real_
        ),
        
        median_motorised_time_to_nearest_level4_for_zone3_population = ifelse(
          z == 3,
          county_zone3_median_time_level4,
          NA_real_
        )
      )
      
      counter <- counter + 1
    }
  }
  
  county_zone_df <- bind_rows(county_zone_list)
  
  write_csv(
    county_zone_df,
    file.path(threshold_dir, paste0(threshold_label, "_county_population_by_access_zone.csv"))
  )
  
  all_county_results[[threshold_label]] <- county_zone_df
  
  # -------------------------------
  # 15.5 Plot national access zones
  # -------------------------------
  
  access_zone_df <- raster_to_plot_df(access_zones, "zone")
  
  p_access <- ggplot() +
    geom_tile(
      data = access_zone_df,
      aes(x = x, y = y, fill = factor(zone))
    ) +
    geom_sf(data = counties_wgs, fill = NA, color = "grey30", linewidth = 0.2) +
    geom_sf_text(
      data = county_labels_wgs,
      aes(label = .data[[counties_name_col]]),
      size = 2.0,
      color = "black",
      check_overlap = TRUE
    ) +
    scale_fill_manual(
      values = access_zone_colors,
      labels = c(
        "1" = "Zone 1",
        "2" = "Zone 2",
        "3" = "Zone 3"
      ),
      name = "Access zone"
    ) +
    annotation_north_arrow(location = "tr", which_north = "true") +
    labs(title = paste0("National Access Zones: ", threshold, "-Minute Threshold")) +
    coord_sf(crs = st_crs(output_epsg)) +
    clean_map_theme
  
  save_plot(
    p_access,
    file.path(threshold_label, paste0(threshold_label, "_access_zones.png"))
  )
  
  # -------------------------------
  # 15.6 County choropleths: Zone 3 %
  # -------------------------------
  
  zone3_county <- county_zone_df %>%
    filter(zone == 3)
  
  counties_zone3_wgs <- counties_wgs %>%
    left_join(zone3_county, by = setNames("county", counties_name_col))
  
  county_labels_zone3 <- make_county_labels(counties_zone3_wgs)
  
  p_zone3_pct <- ggplot() +
    geom_sf(
      data = counties_zone3_wgs,
      aes(fill = pct_county_population),
      color = "grey30",
      linewidth = 0.2
    ) +
    geom_sf_text(
      data = county_labels_zone3,
      aes(label = .data[[counties_name_col]]),
      size = 2.0,
      color = "black",
      check_overlap = TRUE
    ) +
    scale_fill_gradient(
      low = "#fee08b",
      high = "#d7191c",
      labels = function(x) paste0(round(x, 1), "%"),
      name = "% population\nin Zone 3"
    ) +
    annotation_north_arrow(location = "tr", which_north = "true") +
    labs(title = paste0("% Population in Zone 3 by County: ", threshold, "-Minute Threshold")) +
    coord_sf(crs = st_crs(output_epsg)) +
    clean_map_theme
  
  save_plot(
    p_zone3_pct,
    file.path(threshold_label, paste0(threshold_label, "_zone3_percent_by_county.png"))
  )
  
  # -------------------------------
  # 15.7 County choropleths: Zone 3 population number
  # -------------------------------
  
  p_zone3_pop <- ggplot() +
    geom_sf(
      data = counties_zone3_wgs,
      aes(fill = population),
      color = "grey30",
      linewidth = 0.2
    ) +
    geom_sf_text(
      data = county_labels_zone3,
      aes(label = .data[[counties_name_col]]),
      size = 2.0,
      color = "black",
      check_overlap = TRUE
    ) +
    scale_fill_gradient(
      low = "#fee08b",
      high = "#d7191c",
      labels = comma,
      name = "Population\nin Zone 3"
    ) +
    annotation_north_arrow(location = "tr", which_north = "true") +
    labs(title = paste0("Population in Zone 3 by County: ", threshold, "-Minute Threshold")) +
    coord_sf(crs = st_crs(output_epsg)) +
    clean_map_theme
  
  save_plot(
    p_zone3_pop,
    file.path(threshold_label, paste0(threshold_label, "_zone3_population_by_county.png"))
  )
  
  # -------------------------------
  # 15.8 Level2_3 point map over counties
  # -------------------------------
  
  p_points <- ggplot() +
    geom_sf(data = counties_wgs, fill = "white", color = "grey30", linewidth = 0.3) +
    geom_sf_text(
      data = county_labels_wgs,
      aes(label = .data[[counties_name_col]]),
      size = 2.0,
      color = "black",
      check_overlap = TRUE
    ) +
    geom_sf(
      data = level2_3_wgs,
      aes(color = .data[[keph_col]]),
      size = 1,
      alpha = 0.8
    ) +
    annotation_north_arrow(location = "tr", which_north = "true") +
    labs(
      title = "Level 2 and 3 Facility Points by KEPH Level",
      color = "KEPH level"
    ) +
    coord_sf(crs = st_crs(output_epsg)) +
    clean_map_theme
  
  save_plot(
    p_points,
    file.path(threshold_label, paste0(threshold_label, "_level2_3_points_by_keph_level.png"))
  )
  
  # -------------------------------
  # 15.9 Dissolved Zone 3 with Level2_3 points
  # -------------------------------
  
  p_points_zone3 <- ggplot() +
    geom_sf(data = counties_wgs, fill = "white", color = "grey30", linewidth = 0.3) +
    geom_sf(data = zone3_dissolved_wgs, fill = "#d7191c", alpha = 0.45, color = NA) +
    geom_sf_text(
      data = county_labels_wgs,
      aes(label = .data[[counties_name_col]]),
      size = 2.0,
      color = "black",
      check_overlap = TRUE
    ) +
    geom_sf(
      data = level2_3_wgs,
      aes(color = .data[[keph_col]]),
      size = 1,
      alpha = 0.8
    ) +
    annotation_north_arrow(location = "tr", which_north = "true") +
    labs(
      title = paste0("Dissolved Zone 3 and Level 2/3 Facilities: ", threshold, "-Minute Threshold"),
      color = "KEPH level"
    ) +
    coord_sf(crs = st_crs(output_epsg)) +
    clean_map_theme
  
  save_plot(
    p_points_zone3,
    file.path(threshold_label, paste0(threshold_label, "_dissolved_zone3_level2_3_points.png"))
  )
  
  # -------------------------------
  # 15.10 Save dissolved Zone 3 polygon
  # -------------------------------
  
  st_write(
    zone3_dissolved_wgs,
    file.path(threshold_dir, paste0(threshold_label, "_dissolved_zone3.gpkg")),
    delete_dsn = TRUE,
    quiet = TRUE
  )
  
  # -------------------------------
  # 15.11 Line chart: county covered population number
  # -------------------------------
  
  coverage_line_df <- county_zone_df %>%
    filter(zone %in% c(1, 2)) %>%
    group_by(county, threshold_min) %>%
    summarise(
      covered_population = sum(population, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(covered_population) %>%
    mutate(county = factor(county, levels = county))
  
  p_line <- ggplot(
    coverage_line_df,
    aes(x = county, y = covered_population, group = 1)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    scale_y_continuous(labels = comma) +
    labs(
      title = paste0("County Covered Population: ", threshold, "-Minute Service Threshold"),
      x = "County",
      y = "Covered population"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  save_plot(
    p_line,
    file.path(threshold_label, paste0(threshold_label, "_line_county_covered_population.png")),
    width = 14,
    height = 7
  )
}

# ===============================
# 16. EXPORT COMBINED RESULTS
# ===============================

all_national_df <- bind_rows(all_national_results)
all_county_df <- bind_rows(all_county_results)

write_csv(all_national_df, combined_national_csv)
write_csv(all_county_df, combined_county_csv)

# ===============================
# 17. COMBINED LINE CHART: ALL THRESHOLDS
# ===============================

combined_line_df <- all_county_df %>%
  filter(zone %in% c(1, 2)) %>%
  group_by(county, threshold_min) %>%
  summarise(
    covered_population = sum(population, na.rm = TRUE),
    .groups = "drop"
  )

county_order <- combined_line_df %>%
  filter(threshold_min == 120) %>%
  arrange(covered_population) %>%
  pull(county)

combined_line_df$county <- factor(combined_line_df$county, levels = county_order)

threshold_colors <- c(
  "30"  = "#1f78b4",
  "60"  = "#33a02c",
  "120" = "#e31a1c"
)

p_combined_line <- ggplot(
  combined_line_df,
  aes(
    x = county,
    y = covered_population,
    color = factor(threshold_min),
    group = threshold_min
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(
    values = threshold_colors,
    name = "Service threshold (mins)",
    labels = c("30 mins", "60 mins", "120 mins")
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "County Covered Population by Service Threshold",
    x = "County",
    y = "Covered population"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(output_dir, "combined_line_county_coverage_all_thresholds.png"),
  plot = p_combined_line,
  width = 16,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ===============================
# 18. TRAVEL TIME SUMMARY
# Mean removed
# ===============================

walk_vals_national  <- as.vector(terra::values(walk_time))
motor_vals_national <- as.vector(terra::values(motor_time))

walk_vals_national  <- walk_vals_national[!is.na(walk_vals_national)]
motor_vals_national <- motor_vals_national[!is.na(motor_vals_national)]

national_tt_summary <- data.frame(
  metric = c(
    "median_walking_travel_time_min",
    "median_motorised_travel_time_min"
  ),
  value = c(
    median(walk_vals_national),
    median(motor_vals_national)
  )
)

write_csv(
  national_tt_summary,
  file.path(output_dir, "national_travel_time_summary.csv")
)

# ===============================
# 19. FINAL MESSAGE
# ===============================

cat("\nDONE.\n")
cat("Outputs saved in:\n", output_dir, "\n")
cat("\nMain updates included:\n")
cat(" - degree_of_urbanisation.tif imported and used for Zone 3 rural/urban population summaries\n")
cat(" - Zone 3 rural/urban population and percentages added to county and national CSVs\n")
cat(" - level2_3 points used instead of generic points layer name\n")
cat(" - level2_3 points plotted with counties and colored by keph_level\n")
cat(" - county names added to maps at size 2.0\n")
cat(" - Zone 3 areas dissolved and mapped with level2_3 points\n")
cat(" - level 4 facilities imported\n")
cat(" - median motorised travel time to nearest level 4 facility for Zone 3 population added to CSVs\n")