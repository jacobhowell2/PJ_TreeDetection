### Modeling Aboveground Biomass of Individual Trees in PJ-woodlands 

### Purpose: Use CHM's derived from NAIP imagery (Allred et al., 2025) to model
### individual tree biomass using tree detection

### Author: Jake Howell : u1537023@umail.utah.edu

### Last Updated: 1/28/26

### Install Packages------------------------------------------------------------
# Look where R is installing packages
.libPaths()

# Required packages for cloud2trees
print("Checking necessary libraries are installed.")
pckgs <- c("pkgbuild",
           "remotes",
           "lasR",
           "lidR",
           "mapview",
           "units",
           "ggplot2",
           "broom",
           "purrr")

newPckgs <- pckgs[!(pckgs %in% installed.packages()[,"Package"])]

if(length(newPckgs) > 0) install.packages(newPckgs)

# Check if build tools are available 
pkgbuild::check_build_tools(debug = TRUE)

## Check for package installations from github 
# TreeLS for extracting DBH
if (!requireNamespace("TreeLS", quietly = TRUE)) {
  # Install TreeLS  if not already installed
  remotes::install_github(repo = "tiagodc/TreeLS", upgrade = F)
}

# LaddarFuelsR for estimating crown base height
if (!requireNamespace("LaddarFuelsR", quietly = TRUE)) {
  # Install LadderFuelsR 
  remotes::install_github(repo = "olgaviedma/LadderFuelsR", upgrade = F)
}

# leafR 
if (!requireNamespace("leafR", quietly = TRUE)) {
  # Install leafR  if not already installed
  remotes::install_github(repo = "DRAAlmeida/leafR", upgrade = F)
}

# cloud2trees
if (!requireNamespace("cloud2trees", quietly = TRUE)) {
  # Install cloud2trees
  remotes::install_github(repo = "georgewoolsey/cloud2trees", upgrade = F)
}


# Update any packages
update.packages(ask = FALSE)

# Load libraries
library(sf)
library(dplyr)
library(sp)
library(raster)
library(terra)
library(lidR)
library(ForestTools)
library(cloud2trees)
library(mapview)
library(units)
library(ggplot2)
library(broom)
library(purrr)

### Set up data for Cloud2Trees-------------------------------------------------
# Location of cloud2trees data needed to extract features to tree list:  If a 
# custom directory is not set up, data will download to users tempdir() folder
locExtData <- r'(N:\Howell\PJ_Forests\99_ExtData\cloud2trees)'

### Download TreeMap, Forest Type Group, and LANDFIRE data ---------------------
# Download ALL external data to estimate tree DBH, FIA forest type, and LANDFIRE 
# fuel loading. This data must be downloaded first before using
cloud2trees::get_data(savedir = locExtData, force = F)

## OR download individual data sets to external data directory
# Download TreeMap 2016 to derive site-specific allometric equations using data
# from FIA. This is a ~3 GB download that must be run the first time using the
# cloud2trees package
cloud2trees::get_treemap(savedir = locExtData, force = F)

# Download the Forest Type Groups of the Continental US data. This group layer 
# includes a 30-meter resolution of forest type data. This is a ~2GB download
# that must be downloaded the first time cloud2trees is run
cloud2trees::get_foresttype(savedir = locExtData, force = F)

# Download LANDFIRE data to estimate tree crown biomass using the forest canopy
# bulk density to determine the fuel loading of a particular area
cloud2trees::get_landfire(savedir = locExtData, force = F)

### Read Data-------------------------------------------------------------------
# Set parent location
loc <- r'(N:\Howell\PJ_Forests)'

## SITE RASTERS ----------
# Set the name of site folder containing all of the kNN CHMs
foldSite <- "01_Sites"

# Folder of site level CHMs
foldRast <- "NAIP_derived_CHMs"

# Set file path
pathSites <- file.path(loc, foldSite, foldRast)

# List the sites within this folder
siteNames <- list.files(pathSites, pattern = "\\.tif$", full.names = TRUE)

# Now read the rasters as a list
siteRast <- lapply(siteNames, rast)

## PLOT DATA ------------
# Set the name of plot folder
foldPlot <- "02_Plots"

# Set the name of the plot level .csv
namePlotData <- "plot_level_data_uofu.csv"

# Set full path to plot level data
fullPlotPath <- file.path(loc, foldPlot, namePlotData)

# Read plot field data
plotData <- read.csv(fullPlotPath)

## TREE DATA ------------
foldTree <- "03_Trees"

# Set the name of the tree data
nameTreeData <- "tree_level_data.csv"

# Set full path to tree level data
fullTreePath <- file.path(loc, foldTree, nameTreeData)

# Read tree data
treeData <- read.csv(fullTreePath)

## BIOMASS EQUATIONS
pathBiomass <- r'(N:\Howell\PJ_Forests\04_Biomass\all_tree_species_in_pj_fia_plots_new.csv)'

# Read biomass equations
biomassEquations <- read.csv(pathBiomass)

### Output Locations------------------------------------------------------------
# Output directory of clipped plots
clippedPlotsOut <- r'(N:\Howell\PJ_Forests\02_Plots\ClippedNAIPCHMs)'

### Input Switches--------------------------------------------------------------
# Biomass Species Name
specBio <- "plants.code"

# Tree Data Species Name
specTree <- "tree_sp"

### 1) Apply Biomass Equations to Tree Data-------------------------------------
# Edit tree data species column
treeData <- treeData |>
  rename(Species = !!specTree)

# Edit biomass equation species column
biomassEquations <- biomassEquations |>
  rename(Species = !!specBio)

# Join the tree data to the biomass equations by Species column
treeData <- left_join(treeData, biomassEquations, by = "Species")

# Now apply the biomass equation based off coefficients to each tree
treeData <- treeData |>
  mutate(tree_biomass_kg = exp(beta.0 + beta.1 * log(tree_diam)))

# Remove outliers
treeData_filtered <- treeData |>
  filter(tree_diam < 200)

# Plot biomass for each species type by DBH and expected biomass
ggplot2::ggplot(treeData_filtered, aes(x = tree_diam, y = tree_biomass_kg, color = Species)) +
  geom_point() +
  labs(
    title = "Tree Biomass by Diameter at Breast Height (DBH)",
    x = "Diameter at Breast Height (cm)",
    y = "Biomass (kg)") +
  theme_minimal()

### 2) Create plot extents from the plot data spreadsheet-----------------------
# # Create extent of site CHM
# chm_extent <- st_as_sf(st_as_sfc(st_bbox(siteCHM)))

# Create function to make one polygon from a row
createFieldPlotExtents <- function(row) {
  coords <- matrix(c(
    row$plot_ne_x, row$plot_ne_y,
    row$plot_se_x, row$plot_se_y,
    row$plot_sw_x, row$plot_sw_y,
    row$plot_nw_x, row$plot_nw_y,
    row$plot_ne_x, row$plot_ne_y  
  ), ncol = 2, byrow = TRUE)
  st_polygon(list(coords))
}

# Create field plot extent geometry in plot data 
fieldPlots_extents <- plotData |>
  rowwise() |>
  mutate(geometry = st_sfc(createFieldPlotExtents(cur_data())))

# Set the current CRS
fieldPlots_extents <- st_as_sf(fieldPlots_extents, crs = 4326)
st_crs(fieldPlots_extents)

# Reproject
fieldPlots_extents <- st_transform(fieldPlots_extents, crs = 32612)
st_crs(fieldPlots_extents)

# Check plot geometries
# plot(fieldPlots_extents$geometry[[149]])
# crs(fieldPlots_extents)

# fieldPlots_extents <- vect(fieldPlots_extents)
# class(fieldPlots_extents)

# Take a subset for testing
siteRast_subset <- siteRast[c(1:2)]

### Loop through each site CHM and clip to plots------------------------------
allPlotsCHMs <- list()

for (i in seq_along(siteRast)) {
  
  chm <- siteRast[[i]]
  
  # Create ONE polygon from raster extent
  chm_poly <- st_as_sf(as.polygons(ext(chm), crs = crs(chm)))
  
  # sf spatial intersection
  fieldPlots_extents <- st_transform(fieldPlots_extents, crs(chm_poly))
  idx <- st_intersects(fieldPlots_extents, chm_poly, sparse = FALSE)
  plot_rows <- which(idx)                
  plots_in_tile <- fieldPlots_extents[plot_rows, ]
  
  if (nrow(plots_in_tile) == 0) next
  
  for (j in 1:nrow(plots_in_tile)) {
    
    plot_sf <- plots_in_tile[j, ]
    
    # Add a 5 m BUFFER (assumes projected CRS in meters)
    plot_sf_buf <- st_buffer(plot_sf, dist = 2.5)
    
    # Convert sf → SpatVector
    plot_vect <- vect(plot_sf_buf)
    
    # Store plot name 
    plot_name <- plot_sf$plot_id
    
    # Clip the CHM
    clipped <- crop(chm, plot_vect) |> mask(plot_vect)
    
    # Add to list with plot name
    allPlotsCHMs[[plot_name]] <- clipped
    
    # Build output filename
    out_name <- file.path(clippedPlotsOut, paste0("plot_", plot_sf$plot_id, ".tif"))
    
    writeRaster(clipped, out_name, overwrite = TRUE)
  }
}
# Save the names
chm_names <- names(allPlotsCHMs)

### 4) Sum the per plot biomass ------------------------------------------------
# Return just my plot names
plot_ids <- unique(fieldPlots_extents$plot_id)

# Subset the tree data to just the unique plot names 
treeDataInPlots <- treeData[treeData$plot_id %in% plot_ids, ]

# Estimate Biomass Summary for each plot
plotLevelBiomass <- treeDataInPlots %>%
  group_by(plot_id) %>%
  summarise(
    total_biomass = sum(tree_biomass_kg, na.rm = TRUE),
    mean_biomass = mean(tree_biomass_kg, na.rm = TRUE),
    n_trees_GPS = n(),
    .groups = 'drop'
  )

### 5) Run Tree Detection on each plot and compute volume-----------------------
# # Reread the plot level rasters in
# plotNames <- list.files(clippedPlotsOut, pattern = "\\.tif$", full.names = TRUE)

# # Now read the rasters as a list
# allPlotsCHMs <- lapply(plotNames, rast)

# # # Set moving window size
# h <- function(x) {
#   y <- 2.6 * (-(exp(-0.08*(x-2)) - 1)) + 3
#   y[x < 2] <- 3
#   y[x > 20] <- 5
#   return(y)
# }

# Window size 2
h <- function(x) {
  cw <- 1.0827 * x - 0.3993
  ws <- pmin(10, pmax(3, cw/2))
  return(ws)
}

heights <- seq(-5,30,0.5)
ws <- h(heights)
plot(heights, ws, type = "l",  ylim = c(0,12))

# Initialize lists for both points and crowns
treePoints_list <- list()
treeCrowns_list <- list()


# Loop through all plot CHMs and segment trees and crowns
for (plot_id in chm_names) {
  
  message("Processing plot: ", plot_id)
  
  # Subset to a single CHM
  chm_i <- allPlotsCHMs[[plot_id]]
  
  # Run tree detection safely
  treeCrowns_i <- tryCatch(
    raster2trees(
      chm_rast   = chm_i,
      outfolder  = tempdir(),
      ws         = h,
      min_height = 0.25
    ),
    error = function(e) {
      message("Tree detection failed for plot ", plot_id, 
              ": ", e$message)
      return(NULL)
    }
  )
  
  # Skip if detection failed or returned nothing
  if (is.null(treeCrowns_i) || nrow(treeCrowns_i) == 0) {
    message("No trees detected for plot ", plot_id, "skipping.")
    next
  }
  
  # Add plot ID
  treeCrowns_i$plotID <- plot_id
  
  # convert detected trees to points
  treePoints_i <- st_as_sf(
    as.data.frame(treeCrowns_i),
    coords = c("tree_x", "tree_y"),
    crs = 32612   # EPSG:32612
  )
  
  # Reproject tree crowns to match points
  treeCrowns_i <- st_transform(treeCrowns_i, crs = crs(treePoints_i))
  
  # Store outputs
  treePoints_list[[plot_id]] <- treePoints_i
  treeCrowns_list[[plot_id]] <- treeCrowns_i
  
}

# Dataframe for detected tree points
allTreePoints <- bind_rows(treePoints_list)

# Dataframe for tree crowns
allTreeCrowns <- do.call(rbind, treeCrowns_list)

### EXPLORE RESULTS-------------------------------------------------------------
# Select the first plot name
first_plot <- names(treeCrowns_list)[2]
# first_plot <- "UOFU_UTMO_P07"

# Extract CHM and crown polygons
chm_i     <- allPlotsCHMs[[first_plot]]
plot(chm_i)

crowns_i  <- treeCrowns_list[[first_plot]]
points_i <- treePoints_list[[first_plot]]

# Ensure CRS matches
crowns_i <- st_transform(crowns_i, st_crs(chm_i))

# Plot CHM
plot(chm_i, main = paste("NAIP-derived CHM using VWS1:", first_plot))

# Overlay crown segments
plot(st_geometry(crowns_i),
     add = TRUE,
     border = "yellow",
     lwd = 2)

# Overlay tree points
plot(st_geometry(points_i),
     add = TRUE,
     col = "red",
     lwd = 2)

### Create biomass allometries at the individual tree level---------------------
# Create tree and plot df
df_tree <- allTreeCrowns
df_plot <- plotLevelBiomass

# Select the variables I want
df_tree <- df_tree[,c("plotID", "tree_height_m", "crown_area_m2")]
df_plot <- df_plot[,c("plot_id", "total_biomass")]

# Subset data to these variables
colnames(df_tree) <- c("plot_id", "tree_height", "tree_area")
colnames(df_plot) <- c("plot_id", "plot_biomass")

# Merge the data
df <- merge(df_plot, df_tree)

# Split into list by plot
plot_ids <- unique(df$plot_id)
tree_list <- split(df, df$plot_id)

# Extract the observed plot biomass
plot_biomass <- tapply(df$plot_biomass, df$plot_id, FUN = function(x) x[1])

# Define a predictive allometric model
predict_plot_biomass <- function(params) {
  a <- params[1]
  b <- params[2]
  c <- params[3]
  
  pred <- numeric(length(tree_list))
  names(pred) <- names(tree_list)
  
  for(i in seq_along(tree_list)) {
    trees <- tree_list[[i]]
    # apply the allometry to each tree, then sum
    tree_vals <- a * (trees$tree_height^b) * (trees$tree_area^c)
    pred[i] <- sum(tree_vals, na.rm = TRUE)
  }
  
  pred
}

# Define the objective function (sum of squared errors)
objective <- function(params) {
  pred <- predict_plot_biomass(params)
  # ensure ordering matches
  resid <- plot_biomass - pred[names(plot_biomass)]
  sum(resid^2)
}

# Fit the model through numerical optimization
start_vals <- c(a = 0.1, b = 1, c = 1)

fit <- optim(
  par = start_vals,
  fn = objective,
  method = "L-BFGS-B",
  lower = c(1e-6, -5, -5),
  upper = c(100,   5,  5)
)

# Return the predicted biomass
predicted <- predict_plot_biomass(fit$par)

# scatterplot
plot(plot_biomass, predicted,
     xlab = "Observed plot biomass",
     ylab = "Predicted plot biomass")
abline(0, 1, col = "red")

# Observed and predicted vectors
obs <- plot_biomass
pred <- predicted[names(plot_biomass)]   # ensure same order

# 1. Root Mean Squared Error
RMSE <- sqrt(mean((obs - pred)^2))
print(RMSE)

# 2. Mean Absolute Error
MAE <- mean(abs(obs - pred))
print(MAE)

# 3. Bias (signed error)
Bias <- mean(pred - obs)

# 4. Mean Absolute Percent Error
MAPE <- mean(abs((obs - pred) / obs)) * 100

# 5. R-squared (coefficient of determination)
SSE <- sum((obs - pred)^2)
SST <- sum((obs - mean(obs))^2)
R2 <- 1 - SSE/SST
print(R2)

# 6. Relative RMSE (RMSE / mean observed)
rRMSE <- RMSE / mean(obs) * 100

# Save metrics
NAIP_chm_vws1_metrics <- data.frame(
  RMSE = RMSE,
  MAE = MAE,
  Bias = Bias,
  MAPE_percent = MAPE,
  R2 = R2,
  rRMSE_percent = rRMSE
)

plot(obs, pred,
     xlab = "Observed plot biomass",
     ylab = "Predicted plot biomass",
     main = "Observed vs Predicted Plot Biomass")
abline(0, 1, col = "red", lwd = 2)



### Create biomass allometries at the individual tree level (log-log)
# Create tree and plot df
df_tree <- allTreeCrowns
df_plot <- plotLevelBiomass

# Select variables
df_tree <- df_tree[, c("plotID", "tree_height_m", "crown_area_m2")]
df_plot <- df_plot[, c("plot_id", "total_biomass")]

# Rename columns
colnames(df_tree) <- c("plot_id", "tree_height", "tree_area")
colnames(df_plot) <- c("plot_id", "plot_biomass")

# Merge plot and tree data
df <- merge(df_plot, df_tree)

# Log transform
df$log_observed_bio <- log(df$plot_biomass)
df_plot$log_observed_bio <- log(df_plot$plot_biomass) 

# Remove invalid values 
df <- df[
  df$tree_height > 0 &
    df$tree_area  > 0 &
    df$plot_biomass > 0,
]

# Split into list by plot
tree_list <- split(df, df$plot_id)

# Extract the observed plot biomass
# This is returning a numeric values of each of the plots observed biomass
plot_biomass <- tapply(
  df$plot_biomass,
  df$plot_id,
  FUN = function(x) x[1]
)

print(plot_biomass)

# Predict log - log biomass 
predict_plot_biomass_loglog <- function(params) {
  
  # Intercept
  a     <- params[1]  
  
  # Coefficient for tree height
  b     <- params[2]
  
  # Coefficient for crown area
  c     <- params[3]
  
  pred <- numeric(length(tree_list))
  names(pred) <- names(tree_list)
  
  for (i in seq_along(tree_list)) {
    trees <- tree_list[[i]]
    
    # Set the equation for AGB for individual trees
    log_tree_agb <- a +  b * log(trees$tree_height) + c * log(trees$tree_area)
    
    # Back-transform
    tree_agb <- exp(log_tree_agb)
    
    # Get the predicted AGB for each tree
    pred[i] <- sum(tree_agb, na.rm = TRUE)
  }
  
  pred
}


# Set the ojective function (plot-level SSE)
objective_loglog <- function(params) {
  pred <- predict_plot_biomass_loglog(params)
  resid <- log(plot_biomass[names(pred)]) - log(pred)
  sum(resid^2)
}

# Initial values (log-space!)
start_vals <- c(
  alpha = log(0.1),  # log(a)
  b     = 1,
  c     = 1
)

# Fit model
fit_loglog <- optim(
  par     = start_vals,
  fn      = objective_loglog,
  method  = "L-BFGS-B",
  lower   = c(log(1e-6), -5, -5),
  upper   = c(log(100),   5,  5)
)

# Predictions
predicted_plot_agb <- predict_plot_biomass_loglog(fit_loglog$par)

# Back-transform coefficients for reporting
a_hat <- exp(fit_loglog$par[1])
b_hat <- fit_loglog$par[2]
c_hat <- fit_loglog$par[3]

# Observed vs predicted (plot level)
diag_df <- data.frame(
  plot_id   = names(plot_biomass),
  observed  = plot_biomass,
  predicted = predicted_plot_agb[names(plot_biomass)]
)

smear <- mean(exp(diag_df$log_resid))
pred_corrected <- predicted_plot_agb * smear

diag_df$pred_corrected <- pred_corrected
diag_df$resid_corrected <- diag_df$observed - diag_df$pred_corrected

# Residuals
diag_df$residuals <- diag_df$observed - diag_df$predicted
diag_df$log_obs   <- log(diag_df$observed)
diag_df$log_pred  <- log(diag_df$predicted)
diag_df$log_resid <- diag_df$log_obs - diag_df$log_pred


# RMSE
rmse <- sqrt(mean(diag_df$residuals^2))

# MAE
mae <- mean(abs(diag_df$residuals))

# R-squared (plot level)
ss_res <- sum(diag_df$residuals^2)
ss_tot <- sum((diag_df$observed - mean(diag_df$observed))^2)
r2 <- 1 - ss_res / ss_tot

rmse
mae
r2


print(df_tree)

plot_predictors <- df_tree %>%
  group_by(plot_id) %>%
  summarise(
    sum_height = sum(tree_height, na.rm = TRUE),
    sum_area   = sum(tree_area, na.rm = TRUE)
  )

# Merge with observed biomass
df_plot_level <- df_plot %>%
  inner_join(plot_predictors, by = "plot_id")








### Predict Plot Level Biomass using Log-Log function---------------------------
# Function to predict tree level biomass
predict_plot_biomass <- function(params) {
  
  # Intercept
  a <- params[1]
  
  # Exponent of Height
  b <- params[2]
  
  # Exponent of crown area
  c <- params[3]
  
  # Create vector to store predictions
  # Pred = the total predicted biomass for each plot
  pred <- numeric(length(tree_list))
  
  # This will store the names of each plot with it's predicted biomass
  names(pred) <- names(tree_list)
  
  # Loops through each individual tree 
  for(i in seq_along(tree_list)) {
    trees <- tree_list[[i]]
    tree_vals <- a * (trees$tree_height^b) * (trees$tree_area^c)
    pred[i] <- sum(tree_vals, na.rm = TRUE)
  }
  
  pred
}

objective <- function(params) {
  pred <- predict_plot_biomass(params)
  resid <- plot_biomass - pred[names(plot_biomass)]
  sum(resid^2)
}

objective_loglog <- function(params) {
  
  pred <- predict_plot_biomass(params)
  
  # Enforce positivity for log transform
  if (any(pred <= 0)) return(Inf)
  
  log_obs  <- log(plot_biomass)
  log_pred <- log(pred[names(plot_biomass)])
  
  resid <- log_obs - log_pred
  
  sum(resid^2)
}

start_vals <- c(a = 0.05, b = 1, c = 1)

fit <- optim(
  par = start_vals,
  fn = objective_loglog,
  method = "L-BFGS-B",
  lower = c(1e-8, -3, -3),
  upper = c(100,   5,  5)
)

predicted <- predict_plot_biomass(fit$par)

plot(plot_biomass, predicted,
     xlab = "Observed plot biomass",
     ylab = "Predicted plot biomass")
abline(0, 1, col = "red")

# Predictions
pred_linear <- predict_plot_biomass(fit$par)

# Log-space values
log_obs  <- log(plot_biomass)
log_pred <- log(pred_linear[names(plot_biomass)])

# Residuals
resid_log    <- log_obs - log_pred
resid_linear <- plot_biomass - pred_linear[names(plot_biomass)]

# Get R2 for log space
R2_log <- 1 - sum(resid_log^2) / sum((log_obs - mean(log_obs))^2)
print(R2_log)

# R2 in linear space
R2_lin <- 1 - sum(resid_linear^2) / sum((plot_biomass - mean(plot_biomass))^2)

# RMSE
RMSE <- sqrt(mean(resid_linear^2))

print(c(R2_log  = R2_log,
  R2_lin  = R2_lin,
  RMSE    = RMSE))



























































































