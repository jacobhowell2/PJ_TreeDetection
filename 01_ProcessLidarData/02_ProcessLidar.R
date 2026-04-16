### Create Lidar-Based CHM's from PJ Sites using GEDI data

### Purpose: Use GEDI data to create lidar-based CHM's across field sites within
## Pinyon-Juniper Woodlands

### Author: Jake Howell : u1537023@umail.utah.edu

### Last Updated: 03/16/26

### Load Libraries--------------------------------------------------------------
library(sf)
library(raster)
library(terra)
library(lidR)
library(mapview)
library(future)
library(parallel)
library(fs)

### Set up parallel processing--------------------------------------------------
detectCores()
plan(multisession, workers = parallel::detectCores() /5)

### Set data paths------------------------------------------------------------- 
# Set location of height-normalized lidar files
lazDir <- r'(S:\ursa\howell\Research\PJ_forests\01_Sites\02_Lidar_CHMs\Final\reprojected_laz)'

# Set location of site polygons
sitePolygonPath <- r'(S:\ursa2\campbell\pj_als_gedi\data\field_data\site_polygons.shp)'

# Set location of field plot polygons
plotPolygonPath <- r'(S:\ursa\howell\Research\PJ_forests\02_Plots\00_PlotBounds\fieldPlotExtents.shp)'

# Read site polygons
sitePolygons <- st_read(sitePolygonPath)
plot(sitePolygons$geometry)

# Read plot polygons
plotPolygons <- st_read(plotPolygonPath)
plot(plotPolygons$geometry)

### Set export location--------------------------------------------------------
# Output location for CHM's
outputPath <- r'(S:\ursa\howell\Research\PJ_forests\01_Sites\02_Lidar_CHMs\Final\CHMs)'

### Set switches----------------------------------------------------------------
# Target CRS
target_crs <- 32612 # UTM Zone 12N

### Prep data ------------------------------------------------------------------
folders <- dir_ls(lazDir, type = "directory")

# Define f as a test folder
f <- folders[1]

# Loop through each folder and process lidar data
for (f in folders) {
  
  message("Processing folder: ", f)
  
  folder_name <- basename(f)
  print(folder_name)
  
  # Get all lidar files 
  las_files <- list.files(f, pattern = "\\.la[szx]$", full.names = TRUE)

  # Read as a LAS catalog
  ctg <- readLAScatalog(f)
  plot(ctg)
  sitePolygons <- st_transform(sitePolygons, crs(ctg))
  plot(sitePolygons$geometry, add = TRUE)
  plot(plotPolygons$geometry, add = TRUE, col = "red", lwd = 7)
  
  # Set processing options
  opt_progress(ctg)     <- TRUE
  
  # Output CHM filename for this folder
  opt_output_files(ctg) <- file.path(outputPath, paste0(folder_name, "_{*}.tif"))
  
  # CHM function
  rasterize_canopy(ctg, 
                   res = 0.6, 
                   algorithm = p2r(subcircle = 0.20),
                   pkg = "terra")
  
  # Identify the directory where outputs were written
  out_dir <- dirname(outputPath)
  
  # Rename the VRT so it is unique per folder
  vrt_src  <- file.path(out_dir, "rasterize_canopy.vrt")
  vrt_dest <- file.path(out_dir, paste0(folder_name, "_mosaic.vrt"))

  if (file.exists(vrt_src)) {
    file.rename(vrt_src, vrt_dest)
  }
  
  # Remove temporary rasters created during mosaicking
  temp_rasters <- list.files(out_dir, pattern = "chunk_.*\\.tif$", full.names = TRUE)
  if (length(temp_rasters) > 0) {
    file.remove(temp_rasters)
  }
  
}


##########----------------------------------------------------------------------


### Create Lidar-Based CHM's from PJ Sites using GEDI data

### Purpose: Use GEDI data to create lidar-based CHM's across field sites within
## Pinyon-Juniper Woodlands

### Author: Jake Howell : u1537023@umail.utah.edu

### Last Updated: 03/16/26

### Load Libraries --------------------------------------------------------------
library(sf)
library(raster)
library(terra)
library(lidR)
library(mapview)
library(future)
library(parallel)
library(fs)

### Set up parallel processing --------------------------------------------------
plan(multisession, workers = max(1, floor(parallel::detectCores() / 5)))

### Set data paths -------------------------------------------------------------- 
lazDir <- r'(S:\ursa\howell\Research\PJ_forests\01_Sites\02_Lidar_CHMs\Final\reprojected_laz)'
sitePolygonPath <- r'(S:\ursa2\campbell\pj_als_gedi\data\field_data\site_polygons.shp)'
plotPolygonPath <- r'(S:\ursa\howell\Research\PJ_forests\02_Plots\00_PlotBounds\fieldPlotExtents.shp)'

### Read polygons once ----------------------------------------------------------
sitePolygons <- st_read(sitePolygonPath, quiet = TRUE)
plotPolygons <- st_read(plotPolygonPath, quiet = TRUE)

### Set export location ---------------------------------------------------------
outputPath <- r'(S:\ursa\howell\Research\PJ_forests\01_Sites\02_Lidar_CHMs\Final\CHMs)'
dir_create(outputPath)

### Set switches ----------------------------------------------------------------
target_crs <- 32612

### Prep data -------------------------------------------------------------------
folders <- dir_ls(lazDir, type = "directory")

f <- folders[1]
### Loop through folders --------------------------------------------------------
for (f in folders) {
  
  folder_name <- basename(f)
  message("Processing folder: ", folder_name)
  
  # Create a site-specific output folder
  site_out <- file.path(outputPath, folder_name)
  dir_create(site_out)
  
  # Read folder as LAS catalog
  ctg <- readLAScatalog(f)
  
  # # Skip if catalog is empty
  # if (is.empty(ctg)) {
  #   message("  Catalog is empty. Skipping.")
  #   next
  # }
  
  # CRS check
  ctg_crs <- tryCatch(st_crs(ctg)$epsg, error = function(e) NA)
  if (!is.na(ctg_crs) && ctg_crs != target_crs) {
    message("  CRS mismatch (found EPSG:", ctg_crs, "). Skipping.")
    next
  }
  
  # Set catalog options
  opt_progress(ctg) <- TRUE
  opt_chunk_size(ctg) <- 0
  opt_chunk_buffer(ctg) <- 0
  
  # Drop overlap-classified points only
  opt_filter(ctg) <- "-drop_overlap"
  
  # Write chunk rasters into the site folder
  opt_output_files(ctg) <- file.path(site_out, paste0(folder_name, "_chunk_{*}"))
  
  # Optional QC plot
  # plot(ctg)
  # plot(st_geometry(st_transform(sitePolygons, st_crs(ctg))), add = TRUE)
  # plot(st_geometry(st_transform(plotPolygons, st_crs(ctg))), add = TRUE, col = "red", lwd = 2)
  
  # Create CHM tiles
  rasterize_canopy(
    ctg,
    res = 0.6,
    algorithm = p2r(subcircle = 0.20),
    pkg = "terra"
  )
  
  # Get written chunk tiles
  tif_files <- dir_ls(site_out, regexp = paste0("^", folder_name, "_chunk_.*\\.tif$"))
  
  if (length(tif_files) == 0) {
    message("  No CHM tiles written for ", folder_name)
    next
  }
  
  # Build VRT
  vrt_file <- file.path(site_out, paste0(folder_name, "_CHM_tiles.vrt"))
  terra::vrt(tif_files, filename = vrt_file, overwrite = TRUE)
  
  # Mosaic tiles
  ras_list <- lapply(tif_files, terra::rast)
  chm_mosaic <- do.call(terra::mosaic, ras_list)
  
  mosaic_file <- file.path(site_out, paste0(folder_name, "_CHM_mosaic.tif"))
  terra::writeRaster(chm_mosaic, mosaic_file, overwrite = TRUE)
  
  message("  Wrote chunk tiles to: ", site_out)
  message("  Wrote VRT: ", vrt_file)
  message("  Wrote mosaic: ", mosaic_file)
}






# PRocess Lidar Data
library(lidR)
library(terra)
library(fs)

# Input directory: contains one folder per plot
input_dir <- r'(N:\Howell\PJ_Forests\01_Sites\plotLidar)'

# Output directory for CHMs
output_dir <- r'(S:\ursa\howell\Research\PJ_forests\02_Plots\03_LidarDerivedCHMs\Original)'
dir_create(output_dir)

# Find only files ending in _buff_10m.laz
las_files <- list.files(
  input_dir,
  pattern = "_buff_10m\\.laz$",
  full.names = TRUE,
  recursive = TRUE
)

# Check what was found
print(las_files)

# Loop through files
for (f in las_files) {
  
  # Remove _buff_10m so output is just like UOFU_AZCF_P01
  plot_name <- sub("_buff_10m$", "", tools::file_path_sans_ext(basename(f)))
  
  message("Processing: ", plot_name)
  
  # Read lidar
  las <- tryCatch(
    readLAS(f),
    error = function(e) {
      message("Failed to read: ", e$message)
      NULL
    }
  )
  
  # Skip bad or empty files
  if (is.null(las) || is.empty(las)) {
    message("Empty or unreadable file. Skipping.")
    next
  }
  
  # Create CHM
  chm <- tryCatch(
    rasterize_canopy(
      las,
      res = 0.5,
      algorithm = p2r(subcircle = 0.2),
      pkg = "terra"
    ),
    error = function(e) {
      message("CHM failed: ", e$message)
      NULL
    }
  )
  
  if (is.null(chm)) next
  
  # Output file name: UOFU_AZCF_P01.tif
  out_file <- file.path(output_dir, paste0(plot_name, ".tif"))
  
  # Write CHM
  writeRaster(chm, out_file, overwrite = TRUE)
  
  message("Wrote: ", out_file)
}
