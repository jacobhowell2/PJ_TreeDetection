library(terra)
library(progress)
library(curl)

#-------------------------------------------------------------------------------
# step 1. download the index shapefile
#-------------------------------------------------------------------------------
message(date(), " Step 1. Downloading the CHM-NAIP tile index...")
in_file <- "http://rangeland.ntsg.umt.edu/data/naip-chm/index.shp.zip"
out_dir <- "S:/ursa/campbell/grad_students/jake"
out_file <- file.path(out_dir, basename(in_file))
if (!file.exists(out_file)){
  download.file(in_file, out_file)
}

#-------------------------------------------------------------------------------
# step 2. extract the shp from zip
#-------------------------------------------------------------------------------
message(date(), " Step 2. Extracting the SHP from the ZIP...")
zip_file <- out_file
shp_file <- file.path(out_dir, "index_shapefile", "index.shp")
if (!file.exists(shp_file)){
  unzip(zip_file, exdir = dirname(zip_file))
}

#-------------------------------------------------------------------------------
# step 3. read in the index file and pj site polygons
#-------------------------------------------------------------------------------
message(date(), " Step 3. Reading in the spatial data...")
chm_polys <- vect(shp_file)
pj_sites <- vect("S:/ursa/campbell/pj_cnn_naip/data/field/site_polygons.shp")

#-------------------------------------------------------------------------------
# step 4. loop through sites, download, mosaic, and clip
#-------------------------------------------------------------------------------
message(date(), " Step 3. Downloading and processing site-level data...")
chm_dir <- file.path(out_dir, "chms")
dir.create(chm_dir, showWarnings = F)
site_ids <- pj_sites$site_id
for (site_id in site_ids){
  message(date(), "   ", site_id)
  out_file_final <- file.path(chm_dir, paste0("CHM_", site_id, ".tif"))
  if (file.exists(out_file_final)) next
  pj_site <- pj_sites[pj_sites$site_id == site_id,]
  tiles <- chm_polys[relate(chm_polys, pj_site, "intersects"),]
  plot(tiles, border = "red", main = site_id)
  plot(pj_site, border = "blue", lwd = 3, add = T)
  temp_dir <- file.path(chm_dir, site_id)
  dir.create(temp_dir, showWarnings = F)
  tiles$out_file <- file.path(temp_dir, basename(tiles$chm_url))
  pb <- progress_bar$new(
    total = nrow(tiles),
    format = ":current/:total [:bar] :eta"
  )
  message(date(), "     Downloading data...")
  for (i in 1:nrow(tiles)){
    pb$tick()
    tile <- tiles[i,]
    in_file <- tile$chm_url
    out_file <- tile$out_file
    options(timeout = 600)
    if (!file.exists(out_file)) curl_download(in_file, destfile = out_file)
  }
  message(date(), "     Mosaicking and clipping...")
  pj_site_cent_crds <- pj_site |>
    centroids() |>
    crds()
  ref_utm_zone <- floor((pj_site_cent_crds[,1] + 180) / 6) + 1
  ref_epsg <- 26900 + ref_utm_zone
  ref_crs <- paste0("epsg:", ref_epsg) |> crs()
  pj_site_utm <- project(pj_site, ref_crs)
  mosaic_tiles <- tiles$out_file
  scalings <- tiles$scale_fact
  ref_tile <- rast(mosaic_tiles[1]) / scalings[1]
  tile_list <- list()
  pb <- progress_bar$new(
    total = length(mosaic_tiles),
    format = ":current/:total [:bar] :eta"
  )
  for (i in 1:length(mosaic_tiles)){
    pb$tick()
    mosaic_tile <- rast(mosaic_tiles[i]) / scalings[i]
    mosaic_tile <- project(mosaic_tile, ref_crs, method = "bilinear")
    tile_list[[i]] <- mosaic_tile
  }
  tile_sprc <- sprc(tile_list)
  tile_mos <- merge(tile_sprc)
  tile_mos <- crop(tile_mos, pj_site_utm)
  writeRaster(tile_mos, out_file_final, overwrite = T)
  unlink(temp_dir, recursive = T)
  gc()
}