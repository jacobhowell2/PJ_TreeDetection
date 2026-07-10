### PJ Woodlands R-Shiny Application
### app.R

### Load Libraries -------------------------------------------------------------
library(terra)
library(sf)
library(shiny)
library(leaflet)
library(base64enc)

### Setup ----------------------------------------------------------------------
# Images
img_file <- normalizePath(file.path("www", "pj.jpg"), mustWork = TRUE)
pj_img_src <- base64enc::dataURI(file = img_file, mime = "image/jpeg")

pinyon_file <- normalizePath(file.path("www", "pinyon.jpg"), mustWork = TRUE)
pinyon_img_src <- base64enc::dataURI(file = pinyon_file, mime = "image/jpeg")

# Precomputed folder
pre_dir <- "data"

required_files <- c(
  file.path(pre_dir, "dem_overview_leaf.tif"),
  file.path(pre_dir, "ppt_overview_leaf.tif"),
  file.path(pre_dir, "finalStrat_leaf.tif"),
  file.path(pre_dir, "pj_boundary_sf_ll.rds"),
  file.path(pre_dir, "sites_sf_ll.rds"),
  file.path(pre_dir, "plots_sf_ll.rds"),
  file.path(pre_dir, "plot_pts_sf_ll.rds"),
  file.path(pre_dir, "trees_ll.rds"),
  file.path(pre_dir, "app_meta.rds")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "The following precomputed files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

### Load Precomputed Data ------------------------------------------------------
dem_overview_leaf <- terra::rast(file.path(pre_dir, "dem_overview_leaf.tif"))
ppt_overview_leaf <- terra::rast(file.path(pre_dir, "ppt_overview_leaf.tif"))
finalStrat_leaf   <- terra::rast(file.path(pre_dir, "finalStrat_leaf.tif"))

pj_boundary_sf_ll <- readRDS(file.path(pre_dir, "pj_boundary_sf_ll.rds"))
sites_sf_ll       <- readRDS(file.path(pre_dir, "sites_sf_ll.rds"))
plots_sf_ll       <- readRDS(file.path(pre_dir, "plots_sf_ll.rds"))
plot_pts_sf_ll    <- readRDS(file.path(pre_dir, "plot_pts_sf_ll.rds"))
trees_ll          <- readRDS(file.path(pre_dir, "trees_ll.rds"))
app_meta          <- readRDS(file.path(pre_dir, "app_meta.rds"))

# Rebuild palettes from metadata
dem_pal <- colorNumeric(
  palette = terrain.colors(20),
  domain = app_meta$dem_range,
  na.color = "transparent"
)

ppt_pal <- colorNumeric(
  palette = colorRampPalette(c(
    "#d9f0ff",
    "#73bfe2",
    "#2b8cbe",
    "#045a8d",
    "#023858"
  ))(20),
  domain = app_meta$ppt_range,
  na.color = "transparent"
)

tree_pal <- colorFactor(
  palette = app_meta$tree_species_cols,
  domain = names(app_meta$tree_species_cols),
  na.color = "gray70"
)

tree_dot_radius <- app_meta$tree_dot_radius
strata_vals     <- app_meta$strata_vals
strata_cols     <- app_meta$strata_cols

strata_pal <- colorFactor(
  palette = strata_cols,
  domain = strata_vals,
  na.color = "transparent"
)

### Initialize User Interface --------------------------------------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML(
      "body {
        background: linear-gradient(135deg, #2f3e2f 0%, #5c4b3b 45%, #a67c52 100%);
        background-attachment: fixed;
        color: white;
      }

      .container-fluid {
        padding: 20px;
      }

      .header-row-custom {
        display: flex;
        align-items: center;
        justify-content: flex-start;
        gap: 20px;
        background: rgba(0, 0, 0, 0.35);
        padding: 20px;
        border-radius: 12px;
        margin-bottom: 15px;
        box-shadow: 0 4px 12px rgba(0,0,0,0.3);
      }

      .header-title-custom {
        margin: 0;
        color: white;
        font-weight: 700;
        font-size: 32px;
      }

      .sidebar-image-custom {
        margin-bottom: 18px;
        text-align: center;
      }

      .sidebar-image-custom img {
        width: 100%;
        height: auto;
        max-height: 220px;
        object-fit: cover;
        border-radius: 12px;
        box-shadow: 0 4px 10px rgba(0,0,0,0.35);
      }

      .sidebar-panel-custom {
        background: rgba(255, 255, 255, 0.12);
        backdrop-filter: blur(4px);
        -webkit-backdrop-filter: blur(4px);
        padding: 18px;
        border-radius: 12px;
        line-height: 1.6;
        box-shadow: 0 4px 12px rgba(0,0,0,0.25);
        min-height: 600px;
      }

      .map-panel-custom {
        background: rgba(0, 0, 0, 0.18);
        padding: 10px;
        border-radius: 12px;
        box-shadow: 0 6px 18px rgba(0,0,0,0.35);
      }

      .strata-info-box {
        margin-top: 18px;
        padding: 14px;
        background: rgba(255, 255, 255, 0.10);
        border-radius: 10px;
        box-shadow: inset 0 0 0 1px rgba(255,255,255,0.08);
      }

      .strata-table {
        width: 100%;
        border-collapse: separate;
        border-spacing: 4px;
        margin-top: 10px;
      }

      .strata-table th {
        text-align: center;
        font-size: 11px;
        color: rgba(255,255,255,0.75);
        font-weight: 600;
        padding: 2px 4px;
      }

      .strata-table th:first-child {
        text-align: right;
      }

      .strata-row-label {
        font-size: 11px;
        color: rgba(255,255,255,0.75);
        font-weight: 600;
        text-align: right;
        padding-right: 6px;
        white-space: nowrap;
      }

      .strata-cell {
        border-radius: 6px;
        text-align: center;
        padding: 8px 4px;
        font-weight: 700;
        font-size: 15px;
        color: #222;
      }

      #map {
        border-radius: 12px;
        overflow: hidden;
      }

      p {
        color: white;
        font-size: 15px;
      }

      .control-label {
        color: white;
        font-weight: 600;
        margin-top: 10px;
      }

      .selectize-input,
      .selectize-dropdown,
      .form-control {
        color: black;
        border-radius: 8px;
      }

      .nav-tabs {
        border-bottom: none;
        margin-bottom: 15px;
      }

      .nav-tabs > li > a {
        color: white;
        background: rgba(255,255,255,0.08);
        border: none;
        border-radius: 10px 10px 0 0;
        margin-right: 4px;
      }

      .nav-tabs > li > a:hover {
        background: rgba(255,255,255,0.16);
        color: white;
      }

      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:hover,
      .nav-tabs > li.active > a:focus {
        color: white;
        background: rgba(255,255,255,0.20);
        border: none;
      }

      .tab-content {
        padding-top: 8px;
      }

      .main-layout-custom {
        display: flex;
        gap: 16px;
        align-items: flex-start;
      }

      .sidebar-col-custom {
        width: 31%;
      }

      .map-col-custom {
        width: 69%;
      }

      .nav-tabs {
        display: flex;
        flex-wrap: nowrap;
        gap: 4px;
        margin-bottom: 15px;
      }

      .nav-tabs > li {
        float: none;
        flex: 1 1 0;
      }

      .nav-tabs > li > a {
        text-align: center;
        font-size: 12px;
        padding: 8px 6px;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
      }

      @media (max-width: 1100px) {
        .main-layout-custom {
          flex-direction: column;
        }

        .sidebar-col-custom,
        .map-col-custom {
          width: 100%;
        }
      }"
    ))
  ),
  
  div(
    class = "header-row-custom",
    tags$h1("Pinyon-Juniper Woodlands R Shiny Dashboard", class = "header-title-custom")
  ),
  
  fluidRow(
    column(
      width = 4,
      div(
        class = "sidebar-panel-custom",
        
        tabsetPanel(
          id = "left_tab",
          selected = "project_overview",
          
          tabPanel(
            title = "Project Overview",
            value = "project_overview",
            
            div(
              class = "sidebar-image-custom",
              tags$img(src = pj_img_src, alt = "Pinyon-Juniper ecosystem photo")
            ),
            
            tags$h4(tags$b("Research Background:")),
            p("Welcome!"),
            p(
              "This dashboard is meant to display the spatial extent of my research's field campaign data, as well as visualize the climatic conditions and topographic",
              "variability within PJ ecosystems"
            ),
            p(
              "The field campaign for this project took place in 2021 and 2022 included measurements from a total of 18 field sites, 180 field",
              "plots, and 9447 field sampled trees."
            ),
            p("There are three primary tabs for viewing the different data products associated with this project:"),
            p(
              tags$ul(
                style = "margin-left: 22px; margin-top: 8px;",
                tags$li(
                  tags$b("Spatial Extent Tab: "),
                  "View the spatial extent of PJ woodlands throughout the southwestern USA, as well as mean 30-year precipitation and landscape elevation."
                ),
                tags$li(
                  tags$b("Stratification Tab: "),
                  "View the 9-class stratification to understand overall variability in topographic and climatic conditions."
                ),
                tags$li(
                  tags$b("Field Data Tab: "),
                  "Zoom to individual sites, plots, and field-sampled trees visited during the field campaign."
                )
              )
            )
          ),
          
          tabPanel(
            title = "Spatial Extent",
            value = "spatial_extent",
            
            div(
              class = "sidebar-image-custom",
              tags$img(src = pinyon_img_src, alt = "Pinyon-Juniper ecosystem photo")
            ),
            
            tags$h4(tags$b("Spatial Extent")),
            p(
              "Pinyon-juniper woodlands represent the most extensive dryland forest ecosystem in the western United States. ",
              "They span broad climatic and topographic gradients and support important ecological functions related to ",
              "dryland forest structure, biodiversity, and carbon storage."
            ),
            p(
              "Use the layer checkboxes below to visualize elevation and mean 30-year precipitation across the PJ ecosystem boundary."
            ),
            
            checkboxGroupInput(
              inputId = "overview_layers",
              label = "Spatial extent layers:",
              choices = c(
                "Elevation (DEM)" = "dem",
                "Precipitation" = "ppt"
              ),
              selected = character(0)
            ),
            
            p(
              "Elevation and precipitation provide broad environmental context for understanding where PJ woodlands occur. ",
              "Together, these layers help illustrate how topographic and climatic gradients structure the landscape."
            )
          ),
          
          tabPanel(
            title = "Stratification",
            value = "stratification",
            
            tags$h4(tags$b("Stratification")),
            p(
              "This is a nine-class stratification derived from two variables: ",
              tags$b("30-year mean precipitation"),
              " and ",
              tags$b("elevation"),
              ". Together, these variables represent the broad climatic and topographic gradients across the PJ ecosystem."
            ),
            
            div(
              class = "strata-info-box",
              tags$h4("Strata Interpretation"),
              p(
                "This landscape stratification can be interpreted as a bivariate classification, where each class represents a unique combination ",
                "of precipitation and elevation conditions. Lower-numbered classes correspond to lower-elevation and drier settings, while ",
                "higher-numbered classes correspond to higher-elevation and wetter settings. The combinations of these variables are shown below."
              ),
              tags$table(
                class = "strata-table",
                tags$thead(
                  tags$tr(
                    tags$th(""),
                    tags$th("Low Precip"),
                    tags$th("Med Precip"),
                    tags$th("High Precip")
                  )
                ),
                tags$tbody(
                  tags$tr(
                    tags$td(class = "strata-row-label", "High Elev"),
                    tags$td(class = "strata-cell", style = "background:#80FF00", "7"),
                    tags$td(class = "strata-cell", style = "background:#228B22; color:#fff", "8"),
                    tags$td(class = "strata-cell", style = "background:#006400; color:#fff", "9")
                  ),
                  tags$tr(
                    tags$td(class = "strata-row-label", "Med Elev"),
                    tags$td(class = "strata-cell", style = "background:#FFC000", "4"),
                    tags$td(class = "strata-cell", style = "background:#FFFF00", "5"),
                    tags$td(class = "strata-cell", style = "background:#C0FF00", "6")
                  ),
                  tags$tr(
                    tags$td(class = "strata-row-label", "Low Elev"),
                    tags$td(class = "strata-cell", style = "background:#8B0000; color:#fff", "1"),
                    tags$td(class = "strata-cell", style = "background:#FF4000; color:#fff", "2"),
                    tags$td(class = "strata-cell", style = "background:#FF8000", "3")
                  )
                )
              )
            )
          ),
          
          tabPanel(
            title = "Field Data",
            value = "field_data",
            
            tags$h4(tags$b("Field Data")),
            p(
              "Field data were collected across 18 sites spanning the environmental conditions present across the PJ landscape. ",
              "Select a site to view plot extents and individual tree locations, then zoom to a specific plot to explore trees in detail."
            ),
            
            selectInput(
              inputId = "site_select",
              label = "Zoom to a site:",
              choices = c("Full extent", sort(unique(sites_sf_ll$site_label))),
              selected = "Full extent"
            ),
            
            selectInput(
              inputId = "plot_select",
              label = "Zoom to a plot:",
              choices = "All plots",
              selected = "All plots"
            ),
            
            checkboxInput(
              inputId = "show_pj_extent_field",
              label = "Show PJ ecosystem extent",
              value = FALSE
            )
          )
        )
      )
    ),
    
    column(
      width = 8,
      div(
        class = "map-panel-custom",
        leafletOutput("map", height = "800px")
      )
    )
  )
)

### Initialize Server ----------------------------------------------------------
server <- function(input, output, session) {
  
  full_bb  <- st_bbox(pj_boundary_sf_ll)
  sites_bb <- st_bbox(sites_sf_ll)
  
  overview_center_lng <- mean(c(as.numeric(full_bb["xmin"]), as.numeric(full_bb["xmax"])))
  overview_center_lat <- mean(c(as.numeric(full_bb["ymin"]), as.numeric(full_bb["ymax"])))
  
  spatial_extent_zoom <- 6
  stratification_zoom <- 6
  
  add_hybrid_basemap <- function(map_obj) {
    map_obj |>
      clearTiles() |>
      addProviderTiles(
        "Esri.WorldImagery",
        options = providerTileOptions(maxNativeZoom = 19, maxZoom = 23)
      ) |>
      addProviderTiles(
        "CartoDB.PositronOnlyLabels",
        options = providerTileOptions(opacity = 0.9)
      )
  }
  
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(zoomControl = TRUE, maxZoom = 23))
  })
  
  draw_project_overview <- function(reset_view = FALSE) {
    proxy <- leafletProxy("map", session = session) |>
      clearControls() |>
      clearImages() |>
      clearShapes() |>
      clearMarkers() |>
      clearGroup("selected_site")
    
    if (reset_view) {
      proxy <- add_hybrid_basemap(proxy) |>
        fitBounds(
          lng1 = as.numeric(sites_bb["xmin"]),
          lat1 = as.numeric(sites_bb["ymin"]),
          lng2 = as.numeric(sites_bb["xmax"]),
          lat2 = as.numeric(sites_bb["ymax"])
        )
    }
    
    legend_html <- HTML("
      <div style='
        background: rgba(0,0,0,0.75);
        padding: 10px 12px;
        border-radius: 10px;
        color: white;
        font-size: 13px;
        line-height: 1.5;
        box-shadow: 0 2px 8px rgba(0,0,0,0.35);
      '>
        <div style='font-weight: 700; margin-bottom: 6px;'>Project Overview</div>

        <div style='display:flex; align-items:center; margin-bottom:5px;'>
          <span style='
            display:inline-block;
            width:18px;
            height:0;
            border-top:3px solid red;
            margin-right:8px;
          '></span>
          <span>PJ Boundary Extent</span>
        </div>

        <div style='display:flex; align-items:center; margin-bottom:5px;'>
          <span style='
            display:inline-block;
            width:18px;
            height:0;
            border-top:3px solid white;
            margin-right:8px;
          '></span>
          <span>Field Site Boundaries</span>
        </div>

        <div style='display:flex; align-items:center;'>
          <span style='
            display:inline-block;
            width:10px;
            height:10px;
            background:#FFD54F;
            border:1px solid black;
            border-radius:50%;
            margin-right:8px;
          '></span>
          <span>Field Plot Locations</span>
        </div>
      </div>
    ")
    
    proxy |>
      addPolygons(
        data = pj_boundary_sf_ll,
        color = "red",
        weight = 2,
        fill = TRUE,
        fillColor = "red",
        fillOpacity = 0.08,
        label = "PJ Ecosystem Boundary"
      ) |>
      addPolygons(
        data = sites_sf_ll,
        color = "white",
        weight = 2,
        fill = FALSE,
        label = ~site_label
      ) |>
      addCircleMarkers(
        data = plot_pts_sf_ll,
        radius = 5,
        color = "black",
        weight = 1,
        fillColor = "#FFD54F",
        fillOpacity = 0.95,
        stroke = TRUE,
        label = ~plot_label
      ) |>
      addControl(
        html = legend_html,
        position = "bottomright"
      )
  }
  
  draw_spatial_extent <- function(reset_view = FALSE) {
    proxy <- leafletProxy("map", session = session) |>
      clearControls() |>
      clearImages() |>
      clearShapes() |>
      clearMarkers() |>
      clearGroup("selected_site")
    
    if (reset_view) {
      proxy <- add_hybrid_basemap(proxy) |>
        setView(
          lng = overview_center_lng,
          lat = overview_center_lat,
          zoom = spatial_extent_zoom
        )
    }
    
    if ("dem" %in% input$overview_layers) {
      proxy <- proxy |>
        addRasterImage(
          dem_overview_leaf,
          colors = dem_pal,
          opacity = 0.8,
          project = FALSE
        ) |>
        addLegend(
          position = "bottomleft",
          pal = dem_pal,
          values = app_meta$dem_range,
          title = "Elevation (m)"
        )
    }
    
    if ("ppt" %in% input$overview_layers) {
      proxy <- proxy |>
        addRasterImage(
          ppt_overview_leaf,
          colors = ppt_pal,
          opacity = 0.75,
          project = FALSE
        ) |>
        addLegend(
          position = "bottomright",
          pal = ppt_pal,
          values = app_meta$ppt_range,
          title = "Precipitation"
        )
    }
    
    show_env_layer <- any(c("dem", "ppt") %in% input$overview_layers)
    
    proxy |>
      addPolygons(
        data = pj_boundary_sf_ll,
        color = "red",
        weight = 3,
        fill = !show_env_layer,
        fillColor = "red",
        fillOpacity = 0.25,
        label = "PJ Ecosystem Boundary"
      ) |>
      addLegend(
        position = "topright",
        colors = "red",
        labels = "PJ Ecosystem Boundary",
        opacity = 1,
        title = NULL
      )
  }
  
  draw_field_data <- function(reset_view = FALSE) {
    site_sel <- input$site_select
    
    selected_site <- if (!identical(site_sel, "Full extent")) {
      sites_sf_ll[trimws(sites_sf_ll$site_label) == trimws(site_sel), , drop = FALSE]
    } else {
      NULL
    }
    
    plot_sel <- if (
      !is.null(input$plot_select) &&
      !identical(input$plot_select, "All plots")
    ) input$plot_select else "All plots"
    
    target_bb <- if (!is.null(selected_site) && nrow(selected_site) > 0) {
      if (!identical(plot_sel, "All plots")) {
        sel_plot <- plots_sf_ll[plots_sf_ll$plot_label == plot_sel, , drop = FALSE]
        if (nrow(sel_plot) > 0) {
          st_bbox(sel_plot)
        } else {
          st_bbox(st_union(st_make_valid(selected_site)))
        }
      } else {
        st_bbox(st_union(st_make_valid(selected_site)))
      }
    } else {
      sites_bb
    }
    
    proxy <- leafletProxy("map", session = session) |>
      clearControls() |>
      clearImages() |>
      clearShapes() |>
      clearMarkers() |>
      clearGroup("selected_site")
    
    if (reset_view) {
      proxy <- add_hybrid_basemap(proxy)
    }
    
    proxy <- proxy |>
      fitBounds(
        lng1 = as.numeric(target_bb["xmin"]),
        lat1 = as.numeric(target_bb["ymin"]),
        lng2 = as.numeric(target_bb["xmax"]),
        lat2 = as.numeric(target_bb["ymax"])
      )
    
    if (isTRUE(input$show_pj_extent_field)) {
      proxy <- proxy |>
        addPolygons(
          data = pj_boundary_sf_ll,
          color = "red",
          weight = 2,
          fill = FALSE,
          label = "PJ Ecosystem Boundary"
        )
    }
    
    proxy <- proxy |>
      addPolygons(
        data = sites_sf_ll,
        color = "white",
        weight = 2,
        fill = FALSE,
        label = ~site_label
      )
    
    if (!is.null(selected_site) && nrow(selected_site) > 0) {
      site_id_val <- selected_site$site_id[1]
      
      site_plots <- plots_sf_ll[plots_sf_ll$site_id == site_id_val, , drop = FALSE]
      if (nrow(site_plots) > 0) {
        proxy <- proxy |>
          addPolygons(
            data = site_plots,
            color = "#FFD54F",
            weight = 2,
            fill = FALSE,
            label = ~plot_label
          )
      }
      
      site_trees <- trees_ll[
        !is.na(trees_ll$site_id) & trees_ll$site_id == site_id_val,
        ,
        drop = FALSE
      ]
      
      if (!identical(plot_sel, "All plots")) {
        sel_plot <- plots_sf_ll[plots_sf_ll$plot_label == plot_sel, , drop = FALSE]
        if (nrow(sel_plot) > 0) {
          sel_plot_id <- sel_plot$plot_id[1]
          site_trees <- site_trees[site_trees$plot_id == sel_plot_id, , drop = FALSE]
        } else {
          site_trees <- site_trees[0, , drop = FALSE]
        }
      }
      
      if (nrow(site_trees) > 0) {
        proxy <- proxy |>
          addCircleMarkers(
            data = site_trees,
            radius = tree_dot_radius,
            color = "black",
            weight = 0.7,
            fillColor = ~tree_pal(Species),
            fillOpacity = 0.95,
            popup = ~paste0(
              "<b>", tree_id, "</b><br>",
              "Species: ", Species, "<br>",
              "Diameter: ", tree_dm, " cm<br>",
              "Biomass: ", round(tr_bms_, 1), " kg<br>",
              "Condition: ", tre_cnd
            )
          ) |>
          addLegend(
            position = "bottomright",
            pal = tree_pal,
            values = site_trees$Species,
            title = "Tree Species",
            opacity = 1
          )
      }
      
      proxy <- proxy |>
        addPolygons(
          data = selected_site,
          color = "cyan",
          weight = 3,
          fill = FALSE,
          group = "selected_site"
        )
      
      if (!identical(plot_sel, "All plots")) {
        sel_plot <- plots_sf_ll[plots_sf_ll$plot_label == plot_sel, , drop = FALSE]
        if (nrow(sel_plot) > 0) {
          proxy <- proxy |>
            addPolygons(
              data = sel_plot,
              color = "cyan",
              weight = 2,
              fill = FALSE
            )
        }
      }
      
    } else {
      proxy <- proxy |>
        addCircleMarkers(
          data = plot_pts_sf_ll,
          radius = 5,
          color = "black",
          weight = 1,
          fillColor = "#FFD54F",
          fillOpacity = 0.95,
          stroke = TRUE,
          label = ~plot_label
        )
    }
  }
  
  observeEvent(input$left_tab, ignoreInit = FALSE, {
    req(input$left_tab)
    
    if (input$left_tab == "project_overview") {
      
      draw_project_overview(reset_view = TRUE)
      
    } else if (input$left_tab == "spatial_extent") {
      
      draw_spatial_extent(reset_view = TRUE)
      
    } else if (input$left_tab == "stratification") {
      
      proxy <- leafletProxy("map", session = session) |>
        clearControls() |>
        clearImages() |>
        clearShapes() |>
        clearMarkers() |>
        clearGroup("selected_site")
      
      proxy <- add_hybrid_basemap(proxy)
      
      proxy |>
        setView(
          lng = overview_center_lng,
          lat = overview_center_lat,
          zoom = stratification_zoom
        ) |>
        addRasterImage(
          finalStrat_leaf,
          colors = strata_pal,
          opacity = 0.55,
          project = FALSE
        ) |>
        addPolygons(
          data = pj_boundary_sf_ll,
          color = "red",
          weight = 2,
          fill = FALSE
        ) |>
        addLegend(
          position = "bottomright",
          pal = strata_pal,
          values = strata_vals,
          title = "Strata"
        )
      
    } else if (input$left_tab == "field_data") {
      
      draw_field_data(reset_view = TRUE)
    }
  })
  
  observeEvent(input$overview_layers, ignoreInit = TRUE, {
    req(input$left_tab == "spatial_extent")
    draw_spatial_extent(reset_view = FALSE)
  })
  
  observeEvent(input$site_select, ignoreInit = TRUE, {
    req(input$left_tab == "field_data")
    
    if (!identical(input$site_select, "Full extent")) {
      selected_site <- sites_sf_ll[
        trimws(sites_sf_ll$site_label) == trimws(input$site_select),
        ,
        drop = FALSE
      ]
      
      if (nrow(selected_site) > 0) {
        site_id_val <- selected_site$site_id[1]
        site_plots <- plots_sf_ll[plots_sf_ll$site_id == site_id_val, , drop = FALSE]
        
        updateSelectInput(
          session,
          "plot_select",
          choices = c("All plots", sort(unique(site_plots$plot_label))),
          selected = "All plots"
        )
      }
    } else {
      updateSelectInput(
        session,
        "plot_select",
        choices = "All plots",
        selected = "All plots"
      )
    }
    
    draw_field_data(reset_view = FALSE)
  })
  
  observeEvent(input$plot_select, ignoreInit = TRUE, {
    req(input$left_tab == "field_data")
    draw_field_data(reset_view = FALSE)
  })
  
  observeEvent(input$show_pj_extent_field, ignoreInit = TRUE, {
    req(input$left_tab == "field_data")
    draw_field_data(reset_view = FALSE)
  })
}

### Run Application ------------------------------------------------------------
shinyApp(ui, server)