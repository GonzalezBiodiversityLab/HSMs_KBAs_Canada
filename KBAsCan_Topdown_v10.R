


############################################################################
#                           KBAS Canada                                    #
#                                                                          #
#                  "Top-down approach to help
#     identifying Key Biodiversity Areas (KBAs) in Canada.                 #
#                                                                          #
#  This top-down approach has been applied to a list of species triggering #
# A1 (threatened species) and B1 (geographically restricted species)
#               criteria defined in the KBA standard.                      #
#                                                                          #
#  by Juan Zuloaga (Post-doctoral Researcher) & A7ysundrew Gonzalez (Professor)# 
#                     McGill University                                    #
#                                                                          #
#                         August 2020                                    #
############################################################################ 


############################## CONTENT ##############################

# 0. Installing and loading packages
# 1. General settings 
# 2. Loading data
# 2.1 Species
# 2.2 Maps 
# 2.3 Geographic ranges
# 2.4 Predictors
# 3. Looping to create SDMs reports for every single species

###

#options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx12288M"))
options(java.parameters = "-Xmx12g")

# 0. INSTALLING A LOADING PACKAGES ----------

Packages <- c(
  "stringr",
  "knitr",
  "visNetwork", #workflow visualization
  "readxl", # read Excel .xls
  "plyr", # Combine data.frames by row, filling in missing columns.
  "dplyr", # create a new variable using mutate function
  "raster", # for raster
  "rgdal", # for spatial data
  "maptools", # read polygons
  "rgeos", # dissolve and intersect polygons
  "sp",
  "SSDM", # SDMs
  "leaflet", # interactive maps
  "shiny", # interactive maps
  "dismo", #read bioclim and SDMs
  "virtualspecies", # Remove collinearity among variables of a a raster stack
  "sf", # intersect polygonsingles
  "ecospat", # cross validation maxent
  "ggplot2", #plots
  "KernSmooth", # hotspots
  "spatstat", # kernel density
  "forcats", #order numbers
  "stars", # raster to polygon
  "classInt",
  "XML",
  "SSDM", #Ensemble distribution models 
  "maps",
  "arcgisbinding", # connect to NatureServe data set, using R-ArcGIS pro
  "tidyverse",
  "tidyr", # extract raster values
  "spThin", # cleaning data
  "here", # relative paths
  "adehabitatHR", # to calculate the maximum convex polygon (number of grids used as a background point in maxent)
  "exactextractr", # to extract values from raster land cover
  "ENMeval", # v2.0to evaluate optimal model complexity for Maxent models
  "tidyselect", #"all_of function was moved from tinyselect to tidyselect"
  "profvis", # profiling
  "fasterize", # fast rasterize
  "cowplot", # adding one legend to ggplot
  "gridExtra",
  "kableExtra"
) 

# Loading packages
lapply(Packages, library, character.only =TRUE) 

# 1. GENERAL SETTINGS ----------
# Set working directory
setwd("C:/HSMs_KBAsCan")
getwd() # verify
 

# Creating directories to save points/maps
dir.create("./Intermediate/", showWarnings=F)
dir.create("./Data_inputs/", showWarnings=F)
dir.create("./HSMs_html_outputs/", showWarnings=F)
dir.create("./Intermediate/Points_thinned_all_wgs84/", showWarnings=F)
dir.create("./Intermediate/Points_thinned_EBARs_aeac/", showWarnings=F)
dir.create("./Intermediate/hsm_map_best_all/", showWarnings=F)
dir.create("./Intermediate/hsm_map_uncertanty_all/", showWarnings=F)
dir.create("./Intermediate/sdm_polygons_KBAs_Top3/", showWarnings=F)
dir.create("./Intermediate/hsm_summary_performance/", showWarnings=F)
dir.create("./Intermediate/hsm_map_future/", showWarnings=F)
dir.create("./Intermediate/hsm_novel_climate/", showWarnings=F)
dir.create ("./Intermediate/hsm_polygons_KBAs/", showWarnings=F)
dir.create ("./Intermediate/hsm_polygons_KBAs_Top3/", showWarnings=F)
dir.create ("./Intermediate/hsm_polygons_area/", showWarnings=F)
dir.create ("./Intermediate/hsm_hfp_area/", showWarnings=F)


# Creating a directory to save temporary files that will be deleted after each loop
dir.create (file.path("./temp_SDMs_to_REMOVE"), showWarnings = FALSE)

# Setting temp directory
rasterOptions(tmpdir=file.path("./temp_SDMs_to_REMOVE"))
rm(list=ls(all.names = TRUE))

# Setting Projection to preserve area for all files 
# North America Albers Equal Area Conic
#[See:](https://spatialreference.org/ref/esri/102008/proj4/) 
aeac="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

# Lon-Lat (for some analysis)
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "

# Laea
laea <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"

# Increasing memory size for Maxent processing (Java demands more memory)
options(java.parameters = "-Xmx8000m")

# Checking ArcGIS Pro (It requieres having ArcGIS Pro and get autorization from NatureServe to connect with its data sets)
arc.check_product()

# Connecting to the portal URL (NatureServe Data sets)
url <- "https://gis.natureserve.ca/arcgis/rest/services/"


# 2. LOADING DATA ----------
# 2.1 Maps/Templates ----------
#(potential extent of the Analysis is North America)
print("Loading Maps")

# Canada to use in plotting
Canada_cont_prov<-st_read("./Data_inputs/GeopoliticalAreas_p_nolakes.shp")  %>%
st_transform(crs=aeac)   %>% # Projecting it
as('Spatial') #  making it spatial

# Canada lon/lat leaflet map
Canada_WGS84 <- st_read("./Data_inputs/GeopoliticalAreas.shp")  %>%
st_transform(CRS("+proj=longlat +datum=WGS84"))

# North America
US_Canada_Mex_WGS84 <- st_read("./Data_inputs/US_Canada_Mex_WGS84.shp")

# North America
US_Canada_Mex_aeac <- st_read("./Data_inputs/US_Canada_Mex_WGS84.shp")  %>%
st_transform(crs=aeac) %>%
as('Spatial')

# North America Diss
area_thresh <- units::set_units(200000, km^2)
US_Canada_Mex_aeac_diss <- sfheaders::sf_remove_holes(US_Canada_Mex_WGS84)

US_Canada_Mex_aeac_diss <- smoothr::fill_holes(US_Canada_Mex_WGS84, area_thresh)

#st_write(US_Canada_Mex_aeac_diss, "./Intermediate/US_Canada_Mex_WGS84_diss.shp", delete_layer=T)

US_Canada_Mex_aeac_diss2 <- sfheaders::sf_remove_holes(US_Canada_Mex_aeac_diss)

# US states (from https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html)
US_states_wgs84 <- st_read("./Data_inputs/cb_2018_us_state_500k.shp") %>%
st_transform(crs=wgs84)

# Protected Areas (PAs)
# from https://www.iucn.org/theme/protected-areas/our-work/world-database-protected-areas
Protected_Area <- st_read("./Data_inputs/WorldDatabaseOfProtectedAreas_poly.shp") %>%
st_transform(CRS("+proj=longlat +datum=WGS84"))

# Towns (from https://www.nrcan.gc.ca/earth-sciences/geography/download-geographical-names-data/9245)
cities_ca_t <- st_read("./Data_inputs/Can_pop_cities.shp")  %>%
st_transform(CRS(aeac))  %>%
st_transform(CRS("+proj=longlat +datum=WGS84"))

# Human footprint map
hfp_Can_p_r <- raster("./Data_inputs/hfp_Can_p_aeac1000_r.tif")

# Ecoregrions
Ecoregions <- st_read("./Data_inputs/Ecoregion_NA_wgs84.shp")

# 2.2 List of species that might trigger A1 and B1 ----------
print("Loading species data from excel table")

# Species list (from Excel) provided by Chloé Debyser (WCS-KBAs)
species_list <- read_xlsx("./Data_inputs/KBAsamplespecies_modelrun_20200604.xlsx")

# Scientific name
sp_eos_short_names <- unlist(species_list$NATIONAL_SCIENTIFIC_NAME)

# Replacing space between genus and species name with an underscore
sp_eos <- species_list
sp_eos$GLOBAL_SCIENTIFIC_NAME <- chartr(" ", "_", sp_eos$NATIONAL_SCIENTIFIC_NAME)

# Global id
element_id <- unique(species_list$ELEMENT_GLOBAL_ID)

# Some rules for cleaning data
# Removing low accuracy features, If needed (we set here a big value to include all features)
low_accuracy <- 1000000
# Excluding features before ???, if needed (we set here an aold value to include all features)
old_features <- "1400-01-01"

# Information for summary tables
for(i in 1:length(element_id)){
  # Creating a directory to save temporary files that will be deleted after each loop
  dir.create (file.path("./temp_SDMs_to_REMOVE"), showWarnings = FALSE)
  
  # Setting temp directory
  rasterOptions(tmpdir=file.path("./temp_SDMs_to_REMOVE"))
  listObjects <-ls()
  
namesp <-  sp_eos_short_names[i] # Species names
cat(paste0("Species: ", namesp), '\n')
Global_trigger <-  species_list$KBATrigger_G[i] # Global trigger? Y or N
National_trigger <-  species_list$KBATrigger_N[i]  # National trigger? Y or N
National_name <- species_list$NATIONAL_ENGL_NAME[i] #National name
global_name <- species_list$GLOBAL_ENGL_NAME[i] # Global name
global_name_fr <- species_list$GLOB_FR_NAME[i] # Global name French
sp_table_IUCN <- species_list$IUCN_CD[i] # Status IUCN
sp_global_synonyms <- species_list$GLOBAL_SYNONYMS[i] # Globla synonymus
Global_A1 <- species_list$KBATrigger_G_A1[i] # Global criterion A1
Global_B1 <- species_list$KBATrigger_G_B1[i] # Global criterion B1
National_A1 <- species_list$KBATrigger_N_A1[i]
National_B1 <- species_list$KBATrigger_N_B1[i]
cosewic <- species_list$COSEWIC_STATUS[i] # Cosewic status
natureserveG <- species_list$ROUNDED_G_RANK[i] # NAtureServe Status Global
natureserveN <- species_list$ROUNDED_N_RANK[i] # NatureServe Status National
Element_code <- species_list$ELEMENT_CODE[i] # Element code number

# 2.3 Loading data from NAtureServe (Tabular data) ----------
# You required permission to access Nature Serve dataset 9contact (Randal Green)
print("Loading data sets from NatureServe dataset")

# Tabular biotics
biotics <- arc.open(paste0(url, 'EBAR-KBA/Restricted/FeatureServer/4')) %>% # oLD PATH
arc.select(.) # Convert to data frame

# Select species
# biotics_sel <- biotics[biotics$element_global_id %in% i, ]
biotics_sel <- subset(biotics, biotics$element_global_id == element_id[i])

# Species ID
sp_id <-  biotics_sel$speciesid

# where_clause to select species
myWhereClause_species <- paste("speciesid  =", sp_id)

# Tabular Species
species <- arc.open(paste0(url, 'EBAR-KBA/KBA/FeatureServer/8')) %>% # Load
arc.select(.) # Convert to data frame

# Tabular inputDataset
inputDataset <- arc.open(paste0(url, 'EBAR-KBA/Restricted/FeatureServer/7')) %>% # Load
arc.select(.) # Convert to data frame

# Tabular datasetSource
datasetSource <- arc.open(paste0(url, 'EBAR-KBA/Restricted/FeatureServer/5')) %>% # Load
arc.select(.) # Convert to data frame

# Range map
rangeMap <- arc.open(paste0(url, 'EBAR-KBA/Restricted/FeatureServer/10'))  # Load
rangeMap_sel <-  arc.select(rangeMap, where_clause = myWhereClause_species) # Convert to data frame

# range map ID
rangemap_id <- rangeMap_sel$rangemapid

# where_clause to select range map id
myWhereClause_rangemap <- paste("rangemapid  =", rangemap_id)


# 2.4 Loading data from NAtureServe (Spatial data) ----------

# 2.4.1. EcoshapeRangeMap ----------
print("Loading EBARs")

# Loading Ecoshapes from NatureServe
ecoshapeRangeMap <- arc.open(paste0(url, 'EBAR-KBA/EcoshapeRangeMap/FeatureServer/0'))
ecoshapeRangeMap_sel <- list()

# Selecting only species' Ecoshapes
for(ecs in 1:length(myWhereClause_rangemap)){
  ecoshapeRangeMap_sel[[ecs]] <-  arc.select(ecoshapeRangeMap,  where_clause = myWhereClause_rangemap[ecs]) %>%
  arc.data2sp() %>%
  st_as_sf()%>%
  st_union()
  ecoshapeRangeMap_sel_list2 <- do.call("c", ecoshapeRangeMap_sel)
}

# Projecting Ecoshapes to aeac
Ecod_sp_join_id_diss <- st_transform(ecoshapeRangeMap_sel_list2, st_crs(aeac))
Ecod_sp_join_id_diss_s <- as(Ecod_sp_join_id_diss, 'Spatial')
Ecod_sp_join_id_diss_s_wgs84 <- st_transform(st_as_sf(Ecod_sp_join_id_diss_s), crs = wgs84)

# 2.4.2. Polygons (other than EOs, that is ' Critical Habitat') ----------
print("Creating points from polygons")

# Loading Polygons from NatureServe
polygons <- arc.open(paste0(url, 'EBAR-KBA/Restricted/FeatureServer/2'))
polygons_sel <-  arc.select(polygons, where_clause = myWhereClause_species)
polygons_sel_shp <- arc.shape(polygons_sel)

# Checking if data set is empty
  if(length(polygons_sel_shp$shape_buffer) > 0){
    # Convert an ArcGIS arc.data to the equivalent sp data frame type.
      polygons_sel2 <- arc.data2sp(polygons_sel) %>%
      st_as_sf() # Convert foreign object to an sf object
      
      # Projecting them to aeac
        polygons_sel_p <- st_transform(polygons_sel2, crs(aeac))
      # Excluding polygons with low accuracy, If needed
        polygons_sel_p_a <- polygons_sel_p[is.na(polygons_sel_p$accuracy) < low_accuracy, ]
      # Excluding points before ??? or NA If needed 
        polygons_sel_p_amax <- polygons_sel_p_a[(polygons_sel_p_a$maxdate > as.Date(old_features) | is.na(polygons_sel_p_a$maxdate)),]
      # Merging points and tables
        polygons_sel_p_amax_j <- merge(polygons_sel_p_amax, inputDataset, by="inputdatasetid", all.X=T)  # with inputDataset
        polygons_sel_p_amax_j2 <- merge(polygons_sel_p_amax_j, datasetSource, by="datasetsourceid", all.X=T) # with datasetSource
      # Excluding points from these data sets if needed: (replace here for 'GBIF' , 'Ecoengine', 'iDigBio', etc)
        polygons_sel_p_amax_j2_sel <- polygons_sel_p_amax_j2[!(polygons_sel_p_amax_j2$datasetsourcename == 'none'),]
      # Including only Critical Habitat
        polygons_sel_p_amax_j2_sel_p <- polygons_sel_p_amax_j2_sel[(polygons_sel_p_amax_j2_sel$datasettype == 'Critical Habitat'),]
        
        if(!is.empty(polygons_sel_p_amax_j2_sel_p$geometry)){
        
          # Excluding	"Canadian Wildlife Service records"
          polygons_sel_p_amax_j2_sel2 <- polygons_sel_p_amax_j2_sel_p[!(polygons_sel_p_amax_j2_sel_p$datasetname == 'CWS Priority Species'),]

        # remove Unobscured versions of iNaturalist
          #  if(any(unique(polygons_sel_p_amax_j2_sel2$datasetname == "iNaturalist"), na.rm=T) || any(unique(polygons_sel_p_amax_j2_sel2$taxongeoprivacy == "obscure"), na.rm=T)){
          #   polygons_sel_p_amax_j2_sel2 <- polygons_sel_p_amax_j2_sel2[!(polygons_sel_p_amax_j2_sel2$datasetname == "iNaturalist"), ]
          #   polygons_sel_p_amax_j2_sel2 <- polygons_sel_p_amax_j2_sel2[!(polygons_sel_p_amax_j2_sel2$taxongeoprivacy == "obscure"), ]
                 # }else{
              #  }
        
        # remove RECORDS FROM BC DATA_SENSITIVE
          if(any(unique(polygons_sel_p_amax_j2_sel2$subnation == 'BC'), na.rm=T) || any(unique(polygons_sel_p_amax_j2_sel2$datasensitivity =="Yes"), na.rm=T)){
          
          polygons_sel_p_amax_j2_sel2 <- polygons_sel_p_amax_j2_sel2[(polygons_sel_p_amax_j2_sel2$datasensitivity == "No"), ]
          polygons_sel_p_amax_j2_sel2 <- polygons_sel_p_amax_j2_sel2[(polygons_sel_p_amax_j2_sel2$datasensitivity == "NA"), ]

          }else{
          }
        
      # Creating points
        if(!is.empty(polygons_sel_p_amax_j2_sel2$geometry)){
          # Create a raster template to convert polygons to rasters and then to points, using DEMs 30m resolution
            # (it has to be high resolution, otherwise small polygons will be ignored)
            polygons_sel_pe <- st_bbox(polygons_sel_p_amax_j2_sel2) # get extent
            epol <- extent(polygons_sel_pe[c(1,3,2,4)]) # create extent
            epol_raster <- raster(epol, res=50, crs=aeac) # create raster 30m raster
            values(epol_raster)<-1 # add values
            epol_raster_1km <- raster(epol, res=1000, crs=aeac) # create raster 1km
            
            # Rasterizing
            poly_raster <- rasterize(polygons_sel_p_amax_j2_sel2[1], epol_raster, fun = "first", na.rm=TRUE)  #rasterize 30m
            poly_raster_1km <- resample(aggregate(poly_raster, fact=20), epol_raster_1km) # rasterize 1000m
            
            # Converting to points
            poly_points <- rasterToPoints(poly_raster_1km, spatial=T) 
            
        }else{
          
        }
        }else{
         # poly_points <- st_sf(st_sfc())
          poly_points <-  SpatialPoints(data.frame(x = 0, y = 0))[-1,]
        }  
            }else{
              
              #poly_points <- st_sf(st_sfc())
              poly_points <-  SpatialPoints(data.frame(x = 0, y = 0))[-1,]
              
              }
      

# 2.4.3. Points ----------
print("Loading point features")
# loading all points
  points <- arc.open(paste0(url, 'EBAR-KBA/Restricted/FeatureServer/0'))
# selecting points for each species
  points_sel <- arc.select(points, where_clause = myWhereClause_species)
  points_sel_shp <- arc.shape(points_sel)
# Checking if data set is empty
  if(length(points_sel_shp) > 0){
    # Convert an ArcGIS arc.data to the equivalent sp data frame type.
      points_sel2 <- arc.data2sp(points_sel) %>%
        st_as_sf() # Convert foreign object to an sf object
      # Projecting them to aeac
        points_sel_p <- st_transform(points_sel2, crs(aeac))
      # Excluding points with low accuracy, if needed
        points_sel_p_a <- points_sel_p[is.na(points_sel_p$accuracy) < low_accuracy, ]  # Accuracy: Exclude data with an uncertainty distance greater than 32 km.
      # Excluding points before ???? and NAs, If needed
        points_sel_p_amax <- points_sel_p_a[(points_sel_p_a$maxdate > as.Date(old_features) | is.na(points_sel_p_a$maxdate)),] # Exclude points older than 1950
      # Merging points and tables
        points_sel_p_amax_j <- merge(points_sel_p_amax, inputDataset, by="inputdatasetid", all.X=T)  # with inputDataset
        points_sel_p_amax_j2 <- merge(points_sel_p_amax_j, datasetSource, by="datasetsourceid", all.X=T) # with datasetSource
      # Excluding points from these data sets, if needed (e.g., ) 'GBIF' , 'Ecoengine', 'iDigBio', etc)
        points_sel_p_amax_j2_sel <- points_sel_p_amax_j2[!(points_sel_p_amax_j2$datasetsourcename == 'none'),]
      # Excluding some polygons (Area of Occupancy, Other,  Other Observations, Other Range), if needed
        points_sel_p_amax_j2_sel_p <- points_sel_p_amax_j2_sel[!(points_sel_p_amax_j2_sel$datasettype == 'none'),]
      
     #  st_coordinates(points_sel_p_amax_j2_sel2)
        
        if(!is.empty(points_sel_p_amax_j2_sel_p$geometry)){
          
        # Excluding	"Canadian Wildlife Service records"
          points_sel_p_amax_j2_sel2 <- points_sel_p_amax_j2_sel_p[!(points_sel_p_amax_j2_sel_p$datasetname == "Ecoengine"), ]
          
        # remove Unobscured versions of iNaturalist
        if(any(unique(points_sel_p_amax_j2_sel2$datasetname == "iNaturalist"), na.rm=T) || any(unique(points_sel_p_amax_j2_sel2$taxongeoprivacy == "obscure"), na.rm=T)){
          
          points_sel_p_amax_j2_sel2 <- points_sel_p_amax_j2_sel2[!(points_sel_p_amax_j2_sel2$datasetname == "iNaturalist"), ] &&  points_sel_p_amax_j2_sel2[!(points_sel_p_amax_j2_sel2$taxongeoprivacy == "obscure"), ]
          
        }else{
        }

        # remove RECORDS FROM BC DATA_SENSITIVE
                  
        if(any(unique(points_sel_p_amax_j2_sel2$subnation == 'BC'), na.rm=T) && any(unique(points_sel_p_amax_j2_sel2$datasensitivity =="Yes"), na.rm=T)){
        
           points_sel_p_amax_j2_sel2 <- points_sel_p_amax_j2_sel2[!(points_sel_p_amax_j2_sel2$subnation == 'BC'), ] && points_sel_p_amax_j2_sel2[!(points_sel_p_amax_j2_sel2$datasensitivity == 'Yes'), ]
        
      }else{
       
      }
        }
  }else{
        points_sel_p_amax_j2_sel2 <- st_sf(st_sfc())
      }
  
# 2.4.4. Lines ----------
  print("Creating points from lines")
  # loading all points
    lines <- arc.open(paste0(url, 'EBAR-KBA/Restricted/FeatureServer/1'))
  # Convert an ArcGIS arc.data to the equivalent sp data frame type.
    lines_sel <- arc.select(lines,  where_clause = myWhereClause_species)
    lines_sel_shp <- arc.shape(lines_sel)
    
  # Checking if data set is empty
    if(length(lines_sel_shp) > 0){
      lines_sel2 <- arc.data2sp(lines_sel) %>%
        st_as_sf()
      # Projecting them to aeac
        line_test_p <- st_transform(lines_sel2, aeac)
      # Excluding lines with low accuracy. Excluding data with an uncertainty distance greater than 32 km.
        line_test_p_a <- line_test_p[is.na(line_test_p$accuracy) < low_accuracy, ]
      # Excluding lines before ???? and NAs, if needed
        line_test_p_amax <- line_test_p_a[(line_test_p_a$maxdate > as.Date(old_features) | is.na(line_test_p_a$maxdate)),]
      # Merging lines and tables
        line_test_p_amax_j <- merge(line_test_p_amax, inputDataset, by="inputdatasetid", all.X=T)  # with inputDataset
        line_test_p_amax_j2 <- merge(line_test_p_amax_j, datasetSource, by="datasetsourceid", all.X=T) # with datasetSource
      # Excluding lines if needed (e.g., 'GBIF' , 'Ecoengine', 'iDigBio')
        line_test_p_amax_j2_sel <- line_test_p_amax_j2[!(line_test_p_amax_j2$datasetsourcename == 'none'),]
      # Buffering lines (25m)
        line_test_p_buf <- st_buffer(line_test_p_amax_j2_sel, dist=25) # buffer 25 m each side
      # Create a raster template to convert polygons to rasters and then to points, using DEMs 30m resolution
        # (it has to be high resolution, otherwise small polygons will be ignored)
        lines_sel_pe <- st_bbox(line_test_p_buf) # get extent
        elines <- extent(lines_sel_pe[c(1,3,2,4)]) # create extent
        elines_raster <- raster(elines, res=50, crs=aeac) # create raster
        values(elines_raster)<-1 # add values
        elines_raster_1km <- raster(elines, res=1000, crs=aeac) # create raster 1km
      # line_test_p_buf_s <- as(line_test_p_buf, "Spatial")         # convert to spatialdataframe
        lines_raster <- rasterize(line_test_p_buf, elines_raster, fun = "first", na.rm=TRUE)  #rasterize
        lines_raster_1km <- resample(aggregate(lines_raster, fact=20), elines_raster_1km) # rasterize 1000m
        lines_points <- rasterToPoints(lines_raster_1km, spatial=T) # Converting to points
        
        }else{
          lines_points <- SpatialPoints(data.frame(x = 0, y = 0))[-1,]
        }  

    
# 2.4.5. EOs ----------
  print("Creating points from EOs")
  # loading all points
    eos <- arc.open(paste0(url, 'EO_Polygons/FeatureServer/0'))
    eos_sel <- arc.select(eos, where_clause = myWhereClause_species)
    eos_sel_shp <- arc.shape(eos_sel)
  # Checking if data set is empty
    if(length(eos_sel_shp) > 0){
      eos_sel2 <- arc.data2sp(eos_sel) %>%
        st_as_sf()
      # Projection to aeac
        eos_sel_p <- st_transform(eos_sel2, crs(aeac))
      # Excluding EOs with low accuracy, IF NEEDED
        eos_sel_p_a <- eos_sel_p[is.na(eos_sel_p$accuracy) < low_accuracy, ]
      # Excluding points before ???? and NAs, if needed
        eos_sel_p_amax <- eos_sel_p_a[(eos_sel_p_a$maxdate > as.Date(old_features) | is.na(eos_sel_p_a$maxdate)),]
      # Merging points and tables
        eos_sel_p_amax_j <- merge(eos_sel_p_amax, inputDataset, by="inputdatasetid", all.X=T)  # with inputDataset
      eos_sel_p_amax_j2 <- merge(eos_sel_p_amax_j, datasetSource, by="datasetsourceid", all.X=T) # with datasetSource
      # Excluding some polygons (Area of Occupancy, Other,  Other Observations, Other Range)
        eos_sel_p_amax_j2_p <- eos_sel_p_amax_j2[!(eos_sel_p_amax_j2$datasettype == 'Area of Occupancy' |
                                                   eos_sel_p_amax_j2$datasettype == 'Area of Other' |
                                                   eos_sel_p_amax_j2$datasettype == 'Other' |
                                                   eos_sel_p_amax_j2$datasettype == 'Other Range'),]
      # Excluding points if needed (e.g., from 'GBIF' , 'Ecoengine', 'iDigBio.)
        eos_sel_p_amax_j2_sel <- eos_sel_p_amax_j2_p[!(eos_sel_p_amax_j2_p$datasetsourcename == 'none'),]
        eos_sel_p_amax_j2_sel2 <- eos_sel_p_amax_j2_sel

        ##############
 
 # Creating points
 
      if(!is.null(eos_sel_p_amax_j2_sel2$geometry)){
        # Create a raster template to convert polygons to rasters and then to points, using DEMs 30m resolution
          # (it has to be high resolution, otherwise small polygons will be ignored)
          eos_sel_pe <- st_bbox(eos_sel_p_amax_j2_sel2) # get extent
          eeos <- extent(eos_sel_pe[c(1,3,2,4)]) # create extent
          eeos_raster <- raster(eeos, res=50, crs=aeac) # create raster
          values(eeos_raster)<-1 # add values
          eeos_raster_1km <- raster(eeos, res=1000, crs=aeac) # create raster 1km
          # eos_sel_ps <- as(eos_sel_p_amax_j2_sel2,'Spatial') # convert to spatialdataframe
          eos_raster <- rasterize(eos_sel_p_amax_j2_sel2[1], eeos_raster, fun = "first", na.rm=TRUE)  #rasterize
          eos_raster_1km <- resample(aggregate(eos_raster, fact=20), eeos_raster_1km) # rasterize 1000m
          eos_points <- rasterToPoints(eos_raster_1km, spatial=T) # Converting to points
         
        }else{
          
        }
        }else{
          eos_points <-  SpatialPoints(data.frame(x = 0, y = 0))[-1,]
          
        }  

# 2.4.6 Merging points --------- 
  print("Merging points")
  zero = SpatialPoints(data.frame(x = 0, y = 0))[-1,]
  # Converting to data frames
    # Points
    if(!is.empty(points_sel_p_amax_j2_sel2$geometry)){
      points_df <- data.frame(as_Spatial(points_sel_p_amax_j2_sel2 $geometry)) # points
        names(points_df)[1] <- "x"
      names(points_df)[2] <- "y"
      }else{
        points_df <- as.data.frame(if(exists("zero@coords")) zero@coords)
        }
  # EOs
    if(!is.empty(eos_points)){
      eos_points_df <- data.frame(eos_points@coords)
      }else{
        eos_points_df <- as.data.frame(if(exists("zero@coords")) zero@coords)
        }
  #Polygons
    if(!is.empty(poly_points)){
      poly_points_df <- data.frame(poly_points@coords)
      }else{
        poly_points_df <- as.data.frame(if(exists("zero@coords")) zero@coords)
        }
  # Lines
    if(is.empty(lines_points)){
      lines_df <- data.frame(lines_points@coords)
      }else{
        lines_df <- as.data.frame(if(exists("zero@coords")) zero@coords)
        }
  
  # Merging them
    dat_df <- rbind(if(exists("poly_points_df")) poly_points_df, if(exists("eos_points_df")) eos_points_df)
    dat_df_2 <- rbind(if(exists("dat_df")) dat_df, if(exists("points_df")) points_df)
    dat_df_3 <- rbind(if(exists("dat_df_2")) dat_df_2, if(exists("lines_df")) lines_df)
  
  # Removing duplicates and NAs
    merged_points_duplicated_nona <- dat_df_3[!duplicated(c(dat_df_3$x, dat_df_3$y)),] %>%
    na.omit()
    
    if(nrow(merged_points_duplicated_nona) <= 2){
      
      rmarkdown::render(input = "./Output_KBAs_v7_empty_v2.Rmd",
                        output_format = "html_document",
                        output_file =  paste0(sp_eos$GLOBAL_SCIENTIFIC_NAME[i], ".html"),
                        output_dir = "./HSMs_html_outputs")
      
    #############
      print("Unable to run models and any further analyses due to occurrences <=2")
    ##############
        unlink(file.path("./temp_SDMs_to_REMOVE"), recursive = TRUE)
      rm(list=setdiff(ls(), listObjects))
     
     next
   
    }else{
    }
    
  # Creating geometry and projecting them
    merged_points_p_aeac <- st_geometry(st_as_sf(merged_points_duplicated_nona, coords =c("x", "y")))  %>%
    st_set_crs(crs(aeac))
  # Removing points that are closer 4Km to each other (Alternative method)
    merged_points_p_wgs84 <- st_transform(merged_points_p_aeac, crs(wgs84))%>%
    as('Spatial')%>%
    as.data.frame()
    merged_points_p_wgs84$SPEC <- "sp1"
  # Spatial thining
    print("Spatial thining (points)")
    merged_points_p_thin <- thin(merged_points_p_wgs84,
                               lat.col = "coords.x2",
                               long.col = "coords.x1",
                               spec.col = "SPEC",
                               thin.par =  5, reps=1,
                               locs.thinned.list.return = TRUE,
                               write.files = FALSE,
                               write.log.file = FALSE)
  # Number of points
    number_of_points <- data.frame(merged_points_p_thin)
    
    if(nrow(number_of_points) <= 2){
      rmarkdown::render(input = "./Output_KBAs_v7_empty_v2.Rmd",
                        output_format = "html_document",
                        output_file =  paste0(sp_eos$GLOBAL_SCIENTIFIC_NAME[i], ".html"),
                        output_dir = "./HSMs_html_outputs")
      
    ###################
      print("Unable to run models and any further analyses due to occurrences <=2")
    #################
   
      unlink(file.path("./temp_SDMs_to_REMOVE"), recursive = TRUE)
      rm(list=setdiff(ls(), listObjects))
   
      next
   
    }else{
    }
    
    
    
    
  # Total points WGS84
    points_thin_sf <- st_as_sf(data.frame(number_of_points), coords=c("Longitude","Latitude"), crs=wgs84)
    points_thin_sf_c <- st_coordinates(points_thin_sf)
    st_write(points_thin_sf, paste0("./Intermediate/Points_thinned_all_wgs84/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]), sp_eos$GLOBAL_SCIENTIFIC_NAME[i], driver="ESRI Shapefile", append=FALSE)  # create to a shapefile )
    
  # Total points aeac
    points_thin_sf_total <- st_coordinates(st_transform(points_thin_sf, aeac))
    points_thin_sf_total_s <- SpatialPoints(as(st_transform(points_thin_sf, aeac), 'Spatial'))
    
  # Total points within EBARs aeac (if needed)
    print("points within EBARs for model fitting")
    points_thin_aeac <- SpatialPoints(as(st_transform(points_thin_sf, aeac), 'Spatial'), proj4string=CRS(aeac))[as(Ecod_sp_join_id_diss, 'Spatial'), ]
    points_thin_aeac_c <- st_coordinates(st_as_sf(points_thin_aeac))
  # Saving points
    st_write(st_as_sf(points_thin_aeac), paste0("./Intermediate/Points_thinned_EBARs_aeac/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]), sp_eos$GLOBAL_SCIENTIFIC_NAME[i], driver="ESRI Shapefile", append=FALSE)  # create to a shapefile )
    
# 2.5 Predictors --------
  
# 2.5.1 Topographic heterogeneity (based on Digital Elevation Models, DEMs) ----------
  print("Topographic heterogeneity")
  
  # Vector_Ruggedness_Measure
    vrm_NA_aeac_ebar <- raster("./Data_inputs/vrm_1KMmn_GMTEDmd.tif") %>%
    crop(Ecod_sp_join_id_diss_s)
    names(vrm_NA_aeac_ebar) <- "Vector_Ruggedness_Measure"
    
  # Roughness
    roughness_NA_aeac_ebar <- raster("./Data_inputs/roughness_1KMmn_GMTEDmd.tif") %>%
    crop(Ecod_sp_join_id_diss_s)  %>%
    resample(vrm_NA_aeac_ebar, method = "ngb")
    names(roughness_NA_aeac_ebar) <- "Roughness"

  # Slope
    Slope_NA_aeac_ebar <- raster("./Data_inputs/slope_1KMmn_GMTEDmd.tif") %>%
      crop(Ecod_sp_join_id_diss_s)  %>%
      resample(vrm_NA_aeac_ebar, method = "ngb")
      names(Slope_NA_aeac_ebar) <- "Slope"
    
   # Eastness
      Eastness_NA_aeac_ebar <- raster("./Data_inputs/eastness_1KMmn_GMTEDmd.tif") %>%
        crop(Ecod_sp_join_id_diss_s)  %>%
        resample(vrm_NA_aeac_ebar, method = "ngb")
        names(Eastness_NA_aeac_ebar) <- "Eastness"
      
   # Northness
      Northness_NA_aeac_ebar <- raster("./Data_inputs/northness_1KMmn_GMTEDmd.tif") %>%
      crop(Ecod_sp_join_id_diss_s)  %>%
      resample(vrm_NA_aeac_ebar, method = "ngb")
      names(Northness_NA_aeac_ebar) <- "Northness"
    
# 2.5.2 Dynamic Habitat Index (fpar, 3 bands) ----------
  print("Dynamic Habitat Index (fpar)")
  # We need to reproject Ecodistricts to crop bands
    EcoDistCan_wgs84 <- st_transform( Ecod_sp_join_id_diss, wgs84)
  # Loading Dynamic Habitat Index  (all .tiff files)
    fpar_all <- list.files("./Data/fpar_can", pattern=".tif$", full.names = T)
  
    # Function to calculate fpar_mean faster
      rasterstack_meansd_fast <- function(x) {
        s0 <- nlayers(x)
        s1 <- raster(x, layer=1)
        for(ri in 2:s0) {
          r <- raster(x, layer=ri)
          s1 <- s1 + r
          }
        mean=s1/s0
      }
    
  #  Band 1, read rasters, stack them, calculate the mean of all years, projectRaster to aeac and res=1000m
    print("Band_1")
    # Loading rasters and calcualting mean
      fpar_raster_1s_mean <- lapply(fpar_all, raster, band=1) %>%
      stack()  %>%
      rasterstack_meansd_fast()
    # Cropinng, projecting
      Band_1_aeac_r <- crop(fpar_raster_1s_mean, as(EcoDistCan_wgs84, 'Spatial'))%>%
      projectRaster(res=1000, crs=aeac, method = "bilinear")%>%
      resample(roughness_NA_aeac_ebar, method = "ngb")
      names(Band_1_aeac_r) <- "Cummulative_annual_productivity_b1"
    
  #  Band 2, read rasters, stack them, calculate the mean of all years, project Raster to aeac and res=1000m
    print("Band_2")
    # Loading rasters and calcualting mean
      fpar_raster_2s_mean <- lapply(fpar_all, raster, band=2) %>%
      stack() %>%
      rasterstack_meansd_fast()
    # Croping, projecting
      Band_2_aeac_r <- crop(fpar_raster_2s_mean, as(EcoDistCan_wgs84, 'Spatial'))%>%
      projectRaster(res=1000, crs=aeac, method = "bilinear")%>%
      resample(roughness_NA_aeac_ebar, method = "ngb")
      names(Band_2_aeac_r) <- "Minimum_annual_productivity_b2"
    
  #  Band 3, read rasters, stack them, calcualte the mean of all years, projectRaster to aeac and res=1000m
    print("Band_3")
    # Loading rasters and calcualting mean
      fpar_raster_3s_mean <- lapply(fpar_all, raster, band=3) %>%
      stack()%>%
      rasterstack_meansd_fast()
    # Croping and projecting
      Band_3_aeac_r <- crop(fpar_raster_3s_mean, as(EcoDistCan_wgs84, 'Spatial'))%>%
      projectRaster(res=1000, crs=aeac, method = "bilinear")%>%
      resample(roughness_NA_aeac_ebar, method = "ngb")
      names(Band_3_aeac_r) <- "Variation_annual_productivity_b3"
  
# 2.5.4 Soils (WE DID NOT USED THIS DATASETS: LOTS OF NA values (especially in t the coastlines),
        # Observations will be removed when fall in these NAs values. The effect is particularly strong for
        # species with small number of observations 
      
      # from https://files.isric.org/soilgrids/latest/data/
      # various soils grids - spatial resolution 250m
      # libraries needed
       #  library(gdalUtils)
       # library(gdalUtilities)
      # For instance, the .vrt file contains the soils variable global data set.

     #soils_ocs <- gdalwarp(t_srs="EPSG:4326",
                           #multi=TRUE,
                           #wm=200, 
                           #co=c("BIGTIFF=YES", "COMPRESS=DEFLATE", "TILED=TRUE"),
                           #tr=c(0.01,0.01), # Desired output resolution
                           #"/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/ocs/ocs_0-30cm_mean.vrt", # Input VRT
                           #"./test1/ocs_0-30cm_mean_4.tif") # Output

     
# 2.5.5 Geology (WE DID NOT USED THIS DATASET: COARSE RESOLUTION)
      
   # geology <- raster("./Intermediate/Geology/geology_na_p.tif")%>%
   # resample(roughness_NA_aeac_ebar, method = "ngb")
   # names(geology) <- "Geology_Rocktype"
      
# 2.5.6 Lakes
     # We used the HydroLAkes dataset https://hydrosheds.org/page/hydrolakes
      # Lake polygons (including all attributes) in a Shapefile (782 MB zip-file).
      # to create two metrics
  
  # 2.5.6.1   Distance to lakes

    lakesNA_1000m <- raster("./Data_inputs/Lakes_1km2_01_bf250_aeac.tif")%>% # Lakes North America 1000m resolution
    crop(Ecod_sp_join_id_diss_s)
    lakesNA_1000m[lakesNA_1000m ==0]<-NA
    lakesNA_1000m_c_dist <- raster::distance(lakesNA_1000m, doEdge=TRUE)%>%
    resample(roughness_NA_aeac_ebar, method = "ngb")
    names(lakesNA_1000m_c_dist) <- "Distance_to_Lakes"

  # 2.5.6.2 Percentage of lakes
    
    lakes_NA_100m <- raster("./Intermediate/Lakes/Lakes_100m_01_bf250_aeac_NA.tif")%>% # Lakes North America 100m resolution
      crop(Ecod_sp_join_id_diss_s)
    lakes_NA_100m_agg <- aggregate(lakes_NA_100m, fact=10, fun=sum)%>% # aggregating 100m cells into 1000m cells, using 'sum' function
    resample(roughness_NA_aeac_ebar, method = "ngb")
    names(lakes_NA_100m_agg) <- "Lakes_percentage"


# 2.5.7 Climate from BIOCLIM ----------
    print("Bioclim")
    # Getting coordinates from Ecodistricts to select tiles from bioclim
      rdema_cc2 <- as.data.frame(st_coordinates(st_cast(EcoDistCan_wgs84)))
    EcoDistCan_wgs84_bioclim <-   as(EcoDistCan_wgs84, 'Spatial')
    # Creating directory
      if(dir.exists("./Intermediate/bioclim_tiles")){
        }else{
          out_sdm_mapsBIO <-   dir.create("./Intermediate/bioclim_tiles")}
    # Selecting tiles
      if(dir.exists(paste0("./Intermediate/bioclim_tiles", "/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))){
        }else{
          out_sdm_mapsBIO <-   dir.create(paste0("./Intermediate/bioclim_tiles","/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))}
    # tiles_NA <- tile_name(as(US_Canada_Mex_WGS84, 'Spatial'))
      tilex_min <- min(rdema_cc2$X)
      tilex_max <- max(rdema_cc2$X)
      tiley_min <- min(rdema_cc2$Y)
      tiley_max <- max(rdema_cc2$Y)
      devtools::install_github("kapitzas/WorldClimTiles")
      require("WorldClimTiles")
      wc_tiles <- tile_name(EcoDistCan_wgs84_bioclim, "worldclim") #for 0.5 arcmin worldclim tiles of 
      tiles <- tile_get(wc_tiles, name =  "worldclim", var="bio", path = paste0("./Intermediate/bioclim_tiles", "/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i])) #download and load worldclim
      tiles_merge <- tile_merge(tiles) #Reprojects data to 10 times smaller res, i.e. to get srtm data from 90m to approx 1km resolution.
      glo_all_c <- crop(tiles_merge, EcoDistCan_wgs84_bioclim)
    # Transforming it to aeac and 1000m resolution
      bioclim_ecod_p <- projectRaster(glo_all_c, crs=aeac, res = 1000, method = "bilinear")
      # Adding names
      bioclim_ecod_p_names <- bioclim_ecod_p
    names(bioclim_ecod_p_names)  <- c("BIO1_Annual_Mean_Temperature",
                                                    "BIO2_Mean_Diurnal_Range",
                                                    "BIO3_Isothermality",
                                                    "BIO4_Temperature_Seasonality",
                                                    "BIO5_Max_Temperature_Warmest_Month",
                                                    "BIO6_Min_Temperature_Coldest_Month",
                                                    "BIO7_Temperature_Annual_Range",
                                                    "BIO8_Mean_Temperature_Wettest_Quarter",
                                                    "BIO9_Mean_Temperature_Driest_Quarter",
                                                    "BIO10_Mean_Temperature_Warmest_Quarter",
                                                    "BIO11_Mean_Temperature_Coldest_Quarter",
                                                    "BIO12_Annual_Precipitation",
                                                    "BIO13_Precipitation_Wettest_Month",
                                                    "BIO14_Precipitation_Driest_Month",
                                                    "BIO15_Precipitation_Seasonality",
                                                    "BIO16_Precipitation_Wettest_Quarter",
                                                    "BIO17_Precipitation_Driest_Quarter",
                                                    "BIO18_Precipitation_Warmest_Quarter",
                                                    "BIO19_Precipitation_Coldest_Quarter")
    # Calculating correlation
    bioclim_ecod_p_names_all <- resample(bioclim_ecod_p_names, vrm_NA_aeac_ebar)
      rasters.crop.reduced <- removeCollinearity(bioclim_ecod_p_names,  multicollinearity.cutoff = 0.70, plot = F, select.variables = T, sample.points = FALSE)
    # Selecting non-colinear variables
      bioclim_ecod_p_noc <- subset(bioclim_ecod_p_names, rasters.crop.reduced)
    # Resampling with other predictors (same extent)
      bioclim_ecod_p_noc_r<- resample(bioclim_ecod_p_noc, vrm_NA_aeac_ebar)
      bioclim_ecod_p_noc_r_s <- stack(bioclim_ecod_p_noc_r)
      bioclim_ecod_p_noc_r_s_t <- projectRaster(bioclim_ecod_p_noc_r_s, crs = wgs84, res=0.008)
  
      
# 2.5.8  Stacking Predictors ----------
    print("Stacking predictors")
      Env_l <- stack(vrm_NA_aeac_ebar,
                     roughness_NA_aeac_ebar,
                     Slope_NA_aeac_ebar,
                     Eastness_NA_aeac_ebar,
                     Northness_NA_aeac_ebar,
                     lakesNA_1000m_c_dist,
                     lakes_NA_100m_agg,
                     if(maxValue(Band_1_aeac_r) != 0 && minValue(Band_1_aeac_r) !=0){Band_1_aeac_r}else{},
                     if(maxValue(Band_2_aeac_r) != 0 && minValue(Band_2_aeac_r) !=0){Band_2_aeac_r}else{},
                     if(maxValue(Band_3_aeac_r) != 0 && minValue(Band_3_aeac_r) !=0){Band_3_aeac_r}else{},
                     bioclim_ecod_p_names_all, na.rm=T)
      Env_lc <- mask(Env_l, calc(Env_l, fun=sum))
      Env_lc_df <- data.frame(cbind(getValues(Env_lc)))
      
    # Remove colinearity second time
      cor_clim_top_green <- removeCollinearity(Env_lc,  multicollinearity.cutoff = 0.70, plot = F, select.variables = T, sample.points = FALSE)
      cor_clim_top_green_s <- stack(subset(Env_lc, cor_clim_top_green))
      cor_clim_top_green_s_t <- projectRaster(cor_clim_top_green_s, crs = wgs84, res=0.008)
   
    #####
    # checking for intersection of predictors stack and observations
      points_predictors_overlap <- raster::extract(cor_clim_top_green_s, st_transform(points_thin_sf, st_crs(Ecod_sp_join_id_diss_s)), df=T)
      points_predictors_overlap_noNA <-  na.omit(points_predictors_overlap)
      num_occs_initial <- length(points_predictors_overlap_noNA[,1])
      
      points_thin_sf_c_sel <- as.data.frame(points_thin_sf_c)[as.vector(points_predictors_overlap_noNA$ID),]

    if(num_occs_initial <= 2){
      rmarkdown::render(input = "./Output_KBAs_v7_empty_v2.Rmd",
                        output_format = "html_document",
                        output_file =  paste0(sp_eos$GLOBAL_SCIENTIFIC_NAME[i], ".html"),
                        output_dir = "./HSMs_html_outputs")
    ###########
      print("Unable to run models and any further analyses due to occurrences <=2")
    #############
      
      #rm(list=ls())
      #removes entire temp directory without affecting other running processes
      unlink(file.path("./temp_SDMs_to_REMOVE"), recursive = TRUE)
      unlink(paste0("./Intermediate/bioclim_tiles/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]), recursive = T)
      rm(list=setdiff(ls(), listObjects))
      
      next
     
    }else{
    }
    
        
# 3. BACKGROUND POINTS ----------
    # Identify number of background points (bg) for HSMs within the geographic range (EBARs), testing three bg.
    
    # Selecting cell grids within geographic range (EBARs)
      georange_raster <- raster::mask(cor_clim_top_green_s[[1]], Ecod_sp_join_id_diss_s)
    
    # Calculating number of grid cells
      pixel_1000Km2<-(1000*1000)/1000000
      georange_raster_number2<-cellStats((((georange_raster*0)+1)*pixel_1000Km2), 'sum', na.rm=T)
      georange_raster_number2_times <- round(georange_raster_number2*1.2,0)
    
    # Selecting bg (two options: in case small geographic ranges)
      if(georange_raster_number2 < 100000){
      # Taking 20%, 50% and 60% (geographic range)
        background_points <- c(
          round(georange_raster_number2*0.20, 0),
          round(georange_raster_number2*0.50, 0),
          round(georange_raster_number2*0.60, 0)
          )
        }else{
          # Taking 10%, 20% and 40% (geographic range)
            background_points <- c(
              round(georange_raster_number2*0.05, 0),
              round(georange_raster_number2*0.10, 0),
              round(georange_raster_number2*0.20, 0)
              )
            }
   
    
# 4. SAMPLING BIAS LAYER ----------
    # Creating a sampling bias layer to select background point from these areas
     
     # UsIng MASS::kde2d
      points_with_data <- as.data.frame(raster::extract(cor_clim_top_green_s, coordinates(points_thin_aeac)))%>%
      drop_na()

          print("generating density map")
          dens.ras_t <- MASS::kde2d(st_coordinates(st_as_sf(st_transform(points_thin_sf, aeac)))[,1],
                                    st_coordinates(st_as_sf(st_transform(points_thin_sf, aeac)))[,2],
                                    lims = c(range(roughness_NA_aeac_ebar@extent@xmin, roughness_NA_aeac_ebar@extent@xmax),
                                           range(roughness_NA_aeac_ebar@extent@ymin, roughness_NA_aeac_ebar@extent@ymax))) %>%
            raster()
          proj4string(dens.ras_t) <- aeac
          dens.ras_t3 <- raster::disaggregate(dens.ras_t, fact=22, method="bilinear") %>%
            resample(roughness_NA_aeac_ebar)
       
          range_density_map <- mask(dens.ras_t3, Ecod_sp_join_id_diss_s)  
          range_density_map_wgs84 <- projectRaster(range_density_map, crs = crs(cor_clim_top_green_s_t[[1]]), res = 0.008)


# 5. MODEL  EVALUATION --------------
  
# 5.1 Selecting features based on number of occurrences ---------------
      if(nrow(points_thin_sf_c_sel) <=10){
      meth = 'jackknife'
      features <- "L"
      }else if(nrow(points_thin_sf_c_sel) >10 && nrow(points_thin_sf_c_sel) <=15){
        meth = 'jackknife'
        features <- c("L", "Q", "LQ")
        }else if(nrow(points_thin_sf_c_sel) >15 && nrow(points_thin_sf_c_sel) <=25){
          meth = 'jackknife'
          features <- c("L", "Q", "H", "LQ", "QH")
          }else if(nrow(points_thin_sf_c_sel) >25 && nrow(points_thin_sf_c_sel) <=80){
            meth = 'randomkfold'
            features <- c("L", "Q", "H", "LQ", "QH")
            }else if(nrow(points_thin_sf_c_sel) > 80){
              meth = 'randomkfold'
              features <- c("Q",  "H", "LQ", "LQP", "QPT")
              }
   
# 5.2 Running ENMeval
    
# Geographic range (EBArs) raster
    georange_raster_wgs84 <- projectRaster(georange_raster, crs = wgs84, res = 0.008)
    
# List objects to store ENMeval results
    e.mx_nobias <- list()
    e.mx_bias <- list()
    e.mx_nobias_clim_current <- list()
    e.mx_bias_clim_current <- list()
   
# 6. Model evaluation -------------

    ### LOOP STARTS ########
# Testing combinations of: multiple background points, bias/non_bias, and all predictors/only current climate.

  for(bgp in 1:length(background_points)){

    # Creating background points
        bgp_testing <- as.data.frame(randomPoints(georange_raster_wgs84, n=background_points[bgp]))
        colnames(bgp_testing)[1] <- "longitude"
        colnames(bgp_testing)[2] <- "latitude"

    # Occurrences
          occs<- as.data.frame(points_thin_sf_c_sel)
          colnames(occs)[1] <- "longitude"
          colnames(occs)[2] <- "latitude"
    # Removing occurrences that have the same coordinates is good practice to avoid pseudoreplication.
          occs <- occs[!duplicated(occs),]
          
    # Predictors
          envs <- stack(cor_clim_top_green_s_t)
          envs_clim <- stack(bioclim_ecod_p_noc_r_s_t)
          
    # Let's now remove occurrences that are cell duplicates -- these are
          occs.cells <- raster::extract(envs[[1]], occs, cellnumbers = TRUE)
          occs.cellDups <- duplicated(occs.cells[,1])
          occs <- occs[!occs.cellDups,]
          
        # First we extract the climatic variable values at the occurrence points -- these values are 
        # our "reference".

        # Let's make our occs into a sf object -- here we specify the coordinate reference system  
        # (crs) in order to perform spatial functions.
        occs.sf <- sf::st_as_sf(occs, coords = c("longitude","latitude"), 
                                crs = crs(cor_clim_top_green_s_t))
        # Make sure the RasterStack has the same coordinate reference system (CRS) string.
        # Both are the same CRS, but when the strings are different, some spatial functions error.
        crs(envs) <- raster::crs(occs.sf)
        
        occ_predictors_overlap <- raster::extract(envs, occs.sf, df=T)
        occ_predictors_overlap_noNA <-  na.omit(occ_predictors_overlap)
        num_occs_final <- length(occ_predictors_overlap_noNA[,1])
        
        occ_points_thin_sf_c_sel <- as.data.frame(occs)[as.vector(occ_predictors_overlap_noNA$ID),]
        
        occs <- occ_points_thin_sf_c_sel
        
        # Number of occurrences after removing points without environmental data
        num_occs <- length(occs[,1])
        
# 6.1 Evaluating No_bias     
  # Background points
        bg_nobias <-  randomPoints((georange_raster_wgs84*0), n=nrow(bgp_testing)) %>% as.data.frame()
        colnames(bg_nobias) <- colnames(occs)

  # Partition
        if(nrow(points_thin_sf_c_sel) <=25){
          user_partition  <- get.jackknife(occs, bg_nobias)
          
        }else if(nrow(points_thin_sf_c_sel) > 25){
          user_partition  <- get.randomkfold(occs, bg_nobias, 10)
          
        }

  # Running ENMeval

    # All predictors    
      cat(paste0("Evaluating Non_bias= ", background_points[bgp], "...All_predictors"), '\n')
      
      
          e.mx_nobias[[bgp]] <- ENMeval::ENMevaluate(occs = occs, envs = envs, bg = bg_nobias, 
                                algorithm = 'maxent.jar',
                                partitions = 'user',
                                user.grp = list(occs.grp = user_partition$occs.grp,  bg.grp = user_partition$bg.grp),
                                tune.args = list(fc = features, rm = c(0.05,0.5, 1)),
                                parallel =  TRUE,
                                doClamp = T,
                                updateProgress = TRUE,
                                parallelType = "doParallel"
                              )
  
    # Current climate
          cat(paste0("Evaluating Non_bias= ", background_points[bgp], "...Curret climate"), '\n')
          
          e.mx_nobias_clim_current[[bgp]] <- ENMeval::ENMevaluate(occs = occs, envs = envs_clim, bg = bg_nobias, 
                                                     algorithm = 'maxent.jar',
                                                     partitions = 'user',
                                                     user.grp = list(occs.grp = user_partition$occs.grp,  bg.grp = user_partition$bg.grp),
                                                     tune.args = list(fc = features, rm = c(0.05,0.5, 1)),
                                                     parallel =  TRUE,
                                                     doClamp = T,
                                                     updateProgress = TRUE,
                                                     parallelType = "doParallel"
                                                     
          )
          
          
# 6.2 Testing the effect of sampling bias

    # background using bias
        bg_bias <- xyFromCell(!is.na(range_density_map_wgs84),
                              sample(ncell(!is.na(range_density_map_wgs84)),
                                     nrow(bgp_testing),
                                     prob =  values(!is.na(range_density_map_wgs84))))
        colnames(bg_bias) <- colnames(occs)

    # Partition
        if(nrow(points_thin_sf_c_sel) <=25){
          user_partition2  <- get.jackknife(occs, bg_bias)
          
        }else if(nrow(points_thin_sf_c_sel) > 25){
          user_partition2  <- get.randomkfold(occs, bg_bias, 10)
          
        }
        
    # Running ENMevaluate 
        
      # All predictors
        cat(paste0("Evaluating bias= ", background_points[bgp], "...All predictors"), '\n')
        
        e.mx_bias[[bgp]]  <-ENMeval::ENMevaluate(occs = occs, envs = envs, bg = bg_bias, 
                                    algorithm = 'maxent.jar',
                                    partitions = 'user',
                                    user.grp = list(occs.grp = user_partition2$occs.grp,  bg.grp = user_partition2$bg.grp),
                                    tune.args = list(fc = features, rm = c(0.05,0.5, 1)),
                                    parallel =  TRUE,
                                    doClamp = T,
                                    updateProgress = TRUE,
                                    parallelType = "doParallel"
        )
        
        
    # Current climate
          cat(paste0("Evaluating bias= ", background_points[bgp], "...Current climate"), '\n')
          
          e.mx_bias_clim_current[[bgp]] <- ENMeval::ENMevaluate(occs = occs, envs = envs_clim, bg = bg_bias, 
                                                                  algorithm = 'maxent.jar',
                                                                  partitions = 'user',
                                                                  user.grp = list(occs.grp = user_partition2$occs.grp,  bg.grp = user_partition2$bg.grp),
                                                                  tune.args = list(fc = features, rm = c(0.05,0.5, 1)),
                                                                  parallel =  TRUE,
                                                                  doClamp = T,
                                                                  updateProgress = TRUE,
                                                                parallelType = "doParallel"
          )
          
       
     }

###### LOOP ENDS ########   

# 6.3 Overall results ALL PREDICTORS
    
  # No_Minimizing_bias
    res_nobias <- list()
    for (ev in 1:length(e.mx_nobias)){
      res_nobias[[ev]] <- eval.results(e.mx_nobias[[ev]])%>%
        filter(!is.na(or.10p.avg))%>%
        filter(or.10p.avg == min(or.10p.avg)) %>% 
        filter(auc.val.avg == max(auc.val.avg))
   
       
     res_nobias[[ev]]$bg_points <- length(e.mx_nobias[[ev]]@bg[,1])
     res_nobias[[ev]]$occs_points <- length(e.mx_nobias[[ev]]@occs.grp)
      
      res_nobias[[ev]]$bias_layer <- "No_Minimizing_bias"
      
    }
    res_nobias_d <- do.call(rbind, res_nobias)
    

  # Minimizing_bias
    res_bias <- list()
    for (ev in 1:length(e.mx_bias)){
      res_bias[[ev]] <- eval.results(e.mx_bias[[ev]])%>%
        filter(!is.na(or.10p.avg))%>%
        filter(or.10p.avg == min(or.10p.avg)) %>% 
        filter(auc.val.avg == max(auc.val.avg))
      
      res_bias[[ev]]$bg_points <- length(e.mx_bias[[ev]]@bg[,1])
      res_bias[[ev]]$occs_points <- length(e.mx_bias[[ev]]@occs.grp)
      
      res_bias[[ev]]$bias_layer <- "Minimizing_bias"
      
    }
    res_bias_d <- do.call(rbind, res_bias)
    
  # binding all models
    
    AUC_bestmod_full2 <- rbind(res_nobias_d, res_bias_d)
    
    best_e.mx_bias1 <- e.mx_bias[[1]]@models[e.mx_bias[[1]]@results$or.10p.avg == min(AUC_bestmod_full2$or.10p.avg)]
    best_e.mx_bias2 <- e.mx_bias[[2]]@models[e.mx_bias[[2]]@results$or.10p.avg == min(AUC_bestmod_full2$or.10p.avg)]
    best_e.mx_bias3 <- e.mx_bias[[3]]@models[e.mx_bias[[3]]@results$or.10p.avg == min(AUC_bestmod_full2$or.10p.avg)]
     
    best_e.mx_nobias1 <- e.mx_nobias[[1]]@models[e.mx_nobias[[1]]@results$or.10p.avg == min(AUC_bestmod_full2$or.10p.avg)]
    best_e.mx_nobias2 <- e.mx_nobias[[2]]@models[e.mx_nobias[[2]]@results$or.10p.avg == min(AUC_bestmod_full2$or.10p.avg)]
    best_e.mx_nobias3 <- e.mx_nobias[[3]]@models[e.mx_nobias[[3]]@results$or.10p.avg == min(AUC_bestmod_full2$or.10p.avg)]
    
    
  # Selecting best model
    best_model <- list(best_e.mx_nobias1, best_e.mx_nobias2, best_e.mx_nobias3,
                         best_e.mx_bias1, best_e.mx_bias2, best_e.mx_bias3)
    
      if(!is.empty(names(best_e.mx_bias1))){
        best_model <- e.mx_bias[[1]]}else{
          if(!is.empty(names(best_e.mx_bias2))){
            best_model <- e.mx_bias[[2]]}else{
              if(!is.empty(names(best_e.mx_bias3))){
                best_model <- e.mx_bias[[3]]}else{
                  if(!is.empty(names(best_e.mx_nobias1))){
                    best_model <- e.mx_nobias[[1]]}else{
                      if(!is.empty(names(best_e.mx_nobias2))){
                        best_model <- e.mx_nobias[[2]]}else{
                          if(!is.empty(names(best_e.mx_nobias3))){
                            best_model <- e.mx_nobias[[3]]
              }
            }
                      }}}}
                      

      best_model@results
      # Generate a rangeModelMetadata object based on the information stored in the 
      # ENMevaluate object.
      rmm <- eval.rmm(best_model)

  # Model selection
      
    # Sequential method that uses cross-validation results by selecting models
         # with the lowest average test omission rate, and to break ties,
         #with the highest average validation AUC
        
        res <- eval.results(best_model)%>%
        filter(!is.na(or.10p.avg))
        opt.seq <- res %>% 
          filter(or.10p.avg == min(or.10p.avg)) %>% 
          filter(auc.val.avg == max(auc.val.avg))
        
    
    # Let?s now choose the optimal model settings based on the sequential criteria and examine it.
        mod.seq <- eval.models(best_model)[[opt.seq$tune.args]]
        
    # Variable importance
      var_imp <- as.data.frame(best_model@variable.importance[[as.vector(opt.seq$tune.args[1])]])[,-2]
      
      # We can select the model predictions for our optimal model the same way we did for the 
      # model object above.
      pred.seq <- eval.predictions(best_model)[[opt.seq$tune.args]]
      pred.seq_aeac <- projectRaster(pred.seq, crs = aeac, res=1000, method="bilinear")%>%
        crop(bbox(Ecod_sp_join_id_diss_s))
      
      if(dir.exists(paste0("./Intermediate/hsm_map_best_all/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))){
      }else{
         dir.create(paste0("./Intermediate/hsm_map_best_all/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))}
      
      writeRaster(pred.seq_aeac, paste0("./Intermediate/hsm_map_best_all/",
                                        sp_eos$GLOBAL_SCIENTIFIC_NAME[i], "/",
                                        "best_model_all_",
                                        sp_eos$GLOBAL_SCIENTIFIC_NAME[i], ".tif"),
      format="GTiff", overwrite=T)


  # Calculating some stats
    
    pred.seq_aeac_ebar <- mask(pred.seq_aeac, Ecod_sp_join_id_diss_s)
    
    pred.seq_aeac_stats <- cellStats((((pred.seq_aeac_ebar*0)+1)*pixel_1000Km2), 'sum', na.rm=T)
    
    
    pred.seq_aeac_75 <- pred.seq_aeac_ebar
    pred.seq_aeac_75[pred.seq_aeac_75 < 0.75] <- NA
    pred.seq_aeac_75_stats <- cellStats((((pred.seq_aeac_75*0)+1)*pixel_1000Km2), 'sum', na.rm=T)
    
    pred.seq_aeac_50 <- pred.seq_aeac_ebar
    pred.seq_aeac_50[pred.seq_aeac_50 < 0.50] <- NA
    pred.seq_aeac_50[pred.seq_aeac_50 > 0.75] <- NA
    pred.seq_aeac_50_stats <- cellStats((((pred.seq_aeac_50*0)+1)*pixel_1000Km2), 'sum', na.rm=T)
    
    pred.seq_aeac_25 <- pred.seq_aeac_ebar
    pred.seq_aeac_25[pred.seq_aeac_25 < 0.25] <- NA
    pred.seq_aeac_25[pred.seq_aeac_25 > 0.50] <- NA
    pred.seq_aeac_25_stats <- cellStats((((pred.seq_aeac_25*0)+1)*pixel_1000Km2), 'sum', na.rm=T)
    
    pred.seq_aeac_0 <- pred.seq_aeac_ebar
    pred.seq_aeac_0[pred.seq_aeac_0 > 0.25] <- NA
    pred.seq_aeac_0_stats <- cellStats((((pred.seq_aeac_0*0)+1)*pixel_1000Km2), 'sum', na.rm=T)
    
    
    Area_Km2 <- c(pred.seq_aeac_0_stats,
                        pred.seq_aeac_25_stats,
                        pred.seq_aeac_50_stats,
                        pred.seq_aeac_75_stats
                        )

    bin <- c("0.00 - 0.25",
             "0.25 - 0.50",
             "0.50 - 0.75",
             "0.75 - 1.00")

    total_area <- as.data.frame(Area_Km2, bin)
    
  # best model parameters
     best_para <- AUC_bestmod_full2[AUC_bestmod_full2$or.10p.avg == min(AUC_bestmod_full2$or.10p.avg),]
    best_param <- best_para[best_para$auc.val.avg == max(best_para$auc.val.avg),]
    
  # Optimal model parameters
    opt_param <- opt.seq
    
  # Uncertainty map (running 10 time best model and calculate variance fro model predcitions)
   
    e.mx_best <- list()

    for(b in 1:10){
      
      # Creating background points
      bgp_testing_best <- as.data.frame(randomPoints(georange_raster_wgs84, n=best_param$bg_points))
      colnames(bgp_testing_best)[1] <- "longitude"
      colnames(bgp_testing_best)[2] <- "latitude"
      
      #class(best_param$bias_layer)
      
      if(best_param$bias_layer == "No_Minimizing_bias"){
        # Evaluating No_bias     
        # Background points
        bg_best <-  randomPoints((georange_raster_wgs84*0), n=nrow(bgp_testing_best)) %>% as.data.frame()
        colnames(bg_best) <- colnames(occs)
        
      }else{
        
        # background using bias
        bg_best <- xyFromCell(!is.na(range_density_map_wgs84),
                              sample(ncell(!is.na(range_density_map_wgs84)),
                                     nrow(bgp_testing_best),
                                     prob =  values(!is.na(range_density_map_wgs84))))
        colnames(bg_best) <- colnames(occs)
      }
      
      # Partition
      if(nrow(points_thin_sf_c_sel) <=25){
        user_partition_best  <- get.jackknife(occs, bg_best)
        
      }else if(nrow(points_thin_sf_c_sel) > 25){
        user_partition_best  <- get.randomkfold(occs, bg_best, 10)
        
      }
      
    cat(paste0("Running best model...", b), '\n')
      e.mx_best[[b]]  <-ENMeval::ENMevaluate(occs = occs, envs = envs, bg = bg_best, 
                                             algorithm = 'maxent.jar',
                                             partitions = 'user',
                                             user.grp = list(occs.grp = user_partition_best$occs.grp,  bg.grp = user_partition_best$bg.grp),
                                             tune.args = list(fc = best_param$fc, rm = as.numeric(best_param$rm)),
                                             parallel =  TRUE,
                                             doClamp = T,
                                             updateProgress = TRUE,
                                             parallelType = "doParallel"
    )
    }
    
    
    pred.best <- stack(lapply(e.mx_best, eval.predictions))
    uncertanty <- cv(pred.best)
    uncertanty_aeac <- projectRaster(uncertanty, crs = aeac, res=1000, method = 'bilinear')%>%
      crop(bbox(Ecod_sp_join_id_diss_s))
    
    uncertanty_aeac_stand <- uncertanty_aeac/maxValue(uncertanty_aeac)

    
    if(dir.exists(paste0("./Intermediate/hsm_map_uncertanty_all/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))){
    }else{
      dir.create(paste0("./Intermediate/hsm_map_uncertanty_all/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))}
    
    writeRaster(uncertanty_aeac_stand, paste0("./Intermediate/hsm_map_uncertanty_all/",
                                      sp_eos$GLOBAL_SCIENTIFIC_NAME[i], "/",
                                      "Uncertanty_all_",
                                      sp_eos$GLOBAL_SCIENTIFIC_NAME[i], ".tif"),
                format="GTiff", overwrite=T)
    

  # Null models

    # We first run the null simulations with 100 iterations to get a reasonable null distribution 
    # for comparisons with the empirical values
      mod.null <- ENMnulls(best_model,
                           mod.settings = list(fc = as.character(best_param$fc[1]),
                                               rm = as.numeric(as.vector(best_param$rm[1]))),
                           user.eval.type	 = "kspatial",
                           no.iter = 100,
                           parallel =  T,
                           removeMxTemp = T,
                           parallelType = "doParallel" 
      )
      
   
## 7. Kernell density analysis ----------------
    # Identifyin cluster of high suitable areas (>0.75)
      
  if(!is.na(minValue(pred.seq_aeac_75)) &&  !is.na(minValue(pred.seq_aeac_75))){
    
    pred.seq_aeac_75_points <- rasterToPoints(pred.seq_aeac_75, spatial=F)%>%
    SpatialPoints(proj4string=CRS(aeac))%>%
    as("ppp")
    
  # Define window of analysis (Ecodistricts)
    range <- Ecod_sp_join_id_diss_s
    range_diss <- gUnaryUnion(range)
    geo_range <- as(range_diss,"owin")
    Window(pred.seq_aeac_75_points) <- geo_range
    
  # Calculating Kernel for three radious (1 to 5 km)  
    K1 <- density(pred.seq_aeac_75_points, sigma=1000) # Using the default bandwidth
    K1_r<-raster(K1)
    
    K2 <- density(pred.seq_aeac_75_points, sigma=2000) # Using the default bandwidth
    K2_r<-raster(K2)
    
    K3 <- density(pred.seq_aeac_75_points, sigma=3000) # Using the default bandwidth
    K3_r<-raster(K3)
    
    K4 <- density(pred.seq_aeac_75_points, sigma=4000) # Using the default bandwidth
    K4_r<-raster(K4)
    
    K5 <- density(pred.seq_aeac_75_points, sigma=5000) # Using the default bandwidth
    K5_r<-raster(K5)
    
    
  # Creating polygons using Radius 5km and then mask HSA highest quantile
    print("Creating polygons using Radius 5km and then mask HSA highest quantile")
    projection(K5_r) <- CRS(aeac)
    
    K5_r_d <- resample(K5_r, pred.seq_aeac_75, method = "ngb")
    K5_r_q75<-K5_r_d
    K5_r_q75[K5_r_q75 < quantile(K5_r_d, 0.75)]<-NA
    K5_r_q75_1<- (K5_r_q75*0)+1
    
    
    K5_r_q75_pol <- rasterToPolygons(K5_r_q75_1)

    K5_r_q75_pol_s<- st_as_sf(K5_r_q75_pol)
    
    K5_r_q75_pol_s_diss <- st_union(K5_r_q75_pol_s)
    K5_r_q75_pol_s_diss_m<- st_cast(K5_r_q75_pol_s_diss, "POLYGON")
    K5_r_q75_pol_s_diss_m_s <- as( K5_r_q75_pol_s_diss_m, 'Spatial')
    
  # calcualting area
    K5_r_q75_AREA <- raster::extract(K5_r_q75_1, K5_r_q75_pol_s_diss_m_s, fun = sum)
        K5_r_q75_AREA_area <- as.vector(format(round(K5_r_q75_AREA*(1000*1000)/1000000, 0),big.mark=","))
    
  # calculating HSA for each polygon
    print("Stats Maxent model")
    
    K5_r_polygons <- st_cast(K5_r_q75_pol_s_diss_m)
    
  # Create directory to save polygons
    
  # Saving Table
    if(dir.exists(file.path(paste0("./Intermediate/sdm_polygons_KBAs", "/", namesp)))) {
    }else{
      
      dir.create (file.path(paste0("./Intermediate/sdm_polygons_KBAs", "/", namesp)))
    }
    
  # Creating objects to save list of stats
    total_area_sitesHSA <- list()
    master_table_to_plot <- list()
    
  # Looping to get stats for each polygon
    
    for(c in  1:length(K5_r_polygons)){
      polyt <- st_cast(K5_r_q75_pol_s_diss_m[c])
      
      st_write(polyt, paste0("./Intermediate/sdm_polygons_KBAs", "/", namesp), paste0("pol_", c), driver="ESRI Shapefile", delete_layer = TRUE)  # create to a shapefile 
      
  
  # Masking HSA
      polyt_r <- mask(pred.seq_aeac, st_as_sf(polyt))
      plot(pred.seq_aeac)
      plot(Ecod_sp_join_id_diss_s, add=T)
      plot(polyt, add=T)
      polyt_r_075 <- polyt_r
      polyt_r_075[polyt_r_075 <0.75]<-NA
      
      polyt_r_050 <- polyt_r
      polyt_r_050[polyt_r_050 < 0.50]<-NA
      polyt_r_050[polyt_r_050 > 0.75]<-NA
      
      polyt_r_025 <- polyt_r
      polyt_r_025[polyt_r_025 < 0.25]<-NA
      polyt_r_025[polyt_r_025 > 0.50]<-NA
      
      polyt_r_0 <- polyt_r
      polyt_r_0[polyt_r_0 > 0.25]<-NA
      
      
      polyt_r_25 <- polyt_r
      polyt_r_25[polyt_r_25 > 0.25]<-NA
      
    # Calculating areas
      area_polyt_r_0<-cellStats((((polyt_r*0)+1)*pixel_1000Km2), 'sum')
      area_polyt_r_25<-cellStats((((polyt_r_025*0)+1)*pixel_1000Km2), 'sum', na.rm=T)
      area_polyt_r_50<-cellStats((((polyt_r_050*0)+1)*pixel_1000Km2), 'sum', na.rm=T)
      area_polyt_r_75<-cellStats((((polyt_r_075*0)+1)*pixel_1000Km2), 'sum', na.rm=T)
      area_polyt_r_total<-cellStats((((polyt_r*0)+1)*pixel_1000Km2), 'sum', na.rm=T)
      total_polyt_r_area<- data.frame(bin=c("0 - 0.25", "0.26 - 0.50", "0.51 - 0.75", "0.76 - 1"), 
                                      AreaKm2=c(area_polyt_r_0,area_polyt_r_25,area_polyt_r_50,area_polyt_r_75))
    # Total area of HSA within geographic range

    # Proportion
      total_area$Area_proportion<- c((total_area$Area_Km2[1]/pred.seq_aeac_stats)*100,
                                      (total_area$Area_Km2[2]/pred.seq_aeac_stats)*100,
                                      (total_area$Area_Km2[3]/pred.seq_aeac_stats)*100,
                                      (total_area$Area_Km2[4]/pred.seq_aeac_stats)*100
      )  
      
    # Cummulative area
      total_area$cummulative_Percentage <- c(total_area$Area_proportion[4]+total_area$Area_proportion[3]+total_area$Area_proportion[2]+total_area$Area_proportion[1],
                                              total_area$Area_proportion[4]+total_area$Area_proportion[3]+total_area$Area_proportion[2],
                                              total_area$Area_proportion[4]+total_area$Area_proportion[3],
                                              total_area$Area_proportion[4] )
    # Total area
      total_area$range <- c("full_range", "full_range", "full_range", "full_range")
      
    # Total area of HSA within each polgon

      # Proportion
        total_polyt_r_area$Area_proportion<- c((total_polyt_r_area$AreaKm2[1]/area_polyt_r_total)*100,
                                              (total_polyt_r_area$AreaKm2[2]/area_polyt_r_total)*100,
                                              (total_polyt_r_area$AreaKm2[3]/area_polyt_r_total)*100,
                                              (total_polyt_r_area$AreaKm2[4]/area_polyt_r_total)*100)                    
      # Cummulative area
        total_polyt_r_area$cummulative_Percentage <- c(total_polyt_r_area$Area_proportion[4]+total_polyt_r_area$Area_proportion[3]+total_polyt_r_area$Area_proportion[2]+total_polyt_r_area$Area_proportion[1],
                                                      total_polyt_r_area$Area_proportion[4]+total_polyt_r_area$Area_proportion[3]+total_polyt_r_area$Area_proportion[2],
                                                      total_polyt_r_area$Area_proportion[4]+total_polyt_r_area$Area_proportion[3],
                                                      total_polyt_r_area$Area_proportion[4] )
      # Polygons' numbers
        total_polyt_r_area$range <- c(paste0("pol_", as.integer(c)), paste0("pol_", as.integer(c)), paste0("pol_", as.integer(c)), paste0("pol_", as.integer(c)))
      
      # Area within polygons
        total_area_sitesHSA[[c]] <- ((total_polyt_r_area$AreaKm2[1])+
                                     (total_polyt_r_area$AreaKm2[2])+
                                     (total_polyt_r_area$AreaKm2[3])+
                                     (total_polyt_r_area$AreaKm2[4] ))     
        total_area_sitesHSA1<-do.call("rbind", total_area_sitesHSA) # Summary
      
      # Master table for plotting
        master_table_to_plot[[c]] <- total_polyt_r_area
        master_table_to_plot1 <-do.call("rbind", master_table_to_plot)
      
      
    ## For plots (c. size)
      # Subset HSA
        high_suitable <- subset(master_table_to_plot1, bin  == "0.76 - 1")
      
      # Total HSA within the geographic range
        h =total_area$Area_Km2[4]
      
      # Percentage within range
        high_suitable$perc_range <- round(((high_suitable$AreaKm2/h)*100), 1)
      
      #add total hsa in the geographic range
        high_suitable$hsa_total <- total_area$AreaKm2[4]
        
      #Saving table
      write.csv(high_suitable, paste0("./Intermediate/polygons_area/", namesp, "_polygons", ".csv"))
      
      # Top three polygons with high HSA (**d. Top three potential sites)
        top_3_pol <-  high_suitable
        top_3_pol_s <-  top_n(top_3_pol,3, AreaKm2)%>%
          filter(AreaKm2 > 0)
      
        top_3_pol_s_area <- top_3_pol_s$AreaKm2
        top_3_pol_s_perc <- scales::percent(top_3_pol_s$AreaKm2/area_polyt_r_total)
        
        top_3_pol_s__range <- top_3_pol_s$range
        top_3_pol_s__range_short <- as.numeric(str_sub(top_3_pol_s__range, 5, 10))
        
      # Concatenate values
      
        high_suit_pol_selected <- c(high_suitable$AreaKm2[top_3_pol_s__range_short[1]]/total_area_sitesHSA1[top_3_pol_s__range_short[1]],high_suitable$AreaKm2[top_3_pol_s__range_short[2]]/total_area_sitesHSA1[top_3_pol_s__range_short[2]],high_suitable$AreaKm2[top_3_pol_s__range_short[3]]/total_area_sitesHSA1[top_3_pol_s__range_short[3]])
      
      # Saving top thre polygons
      
      st_write(K5_r_q75_pol_s_diss_m[top_3_pol_s__range_short], paste0("./Intermediate/sdm_polygons_KBAs_Top3", "/", namesp), paste0("pol_"), driver="ESRI Shapefile", delete_layer = TRUE)  # create to a shapefile 
     
    }
  }else{
      
    }
    
    
  # summary table
    
    summary_table <- c(National_name,
                       namesp,
                       as.character(best_param$fc[1]),
                       as.numeric(as.vector(best_param$rm[1])),
                       round(as.numeric(as.vector(best_param$or.10p.avg[1])),2),
                       round(as.numeric(as.vector(best_param$auc.val.avg[1])),2)
    )
         
    names(summary_table) <- c("national_name", "Scientific_name", "features", "rm", "or.10p.avg", "auc.val.avg")
    
#####
    
    
# 8. Overall results CLIMATE
    
  # No_Minimizing_bias
    res_nobias_clim <- list()
    
    for (ev in 1:length(e.mx_nobias_clim_current)){
      res_nobias_clim[[ev]] <- eval.results(e.mx_nobias_clim_current[[ev]])%>%
        
        filter(or.10p.avg == min(or.10p.avg)) %>% 
        filter(auc.val.avg == max(auc.val.avg))
      
      
      res_nobias_clim[[ev]]$bg_points <- length(e.mx_nobias_clim_current[[ev]]@bg[,1])
      res_nobias_clim[[ev]]$occs_points <- length(e.mx_nobias_clim_current[[ev]]@occs.grp)
      
      res_nobias_clim[[ev]]$bias_layer <- "No_Minimizing_bias"
      
    }
    
    res_nobias_d_clim <- do.call(rbind, res_nobias_clim)
    
  # Minimizing_bias
    res_bias_clim <- list()
    for (ev in 1:length(e.mx_bias_clim_current)){
      res_bias_clim[[ev]] <- eval.results(e.mx_bias_clim_current[[ev]])%>%
        filter(or.10p.avg == min(or.10p.avg)) %>% 
        filter(auc.val.avg == max(auc.val.avg))
      
      res_bias_clim[[ev]]$bg_points <- length(e.mx_bias_clim_current[[ev]]@bg[,1])
      res_bias_clim[[ev]]$occs_points <- length(e.mx_bias_clim_current[[ev]]@occs.grp)
      
      res_bias_clim[[ev]]$bias_layer <- "Minimizing_bias"
      
    }
    res_bias_d_clim <- do.call(rbind, res_bias_clim)
    
  # binding all models
    
    AUC_bestmod_full2_clim <- rbind(res_nobias_d_clim, res_bias_d_clim)
    
    best_e.mx_bias1_clim <- e.mx_bias_clim_current[[1]]@models[e.mx_bias_clim_current[[1]]@results$or.10p.avg == min(AUC_bestmod_full2_clim$or.10p.avg)]
    best_e.mx_bias2_clim <- e.mx_bias_clim_current[[2]]@models[e.mx_bias_clim_current[[2]]@results$or.10p.avg == min(AUC_bestmod_full2_clim$or.10p.avg)]
    best_e.mx_bias3_clim <- e.mx_bias_clim_current[[3]]@models[e.mx_bias_clim_current[[3]]@results$or.10p.avg == min(AUC_bestmod_full2_clim$or.10p.avg)]
    
    best_e.mx_nobias1_clim <- e.mx_nobias_clim_current[[1]]@models[e.mx_nobias_clim_current[[1]]@results$or.10p.avg == min(AUC_bestmod_full2_clim$or.10p.avg)]
    best_e.mx_nobias2_clim <- e.mx_nobias_clim_current[[2]]@models[e.mx_nobias_clim_current[[2]]@results$or.10p.avg == min(AUC_bestmod_full2_clim$or.10p.avg)]
    best_e.mx_nobias3_clim <- e.mx_nobias_clim_current[[3]]@models[e.mx_nobias_clim_current[[3]]@results$or.10p.avg == min(AUC_bestmod_full2_clim$or.10p.avg)]
    
    
  # Selecting best model
    best_model_clim <- list(best_e.mx_nobias1_clim, best_e.mx_nobias2_clim, best_e.mx_nobias3_clim,
                       best_e.mx_bias1_clim, best_e.mx_bias2_clim, best_e.mx_bias3_clim)
    
    if(!is.empty(names(best_e.mx_bias1_clim))){
      best_model_clim <- e.mx_bias_clim_current[[1]]}else{
        if(!is.empty(names(best_e.mx_bias2_clim))){
          best_model_clim <- e.mx_bias_clim_current[[2]]}else{
            if(!is.empty(names(best_e.mx_bias3_clim))){
              best_model_clim <- e.mx_bias_clim_current[[3]]}else{
                if(!is.empty(names(best_e.mx_nobias1_clim))){
                  best_model_clim <- e.mx_nobias_clim_current[[1]]}else{
                    if(!is.empty(names(best_e.mx_nobias2_clim))){
                      best_model_clim <- e.mx_nobias_clim_current[[2]]}else{
                        if(!is.empty(names(best_e.mx_nobias3_clim))){
                          best_model_clim <- e.mx_nobias_clim_current[[3]]
                        }
                      }
                  }}}}
    
    
    
      num_occs_clim <-length(best_model_clim@occs.grp)    
    
    
  # Generate a rangeModelMetadata object based on the information stored in the 
  # ENMevaluate object.
      rmm_clim <- eval.rmm(best_model_clim)
    
  # Model selection
    
    # Sequential method that uses cross-validation results by selecting models
    # with the lowest average test omission rate, and to break ties,
    #with the highest average validation AUC
    
      res_clim <- eval.results(best_model_clim)
      opt.seq_clim <- res %>% 
        filter(or.10p.avg == min(or.10p.avg)) %>% 
        filter(auc.val.avg == max(auc.val.avg))
      opt.seq_clim
    
  # Lets now choose the optimal model settings based on the sequential criteria and examine it.
    mod.seq_clim <- eval.models(best_model_clim)[[opt.seq_clim$tune.args]]

  # Variable importance
    var_imp_clim <- as.data.frame(best_model_clim@variable.importance[[as.vector(opt.seq_clim$tune.args[1])]])[,-2]
    
  # We can select the model predictions for our optimal model the same way we did for the 
  # model object above.
    pred.seq_clim <- eval.predictions(best_model_clim)[[opt.seq_clim$tune.args]]
    pred.seq_aeac_clim <- projectRaster(pred.seq_clim, crs = aeac, res=1000, method="bilinear")%>%
      crop(bbox(Ecod_sp_join_id_diss_s))
    

    if(dir.exists(paste0("./Intermediate/hsm_map_best_clim/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))){
    }else{
      dir.create(paste0("./Intermediate/hsm_map_best_clim/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))}
    
    writeRaster(pred.seq_aeac_clim, paste0("./Intermediate/hsm_map_best_clim/",
                                      sp_eos$GLOBAL_SCIENTIFIC_NAME[i], "/",
                                      "best_model_clim_",
                                      sp_eos$GLOBAL_SCIENTIFIC_NAME[i], ".tif"),
                format="GTiff", overwrite=T)
 
  # best model parameters
    best_param_cli <- AUC_bestmod_full2_clim[AUC_bestmod_full2_clim$or.10p.avg == min(AUC_bestmod_full2_clim$or.10p.avg), ]
    best_param_clim <- best_param_cli[best_param_cli$auc.val.avg == max(best_param_cli$auc.val.avg), ]
 
    # Uncertainty map
    
      e.mx_best_clim <- list()
      
      for(b in 1:10){
        
        # Creating background points
        bgp_testing_best_clim <- as.data.frame(randomPoints(georange_raster_wgs84, n=best_param_clim$bg_points))
        colnames(bgp_testing_best_clim)[1] <- "longitude"
        colnames(bgp_testing_best_clim)[2] <- "latitude"
        
  
        if(best_param_clim$bias_layer == "No_Minimizing_bias"){
          # Evaluating No_bias     
          # Background points
          bg_best_clim <-  randomPoints((georange_raster_wgs84*0), n=nrow(bgp_testing_best_clim)) %>% as.data.frame()
          colnames(bg_best_clim) <- colnames(occs)
          
        }else{
          
          # background using bias
          bg_best_clim <- xyFromCell(!is.na(range_density_map_wgs84),
                                sample(ncell(!is.na(range_density_map_wgs84)),
                                       nrow(bgp_testing_best_clim),
                                       prob =  values(!is.na(range_density_map_wgs84))))
          colnames(bg_best_clim) <- colnames(occs)
        }
        
        # Partition
        if(nrow(points_thin_sf_c_sel) <=25){
          user_partition_best_clim  <- get.jackknife(occs, bg_best_clim)
          
        }else if(nrow(points_thin_sf_c_sel) > 25){
          user_partition_best_clim  <- get.randomkfold(occs, bg_best_clim, 10)
          
        }
        
      cat(paste0("Running best model CLIM...", b), '\n')
        e.mx_best_clim[[b]]  <-ENMeval::ENMevaluate(occs = occs, envs = envs_clim, bg = bg_best_clim, 
                                               algorithm = 'maxent.jar',
                                               partitions = 'user',
                                               user.grp = list(occs.grp = user_partition_best_clim$occs.grp,  bg.grp = user_partition_best_clim$bg.grp),
                                               tune.args = list(fc = best_param_clim$fc, rm = as.numeric(best_param_clim$rm)),
                                               parallel =  TRUE,
                                               doClamp = T,
                                               updateProgress = TRUE
        )
      }
      
    
      pred.best_clim <- stack(lapply(e.mx_best_clim, eval.predictions))
      uncertanty_clim <- cv(pred.best_clim)
      uncertanty_aeac_clim <- projectRaster(uncertanty_clim, crs = aeac, res=1000, method = 'bilinear')%>%
      crop(bbox(Ecod_sp_join_id_diss_s))
      
      uncertanty_aeac_clim_stand <- uncertanty_aeac_clim/maxValue(uncertanty_aeac_clim)
      
   
    
    
    if(dir.exists(paste0("./Intermediate/hsm_map_uncertanty_clim/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))){
    }else{
      dir.create(paste0("./Intermediate/hsm_map_uncertanty_clim/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))}
    
    writeRaster(uncertanty_aeac_clim_stand, paste0("./Intermediate/hsm_map_uncertanty_clim/",
                                              sp_eos$GLOBAL_SCIENTIFIC_NAME[i], "/",
                                              "Uncertanty_clim_",
                                              sp_eos$GLOBAL_SCIENTIFIC_NAME[i], ".tif"),
                format="GTiff", overwrite=T)
    
    
    
    
# 9. Projected climate  ----------
    print("Running HSM: Maxent model - Projected climate")
    climate_scenarios <- read_xlsx("./Data/Climate_scenarios.xlsx")
    
  # Loading files (two scenarios rcp45 and rcp85)
    future_bio <- list.dirs("./Data_inputs/bioclim_future_mean", recursive=F)

  # Objects to save lists       
    CMIP5_predict <- list()
    CMIP5_predict_clamped <-list()
    hsm_change <- list() 
    hsm_change_clamped <- list()
    future_bio_name <- list()
    
    for (f in 1:length(future_bio)){
      cat(paste0("model_projected_clim", f), '\n')
      
      # Loading files
        CMIP5_file <- list.files(future_bio[f], pattern='.tif$', full.names=TRUE)
      
      # Stacking them
        CMIP5_files_r_stack <- stack(CMIP5_file)
      
      # fiel short name
        future_bio_name[[f]] <- str_sub(future_bio[f], -5, -1)
      
      # Selecting tiles
        # We can download directly from the URL (https://worldclim.org/CMIP5_30s), but no connection using R.
        #globalbioclim_projected <- raster::getData('CMIP5', var='bio', res=0.5, rcp=85, model='CC', year=50, download = TRUE, path="C:/KBAsCan/Data/rcp85_bioclim_tiles_projectionsrcp85")
        #globalbioclim_projected <- raster::getData('CMIP5', var='bio', res=0.5, rcp=60, model='CC', year=50, download = TRUE, path="C:/KBAsCan/Data/rcp60_bioclim_tiles_projections")
        
      #list rasters
        #CMIP5_files <- list.files("C:/KBAsCan/Data/rcp85_bioclim_tiles_projections", pattern='tif$', full.names=TRUE )
        #stack them
        #CMIP5_files_stack <- stack(CMIP5_files)
        
      # Cropping it
        bioclim_ecod2 <- crop(CMIP5_files_r_stack, as(EcoDistCan_wgs84, 'Spatial'))
      
      # Transforming it to aeac and 1000m resolution
        bioclim_ecod_p2 <- projectRaster(bioclim_ecod2, crs=aeac, res = 1000)
      
      # Adding names
        names(bioclim_ecod_p2)  <- c("BIO1_Annual_Mean_Temperature",
                                     "BIO10_Mean_Temperature_Warmest_Quarter",
                                     "BIO11_Mean_Temperature_Coldest_Quarter",
                                     "BIO12_Annual_Precipitation",
                                     "BIO13_Precipitation_Wettest_Month",
                                     "BIO14_Precipitation_Driest_Month",
                                     "BIO15_Precipitation_Seasonality",
                                     "BIO16_Precipitation_Wettest_Quarter",
                                     "BIO17_Precipitation_Driest_Quarter",
                                     "BIO18_Precipitation_Warmest_Quarter",
                                     "BIO19_Precipitation_Coldest_Quarter",
                                     "BIO2_Mean_Diurnal_Range",
                                     "BIO3_Isothermality",
                                     "BIO4_Temperature_Seasonality",
                                     "BIO5_Max_Temperature_Warmest_Month",
                                     "BIO6_Min_Temperature_Coldest_Month",
                                     "BIO7_Temperature_Annual_Range",
                                     "BIO8_Mean_Temperature_Wettest_Quarter",
                                     "BIO9_Mean_Temperature_Driest_Quarter"
                                     )
      
      
      # Subsetting
        bioclim_ecod_p2_drop <- subset(bioclim_ecod_p2, c(names(bioclim_ecod_p_noc_r)))
      
      # Resample (to maintain same extent) 
        bioclim_ecod_p_noc2_r <- resample(bioclim_ecod_p2_drop, Band_3_aeac_r)
      
      # staking them
      #   Env_lc_future <- stack(vrm_NA_aeac_ebar, roughness_NA_aeac_ebar, Band_1_aeac_r, Band_2_aeac_r, Band_3_aeac_r,  bioclim_ecod_p_noc2_r)
      
      # Predicting HSA with projcted climate data
        CMIP5_predict[[f]] = dismo::predict(mod.seq_clim, stack(bioclim_ecod_p_noc2_r))%>%
        resample(pred.seq_aeac_clim)%>%
          crop(bbox(Ecod_sp_join_id_diss_s))
        
        CMIP5_predict_clamped[[f]] = dismo::predict(mod.seq_clim, bioclim_ecod_p_noc2_r, args=c('doclamp=true', 'fadebyclamping=false'))%>%
          resample(pred.seq_aeac_clim)%>%
          crop(bbox(Ecod_sp_join_id_diss_s))

      # Calculating change in HSA
        hsm_change[[f]] <- CMIP5_predict[[f]] - pred.seq_aeac_clim%>%
        resample(pred.seq_aeac_clim)%>%
          crop(bbox(Ecod_sp_join_id_diss_s))
        
        hsm_change_clamped[[f]] <- CMIP5_predict_clamped[[f]] - pred.seq_aeac_clim%>%
          resample(pred.seq_aeac_clim)%>%
          crop(bbox(Ecod_sp_join_id_diss_s))
      

      if(dir.exists(paste0("./Intermediate/hsm_map_projected_clim/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))){
      }else{
        dir.create(paste0("./Intermediate/hsm_map_projected_clim/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))}
      
        writeRaster(stack(CMIP5_predict[f]), by.layer=T,  paste0("./Intermediate/hsm_map_projected_clim/",
                                                       sp_eos$GLOBAL_SCIENTIFIC_NAME[i], "/",
                                                       "projected_clim_",
                                                       sp_eos$GLOBAL_SCIENTIFIC_NAME[i], "_", future_bio_name[f], ".tif"),
                    format="GTiff", overwrite=T)
        
      
      
      
        
        if(dir.exists(paste0("./Intermediate/hsm_map_projected_clim_change/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))){
        }else{
          dir.create(paste0("./Intermediate/hsm_map_projected_clim_change/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]))}
        
        writeRaster(stack(hsm_change[f]), by.layer=T,  paste0("./Intermediate/hsm_map_projected_clim_change/",
                                                                 sp_eos$GLOBAL_SCIENTIFIC_NAME[i], "/",
                                                                 "hsm_projected_clim_change_",
                                                                 sp_eos$GLOBAL_SCIENTIFIC_NAME[i], "_", future_bio_name[f], ".tif"),
                    format="GTiff", overwrite=T)
        
      
    }
    
  
    #####
    
# 10. Protected areas -------
    PA_ecod <- st_crop(Protected_Area,Ecod_sp_join_id_diss_s_wgs84)

# 11. Cities ---------
      cities_ca_t_w <- st_crop(cities_ca_t,  Ecod_sp_join_id_diss_s_wgs84)
    
# 12. Species-KBAs delineated by experts
    
    # bottom_up_sites <- list.files(path =  paste0(getwd(), "/Data/KBA_bottom_up/"), pattern = ".shp$", full.names=T)
    
    
    if(file.exists(paste0("./Data_inputs/KBA_bottom_up/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i], ".shp"))){
      
      
      kba_site_bottom_up <- st_read(paste0("./Data_inputs/KBA_bottom_up/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i], ".shp"))
      kba_site_bottom_up_t <- st_transform(kba_site_bottom_up, CRS("+proj=longlat +datum=WGS84"))
      
    }else{
      
    }
    

#13. Summary performance
    
    Number_observations <- as.data.frame(c(paste0("> ", 10*length(names(cor_clim_top_green_s))),
                             paste0("Between 15 ", "and ", 10*length(names(cor_clim_top_green_s))),
                             "<= 14"))
    
   
    Omission_rate <- as.data.frame(c("<= 0.25",
                       "Between 0.25-0.50",
                       "> 0.50"))
    
  
    AUC_values <- as.data.frame(c(">= 0.75",
                    "< 0.75", "-"))
    
  
    Summary_performance <- cbind(Number_observations,
                                 Omission_rate,
                                 AUC_values)
      
    
    names(Summary_performance) <- c("Number_observations",
                                    "Omission_rate",
                                    "AUC_values")
    
    summary_table <- c(National_name, namesp, "14", num_occs, length(names(cor_clim_top_green_s)), 10*length(names(cor_clim_top_green_s)), round(best_param$or.10p.avg,2), round(best_param$auc.val.avg, 2))
    names_table <- c("NATL_ENGL_NAME", "SCIENTIFIC", "minimum small range", "Observations", "Number of predictors", "minimum rule of 10", "Omission rate", "AUC")
    summary_table_f <- data.frame(rbind(names_table, summary_table))
  
    write.csv(summary_table_f,
              paste0("./Intermediate/hsm_summary_performance/",  sp_eos$GLOBAL_SCIENTIFIC_NAME[i], ".csv"))
    

#14.  Creating HTML files

      
      rmarkdown::render(input = "./Output_KBAs_v7.Rmd",
                                      output_format = "html_document",
                                      output_file =  paste0(sp_eos$GLOBAL_SCIENTIFIC_NAME[i], ".html"),
                                      output_dir = "./HSMs_html_outputs")
      #removes entire temp directory without affecting other running processes
      unlink(file.path("./temp_SDMs_to_REMOVE"), recursive = TRUE)
      unlink(paste0("./Intermediate/bioclim_tiles/", sp_eos$GLOBAL_SCIENTIFIC_NAME[i]), recursive = T)
      rm(list=setdiff(ls(), listObjects))
      
     
}

      
########## End of the script  ##########

     
    #######################################
  
    
