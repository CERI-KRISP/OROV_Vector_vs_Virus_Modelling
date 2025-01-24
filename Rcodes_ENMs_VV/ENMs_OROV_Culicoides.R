# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Paper: Mapping Vector Versus Viral Suitability for Oropouche Transmission in the Americas: 
#         Evaluating Environmental Drivers Using Random Forest

# Author: Jenicca Poongavanan, Graeme D’or, Marcel D. Dunaiski, Moritz Kraemer, Cheryl Baxter, Marta Giovanetti,
#         Tulio de Oliveira, Houriiyah Tegally
# Date: January 2025
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
packages <- c("corrr", "moments", "maps", "geosphere", "tidyverse", "dplyr", "ggplot2", 
              "rworldmap", "data.table", "ggthemes", "airportr", "hrbrthemes", "lubridate", 
              "paletteer", "RColorBrewer", "readr", "MMWRweek", "scales", "showtext", 
              "sysfonts", "sf", "rnaturalearth", "janitor", "raster", "viridis", 
              "ggnewscale", "ggrepel", "terra", "party", "spatstat", "dismo", 
              "rgbif", "biomod2", "gridExtra", "knitr", "ade4", "cleangeo", "rasterVis")

# Load all packages
for(i in packages) { library(i, character.only=T) }


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#  \ Modelling Oropouche Virus (OROV) /  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Load required data 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### OROV Data
# Make sure you specify the data's location
## Presence Data
Presence_only_virus_brazil <- fread("Data/../Orov_presence_cleaned.csv")
# Predictor variables
Env.virus <- rast("Data/../Envi_37layers.tif")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Random Sampling of Pseudo Absences 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SPC_PresAbs_Random_virus <-  BIOMOD_FormatingData(resp.var = Presence_only_virus_brazil$occurrence ,
                                                  expl.var = Env.virus,
                                                  resp.xy = Presence_only_virus_brazil[,c('Longitude', 'Latitude')],
                                                  resp.name = "Culicoides",
                                                  PA.nb.rep = 3,
                                                  PA.nb.absences = 10000,
                                                  PA.strategy = 'random',
                                                  na.rm=T)

SPC_PresAbs_Random_virus
plot(SPC_PresAbs_Random_virus) ## Quick integrated plot 

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Running Biomod Models
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MySpc_models_random_virus <- BIOMOD_Modeling (bm.format = SPC_PresAbs_Random_virus,
                                              models = c("RF"),
                                              OPT.user = MySpc_options_2_virus,
                                              #OPT.strategy	= 'bigboss',
                                              CV.strategy = 'block',  
                                              #var.import = 5,
                                              metric.eval =c('TSS','ROC'),
                                              #CV.do.full.models = FALSE,
                                              modeling.id = "V1")

MySpc_models_random_virus


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Geographic of Pseudo Absence Sampling of Pseudo Absences 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SPC_PresAbs_Geographic_virus <-  BIOMOD_FormatingData(resp.var = Presence_only_virus_brazil$occurrence ,
                                                      expl.var = Env.virus,
                                                      resp.xy = Presence_only_virus_brazil[,c('Longitude', 'Latitude')],
                                                      resp.name = "Culicoides",
                                                      PA.nb.rep = 3,
                                                      PA.nb.absences = 450,
                                                      PA.strategy = 'disk',
                                                      PA.dist.min = 50000,
                                                      PA.dist.max = 500000,
                                                      na.rm=T)

SPC_PresAbs_Geographic_virus
plot(SPC_PresAbs_Geographic_virus) ## Quick integrated plot 

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BIOMOD MODELLING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MySpc_models_geographic_virus <- BIOMOD_Modeling (bm.format = SPC_PresAbs_Geographic_virus,
                                                  models = c("RF"),
                                                  OPT.user = MySpc_options_2_virus,
                                                  #OPT.strategy	= 'bigboss',
                                                  CV.strategy = 'block',  
                                                  #CV.perc = 0.7,
                                                  #var.import = 5,
                                                  metric.eval =c('TSS','ROC'),
                                                  #CV.do.full.models = FALSE,
                                                  modeling.id = "V1")


### get models evaluation scores
MyModels_scores_geographic_virus <- get_evaluations(MySpc_models_geographic_virus)
MyModels_scores_geographic_virus


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Pseudo - Absence sampling based on Sampling rate of OROV (only) and with population density 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Convert your presence points to an `ppp` object
presence_points_ppp_virus <- ppp(PP_virus_nodups$Longitude, PP_virus_nodups$Latitude, 
                                 window = owin(c(min(PP_virus_nodups$Longitude), max(PP_virus_nodups$Longitude)), 
                                               c(min( PP_virus_nodups$Latitude), max( PP_virus_nodups$Latitude))))

# Perform kernel density estimation
kde_virus <- density.ppp(presence_points_ppp_virus, sigma = 3)
# Assuming population_density is a raster layer
kde_virus_raster <- raster(kde_virus, crs=crs(pop_density_raster))
extent(kde_virus_raster) <- extent(pop_density_raster)
kde_virus_raster_masked <- mask(kde_virus_raster, brazil.shp)

## only run these lines if you are including population density 
#pop_density_raster2 <- resample(pop_density_raster,kde_raster_masked)
#combined_density <- overlay(pop_density_raster2, kde_raster_masked, fun=function(x, y) {x * y})

## If you dont want to include population density as a variable into the pseudo absence sampling 
combined_density <- kde_virus_raster_masked
### Raster-Based Approach 
# Use a distance raster approach (optional step for sampling purposes)
distance_raster <- distanceFromPoints(kde_raster_masked, presence_points_sf)

# Create a mask for distances between xx and xx  ## Change the distance as required: 1000 = 1 km 
mask_within_range <- (distance_raster > 50000) & (distance_raster < 300000)
plot(mask_within_range)

#mask_within_range <- resample(mask_within_range, combined_density, method="bilinear")
areaToSample <- combined_density*mask_within_range
plot(areaToSample)


# Number of replicates
n_replicates <- 3

# List to store pseudo-absence replicates
pseudo_absence_replicates <- list()

for (i in 1:n_replicates) {
  # Step 1: Perform kernel density estimation (already defined)
  kde_values <- getValues(areaToSample)
  kde_values[is.na(kde_values)] <- 0
  kde_probabilities <- kde_values / sum(kde_values)
  
  # Step 2: Sample pseudo-absence points
  n_points <- 450
  valid_cells <- which(kde_probabilities > 0)
  sampled_cells <- sample(valid_cells, size = n_points, prob = kde_probabilities[valid_cells], replace = TRUE)
  sampled_coords <- xyFromCell(areaToSample, sampled_cells)
  
  # Step 3: Store the sampled points in a data frame
  pseudo_absences <- data.frame(
    Longitude = sampled_coords[, 1],
    Latitude = sampled_coords[, 2],
    occurrence = 0
  )
  
  pseudo_absence_replicates[[i]] <- pseudo_absences
}

PA1 <- pseudo_absence_replicates[[1]]
PA2 <- pseudo_absence_replicates[[2]]
PA3 <- pseudo_absence_replicates[[3]]

# Combine all pseudo-absence replicates with the original presence points
Presence_only_virus_brazil_PA <- rbind(
  Presence_only_virus_brazil[, c(2:4)],  # Presence points
  PA1,  # Pseudo-absence replicate 1
  PA2,  # Pseudo-absence replicate 2
  PA3   # Pseudo-absence replicate 3
)

# Create the PA.user.table
# Total rows = all presence + all pseudo-absence points
n_rows <- nrow(Presence_only_virus_brazil_PA)

# Create a matrix with zeros
PA.table <- data.frame(
  PA1 = rep(FALSE, nrow(Presence_only_virus_brazil_PA)),
  PA2 = rep(FALSE, nrow(Presence_only_virus_brazil_PA)),
  PA3 = rep(FALSE, nrow(Presence_only_virus_brazil_PA))
)

# Assign presence points as TRUE for all replicates
PA.table[1:nrow(Presence_only_virus_brazil), ] <- TRUE


# Assign pseudo-absence points for each replicate
PA.table[(nrow(Presence_only_virus_brazil) + 1):(nrow(Presence_only_virus_brazil) + nrow(PA1)), 1] <- TRUE
PA.table[(nrow(Presence_only_virus_brazil) + nrow(PA1) + 1):(nrow(Presence_only_virus_brazil) + nrow(PA1) + nrow(PA2)), 2] <- TRUE
PA.table[(nrow(Presence_only_virus_brazil) + nrow(PA1) + nrow(PA2) + 1):(nrow(Presence_only_virus_brazil) + nrow(PA1) + nrow(PA2) + nrow(PA3)), 3] <- TRUE



SPC_PresAbs_SamplingKDE_virus <-  BIOMOD_FormatingData(resp.var = Presence_only_virus_brazil_PA$occurrence ,
                                                       expl.var = Env.vector.proj, 
                                                       #expl.var = Env.virus[[-11]],
                                                       resp.xy = Presence_only_virus_brazil_PA[,c('Longitude', 'Latitude')],
                                                       PA.strategy = 'user.defined',
                                                       resp.name = "Culicoides",
                                                       #PA.nb.rep = 3,               # Number of replicates
                                                       PA.nb.absences = nrow(PA1),  # Absences per replicate
                                                       PA.user.table = PA.table,     
                                                       na.rm=T)



SPC_PresAbs_SamplingKDE_virus
plot(SPC_PresAbs_SamplingKDE_virus) ## Quick integrated plot 



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BIOMOD MODELLING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MySpc_models_SamplingKDE_virus <- BIOMOD_Modeling (bm.format = SPC_PresAbs_SamplingKDE_virus,
                                                   models = c("RF"),
                                                   OPT.user = MySpc_options_2_virus,
                                                   #OPT.strategy	= 'bigboss',
                                                   CV.strategy = 'block',  
                                                   #CV.perc = 0.7,
                                                   #var.import = 3,
                                                   metric.eval =c('TSS','ROC'),
                                                   #CV.do.full.models = FALSE,
                                                   modeling.id = "V1")


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Pseudo - Absence sampling based on Sampling of other Arboviruses 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Surveillance map 

disease.list <- list.files("Disease_Occurence_geo/", pattern = ".csv$", full.names = T)

disease.geo <- disease.list %>% 
  lapply(fread) %>% 
  bind_rows()

##Duplicated data for Virus presence
disease.geo_dups <- duplicated(disease.geo[, c('Longitude', 'Latitude')])
sum(disease.geo_dups)
disease.geo_nodups <- disease.geo[!disease.geo_dups, ]  # Keep the records that are not duplicate
disease.geo_nodups$Latitude <- as.numeric(disease.geo_nodups$Latitude)


virus_kde_prep <- st_as_sf(na.omit(disease.geo[,c(1,2)]), coords = c("Longitude", "Latitude"), crs= crs(Brazil))
plot(virus_kde_prep)

points_in_brazil <- virus_kde_prep[st_intersects(virus_kde_prep, brazil.shp, sparse = FALSE), ]
points_in_brazil_df <- as.data.frame(st_coordinates(points_in_brazil))
names(points_in_brazil_df) <- c("Longitude", "Latitude")

## To get unique points for BRazil
a <- duplicated(points_in_brazil_df[, c('Longitude', 'Latitude')])
sum(a)
b <- points_in_brazil_df[!a, ]  # Keep the records that are not duplicate

ggplot()+
  geom_sf(data= world, fill= NA) + 
  geom_sf(data=points_in_brazil)

# Step 1: Convert the extent of the population density raster to a SpatialPolygons object
# Convert your presence points to an `ppp` object
library(sf)
library(raster)
library(spatstat)
library(dismo)

## Count duplicates and assign weights 

viruses_kde1 <- ppp(points_in_brazil_df$Longitude, points_in_brazil_df$Latitude, 
                    window = owin(c(min(points_in_brazil_df$Longitude), max(points_in_brazil_df$Longitude)), 
                                  c(min( points_in_brazil_df$Latitude), max( points_in_brazil_df$Latitude))))


# Perform kernel density estimation
viruses_kde2 <- density.ppp(viruses_kde1, sigma = 3)

viruses_kde2_raster <- raster(viruses_kde2, crs=crs(Env.virus))
extent(viruses_kde2_raster) <- extent(stack(Env.virus))
kde_virus_raster_masked2 <- mask(viruses_kde2_raster, brazil.shp)
plot(kde_virus_raster_masked2)

## If you dont want to include population density as a variable into the pseudo absence sampling 
#combined_density <- kde_virus_raster_masked
### Raster-Based Approach 
# Use a distance raster approach (optional step for sampling purposes)
distance_raster <- distanceFromPoints(kde_virus_raster_masked2, presence_points_sf)

# Create a mask for distances between 25,000 and 50,000 
mask_within_range <- (distance_raster > 50000) & (distance_raster < 500000)
plot(mask_within_range)

#mask_within_range <- resample(mask_within_range, combined_density, method="bilinear")
areaToSample <- kde_virus_raster_masked2*mask_within_range
plot(areaToSample)


# Number of replicates
n_replicates <- 3

# List to store pseudo-absence replicates
pseudo_absence_replicates <- list()

for (i in 1:n_replicates) {
  # Step 1: Perform kernel density estimation (already defined)
  kde_values <- getValues(areaToSample)
  kde_values[is.na(kde_values)] <- 0
  kde_probabilities <- kde_values / sum(kde_values)
  
  # Step 2: Sample pseudo-absence points
  n_points <- 10000
  valid_cells <- which(kde_probabilities > 0)
  sampled_cells <- sample(valid_cells, size = n_points, prob = kde_probabilities[valid_cells], replace = TRUE)
  sampled_coords <- xyFromCell(areaToSample, sampled_cells)
  
  # Step 3: Store the sampled points in a data frame
  pseudo_absences <- data.frame(
    Longitude = sampled_coords[, 1],
    Latitude = sampled_coords[, 2],
    occurrence = 0
  )
  
  pseudo_absence_replicates[[i]] <- pseudo_absences
}

PA1 <- pseudo_absence_replicates[[1]]
PA2 <- pseudo_absence_replicates[[2]]
PA3 <- pseudo_absence_replicates[[3]]

# Combine all pseudo-absence replicates with the original presence points
Presence_only_virus_brazil_PA <- rbind(
  Presence_only_virus_brazil[, c(2:4)],  # Presence points
  PA1,  # Pseudo-absence replicate 1
  PA2,  # Pseudo-absence replicate 2
  PA3   # Pseudo-absence replicate 3
)

# Create the PA.user.table
# Total rows = all presence + all pseudo-absence points
n_rows <- nrow(Presence_only_virus_brazil_PA)

# Create a matrix with zeros
PA.table <- data.frame(
  PA1 = rep(FALSE, nrow(Presence_only_virus_brazil_PA)),
  PA2 = rep(FALSE, nrow(Presence_only_virus_brazil_PA)),
  PA3 = rep(FALSE, nrow(Presence_only_virus_brazil_PA))
)

# Assign presence points as TRUE for all replicates
PA.table[1:nrow(Presence_only_virus_brazil), ] <- TRUE


# Assign pseudo-absence points for each replicate
PA.table[(nrow(Presence_only_virus_brazil) + 1):(nrow(Presence_only_virus_brazil) + nrow(PA1)), 1] <- TRUE
PA.table[(nrow(Presence_only_virus_brazil) + nrow(PA1) + 1):(nrow(Presence_only_virus_brazil) + nrow(PA1) + nrow(PA2)), 2] <- TRUE
PA.table[(nrow(Presence_only_virus_brazil) + nrow(PA1) + nrow(PA2) + 1):(nrow(Presence_only_virus_brazil) + nrow(PA1) + nrow(PA2) + nrow(PA3)), 3] <- TRUE



SPC_PresAbs_KDE_Arbo_virus <-  BIOMOD_FormatingData(resp.var = Presence_only_virus_brazil_PA$occurrence ,
                                                    expl.var = Env.virus,
                                                    resp.xy = Presence_only_virus_brazil_PA[,c('Longitude', 'Latitude')],
                                                    PA.strategy = 'user.defined',
                                                    resp.name = "Culicoides",
                                                    PA.nb.rep = 3,               # Number of replicates
                                                    PA.nb.absences = nrow(PA1),  # Absences per replicate
                                                    PA.user.table = PA.table,     
                                                    na.rm=T)



SPC_PresAbs_KDE_Arbo_virus
plot(SPC_PresAbs_KDE_Arbo_virus) ## Quick integrated plot 



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BIOMOD MODELLING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MySpc_models_KDE_Arbo_virus <- BIOMOD_Modeling (bm.format = SPC_PresAbs_KDE_Arbo_virus,
                                                models = c("RF"),
                                                OPT.user = MySpc_options_2_virus,
                                                #OPT.strategy	= 'bigboss',
                                                CV.strategy = 'block',  
                                                #CV.perc = 0.7,
                                                #var.import = 5,
                                                metric.eval =c('TSS','ROC'),
                                                #CV.do.full.models = FALSE,
                                                modeling.id = "V1")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## PROJECTIONs - once you have identified the model you want to use 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MySpc_models_proj_current_virus <- BIOMOD_Projection( bm.mod = XX, ## replace XX with the correct model 
                                                      new.env = Env.vector.proj, ## read in the raster data onto which you want to predict, it needs to have the same name as the variables at the modeliing stage
                                                      proj.name = "current",
                                                      selected.models = "all",
                                                      binary.meth = "TSS",
                                                      output.format = ".img",
                                                      do.stack = F,
                                                      build.clamping.mask = F)

## BIOMOD displays output in 1000, you need to divide by 1000

Preds_Ind_random_virus <- get_predictions(MySpc_models_proj_current_virus)
plot(Preds_Ind_random_virus/1000, col= custom_palette_vector)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Species Distribution Model = Ensemble of model
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


MySpc_ensemble_models_virus <- BIOMOD_EnsembleModeling(bm.mod = MySpc_models_SamplingKDE_virus,
                                                       models.chosen = "all", ## Choose model you want to combine 
                                                       em.by = 'all', # combine all single models
                                                       #em.algo = "EMcv",
                                                       em.algo = "EMwmean",
                                                       #EMwmean.decay = 2,
                                                       metric.select = 'TSS',
                                                       metric.select.thresh = 0.5,
                                                       metric.eval = c('TSS','ROC')) 

#committee.averaging = TRUE) #Compute the weighted sum of probabilities across predictions



bm_PlotEvalMean(bm.out = MySpc_ensemble_models_random, dataset = 'validation')

# check the scores of the models.
MySpc_ensemble_models_scores_randomsampling <- get_evaluations(MySpc_ensemble_models_virus)
MySpc_ensemble_models_scores_randomsampling

# Ensemble model forecasts

MySpc_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting(bm.em = MySpc_ensemble_models_virus,
                                                                 bm.proj = MySpc_models_proj_current_virus,
                                                                 models.chosen = 'all',
                                                                 metric.binary = 'TSS',
                                                                 metric.filter = 'TSS')

preds_Ens_Virus_RS <- get_predictions(MySpc_ensemble_models_proj_current)
Virus_Projection_RS <- preds_Ens_Virus_RS$Culicoides_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo/1000



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#  \ Modelling Culicoides paraensis distribution /  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Read in Presence and predictor variables 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Make sure you specify the data's location
Presence_only_vector_brazil <- fread("Data/.. /Culi_Pres_Cleaned.csv")
names(Presence_only_vector_brazil)[3] <- "occurrence"
Env.vector <- rast("Data/../Envi_37layers.tif")


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Geographic of Pseudo Absence Sampling of Pseudo Absences 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Presence_only_vector_brazil <- fread("Data/Cleaned_D_to_Model/Orov_presence_cleaned.csv")

Env.vector <- rast("Data/Cleaned_D_to_Model/Envi_37layers.tif")

SPC_PresAbs_Geographic_vector <-  BIOMOD_FormatingData(resp.var = Presence_only_vector_brazil$occurrence ,
                                                       expl.var = Env.vector,
                                                       resp.xy = Presence_only_vector_brazil[,c('Longitude', 'Latitude')],
                                                       resp.name = "Culicoides",
                                                       PA.nb.rep = 3,
                                                       PA.nb.absences = 10000,
                                                       PA.strategy = 'disk',
                                                       PA.dist.min = 50000,
                                                       PA.dist.max = 500000,
                                                       na.rm=T)

SPC_PresAbs_Geographic_vector
plot(SPC_PresAbs_Geographic_vector) ## Quick integrated plot 

### Plot For Pseudo Absences 

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BIOMOD MODELLING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MySpc_models_geographic_vector <- BIOMOD_Modeling (bm.format = SPC_PresAbs_Geographic_vector,
                                                   models = c("RF"),
                                                   OPT.user = MySpc_options_2_vector,
                                                   #OPT.strategy	= 'bigboss',
                                                   CV.strategy = 'block',  
                                                   #CV.perc = 0.7,
                                                   #var.import = 5,
                                                   metric.eval =c('TSS','ROC'),
                                                   #CV.do.full.models = FALSE,
                                                   modeling.id = "V1")



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Pseudo - Absence sampling based on Sampling rate of presence points (only) and with population density 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Convert your presence points to an `ppp` object
presence_points_ppp_vector <- ppp(PP_vector_df$Longitude, PP_vector_df$Latitude, 
                                  window = owin(c(min(PP_vector_df$Longitude), max(PP_vector_df$Longitude)), 
                                                c(min( PP_vector_df$Latitude), max( PP_vector_df$Latitude))))

# Perform kernel density estimation
kde_vector <- density.ppp(presence_points_ppp_vector, sigma = 3)
# Assuming population_density is a raster layer
kde_vector_raster <- raster(kde_vector, crs=crs(pop_density_raster))
extent(kde_vector_raster) <- extent(pop_density_raster)
kde_vector_raster_masked <- mask(kde_vector_raster, brazil.shp)

## only run these lines if you are including population density 
pop_density_raster2 <- resample(pop_density_raster,kde_raster_masked)
combined_density <- overlay(pop_density_raster2, kde_raster_masked, fun=function(x, y) {x * y})

## If you DONT want to include population density as a variable into the pseudo absence sampling 
#combined_density <- kde_vector_raster_masked
### Raster-Based Approach 
# Use a distance raster approach (optional step for sampling purposes)
distance_raster <- distanceFromPoints(kde_raster_masked, presence_points_sf)

# Create a mask for distances between 25,000 and 50,000 
mask_within_range <- (distance_raster > 50000)  & (distance_raster < 500000)
plot(mask_within_range)

#mask_within_range <- resample(mask_within_range, combined_density, method="bilinear")
areaToSample <- combined_density*mask_within_range
plot(areaToSample)


# Number of replicates
n_replicates <- 3

# List to store pseudo-absence replicates
pseudo_absence_replicates <- list()

for (i in 1:n_replicates) {
  # Step 1: Perform kernel density estimation (already defined)
  kde_values <- getValues(areaToSample)
  kde_values[is.na(kde_values)] <- 0
  kde_probabilities <- kde_values / sum(kde_values)
  
  # Step 2: Sample pseudo-absence points
  n_points <- 390
  valid_cells <- which(kde_probabilities > 0)
  sampled_cells <- sample(valid_cells, size = n_points, prob = kde_probabilities[valid_cells], replace = TRUE)
  sampled_coords <- xyFromCell(areaToSample, sampled_cells)
  
  # Step 3: Store the sampled points in a data frame
  pseudo_absences <- data.frame(
    Longitude = sampled_coords[, 1],
    Latitude = sampled_coords[, 2],
    occurrence = 0
  )
  
  pseudo_absence_replicates[[i]] <- pseudo_absences
}

PA1 <- pseudo_absence_replicates[[1]]
PA2 <- pseudo_absence_replicates[[2]]
PA3 <- pseudo_absence_replicates[[3]]

# Combine all pseudo-absence replicates with the original presence points
Presence_only_vector_brazil_PA <- rbind(
  Presence_only_vector_brazil[, c(1:3)],  # Presence points
  PA1,  # Pseudo-absence replicate 1
  PA2,  # Pseudo-absence replicate 2
  PA3   # Pseudo-absence replicate 3
)

# Create the PA.user.table
# Total rows = all presence + all pseudo-absence points
n_rows <- nrow(Presence_only_vector_brazil_PA)

# Create a matrix with zeros
PA.table <- data.frame(
  PA1 = rep(FALSE, nrow(Presence_only_vector_brazil_PA)),
  PA2 = rep(FALSE, nrow(Presence_only_vector_brazil_PA)),
  PA3 = rep(FALSE, nrow(Presence_only_vector_brazil_PA))
)

# Assign presence points as TRUE for all replicates
PA.table[1:nrow(Presence_only_vector_brazil), ] <- TRUE


# Assign pseudo-absence points for each replicate
PA.table[(nrow(Presence_only_vector_brazil) + 1):(nrow(Presence_only_vector_brazil) + nrow(PA1)), 1] <- TRUE
PA.table[(nrow(Presence_only_vector_brazil) + nrow(PA1) + 1):(nrow(Presence_only_vector_brazil) + nrow(PA1) + nrow(PA2)), 2] <- TRUE
PA.table[(nrow(Presence_only_vector_brazil) + nrow(PA1) + nrow(PA2) + 1):(nrow(Presence_only_vector_brazil) + nrow(PA1) + nrow(PA2) + nrow(PA3)), 3] <- TRUE



SPC_PresAbs_SamplingKDE_vector <-  BIOMOD_FormatingData(resp.var = Presence_only_vector_brazil_PA$occurrence ,
                                                        expl.var = Env.vector,
                                                        resp.xy = Presence_only_vector_brazil_PA[,c('Longitude', 'Latitude')],
                                                        PA.strategy = 'user.defined',
                                                        resp.name = "Culicoides",
                                                        PA.nb.rep = 3,               # Number of replicates
                                                        PA.nb.absences = nrow(PA1),  # Absences per replicate
                                                        PA.user.table = PA.table,     
                                                        na.rm=T)



SPC_PresAbs_SamplingKDE_vector
plot(SPC_PresAbs_SamplingKDE_vector) ## Quick integrated plot 



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BIOMOD MODELLING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MySpc_models_SamplingKDE_vector <- BIOMOD_Modeling (bm.format = SPC_PresAbs_SamplingKDE_vector,
                                                    models = c("RF"),
                                                    OPT.user = MySpc_options_2_vector,
                                                    #OPT.strategy	= 'bigboss',
                                                    CV.strategy = 'block',  
                                                    #CV.perc = 0.7,
                                                    var.import = 5,
                                                    metric.eval =c('TSS','ROC'),
                                                    #CV.do.full.models = FALSE,
                                                    modeling.id = "V1")





# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Pseudo - Absence sampling based on Sampling of other Arbovectores 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cera_culi_fam_df <- fread("Data/Cleaned_D_to_Model/cera_culi_fam.csv")

# Step 1: Convert the extent of the population density raster to a SpatialPolygons object
# Convert your presence points to an `ppp` object
library(sf)
library(raster)
library(spatstat)
library(dismo)

## Count duplicates and assign weights 

vectores_kde1 <- ppp(cera_culi_fam_df$Longitude, cera_culi_fam_df$Latitude, 
                     window = owin(c(min(cera_culi_fam_df$Longitude), max(cera_culi_fam_df$Longitude)), 
                                   c(min( cera_culi_fam_df$Latitude), max( cera_culi_fam_df$Latitude))))


# Perform kernel density estimation
vectores_kde2 <- density.ppp(vectores_kde1, sigma = 3)

vectores_kde2_raster <- raster(vectores_kde2, crs=crs(Env.vector))
extent(vectores_kde2_raster) <- extent(stack(Env.vector))
kde_vector_raster_masked2 <- mask(vectores_kde2_raster, brazil.shp)
plot(kde_vector_raster_masked2)

## If you dont want to include population density as a variable into the pseudo absence sampling 
#combined_density <- kde_vector_raster_masked
### Raster-Based Approach 
# Use a distance raster approach (optional step for sampling purposes)
distance_raster <- distanceFromPoints(kde_vector_raster_masked2, presence_points_sf)

# Create a mask for distances between 25,000 and 50,000 
mask_within_range <- (distance_raster > 50000) & (distance_raster < 500000)
plot(mask_within_range)

#mask_within_range <- resample(mask_within_range, combined_density, method="bilinear")
areaToSample <- kde_vector_raster_masked2*mask_within_range
plot(areaToSample)


# Number of replicates
n_replicates <- 3

# List to store pseudo-absence replicates
pseudo_absence_replicates <- list()

for (i in 1:n_replicates) {
  # Step 1: Perform kernel density estimation (already defined)
  kde_values <- getValues(areaToSample)
  kde_values[is.na(kde_values)] <- 0
  kde_probabilities <- kde_values / sum(kde_values)
  
  # Step 2: Sample pseudo-absence points
  n_points <- 390
  valid_cells <- which(kde_probabilities > 0)
  sampled_cells <- sample(valid_cells, size = n_points, prob = kde_probabilities[valid_cells], replace = TRUE)
  sampled_coords <- xyFromCell(areaToSample, sampled_cells)
  
  # Step 3: Store the sampled points in a data frame
  pseudo_absences <- data.frame(
    Longitude = sampled_coords[, 1],
    Latitude = sampled_coords[, 2],
    occurrence = 0
  )
  
  pseudo_absence_replicates[[i]] <- pseudo_absences
}

PA1 <- pseudo_absence_replicates[[1]]
PA2 <- pseudo_absence_replicates[[2]]
PA3 <- pseudo_absence_replicates[[3]]

# Combine all pseudo-absence replicates with the original presence points
Presence_only_vector_brazil_PA <- rbind(
  Presence_only_vector_brazil[, c(1:3)],  # Presence points
  PA1,  # Pseudo-absence replicate 1
  PA2,  # Pseudo-absence replicate 2
  PA3   # Pseudo-absence replicate 3
)

# Create the PA.user.table
# Total rows = all presence + all pseudo-absence points
n_rows <- nrow(Presence_only_vector_brazil_PA)

# Create a matrix with zeros
PA.table <- data.frame(
  PA1 = rep(FALSE, nrow(Presence_only_vector_brazil_PA)),
  PA2 = rep(FALSE, nrow(Presence_only_vector_brazil_PA)),
  PA3 = rep(FALSE, nrow(Presence_only_vector_brazil_PA))
)

# Assign presence points as TRUE for all replicates
PA.table[1:nrow(Presence_only_vector_brazil), ] <- TRUE


# Assign pseudo-absence points for each replicate
PA.table[(nrow(Presence_only_vector_brazil) + 1):(nrow(Presence_only_vector_brazil) + nrow(PA1)), 1] <- TRUE
PA.table[(nrow(Presence_only_vector_brazil) + nrow(PA1) + 1):(nrow(Presence_only_vector_brazil) + nrow(PA1) + nrow(PA2)), 2] <- TRUE
PA.table[(nrow(Presence_only_vector_brazil) + nrow(PA1) + nrow(PA2) + 1):(nrow(Presence_only_vector_brazil) + nrow(PA1) + nrow(PA2) + nrow(PA3)), 3] <- TRUE



SPC_PresAbs_KDE_Arbo_vector <-  BIOMOD_FormatingData(resp.var = Presence_only_vector_brazil_PA$occurrence ,
                                                     expl.var = Env.vector,
                                                     resp.xy = Presence_only_vector_brazil_PA[,c('Longitude', 'Latitude')],
                                                     PA.strategy = 'user.defined',
                                                     resp.name = "Culicoides",
                                                     PA.nb.rep = 3,               # Number of replicates
                                                     PA.nb.absences = nrow(PA1),  # Absences per replicate
                                                     PA.user.table = PA.table,     
                                                     na.rm=T)



SPC_PresAbs_KDE_Arbo_vector
plot(SPC_PresAbs_KDE_Arbo_vector) ## Quick integrated plot 



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BIOMOD MODELLING
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MySpc_models_KDE_Arbo_vector <- BIOMOD_Modeling (bm.format = SPC_PresAbs_KDE_Arbo_vector,
                                                 models = c("RF"),
                                                 OPT.user = MySpc_options_2_vector,
                                                 #OPT.strategy	= 'bigboss',
                                                 CV.strategy = 'block',  
                                                 #CV.perc = 0.7,
                                                 #var.import = 5,
                                                 metric.eval =c('TSS','ROC'),
                                                 #CV.do.full.models = FALSE,
                                                 modeling.id = "V1")


### get models evaluation scores
MyModels_scores_KDE_Arbo_vector <- get_evaluations(MySpc_models_KDE_Arbo_vector)
MyModels_scores_KDE_Arbo_vector












