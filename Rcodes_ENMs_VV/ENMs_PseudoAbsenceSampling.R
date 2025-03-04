
### Exploratory 
library(corrr)
library(moments)
library(maps)
library(geosphere)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(rworldmap)
#library(plyr)
library(data.table)
library(ggthemes)
library(airportr)
library(hrbrthemes)
library(ggplot2)
library(lubridate)
library(paletteer)
library(RColorBrewer)
library(data.table)
library(readr)
library(MMWRweek)
library(scales)
library(showtext)
library(sysfonts)
#library(rgeos)
library(sf)
library(rnaturalearth)
library(janitor)
library(raster)
library(viridis)
library(ggnewscale)
library(ggrepel)
library(terra)
library(viridis)
library(party)
library(sf)
library(raster)
library(spatstat)
library(dismo)
library(factoextra)

LIB <- c("rgbif", "biomod2", "ggplot2", "gridExtra", "knitr", "raster", 
         "ade4", "rworldmap", "cleangeo", "rasterVis")
#for(i in LIB) { install.packages(i) ; library(i, character.only=T) }
for(i in LIB) { library(i, character.only=T) }

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Setting up covariates
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(ade4)
library(factoextra)
Predictors.files <- list.files("Revised_Analysis/Covariates_Brazil/", pattern= "\\.(asc|tif)$", full.names = T)
Predictors.raster <- lapply(Predictors.files,rast)
names(Predictors.raster[[1]]) <- "Humidity"
names(Predictors.raster[[2]]) <- "Precipitation"
names(Predictors.raster[[3]]) <- "Temperature"
names(Predictors.raster[[4]]) <- "Maximum.Temperature"
names(Predictors.raster[[5]]) <- "Minimum.Temperature"
names(Predictors.raster[[6]]) <- "Soil.Moisture"



reference_raster <- Predictors.raster[[24]]
# Align all rasters to the reference raster
Brazil.Env.variables.aligned <- lapply(Predictors.raster, function(r) {
  aligned.raster = resample(r, reference_raster, method = "bilinear")
  #crop raster: 
  crop.aligned.predictors <- terra::crop(aligned.raster, brazil.shp, mask=T)
  crs(crop.aligned.predictors) <- crs(Predictors.raster[[28]])
  return(crop.aligned.predictors)
  
})

Brazil.predictors.model <- rast(Brazil.Env.variables.aligned)
Brazil.predictors.model <- Brazil.predictors.model[[-c(2,14,22,23)]]
names(Brazil.predictors.model) <- make.names(names(Brazil.predictors.model))
#writeRaster(Brazil.predictors.model, "Revised_Analysis/Model_Predictors_28.tif")
Brazil.predictors.model <- rast("Revised_Analysis/Model_Predictors_28.tif")
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Colinearity / PCA analysis 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Predictors.df <-  na.omit(as.data.frame(Brazil.predictors.model, xy=F))
pca <- ade4::dudi.pca(Predictors.df, scannf = F, nf = 2)


#1. Visualise eigenvalues (Scree plot). Show percentage of variances explained by each principal component

fviz_eig(pca)

#Graph of variables. Variables with a similar profile are grouped together.
pca_plot <- fviz_pca_var(pca,
                         col.var = "contrib", # Color by the quality of representation
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         ggtheme = theme_minimal(),
                         repel = TRUE,
                         labelsize=6)   # Avoid text overlapping)

pca_plot2 <- pca_plot +
  labs(x="PCA axis 1", y= "PCA axis 2", colour= "Contribution",
       title= "Principal Component Analysis")+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        plot.title = element_text(size=18),
        legend.title = element_text(size=17),
        legend.text = element_text(size=15))+
  guides(color=guide_colorbar(barwidth = 1.25, barheight= 14))


### Correlation 
cor_matrix <- cor(Predictors.df, method = "pearson")
abs_cor_matrix <- abs(cor_matrix)
high_correlations <- which(abs_cor_matrix > 0.8, arr.ind = TRUE)

# If there are high correlations, print the variable names and values
if (length(high_correlations) > 0) {
  high_correlations_with_names <- data.frame(
    Predictor1 = rownames(cor_matrix)[high_correlations[, 1]],
    Predictor2 = colnames(cor_matrix)[high_correlations[, 2]],
    Correlation = cor_matrix[high_correlations]
  )
  
  # View the results
  print(high_correlations_with_names)
} else {
  print("No correlations greater than 0.8")
}
plot(Env.virus$Brazil_mixed_trees_earthenv)
plot(Americas.variables.aligned.stack)

heatmap(cor_matrix, 
        main = "Correlation Matrix of Raster Layers", 
        xlab = "Raster Layers", 
        ylab = "Raster Layers")


# We now investigate how our species is distributed across the environmental space 
# ...............................................
## Discriminate species presences from the entire African environmental space. 
# ...............................................
SPC_cell_id <- cellFromXY(Brazil.predictors.model, Presence_only_virus_brazil[,c("Longitude", "Latitude")])
Brazil.rast.df <-  na.omit(as.data.frame(Brazil.predictors.model, xy=F))

s.class(pca$li[,1:2], col=c("grey","red"), csta = 0, cellipse = 2, cpoint = .3, pch = 16,
        fac= factor(rownames(Brazil.rast.df) %in% SPC_cell_id, 
                    levels = c("FALSE", "TRUE" ), labels = c("backgroud", "mySpc")))

s.corcircle(pca$co)
library(factoextra)
library(magrittr)
s.corcircle(pca$co,lab = as.character(1:28),clab = 0.5)

# First: Background & SPC points
s.class(pca$li[,1:2], 
        col=c("grey","red"), 
        csta = 0, 
        cellipse = 2, 
        cpoint = 0.3, 
        pch = 16,
        fac = factor(rownames(Brazil.rast.df) %in% SPC_cell_id, 
                     levels = c("FALSE", "TRUE"), 
                     labels = c("Background", "MySpc")))


Prediction_cell_id <-cellFromXY(Americas.variables.aligned.stack2, Presence_only_virus_brazil[,c("Longitude", "Latitude")]) 
Brazil.pred.df <-  na.omit(as.data.frame(Americas.variables.aligned.stack2, xy=F))
pca2 <- ade4::dudi.pca(Brazil.pred.df, scannf = F, nf = 2)


# Plot PCA with ellipses for all three groups
s.class(pca2$li[,1:2], 
        col = c("grey", "blue"),  # Background = grey, MySpc = red, Prediction = blue
        csta = 0, 
        cellipse = 2,  # Control size of ellipses
        cpoint = 0.3, 
        pch = 16,
        fac  = factor(rownames(Brazil.pred.df) %in% Prediction_cell_id, 
                            levels = c("FALSE", "TRUE"), 
                            labels = c("Prediction Space", "Occurrence points")))
par(bty = "n", mgp=c(2,1,0))  

s.class(
  pca2$li[,2:1], 
  col = c("grey", "blue"),  # Background = grey, Prediction = blue
  csta = 0, 
  cellipse = 2,  # Control size of ellipses
  cpoint = 0.7,  # Increase point size for clarity
  pch = 16,  # Solid points
  grid = F,  # Add grid for better readability
  clabel = 1.6,
  addaxes = F,
  fac  = factor(
    rownames(Brazil.pred.df) %in% Prediction_cell_id, 
    levels = c("FALSE", "TRUE"), 
    labels = c("Prediction Space", "Occurrence points")
  )
)

# Add title and axis labels manually
title(main = "PCA Analysis of Occurrence and Prediction Data", cex.main = 2, font.main = 4)
title(xlab = "PC1 (Variance Explained %)", ylab = "PC2 (Variance Explained %)", cex.lab = 2)




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Presence data - seperate from modeling and testing  
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Presence_only_virus_brazil <- fread("Data/Cleaned_D_to_Model/Orov_presence_cleaned.csv")

## Randomly selected 50 points to set appart for testing 
test_indices <- sample(1:nrow(Presence_only_virus_brazil), 50)

# Create the test set using the selected indices
Presence_test_set <- Presence_only_virus_brazil[test_indices, ]

# Create the training set by removing the selected test points
Presence_train_set <- Presence_only_virus_brazil[-test_indices, ]
plot(st_as_sf(Presence_train_set, coords =c("Longitude", "Latitude"), crs= 4326 ))


## Spatial objects 
Presence_test_set_sf <- st_as_sf(Presence_test_set, coords =c("Longitude", "Latitude"), crs= 4326 )
plot(Presence_test_set_sf)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Sampling Techniques
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Random Sampling
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load required libraries
library(dplyr)
library(BIOMOD2)
library(raster)

# Define the set of PA.nb.absences values to loop over
pa_absence_values <- c(400, 800, 2000, 4000, 10000)

# Create a list to store model results
results_list <- list()

# Create an empty dataframe to store validation summaries
validation_results_df_random <- data.frame(
  PA.nb.absences = integer(),
  metric.eval = character(),
  MeanVal = numeric(),
  min = numeric(),
  max = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each PA.nb.absences value
for (pa in pa_absence_values) {
  cat("Running for PA.nb.absences =", pa, "\n")
  
  # Step 1: Format the data
  SPC_PresAbs_Random_virus <- BIOMOD_FormatingData(
    resp.var = Presence_train_set$occurrence ,
    expl.var = Brazil.predictors.model.without,
    resp.xy = Presence_train_set[,c('Longitude', 'Latitude')],
    resp.name = "Culicoides",
    PA.nb.rep = 3,
    PA.nb.absences = pa,
    PA.strategy = 'random',
    na.rm = TRUE
  )
  
  # Create a unique modeling ID
  model_id <- paste0("V_", pa)
  
  # Step 2: Train the model
  MySpc_models_random_virus <- BIOMOD_Modeling(
    bm.format = SPC_PresAbs_Random_virus,
    models = c("RF"),
    OPT.user = MySpc_options_2_virus,
    CV.strategy = 'block',
    metric.eval = c('TSS', 'ROC'),
    modeling.id = model_id
  )
  
  # Step 3: Retrieve evaluation scores
  eval_scores <- get_evaluations(MySpc_models_random_virus)
  
  # Step 4: Summarize validation scores
  mean_of_validation <- eval_scores %>% 
    group_by(metric.eval) %>% 
    summarise(
      MeanVal = mean(validation, na.rm = TRUE),
      min = min(validation, na.rm = TRUE),
      max = max(validation, na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    mutate(PA.nb.absences = pa)  # Add the PA value for tracking
  
  # Step 5: Model Projection (Single Model)
  MySpc_models_proj_current_virus <- BIOMOD_Projection(
    bm.mod = MySpc_models_random_virus,
    new.env = Brazil.predictors.model.without,
    proj.name = "current",
    selected.models = "all",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE,
    build.clamping.mask = FALSE
  )
  
  # Extract single model projection raster
  #proj_raster_single <- raster(MySpc_models_proj_current_virus@proj.out@val)
  
  # Step 6: Ensemble Modeling
  MySpc_ensemble_models_virus <- BIOMOD_EnsembleModeling(
    bm.mod = MySpc_models_random_virus,
    models.chosen = "all",
    em.by = 'all', # Combine all single models
    em.algo = "EMwmean",
    metric.select = 'TSS',
    metric.select.thresh = 0.5,
    metric.eval = c('TSS','ROC')
  )
  
  # Step 7: Ensemble Projection
  MySpc_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting(
    bm.em = MySpc_ensemble_models_virus,
    bm.proj = MySpc_models_proj_current_virus,
    models.chosen = 'all',
    metric.binary = 'TSS',
    metric.filter = 'TSS'
  )
  
  # Extract final ensemble projection raster
  proj_raster_ensemble.1 <- get_predictions(MySpc_ensemble_models_proj_current)
  proj_raster_ensemble.2 <- get_predictions(MySpc_ensemble_models_proj_current)/1000
  
  # Store results in the list
  results_list[[as.character(pa)]] <- list(
    model = MySpc_models_random_virus,
    evaluations = eval_scores,
    summary = mean_of_validation,
    #projection_single = proj_raster_single,
    ensemble_model = MySpc_ensemble_models_virus
    #ensemble_projection = proj_raster_ensemble
  )
  
  # Append the summary to the validation results dataframe
  validation_results_df_random <- bind_rows(validation_results_df_random, mean_of_validation)
  
  
  # Save the final raster outputs (optional)
  writeRaster(proj_raster_ensemble.2, filename = paste0("Revised_Analysis_2/Random_Sampling/Random_Sampling_", pa, ".tif"), overwrite = TRUE)
  #writeRaster(proj_raster_single, filename = paste0("Single_Projection_", pa, ".tif"), format = "GTiff", overwrite = TRUE)
}

# Display the final validation results
print(validation_results_df_random)
write_csv(validation_results_df_random, "Revised_Analysis_2/Random_sampling_models_results.csv")


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Geographic Sampling
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load required libraries
library(dplyr)
library(BIOMOD2)
library(raster)

# Define the set of PA.nb.absences values to loop over
pa_absence_values <- c(400, 800, 2000, 4000, 10000)

# Define the PA.dist.max values for the disk method
pa_distance_values <- c(150000, 300000, 500000)

# Create a list to store model results
results_list <- list()

# Create an empty dataframe to store validation summaries
validation_results_df_geogra <- data.frame(
  PA.nb.absences = integer(),
  PA.dist.max = integer(),
  metric.eval = character(),
  MeanVal = numeric(),
  min = numeric(),
  max = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each PA.nb.absences value
for (pa in pa_absence_values) {
  
  # Loop over each PA.dist.max value
  for (pa_dist in pa_distance_values) {
    cat("\nRunning DISK PA for PA.nb.absences =", pa, "with PA.dist.max =", pa_dist, "\n")
    
    # Step 1: Format the data using DISK strategy
    SPC_PresAbs_Disk_virus <- BIOMOD_FormatingData(
      resp.var = Presence_train_set$occurrence,
      expl.var = Brazil.predictors.model.without,
      resp.xy = Presence_train_set[, c('Longitude', 'Latitude')],
      resp.name = "Culicoides",
      PA.nb.rep = 3,
      PA.nb.absences = pa,
      PA.strategy = 'disk',
      PA.dist.min = 50000,
      PA.dist.max = pa_dist,
      na.rm = TRUE
    )
    
    # Create a unique modeling ID
    model_id <- paste0("Disk_V", pa, "_D", pa_dist)
    
    # Step 2: Train the model
    MySpc_models_disk_virus <- BIOMOD_Modeling(
      bm.format = SPC_PresAbs_Disk_virus,
      models = c("RF"),
      OPT.user = MySpc_options_2_virus,
      CV.strategy = 'block',
      metric.eval = c('TSS', 'ROC'),
      modeling.id = model_id
    )
    
    # Step 3: Retrieve evaluation scores
    eval_scores_disk <- get_evaluations(MySpc_models_disk_virus)
    
    # Step 4: Summarize validation scores
    mean_of_validation_disk <- eval_scores_disk %>% 
      group_by(metric.eval) %>% 
      summarise(
        MeanVal = mean(validation, na.rm = TRUE),
        min = min(validation, na.rm = TRUE),
        max = max(validation, na.rm = TRUE),
        .groups = "drop"
      ) %>% 
      mutate(PA.nb.absences = pa, PA.dist.max = pa_dist)
    
    # Step 5: Model Projection (Single Model)
    MySpc_models_proj_current_virus <- BIOMOD_Projection(
      bm.mod = MySpc_models_disk_virus,
      new.env = Brazil.predictors.model.without,
      proj.name = "current",
      selected.models = "all",
      binary.meth = "TSS",
      output.format = ".img",
      do.stack = FALSE,
      build.clamping.mask = FALSE
    )
    
    # Step 6: Ensemble Modeling
    MySpc_ensemble_models_virus <- BIOMOD_EnsembleModeling(
      bm.mod = MySpc_models_disk_virus,
      models.chosen = "all",
      em.by = 'all',
      em.algo = "EMwmean",
      metric.select = 'TSS',
      metric.select.thresh = 0.5,
      metric.eval = c('TSS','ROC')
    )
    
    # Step 7: Ensemble Projection
    MySpc_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting(
      bm.em = MySpc_ensemble_models_virus,
      bm.proj = MySpc_models_proj_current_virus,
      models.chosen = 'all',
      metric.binary = 'TSS',
      metric.filter = 'TSS'
    )
    
    # Extract final ensemble projection raster
    proj_raster_ensemble <- get_predictions(MySpc_ensemble_models_proj_current)
    proj_raster_ensemble <- proj_raster_ensemble$Culicoides_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo / 1000
    
    # Store results in the list
    results_list[[paste0("Disk_", pa, "_D", pa_dist)]] <- list(
      model = MySpc_models_disk_virus,
      evaluations = eval_scores_disk,
      summary = mean_of_validation_disk,
      ensemble_model = MySpc_ensemble_models_virus
    )
    
    # Append validation results
    validation_results_df_geogra <- bind_rows(validation_results_df_geogra, mean_of_validation_disk)
    
    # Save ensemble raster
    writeRaster(proj_raster_ensemble, filename = paste0("Revised_Analysis_2/Geographic_Sampling/Geographic_Sampling_", pa, "_D", pa_dist, ".tif"), overwrite = TRUE)
  }
}

# Display validation results
print(validation_results_df_geogra)
write_csv(validation_results_df_geogra, "Revised_Analysis_2/Geographic_sampling_models_results.csv")


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Pseudo - Absence sampling based on Sampling rate of OROV (only) 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load required libraries
#library(spatstat)
#library(terra)
#library(dplyr)
#library(BIOMOD2)

# Define parameter values
sigma_values <- c(1, 3, 5)
buffer_ranges <- list(
  "0" = c(0, 0),  # No buffer
  "50k_150k" = c(50000, 150000),
  "50k_300k" = c(50000, 300000),
  "50k_500k" = c(50000, 500000)
)
n_values <- c(400, 800, 2000, 4000, 10000)

# Create a list to store results
results_list <- list()

# Create an empty dataframe to store validation summaries
validation_results_df_den_weight <- data.frame(
  sigma = integer(),
  PA.dist.max = integer(),
  PA.nb.absences = integer(),
  metric.eval = character(),
  MeanVal = numeric(),
  min = numeric(),
  max = numeric(),
  stringsAsFactors = FALSE
)

# Loop over all sigma values
for (sigma in sigma_values) {
  
  # Compute kernel density estimation for the given sigma
  presence_points_ppp_virus <- ppp(
    Presence_train_set$Longitude, Presence_train_set$Latitude, 
    window = owin(
      c(min(Presence_train_set$Longitude), max(Presence_train_set$Longitude)), 
      c(min(Presence_train_set$Latitude), max(Presence_train_set$Latitude))
    )
  )
  
  kde_virus <- density.ppp(presence_points_ppp_virus, sigma = sigma)
  kde_virus_raster <- rast(kde_virus)
  kde_virus_raster <- resample(kde_virus_raster, Brazil.predictors.model.without, method="bilinear")
  kde_virus_raster <- crop(kde_virus_raster, brazil.shp, mask=T)
  kde_virus_raster <- raster(kde_virus_raster)

    
  # Loop over each buffer distance
  for (buffer_name in names(buffer_ranges)) {
    buffer_min <- buffer_ranges[[buffer_name]][1]
    buffer_max <- buffer_ranges[[buffer_name]][2]
    
    Presence_train_set_sf <- st_as_sf(Presence_train_set, coords =c("Longitude", "Latitude"), crs= 4326 )
    
    # Compute distance raster (only if buffer > 0)
    if (buffer_max > 0) {
      distance_raster <- distanceFromPoints(kde_virus_raster, Presence_train_set_sf)
      mask_within_range <- (distance_raster > buffer_min) & (distance_raster <= buffer_max)
      areaToSample <- kde_virus_raster* mask_within_range
    } else {
      areaToSample <- kde_virus_raster  # No buffer
    }
    
    # Loop over each pseudo-absence ratio
    for (n in n_values) {
      
      cat("\nRunning sigma =", sigma, ", buffer =", buffer_name, ", n =", n, "\n")
      
      # Number of replicates
      n_replicates <- 3
      pseudo_absence_replicates <- list()
      
      for (i in 1:n_replicates) {
        # Extract raster values and normalize
        kde_values <- getValues(areaToSample)
        kde_values[is.na(kde_values)] <- 0
        kde_probabilities <- kde_values / sum(kde_values)
        
        # Identify the presence cells and exclude them
        presence_cells <- cellFromXY(areaToSample, Presence_train_set[, c("Longitude", "Latitude")])
        valid_cells <- which(kde_probabilities > 0)
        valid_cells <- setdiff(valid_cells, presence_cells)
        
        # Sample pseudo-absence points
        sampled_cells <- sample(valid_cells, size = n, prob = kde_probabilities[valid_cells], replace = TRUE)
        sampled_coords <- xyFromCell(areaToSample, sampled_cells)
        
        # Store pseudo-absence points
        pseudo_absences <- data.frame(
          Longitude = sampled_coords[, 1],
          Latitude = sampled_coords[, 2],
          occurrence = 0
        )
        
        pseudo_absence_replicates[[i]] <- pseudo_absences
      }
      
      # Extract pseudo-absence replicates
      PA1 <- pseudo_absence_replicates[[1]]
      PA2 <- pseudo_absence_replicates[[2]]
      PA3 <- pseudo_absence_replicates[[3]]
      
      # Combine presence and pseudo-absence points
      Presence_only_virus_brazil_PA <- rbind(
        Presence_train_set[, c("Longitude", "Latitude", "occurrence")], 
        PA1, PA2, PA3
      )
      
      # Create PA.user.table
      n_rows <- nrow(Presence_only_virus_brazil_PA)
      PA.table <- data.frame(
        PA1 = rep(FALSE, n_rows),
        PA2 = rep(FALSE, n_rows),
        PA3 = rep(FALSE, n_rows)
      )
      
      # Assign presence points as TRUE
      PA.table[1:nrow(Presence_train_set), ] <- TRUE
      
      # Assign pseudo-absence replicates
      PA.table[(nrow(Presence_train_set) + 1):(nrow(Presence_train_set) + nrow(PA1)), 1] <- TRUE
      PA.table[(nrow(Presence_train_set) + nrow(PA1) + 1):(nrow(Presence_train_set) + nrow(PA1) + nrow(PA2)), 2] <- TRUE
      PA.table[(nrow(Presence_train_set) + nrow(PA1) + nrow(PA2) + 1):(nrow(Presence_train_set) + nrow(PA1) + nrow(PA2) + nrow(PA3)), 3] <- TRUE
      
      # Create unique model ID
      model_id <- paste0("Sigma", sigma, "_Buffer", buffer_name, "_n", n)
      
      # Step 2: Train the model
      SPC_PresAbs_Virus <- BIOMOD_FormatingData(
        resp.var = Presence_only_virus_brazil_PA$occurrence,
        expl.var = Brazil.predictors.model.without,
        resp.xy = Presence_only_virus_brazil_PA[, c("Longitude", "Latitude")],
        PA.strategy = 'user.defined',
        resp.name = "Culicoides",
        PA.nb.absences = nrow(PA1),
        PA.user.table = PA.table,
        na.rm = TRUE
      )
      
      MySpc_models_random_virus <- BIOMOD_Modeling(
        bm.format = SPC_PresAbs_Virus,
        models = c("RF"),
        OPT.user = MySpc_options_2_virus,
        CV.strategy = 'block',
        metric.eval = c('TSS', 'ROC'),
        modeling.id = model_id
      )
      
      eval_scores <- get_evaluations(MySpc_models_random_virus)
      
      mean_of_validation <- eval_scores %>% 
        group_by(metric.eval) %>% 
        summarise(
          MeanVal = mean(validation, na.rm = TRUE),
          min = min(validation, na.rm = TRUE),
          max = max(validation, na.rm = TRUE),
          .groups = "drop"
        ) %>% 
        mutate(sigma = sigma, PA.dist.max = buffer_max, PA.nb.absences = n)
      
      validation_results_df_den_weight <- bind_rows(validation_results_df_den_weight, mean_of_validation)
      
      # Step 5: Model Projection (Single Model)
      MySpc_models_proj_current_virus <- BIOMOD_Projection(
        bm.mod = MySpc_models_random_virus,
        new.env = Brazil.predictors.model.without,
        proj.name = "current",
        selected.models = "all",
        binary.meth = "TSS",
        output.format = ".img",
        do.stack = FALSE,
        build.clamping.mask = FALSE
      )
      
      # Step 6: Ensemble Modeling
      MySpc_ensemble_models_virus <- BIOMOD_EnsembleModeling(
        bm.mod = MySpc_models_random_virus,
        models.chosen = "all",
        em.by = 'all',
        em.algo = "EMwmean",
        metric.select = 'TSS',
        metric.select.thresh = 0.3,
        metric.eval = c('TSS','ROC')
      )
      
      MySpc_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting(
        bm.em = MySpc_ensemble_models_virus,
        bm.proj = MySpc_models_proj_current_virus,
        models.chosen = 'all',
        metric.binary = 'TSS',
        metric.filter = 'TSS'
      )
      
      proj_raster_ensemble <- get_predictions(MySpc_ensemble_models_proj_current)
      proj_raster_ensemble.2 <- proj_raster_ensemble$Culicoides_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo / 1000
      writeRaster(proj_raster_ensemble.2, filename = paste0("Revised_Analysis_2/Density_Weighted/Density_Weighted_", model_id, ".tif"), overwrite = TRUE)
    }
  }
}

print(validation_results_df_den_weight)
write_csv(validation_results_df_den_weight, "Revised_Analysis_2/Density_weighted_models_results.csv")


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Pseudo - Absence sampling based on Sampling rate of OROV and with population density 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pop_density_raster.loop  <- pop_density_raster 

# Load required libraries
#library(spatstat)
#library(terra)
#library(dplyr)
#library(BIOMOD2)

# Define parameter values
sigma_values <- c(1, 3, 5)
buffer_ranges <- list(
  "0" = c(0, 0),  # No buffer
  "50k_150k" = c(50000, 150000),
  "50k_300k" = c(50000, 300000),
  "50k_500k" = c(50000, 500000)
)
n_values <- c(400, 800, 2000, 4000, 10000)

# Create a list to store results
results_list <- list()

# Create an empty dataframe to store validation summaries
validation_results_df_den_popu <- data.frame(
  sigma = integer(),
  PA.dist.max = integer(),
  PA.nb.absences = integer(),
  metric.eval = character(),
  MeanVal = numeric(),
  min = numeric(),
  max = numeric(),
  stringsAsFactors = FALSE
)

# Loop over all sigma values
for (sigma in sigma_values) {
  
  # Compute kernel density estimation for the given sigma
  presence_points_ppp_virus <- ppp(
    Presence_train_set$Longitude, Presence_train_set$Latitude, 
    window = owin(
      c(min(Presence_train_set$Longitude), max(Presence_train_set$Longitude)), 
      c(min(Presence_train_set$Latitude), max(Presence_train_set$Latitude))
    )
  )
  
  kde_virus <- density.ppp(presence_points_ppp_virus, sigma = sigma)
  kde_virus_raster <- rast(kde_virus)
  kde_virus_raster <- resample(kde_virus_raster, Brazil.predictors.model.without, method="bilinear")
  kde_virus_raster <- crop(kde_virus_raster, brazil.shp, mask=T)
  kde_virus_raster <- raster(kde_virus_raster)
  
  ## only run these lines if you are including population density 
  pop_density_raster2 <- resample(pop_density_raster,kde_virus_raster, method = "bilinear")
  extent(kde_virus_raster) <- extent(pop_density_raster2)
  combined_density <- overlay(pop_density_raster2, kde_virus_raster, fun=function(x, y) {x * y})
  
  
  # Loop over each buffer distance
  for (buffer_name in names(buffer_ranges)) {
    buffer_min <- buffer_ranges[[buffer_name]][1]
    buffer_max <- buffer_ranges[[buffer_name]][2]
    
    Presence_train_set_sf <- st_as_sf(Presence_train_set, coords =c("Longitude", "Latitude"), crs= 4326 )
    
    # Compute distance raster (only if buffer > 0)
    if (buffer_max > 0) {
      distance_raster <- distanceFromPoints(combined_density, Presence_train_set_sf)
      mask_within_range <- (distance_raster > buffer_min) & (distance_raster <= buffer_max)
      areaToSample <- combined_density* mask_within_range
    } else {
      areaToSample <- combined_density  # No buffer
    }
    
    # Loop over each pseudo-absence ratio
    for (n in n_values) {
      
      cat("\nRunning sigma =", sigma, ", buffer =", buffer_name, ", n =", n, "\n")
      
      # Number of replicates
      n_replicates <- 3
      pseudo_absence_replicates <- list()
      
      for (i in 1:n_replicates) {
        # Extract raster values and normalize
        kde_values <- getValues(areaToSample)
        kde_values[is.na(kde_values)] <- 0
        kde_probabilities <- kde_values / sum(kde_values)
        
        # Identify the presence cells and exclude them
        presence_cells <- cellFromXY(areaToSample, Presence_train_set[, c("Longitude", "Latitude")])
        valid_cells <- which(kde_probabilities > 0)
        valid_cells <- setdiff(valid_cells, presence_cells)
        
        # Sample pseudo-absence points
        sampled_cells <- sample(valid_cells, size = n, prob = kde_probabilities[valid_cells], replace = TRUE)
        sampled_coords <- xyFromCell(areaToSample, sampled_cells)
        
        # Store pseudo-absence points
        pseudo_absences <- data.frame(
          Longitude = sampled_coords[, 1],
          Latitude = sampled_coords[, 2],
          occurrence = 0
        )
        
        pseudo_absence_replicates[[i]] <- pseudo_absences
      }
      
      # Extract pseudo-absence replicates
      PA1 <- pseudo_absence_replicates[[1]]
      PA2 <- pseudo_absence_replicates[[2]]
      PA3 <- pseudo_absence_replicates[[3]]
      
      # Combine presence and pseudo-absence points
      Presence_only_virus_brazil_PA <- rbind(
        Presence_train_set[, c("Longitude", "Latitude", "occurrence")], 
        PA1, PA2, PA3
      )
      
      # Create PA.user.table
      n_rows <- nrow(Presence_only_virus_brazil_PA)
      PA.table <- data.frame(
        PA1 = rep(FALSE, n_rows),
        PA2 = rep(FALSE, n_rows),
        PA3 = rep(FALSE, n_rows)
      )
      
      # Assign presence points as TRUE
      PA.table[1:nrow(Presence_train_set), ] <- TRUE
      
      # Assign pseudo-absence replicates
      PA.table[(nrow(Presence_train_set) + 1):(nrow(Presence_train_set) + nrow(PA1)), 1] <- TRUE
      PA.table[(nrow(Presence_train_set) + nrow(PA1) + 1):(nrow(Presence_train_set) + nrow(PA1) + nrow(PA2)), 2] <- TRUE
      PA.table[(nrow(Presence_train_set) + nrow(PA1) + nrow(PA2) + 1):(nrow(Presence_train_set) + nrow(PA1) + nrow(PA2) + nrow(PA3)), 3] <- TRUE
      
      # Create unique model ID
      model_id <- paste0("Sigma", sigma, "_Buffer", buffer_name, "_n", n)
      
      # Step 2: Train the model
      SPC_PresAbs_Virus <- BIOMOD_FormatingData(
        resp.var = Presence_only_virus_brazil_PA$occurrence,
        expl.var = Brazil.predictors.model.without,
        resp.xy = Presence_only_virus_brazil_PA[, c("Longitude", "Latitude")],
        PA.strategy = 'user.defined',
        resp.name = "Culicoides",
        PA.nb.absences = nrow(PA1),
        PA.user.table = PA.table,
        na.rm = TRUE
      )
      
      MySpc_models_random_virus <- BIOMOD_Modeling(
        bm.format = SPC_PresAbs_Virus,
        models = c("RF"),
        OPT.user = MySpc_options_2_virus,
        CV.strategy = 'block',
        metric.eval = c('TSS', 'ROC'),
        modeling.id = model_id
      )
      
      eval_scores <- get_evaluations(MySpc_models_random_virus)
      
      mean_of_validation <- eval_scores %>% 
        group_by(metric.eval) %>% 
        summarise(
          MeanVal = mean(validation, na.rm = TRUE),
          min = min(validation, na.rm = TRUE),
          max = max(validation, na.rm = TRUE),
          .groups = "drop"
        ) %>% 
        mutate(sigma = sigma, PA.dist.max = buffer_max, PA.nb.absences = n)
      
      validation_results_df_den_popu <- bind_rows(validation_results_df_den_popu, mean_of_validation)
      
      # Step 5: Model Projection (Single Model)
      MySpc_models_proj_current_virus <- BIOMOD_Projection(
        bm.mod = MySpc_models_random_virus,
        new.env = Brazil.predictors.model.without,
        proj.name = "current",
        selected.models = "all",
        binary.meth = "TSS",
        output.format = ".img",
        do.stack = FALSE,
        build.clamping.mask = FALSE
      )
      
      # Step 6: Ensemble Modeling
      MySpc_ensemble_models_virus <- BIOMOD_EnsembleModeling(
        bm.mod = MySpc_models_random_virus,
        models.chosen = "all",
        em.by = 'all',
        em.algo = "EMwmean",
        metric.select = 'ROC',
        metric.select.thresh = 0.1,
        metric.eval = c('TSS','ROC')
      )
      
      MySpc_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting(
        bm.em = MySpc_ensemble_models_virus,
        bm.proj = MySpc_models_proj_current_virus,
        models.chosen = 'all',
        metric.binary = 'TSS',
        metric.filter = 'TSS'
      )
      
      proj_raster_ensemble <- get_predictions(MySpc_ensemble_models_proj_current)
      proj_raster_ensemble.2 <- get_predictions(MySpc_ensemble_models_proj_current)/1000
      writeRaster(proj_raster_ensemble.2, filename = paste0("Revised_Analysis_2/Density_weighted_Population/Density_PopulationD", model_id, ".tif"), overwrite = TRUE)
     }
  }
}

print(validation_results_df_den_popu)
write_csv(validation_results_df_den_popu, "Revised_Analysis_2/Density_weighted_Population_models_results.csv")

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Pseudo - Absence sampling based on Sampling rate of other Arboviruses
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load required libraries
#library(spatstat)
#library(terra)
#library(dplyr)
#library(BIOMOD2)

# Define parameter values
sigma_values <- c(1, 3, 5)
buffer_ranges <- list(
  "0" = c(0, 0),  # No buffer
  "50k_150k" = c(50000, 150000),
  "50k_300k" = c(50000, 300000),
  "50k_500k" = c(50000, 500000)
)
n_values <- c(400, 800, 2000, 4000, 10000)

# Create a list to store results
results_list <- list()

# Create an empty dataframe to store validation summaries
validation_results_df_arbo <- data.frame(
  sigma = integer(),
  PA.dist.max = integer(),
  PA.nb.absences = integer(),
  metric.eval = character(),
  MeanVal = numeric(),
  min = numeric(),
  max = numeric(),
  stringsAsFactors = FALSE
)

# Loop over all sigma values
for (sigma in sigma_values) {
  
  # Compute kernel density estimation for the given sigma
  presence_points_ppp_ARBOVIRUS <- ppp(
    points_in_brazil_df$Longitude, points_in_brazil_df$Latitude, 
    window = owin(
      c(min(points_in_brazil_df$Longitude), max(points_in_brazil_df$Longitude)), 
      c(min(points_in_brazil_df$Latitude), max(points_in_brazil_df$Latitude))
    )
  )
  
  kde_virus <- density.ppp(presence_points_ppp_ARBOVIRUS, sigma = sigma)
  kde_virus_raster <- rast(kde_virus)
  kde_virus_raster <- resample(kde_virus_raster, Brazil.predictors.model.without, method="bilinear")
  kde_virus_raster <- crop(kde_virus_raster, brazil.shp, mask=T)
  kde_virus_raster <- raster(kde_virus_raster)

  
  # Loop over each buffer distance
  for (buffer_name in names(buffer_ranges)) {
    buffer_min <- buffer_ranges[[buffer_name]][1]
    buffer_max <- buffer_ranges[[buffer_name]][2]
    
    Presence_train_set_sf <- st_as_sf(Presence_train_set, coords =c("Longitude", "Latitude"), crs= 4326 )
    
    # Compute distance raster (only if buffer > 0)
    if (buffer_max > 0) {
      distance_raster <- distanceFromPoints(kde_virus_raster, Presence_train_set_sf)
      mask_within_range <- (distance_raster > buffer_min) & (distance_raster <= buffer_max)
      areaToSample <- kde_virus_raster* mask_within_range
    } else {
      areaToSample <- kde_virus_raster  # No buffer
    }
    
    # Loop over each pseudo-absence ratio
    for (n in n_values) {
      
      cat("\nRunning sigma =", sigma, ", buffer =", buffer_name, ", n =", n, "\n")
      
      # Number of replicates
      n_replicates <- 3
      pseudo_absence_replicates <- list()
      
      for (i in 1:n_replicates) {
        # Extract raster values and normalize
        kde_values <- getValues(areaToSample)
        kde_values[is.na(kde_values)] <- 0
        kde_probabilities <- kde_values / sum(kde_values)
        
        # Identify the presence cells and exclude them
        presence_cells <- cellFromXY(areaToSample, Presence_train_set[, c("Longitude", "Latitude")])
        valid_cells <- which(kde_probabilities > 0)
        valid_cells <- setdiff(valid_cells, presence_cells)
        
        # Sample pseudo-absence points
        sampled_cells <- sample(valid_cells, size = n, prob = kde_probabilities[valid_cells], replace = TRUE)
        sampled_coords <- xyFromCell(areaToSample, sampled_cells)
        
        # Store pseudo-absence points
        pseudo_absences <- data.frame(
          Longitude = sampled_coords[, 1],
          Latitude = sampled_coords[, 2],
          occurrence = 0
        )
        
        pseudo_absence_replicates[[i]] <- pseudo_absences
      }
      
      # Extract pseudo-absence replicates
      PA1 <- pseudo_absence_replicates[[1]]
      PA2 <- pseudo_absence_replicates[[2]]
      PA3 <- pseudo_absence_replicates[[3]]
      
      # Combine presence and pseudo-absence points
      Presence_only_virus_brazil_PA <- rbind(
        Presence_train_set[, c("Longitude", "Latitude", "occurrence")], 
        PA1, PA2, PA3
      )
      
      # Create PA.user.table
      n_rows <- nrow(Presence_only_virus_brazil_PA)
      PA.table <- data.frame(
        PA1 = rep(FALSE, n_rows),
        PA2 = rep(FALSE, n_rows),
        PA3 = rep(FALSE, n_rows)
      )
      
      # Assign presence points as TRUE
      PA.table[1:nrow(Presence_train_set), ] <- TRUE
      
      # Assign pseudo-absence replicates
      PA.table[(nrow(Presence_train_set) + 1):(nrow(Presence_train_set) + nrow(PA1)), 1] <- TRUE
      PA.table[(nrow(Presence_train_set) + nrow(PA1) + 1):(nrow(Presence_train_set) + nrow(PA1) + nrow(PA2)), 2] <- TRUE
      PA.table[(nrow(Presence_train_set) + nrow(PA1) + nrow(PA2) + 1):(nrow(Presence_train_set) + nrow(PA1) + nrow(PA2) + nrow(PA3)), 3] <- TRUE
      
      # Create unique model ID
      model_id <- paste0("Sigma", sigma, "_Buffer", buffer_name, "_n", n)
      
      # Step 2: Train the model
      SPC_PresAbs_Virus <- BIOMOD_FormatingData(
        resp.var = Presence_only_virus_brazil_PA$occurrence,
        expl.var = Brazil.predictors.model.without,
        resp.xy = Presence_only_virus_brazil_PA[, c("Longitude", "Latitude")],
        PA.strategy = 'user.defined',
        resp.name = "Culicoides",
        PA.nb.absences = nrow(PA1),
        PA.user.table = PA.table,
        na.rm = TRUE
      )
      
      MySpc_models_random_virus <- BIOMOD_Modeling(
        bm.format = SPC_PresAbs_Virus,
        models = c("RF"),
        OPT.user = MySpc_options_2_virus,
        CV.strategy = 'block',
        metric.eval = c('TSS', 'ROC'),
        modeling.id = model_id
      )
      
      eval_scores <- get_evaluations(MySpc_models_random_virus)
      
      mean_of_validation <- eval_scores %>% 
        group_by(metric.eval) %>% 
        summarise(
          MeanVal = mean(validation, na.rm = TRUE),
          min = min(validation, na.rm = TRUE),
          max = max(validation, na.rm = TRUE),
          .groups = "drop"
        ) %>% 
        mutate(sigma = sigma, PA.dist.max = buffer_max, PA.nb.absences = n)
      
      validation_results_df_arbo <- bind_rows(validation_results_df_arbo, mean_of_validation)
      
      # Step 5: Model Projection (Single Model)
      MySpc_models_proj_current_virus <- BIOMOD_Projection(
        bm.mod = MySpc_models_random_virus,
        new.env = Brazil.predictors.model.without,
        proj.name = "current",
        selected.models = "all",
        binary.meth = "TSS",
        output.format = ".img",
        do.stack = FALSE,
        build.clamping.mask = FALSE
      )
      
      # Step 6: Ensemble Modeling
      MySpc_ensemble_models_virus <- BIOMOD_EnsembleModeling(
        bm.mod = MySpc_models_random_virus,
        models.chosen = "all",
        em.by = 'all',
        em.algo = "EMwmean",
        metric.select = 'TSS',
        metric.select.thresh = 0.1,
        metric.eval = c('TSS','ROC')
      )
      
      MySpc_ensemble_models_proj_current <- BIOMOD_EnsembleForecasting(
        bm.em = MySpc_ensemble_models_virus,
        bm.proj = MySpc_models_proj_current_virus,
        models.chosen = 'all',
        metric.binary = 'TSS',
        metric.filter = 'TSS'
      )
      
      proj_raster_ensemble <- get_predictions(MySpc_ensemble_models_proj_current)
      proj_raster_ensemble.2 <- get_predictions(MySpc_ensemble_models_proj_current)/1000
      writeRaster(proj_raster_ensemble.2, filename = paste0("Revised_Analysis_2/Target_Based_Approach//TargetBased_", model_id, ".tif"), overwrite = TRUE)
    }
  }
}

print(validation_results_df_arbo)
write_csv(validation_results_df_arbo, "Revised_Analysis_2/Target_Group_Based_Sampling_Result.csv")

######## 


