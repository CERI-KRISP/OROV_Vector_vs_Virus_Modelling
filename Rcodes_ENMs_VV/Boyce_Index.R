# Load required libraries
library(ecospat)
library(raster)

## Reading in all the raster predictions from the different models 

all.files.raster <- list.files("Revised_Analysis/Results_tiff/", recursive = TRUE, full.names = T)

Random.s.files <-list.files("Revised_Analysis/Results_tiff/Random_Sampling/", recursive = TRUE, full.names = T)
Random_sampling <- stack(Random.s.files)
names(Random_sampling) <- basename(Random.s.files)

Geographic.s.files <- list.files("Revised_Analysis/Results_tiff/Geographic_Sampling/", recursive = TRUE, full.names = T)
Geographic_sampling <- stack(Geographic.s.files)
names(Geographic_sampling) <- basename(Geographic.s.files)

Density_weighted_Sampling.s.files <-list.files("Revised_Analysis/Results_tiff/Density_weighted_Sampling/", recursive = TRUE, full.names = T)
Density_weighted_Sampling <- stack(Density_weighted_Sampling.s.files)
names(Density_weighted_Sampling) <- basename(Density_weighted_Sampling.s.files)
Density_weighted_Sampling <- Density_weighted_Sampling/1000


Density_weighted_Population <-list.files("Revised_Analysis/Results_tiff/Density_weighted_Population/", recursive = TRUE, full.names = T)
Density_weighted_Population_sampling <- stack(Density_weighted_Population)
names(Density_weighted_Population_sampling) <- basename(Density_weighted_Population)
Density_weighted_Population_sampling <- Density_weighted_Population_sampling/1000

Target_based_Approach.files <-list.files("Revised_Analysis/Results_tiff/Target_Based_Geographic/", recursive = TRUE, full.names = T)
Target_based_Approach <- stack(Target_based_Approach.files)
names(Target_based_Approach) <- basename(Target_based_Approach.files)
Target_based_Approach <- Target_based_Approach/1000

## Given that they are all of the same extent and resolution etc, you can stack them directly 
ensemble_projection_stack <- stack(Random_sampling,Geographic_sampling,Density_weighted_Sampling,
                               Density_weighted_Population_sampling,Target_based_Approach)
names(ensemble_projection_stack)

# Independent presence points (assumed format: a data frame with Longitude & Latitude)
presence_points <- Presence_test_set[, c("Longitude", "Latitude")]

######
# Function to compute Boyce Index for each raster in the stack
compute_manual_boyce <- function(suitability_raster, presence_points, num_bins = 20) {
  
  # Extract suitability values at independent presence points
  suitability_values <- terra::extract(suitability_raster, presence_points)
  
  #Generate random background points from the raster extent
  background_points <- randomPoints(suitability_raster, n = 1000)  # Adjust n if needed
  
  # Extract suitability values for the background points
  background_values <- terra::extract(suitability_raster, background_points)
  # Remove NAs
  suitability_values <- suitability_values[!is.na(suitability_values)]
  background_values <- background_values[!is.na(background_values)]
  
  # Ensure values exist
  if (length(suitability_values) == 0 || length(background_values) == 0) {
    warning("Skipping Boyce Index computation due to empty values.")
    return(NA)
  }
  
  # Define bin range to span both presence and background values
  bin_min <- min(c(suitability_values, background_values))
  bin_max <- max(c(suitability_values, background_values))
  bins <- seq(bin_min, bin_max, length.out = num_bins + 1)
  
  # Compute frequency of presence points in each bin
  P <- hist(suitability_values, breaks = bins, plot = FALSE)$counts
  # Compute frequency of background points in each bin
  E <- hist(background_values, breaks = bins, plot = FALSE)$counts
  
  # Avoid division by zero
  E[E == 0] <- 1  
  P_E_ratio <- P / E  # Compute P/E ratio
  
  # Compute Spearman correlation (Boyce Index)
  boyce_index <- cor(P_E_ratio, bins[-1], method = "spearman", use = "complete.obs")
  return(boyce_index)
}

# Initialize a dataframe to store results
boyce_results_df <- data.frame(
  PA.nb.absences = integer(),
  Boyce_Index = numeric(),
  stringsAsFactors = FALSE
)


# Compute Boyce Index for each layer in the ensemble raster stack
for (i in 1:nlayers(ensemble_projection_stack)) {
  pa_value <- gsub("Ensemble_", "", names(ensemble_projection_stack)[i])  # Extract PA value
  raster_layer <- ensemble_projection_stack[[i]]
  
  # Compute Boyce Index
  boyce_value <- compute_manual_boyce(raster_layer, presence_points)
  
  # Store results
  boyce_results_df <- rbind(boyce_results_df, data.frame(PA.nb.absences = pa_value, Boyce_Index = boyce_value))
}

# Display results
print(boyce_results_df)

write_csv(boyce_results_df, "Revised_Analysis/Boyce_Index.csv")

## Boyce index values (after reshaping it in excel)
Boyce.index <- fread("Revised_Analysis/Boyce_Index.csv")


Boyce.index <- Boyce.index %>% 
  mutate(
    Sampling_Technique = factor(Sampling_Technique),
    PA_Ratio = factor(PA_Ratio),
    Radii = factor(Radii),
    `Sigma_(KDE)` = factor(`Sigma_(KDE)`)
  )
Boyce.index <- Boyce.index %>% 
  mutate(Sampling_Technique = factor(Sampling_Technique, 
                                  levels = c("Random Sampling", "Geographic Sampling",
                                             "Density-Weighted Geographic Sampling", 
                                             "Density-Weighted Population-Based Sampling",
                                             "Target Group-Based Sampling")))  # Replace with the desired order


custom_labels <- c(
  "Random Sampling" = "Random Sampling",
  "Geographic Sampling" = "Geographic Sampling",
  "Density-Weighted Geographic Sampling" = "Density-Weighted\nGeographic Sampling",
  "Density-Weighted Population-Based Sampling" = "Density-Weighted\nPopulation-Based Sampling",
  "Target Group-Based Sampling" = "Target Group\nBased Sampling"
)



### Figure 1.Model Performance 
Figure1.A <- Boyce.index %>% 
  ggplot(aes(x = Radii, y = Boyce_Index, color = `Sigma_(KDE)`)) +
  #geom_boxplot(aes(group = interaction(PA_Ratio, `Sigma_(KDE)`)), alpha = 0.8,size=0.7) +  # Boxplot overlay
  geom_boxplot(aes(group = interaction(Radii, PA_Ratio,`Sigma_(KDE)`)), alpha = 0.8,size=0.7) +  # Boxplot overlay
  scale_colour_manual(values= c("#824d5c",'#f55654',"#f99551", "#427f80"))+
  facet_grid(PA_Ratio~Sampling_Technique,labeller = labeller(Sampling_Technique = custom_labels)) +  # Facet by Sampling Technique (rows) and Radii (columns)
  theme_minimal(base_size = 14) +  # Professional minimal theme
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1,size=22),
    strip.text = element_text(size = 25, face = "bold"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size=22),
    axis.title.y = element_text(size=25, face="bold"),
    axis.title.x = element_text(size=25,face="bold"),
    legend.title = element_text(size=24),
    legend.text = element_text(size=22), 
    plot.title = element_text(size=30),
    plot.subtitle = element_text(size=25)) +
  labs(
    title = "Comparison of Boyce Index Across Sampling Techniques, PA Ratios, Radii and KDE Smoothing Factor",
    x = "PA Ratio",
    y = "Boyce Index",
    color = "KDE Smoothing Factor") 
  #scale_color_brewer(palette = "RdYlBu")  # Use a professional color scheme
ggsave("Revised_Analysis/Figures/Figure1A.pdf",width=28, height = 15)


## Zoom in on the best perfromaing strategies etc
Figure1.B <- Boyce.index %>% 
  filter(PA_Ratio %in% c("400", "800")) %>% 
  filter(Sampling_Technique %in% c("Density-Weighted Population-Based Sampling", "Target Group-Based Sampling")) %>% 
  ggplot(aes(x = Radii, y = Boyce_Index, color = `Sigma_(KDE)`)) +
  #geom_boxplot(aes(group = interaction(PA_Ratio, `Sigma_(KDE)`)), alpha = 0.8,size=0.7) +  # Boxplot overlay
  geom_boxplot(aes(group = interaction(Radii, PA_Ratio,`Sigma_(KDE)`)), alpha = 0.8,size=0.7) +  # Boxplot overlay
  scale_colour_manual(values= c('#f55654',"#f99551", "#427f80"))+
  facet_grid(PA_Ratio~Sampling_Technique,labeller = labeller(Sampling_Technique = custom_labels)) +  # Facet by Sampling Technique (rows) and Radii (columns)
  theme_minimal(base_size = 14) +  # Professional minimal theme
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, size=22),
    strip.text = element_text(size = 25, face = "bold"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size=22),
    axis.title.y = element_text(size=25, face="bold"),
    axis.title.x = element_text(size=25,face="bold"),
    legend.title = element_text(size=24),
    legend.text = element_text(size=22), 
    plot.title = element_text(size=30),
    plot.subtitle = element_text(size=25)) +
  labs(
    title = "Zoomed from A",
    x = "PA Ratio",
    y = "Boyce Index",
    color = "KDE Smoothing Factor") 

ggsave("Revised_Analysis/Figures/Figure1B.pdf",Figure1.B)


#### 
##smaller plots to differentiate only across the different categories individually 
## BY Pseudo absence ratio 
Figure1.C <- Boyce.index %>% 
  ggplot(aes(x=PA_Ratio, y=Boyce_Index, colour=Sampling_Technique))+
  geom_boxplot(show.legend = T) +
  scale_colour_manual(values= c("#9d41ff",'#F04923',"#0067A5","#ffbf00", "#00A86B"))+
  theme_minimal(base_size = 14) +  # Professional minimal theme
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size=22),
    strip.text = element_text(size = 25, face = "bold"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size=22),
    axis.title.y = element_text(size=25, face="bold"),
    axis.title.x = element_text(size=25,face="bold"),
    legend.title = element_text(size=24),
    legend.text = element_text(size=22), 
    plot.title = element_text(size=30),
    plot.subtitle = element_text(size=25)) +
    #scale_color_brewer(palette= 'Dark2') + 
  labs( x= "PA Ratio", y= "Boyce Index", colour= "Sampling Technique", title = "Comparison of Boyce Index across Pseudo Absence Ratios")
ggsave("Revised_Analysis/Figures/Figure1C.pdf",Figure1.C)



## 
## BY Geographic buffers ratio 
Figure1.D <- Boyce.index %>% 
  ggplot(aes(x=Radii, y=Boyce_Index, colour=Sampling_Technique))+
  geom_boxplot(show.legend = F) +
  scale_colour_manual(values= c("#9d41ff",'#F04923',"#0067A5","#ffbf00", "#00A86B"))+
  theme_minimal(base_size = 14) +  # Professional minimal theme
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size=22),
    strip.text = element_text(size = 25, face = "bold"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size=22),
    axis.title.y = element_text(size=25, face="bold"),
    axis.title.x = element_text(size=25,face="bold"),
    legend.title = element_text(size=24),
    legend.text = element_text(size=22), 
    plot.title = element_text(size=30),
    plot.subtitle = element_text(size=25)) +
  labs( x= "Radii", y= "Boyce Index", colour= "Sampling Technique", title = "Comparison of Boyce Index across Radii")

ggsave("Revised_Analysis/Figures/Figure1D.pdf",Figure1.D)




## BY KDE soomthing factor ratio 
Figure1.E <- Boyce.index %>% 
  ggplot(aes(x=`Sigma_(KDE)`, y=Boyce_Index, colour=Sampling_Technique))+
  geom_boxplot(show.legend=F) +
  scale_colour_manual(values= c("#9d41ff",'#F04923',"#0067A5","#ffbf00", "#00A86B"))+
  theme_minimal(base_size = 14) +  # Professional minimal theme
  theme(
    legend.position = "right",
    axis.text.x = element_text(hjust = 1, size=22),
    strip.text = element_text(size = 25, face = "bold"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size=22),
    axis.title.y = element_text(size=25, face="bold"),
    axis.title.x = element_text(size=25,face="bold"),
    legend.title = element_text(size=24),
    legend.text = element_text(size=22), 
    plot.title = element_text(size=30),
    plot.subtitle = element_text(size=25)) +
  labs( x= "KDE Smoothing Factor", y= "Boyce Index", colour= "Sampling Technique", title = "Comparison of Boyce Index across KDE Smoothing Factor")

ggsave("Revised_Analysis/Figures/Figure1E_Supp.pdf",Figure1.E)







