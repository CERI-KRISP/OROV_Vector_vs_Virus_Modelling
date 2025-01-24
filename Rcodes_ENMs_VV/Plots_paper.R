
### ### ### ### ### ### ### 
### Plot for Vector 
### ### ### ### ### ### ### 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### Variable Importance Plot ##### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

#####
## Reading it in again (saved after the above)
data_long <- fread("Data/Variable_importance_values.csv")

# Reorder the dataframe by Group and expl.var
data_long <- data_long %>%
  group_by(Group) %>% # Keep grouping intact
  arrange(Group, expl.var,var.imp, .by_group = TRUE) %>% # Sort by group and variable
  ungroup() # Ungroup after arranging


# Convert expl.var to a factor to retain the new order
data_long$expl.var <- factor(data_long$expl.var, levels = unique(data_long$expl.var))

# Plot
VI <- ggplot(data_long, aes(x = var.imp, y = expl.var, fill = Type)) +
  geom_col(position = "identity", width = 0.6) + # Bars for A and B
  geom_text(aes(label = scales::percent(abs(var.imp), accuracy = 0.1)), 
            position = position_dodge(width = 0.6), # Adjust position relative to bars
            hjust = ifelse(data_long$var.imp > 0, -0.4, 1.4), # Position text to the right or left of the bars
            size = 11, family = "Times", colour = "black") + # Adjust text size, family, and color
  scale_fill_manual(
    values = c("var.imp.x" = "#ffbe0b", "var.imp.y" = "#fb5607"), 
    labels = c("Culicoides paraensis", "OROV")
  ) +
  scale_x_continuous(
    limits = c(-0.4 * max(abs(data_long$var.imp)), 1 * max(abs(data_long$var.imp))),
    #breaks = seq(-0.4, 1, by = 0.2),
    #labels = abs(seq(-0.4, 1, by = 0.2)) # Removes negative sign for the left side
  ) +
  geom_vline(xintercept = 0.008, color = "darkred", size = 0.5, linetype = "longdash") + 
  geom_vline(xintercept = -0.007, color = "darkred", size = 0.5, linetype = "longdash") + 
  theme_ipsum() +
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_text(size = 20, family= "Times"), 
    axis.text.y = element_text(size = 28,family= "Times",face = "bold"), # Bold variable names in the center
    axis.text.x = element_text(size = 25, colour="black",family= "Times"),
    legend.title = element_blank(),
    panel.grid.minor.y = element_blank(), # Add x-axis gridlines for clarity
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(), # Add x-axis gridlines for clarity
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.spacing =  unit(5.0, 'cm'),
    #legend.key.h.xeight = unit(0.5, "cm"),
    #legend.justification = c(0.35, 0.5),
    legend.text = element_text(family= "Times",size=28),
    plot.title = element_text(family= "Times", hjust=-1.5, vjust=1, size=40)
  ) +
  guides(fill = guide_legend(byrow = TRUE))+
  labs(
    x = " ",
    title = "Variable Importance for Vector and Virus"
  ) +
  coord_cartesian(clip = "off")
VI


ggsave("Plots/Variable_importance.pdf", VI, height = 30, width = 25)

###########################################################################################
## Projections Across the Americas over current periods - VIRUS 
###########################################################################################
## brazil shapefile
Brazil.admin1 <- st_read("Data/Shapefiles/BRAZIL_shp/gadm41_BRA_1.shp")

custom_palette_vector <- c("#ccd1db","#698a9c","#d2453c","#9b0306")

Virus_IS_v2 <- ggplot()+
  #geom_raster(data = as.data.frame(Preds_Ind_random_virus[[14]]/1000, xy=T), aes(x=x, y=y, fill = Culicoides_PA3_RUN4_RF )) +
  #geom_raster(data = as.data.frame(final_raster_stack_virus$Proj.360.2080, xy=T), aes(x=x, y=y, fill = Proj.360.2080 )) +
  geom_raster(data = as.data.frame(Virus_Projection_RS, xy=T), aes(x=x, y=y, fill = Culicoides_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo )) +
  scale_fill_gradientn(name = "Suitability Index", na.value = "transparent", colours=custom_palette_vector,
                       guide = guide_colorbar(barwidth = 0.7, barheight = 15, title.position = "top"))+
  #geom_sf(data=re_B, fill=NA, colour= "#bc8853",size=0.1)+
  geom_sf(data=America_sf %>% filter(admin !="Canada"), fill=NA, colour= "white",linewidth=0.5)+
  geom_sf(data=re_A, fill=NA, colour="white") +
  coord_sf(xlim = c(-100, -75), ylim = c(2, 22)) +  # adjust these limits as needed to zoom in on Central America
  #coord_sf()+
  #theme_ipsum()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.box= "horizontal")+
  #coord_fixed()+
  xlab(" ") + ylab(" ")+
  labs(title = "Oropouche Virus Suitability Map" , subtitle = " ")
#Virus_IS
Virus_IS_v2


ggsave("Plots/Virus_Current_Projection_Zoom.pdf",Virus_IS_v2)
#writeRaster(preds_Ens/1000,"Vector_Suitability_v2_IS.tiff", overwrite = TRUE)


##########


###########################################################################################
## Projections Across the Americas over current periods - VECTOR 
###########################################################################################

custom_palette_vector <- c("#ccd1db","#698a9c","#d2453c","#9b0306")
final_raster_stack3
Vector_proj <- ggplot()+
  #geom_raster(data = as.data.frame(vector_Projection_RS, xy=T), aes(x=x, y=y, fill = Culicoides_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo )) +
  geom_raster(data = as.data.frame(final_raster_stack3$Proj.360.2080, xy=T), aes(x=x, y=y, fill = Proj.360.2080 )) +
  scale_fill_gradientn(name = "Suitability Index", na.value = "transparent", colours=custom_palette_vector,
                       guide = guide_colorbar(barwidth = 1, barheight = 15, title.position = "top"))+
  #geom_sf(data=re_B, fill=NA, colour= "#bc8853",size=0.1)+
  geom_sf(data=America_sf %>% filter(admin !="Canada"), fill=NA, colour= "white",linewidth=0.5)+
  geom_sf(data=re_A, fill=NA, colour="white") +
  coord_sf(xlim = c(-100, -75), ylim = c(2, 22)) +  # adjust these limits as needed to zoom in on Central America
  #coord_sf()+
  #theme_ipsum()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size=20),
        plot.subtitle = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.box= "horizontal")+
  #coord_fixed()+
  xlab(" ") + ylab(" ")+
  labs(title = "Culicoides Paraensis Suitability Map" , subtitle = "  ")
Vector_proj

ggsave("Plots/Vector_Future_Projection_360_Zoom.pdf",Vector_proj)

###########################################################################################
## Boxplot for TSS and ROC values for models 
###########################################################################################
TSS_ROC <- fread("Data/TSS_ROC_boxplot_data.csv", header = T)

# Split "TSS / ROC" into separate columns
TSS_ROC_long <- TSS_ROC %>%
  pivot_longer(cols = 3:7, names_to = "Pseudo_Absence_Ratio", values_to = "TSS_ROC") %>%
  separate(TSS_ROC, into = c("TSS", "ROC"), sep = " / ", convert = TRUE) %>% 
  separate(TSS, into = c("TSS_min", "TSS_max"), sep = " - ", convert = TRUE) %>% 
  separate(ROC, into = c("ROC_min", "ROC_max"), sep = " - ", convert = TRUE)

# Combine TSS_min and TSS_max into a single column
TSS_ROC_long2 <- TSS_ROC_long %>%
  pivot_longer(cols = c(TSS_min, TSS_max, ROC_min, ROC_max),
               names_to = c("Metric", "Type"),
               names_sep = "_",
               values_to = "Value")


TSS_ROC_long2 <- TSS_ROC_long2 %>%
  mutate(
    Pseudo_Absence_Ratio = factor(
      Pseudo_Absence_Ratio,
      levels = c("1:1 (450)", "1:2 (900)", "1:5 (2250)", "1:10 (4500)", "10 000")
    )
  )

custom_labels <- c(
  "Random Sampling" = "Random Sampling",
  "Geographic Sampling" = "Geographic Sampling",
  "Density-Weighted Geographic Sampling" = "Density-Weighted\nGeographic Sampling",
  "Density-Weighted Population-Based Sampling" = "Density-Weighted\nPopulation-Based Sampling",
  "Target Group Based Sampling" = "Target Group\nBased Sampling"
)

TSS_ROC_long2 <- TSS_ROC_long2 %>%
  mutate(`Sampling Method` = factor(`Sampling Method`, 
                                    levels = c("Random Sampling", "Geographic Sampling",
                                               "Density-Weighted Geographic Sampling", 
                                               "Density-Weighted Population-Based Sampling",
                                               "Target Group Based Sampling")))  # Replace with the desired order

TSS_ROC_long2 %>% 
  filter(Metric == "TSS") %>% 
  group_by(`Sampling Method`) %>% 
  summarise(mmean = median(Value, na.rm=T))


TSS_boxplot <- TSS_ROC_long2 %>% 
  filter(Metric == "TSS") %>% 
  ggplot(aes(x = Pseudo_Absence_Ratio, y = Value, fill = Radius)) +
  geom_boxplot(position= position_dodge(width = 1), width = 0.6) +
  scale_fill_manual(values= c("#c5d86d",'#6d597a',"#82c0cc", "#e56b6f", "#eaac8b"))+
  geom_hline(yintercept=0.66, color= "darkred", linewidth=1, linetype= "longdash")+
  facet_wrap(~`Sampling Method`, nrow=1, labeller = labeller(`Sampling Method`= custom_labels)) +  # Separate plots for TSS and ROC
  labs(
    title = "Boxplots of TSS by Sampling Method",
    subtitle = "Oropouche Virus ENM",
    x = "Pseudo-absence Ratio",
    y = "TSS"
  ) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=19),
        strip.text = element_text (size= 25),
        axis.text.y = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
        legend.title = element_text(size=24),
        legend.text = element_text(size=22), 
        plot.title = element_text(size=30),
        plot.subtitle = element_text(size=25))

TSS_boxplot
ggsave("Plots/TSS_Boxplot.pdf",TSS_boxplot, width=30)

### For Vector 

TSS_ROC_vector <- fread("Data/TSS_ROC_boxplot_Vector_data.csv", header = T)

# Split "TSS / ROC" into separate columns
TSS_ROC_vector_long <- TSS_ROC_vector %>%
  pivot_longer(cols = 3:7, names_to = "Pseudo_Absence_Ratio", values_to = "TSS_ROC") %>%
  separate(TSS_ROC, into = c("TSS", "ROC"), sep = " / ", convert = TRUE) %>% 
  separate(TSS, into = c("TSS_min", "TSS_max"), sep = " - ", convert = TRUE) %>% 
  separate(ROC, into = c("ROC_min", "ROC_max"), sep = " - ", convert = TRUE)

# Combine TSS_min and TSS_max into a single column
TSS_ROC_vector_long2 <- TSS_ROC_vector_long %>%
  pivot_longer(cols = c(TSS_min, TSS_max, ROC_min, ROC_max),
               names_to = c("Metric", "Type"),
               names_sep = "_",
               values_to = "Value")


TSS_ROC_vector_long2 <- TSS_ROC_vector_long2 %>%
  mutate(
    Pseudo_Absence_Ratio = factor(
      Pseudo_Absence_Ratio,
      levels = c("1:1 (78)", "1:2 (156)", "1:5 (390)", "1:10 (780)", "10 000")
    )
  )

custom_labels <- c(
  "Random Sampling" = "Random Sampling",
  "Geographic Sampling" = "Geographic Sampling",
  "Density-Weighted Geographic Sampling" = "Density-Weighted\nGeographic Sampling",
  "Density-Weighted Population-Based Sampling" = "Density-Weighted\nPopulation-Based Sampling",
  "Target Group Based Sampling" = "Target Group\nBased Sampling"
)

TSS_ROC_vector_long2 <- TSS_ROC_vector_long2 %>%
  mutate(`Sampling Method` = factor(`Sampling Method`, 
                                    levels = c("Random Sampling", "Geographic Sampling",
                                               "Density-Weighted Geographic Sampling", 
                                               "Density-Weighted Population-Based Sampling",
                                               "Target Group Based Sampling")))  # Replace with the desired order


TSS_ROC_vector_long2 %>% 
  filter(Metric == "TSS") %>% 
  group_by(`Sampling Method`) %>% 
  summarise(mmean = median(Value, na.rm=T))

TSS_boxplot_vector <- TSS_ROC_vector_long2 %>% 
  filter(Metric == "TSS") %>% 
  ggplot(aes(x = Pseudo_Absence_Ratio, y = Value, fill = Radius)) +
  geom_boxplot(position= position_dodge(width = 1), width = 0.6) +
  geom_hline(yintercept=0.55, color= "darkred", linewidth=1, linetype= "longdash")+
  scale_fill_manual(values= c("#c5d86d",'#6d597a',"#82c0cc", "#e56b6f", "#eaac8b"))+
  facet_wrap(~`Sampling Method`, nrow=1, labeller = labeller(`Sampling Method`= custom_labels,
                                                             scales = "fixed")) +  # Separate plots for TSS and ROC
  labs(
    title = "Boxplots of TSS by Sampling Method ",
    subtitle = "Culicoides paraensis ENM",
    x = "Pseudo-absence Ratio",
    y = "TSS"
  ) +
  theme_clean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=19),
        strip.text = element_text (size= 25),
        axis.text.y = element_text(size=22),
        axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
        legend.title = element_text(size=24),
        legend.text = element_text(size=22), 
        plot.title = element_text(size=30),
        plot.subtitle = element_text(size=25))


TSS_boxplot_vector
ggsave("Plots/TSS_Boxplot_Vector.pdf",TSS_boxplot_vector, width=30)













