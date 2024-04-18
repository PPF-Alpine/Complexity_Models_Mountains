# 0 setup ----------------------------------------------------------------------

rm(list = ls()) #full reset

#install.packages("dplyr") # installing packages
#install.packages( c( "dplyr", "ggplot2", "plyr" ) )

## general
library("tidyverse")    #activating the package for running code for it.
library("plyr")

library("ReacTran")
library("fields")
library("animation")

## project specific

library("terra")
library("sf")
library("raster")

wd <- ("C:/Users/repap5991/Downloads/Data/R/Rdocs/")
setwd(wd) #setting the working directory

# 1.0 loading data -----------------------------------------------------------------

#DEM <- raster(paste0(wd, "DataStorage/DEM/dem_geo.sdat"))  ## adding dem_geo chelsa dem

GMBApoly <- st_read(paste0(wd, "DataStorage/ArcGIs Data/GMBA_Level4_ExportTable_TableToExcel.csv"))  ## arcgis dataset which contains all the variables according to the GMBA_level4.


Ama_dx <- raster(paste0(wd, "DataStorage/Amatulli/dtm_dx_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_dxx <- raster(paste0(wd, "DataStorage/Amatulli/dtm_dxx_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_dy <- raster(paste0(wd, "DataStorage/Amatulli/dtm_dy_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_dyy <- raster(paste0(wd, "DataStorage/Amatulli/dtm_dyy_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_slope <- raster(paste0(wd, "DataStorage/Amatulli/dtm_slope_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_vrm <- raster(paste0(wd, "DataStorage/Amatulli/dtm_vrm_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))

Ham_Landfrom <- raster(paste0(wd, "DataStorage/NewHammond/World Ecological Facets L.tif"))
Saph_litho <- raster(paste0(wd, "DataStorage/Saphre/dtm_lithology_usgs.ecotapestry_c_250m_s0..0cm_2014_v1.0.tif"))



# 2.0 data type fixing --------------------------------------------------------------

str(GMBApoly) ## overview of data

numeric_dataconversion <- c("Area..km2.", "Uplift.start..MYR.","CO_slo","CO_VRM", "CO_Dy", "CO_Dx", "CO_Dxx", "CO_Dyy", "CO_Lith", "MAX_slo", "MAX_VRM", "MAX_Dy", "MAX_Dx", "MAX_Dxx", "MAX_Dyy", "MEA_slo", "MEA_VRM", "MEA_Dy","MEA_Dx", "MEA_Dxx", "MEA_Dyy")

 

Catagorical_dataconversion <- c("VAR_Lith")


GMBApoly[numeric_dataconversion] <- lapply(GMBApoly[numeric_dataconversion], as.numeric)

GMBApoly[Catagorical_dataconversion] <- lapply(GMBApoly[Catagorical_dataconversion], as.factor)




str(GMBApoly) ## overview of data

# 3.0 Basic statistics----------------------------------------------------------------------

ggplot(GMBApoly, aes(x = GMBApoly$VAR_LaFor)) +
  geom_bar(stat = "count")


## 3.1 Set up a multi-panel plot for orginal data--------------------------------------------------------


par(mfrow=c(4, 2), mar=c(3,3,1,1))


hist(Ama_dx, main = "Histogram of dx", xlab = "dx")
hist(Ama_dxx, main = "Histogram of dxx", xlab = "dxx")
hist(Ama_dy, main = "Histogram of dy", xlab = "dy")
hist(Ama_dyy, main = "Histogram of dyy", xlab = "dyy")
hist(Ama_slope, main = "Histogram of slope", xlab = "slope")
hist(Ama_vrm, main = "Histogram of vrm", xlab = "vrm")
hist(Ham_Landfrom, main = "Histogram of Landforms", xlab = "Landforms")
hist(Saph_litho, main = "Histogram of LithoClass", xlab = "Lithclasses")



## 3.2  Plot histograms of GMBApoly data--------------------------------------------------------

# Variables and corresponding titles
variablesHis <- c("Area..km2.", "MAX_VRM", "MEA_VRM", "MAX_slo", "MEA_slo", "MAX_Dy", 
               "MEA_Dy", "Uplift.start..MYR.", "MAX_Dx", "MEA_Dx", "MAX_Dxx", 
               "MAX_Dyy", "MEA_Dxx", "MEA_Dyy")
titlesHis <- c("Area", "MAX_VRM", "MEAN_VRM", "MAX_Slope", "MEAN_Slope", "MAX_Dy", 
            "MEAN_Dy", "Uplift_Sta", "MAX_dx", "MEAN_dx", "MAX_Dxx", "MAX_Dyy", 
            "MEAN_Dxx", "MEAN_Dyy")

# Create histograms for each variable
par(mfrow = c(4, 4))  # Set up 4x4 layout for multiple plots

for (i in seq_along(variablesHis)) {
  hist(GMBApoly[[variablesHis[i]]], 
       main = paste("Histogram of", titlesHis[i]), 
       xlab = titlesHis[i],
       ylab = "Frequency",  # Add y-axis label
       col = "#3366FF",  # Change histogram color
       border = "black",  # Add border to bars
       breaks = 20,  # Adjust number of bins
       xlim = range(GMBApoly[[variablesHis[i]]]),  # Set x-axis limits
       xaxt = "n",  # Remove x-axis labels
       yaxt = "s"  # Show y-axis labels
  )
  box()  # Add box around each plot
  abline(v = mean(GMBApoly[[variablesHis[i]]]), col = "red", lwd = 2)  # Add vertical line for mean
  abline(v = median(GMBApoly[[variablesHis[i]]]), col = "blue", lwd = 2, lty = 2)  # Add vertical line for median
}


## 3,3 Violin plots of GMBApoly data --------------------------------------------------------


variables_Violi <- c("Area..km2.", "Uplift.start..MYR.","MAX_slo", "MAX_VRM",
               "MAX_Dy", "MAX_Dx", "MAX_Dxx", "MAX_Dyy", "MEA_slo",
               "MEA_VRM", "MEA_Dy","MEA_Dx", "MEA_Dxx", "MEA_Dyy")

# Reshape the data to long format
melted_data_viol <- gather(GMBApoly, key = "variable", value = "value", all_of(variables))


melted_data_viol$variable <- factor(melted_data$variable, levels = variables)

# Create violin plots with facets and add a boxplot
ggplot(melted_data, aes(x = variable, y = value, fill = variable)) +
  geom_violin() +
  geom_boxplot(width = 0.3, color = "black", outlier.shape = NA) + # Add a boxplot
  facet_wrap(~ variable, scales = "free", ncol = 4) +
  xlab("Variable") +
  ylab("Value") +
  ggtitle("Violin Plots of All Variables with Boxplot") +
  scale_fill_manual(values = rep(rainbow(length(variables)), each = 2)) +
  theme_minimal() +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.key.size = unit(0.5, "cm"))






