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

desired_order <- c(1, 5, 2, 3, 4, 6, 7) # Define the desired order for the first columns
GMBApoly <- GMBApoly[, c(desired_order, setdiff(order(names(GMBApoly)), desired_order))] # Reorder the columns in the GMBApoly dataframe


Ama_dx <- raster(paste0(wd, "DataStorage/Amatulli/dtm_dx_merit.dem_m_250m_s0..0cm_2018_v1.0.tif")) #Original datasets
Ama_dxx <- raster(paste0(wd, "DataStorage/Amatulli/dtm_dxx_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_dy <- raster(paste0(wd, "DataStorage/Amatulli/dtm_dy_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_dyy <- raster(paste0(wd, "DataStorage/Amatulli/dtm_dyy_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_slope <- raster(paste0(wd, "DataStorage/Amatulli/dtm_slope_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_vrm <- raster(paste0(wd, "DataStorage/Amatulli/dtm_vrm_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))

Ham_Landfrom <- raster(paste0(wd, "DataStorage/NewHammond/World Ecological Facets L.tif"))
Saph_litho <- raster(paste0(wd, "DataStorage/Saphre/dtm_lithology_usgs.ecotapestry_c_250m_s0..0cm_2014_v1.0.tif"))


# 2.0 data type fixing --------------------------------------------------------------

str(GMBApoly) ## overview of data

numeric_dataconversion <- c(4:30, 32, 34:36)
Catagorical_dataconversion <- c(31, 33, 3, 2)


GMBApoly[numeric_dataconversion] <- lapply(GMBApoly[numeric_dataconversion], as.numeric)
GMBApoly[Catagorical_dataconversion] <- lapply(GMBApoly[Catagorical_dataconversion], as.factor)


# 3.0 Basic statistics and overview---------------------------------------------------


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


# Create a list of raster datasets
raster_datasets <- list(
  Ama_dx = list(title = "Histogram of dx", x_label = "dx"),
  Ama_dxx = list(title = "Histogram of dxx", x_label = "dxx"),
  Ama_dy = list(title = "Histogram of dy", x_label = "dy"),
  Ama_dyy = list(title = "Histogram of dyy", x_label = "dyy"),
  Ama_slope = list(title = "Histogram of slope", x_label = "slope"),
  Ama_vrm = list(title = "Histogram of vrm", x_label = "vrm"),
  Ham_Landfrom = list(title = "Histogram of Landforms", x_label = "Landforms"),
  Saph_litho = list(title = "Histogram of LithoClass", x_label = "Lithclasses")
)

# Set up the plotting layout
par(mfrow=c(4, 2), mar=c(3,3,1,1))

# Generate histograms for each dataset using ggplot2
for (name in names(raster_datasets)) {
  raster_data <- get(name)
  values <- raster_data[]
  ggplot(data.frame(Data = values), aes(x = Data)) +
    geom_histogram(binwidth = 1, color = "black") +
    labs(x = raster_datasets[[name]]$x_label, y = "Frequency", title = raster_datasets[[name]]$title)
}
## 3.2  Plot histograms of GMBApoly data--------------------------------------------------------

print(names(GMBApoly)[variablesHis])



Vari_Num <- c(4, 6, 19:30, 35, 36) # Define the column indices of interest Max, mean, uplift, area
Title_var_Num <- names(GMBApoly)[Vari_Num] # Extract the column names of interest
# Create histograms for each variable
par(mfrow = c(4, 4))  # Set up 4x4 layout for multiple plots

for (i in seq_along(Vari_Num)) {
  hist(GMBApoly[[Vari_Num[i]]], 
       main = paste("Histogram of", Title_var_Num[i]), 
       xlab = Title_var_Num[i],
       ylab = "Frequency",  # Add y-axis label
       col = "#3366FF",  # Change histogram color
       border = "black",  # Add border to bars
       breaks = 20,  # Adjust number of bins
       xlim = range(GMBApoly[[Vari_Num[i]]]),  # Set x-axis limits
       xaxt = "n",  # Remove x-axis labels
       yaxt = "s"  # Show y-axis labels
  )
  box()  # Add box around each plot
  abline(v = mean(GMBApoly[[variablesHis[i]]]), col = "red", lwd = 2)  # Add vertical line for mean
  abline(v = median(GMBApoly[[variablesHis[i]]]), col = "blue", lwd = 2, lty = 2)  # Add vertical line for median
}


## 3,3 Violin plots of GMBApoly data --------------------------------------------------------


# Reshape the data to long format
melted_data <- gather(GMBApoly, key = "Variable", value = "Value", all_of(Vari_Num))


melted_data$variable <- factor(melted_data$variable, levels = Vari_Num)

# Create violin plots with facets and add a boxplot
ggplot(melted_data, aes(x = Variable, y = Value, fill = Variable)) +
  geom_violin() +
  geom_boxplot(width = 0.3, color = "black", outlier.shape = NA) + # Add a boxplot
  facet_wrap(~ Variable, scales = "free", ncol = 4) +
  xlab("Variable") +
  ylab("Value") +
  ggtitle("Violin Plots of All Numeric Variables with Boxplot") +
  scale_fill_manual(values = rep(rainbow(length(Title_var_Num)), each = 2)) +
  theme_minimal() +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.key.size = unit(0.5, "cm"))






