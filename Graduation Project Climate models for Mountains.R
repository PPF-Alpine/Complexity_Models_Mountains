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
library("terra") #for working with rasters
library("sf") #for working with rasters
library("raster") #for working with rasters
library("factoextra")  # For PCA visualization


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
DEM <- raster(paste0(wd, "DataStorage/DEM/Chelsa V2/dem_latlong.sdat"))
Ham_Landfrom <- raster(paste0(wd, "DataStorage/NewHammond/World Ecological Facets L.tif"))
Saph_litho <- raster(paste0(wd, "DataStorage/Saphre/dtm_lithology_usgs.ecotapestry_c_250m_s0..0cm_2014_v1.0.tif"))


# 2.0 data type fixing --------------------------------------------------------------

str(GMBApoly) ## overview of data

numeric_Data <- c(4:19, 23, 25:41) 
GMBApoly[numeric_Data] <- lapply(GMBApoly[numeric_Data], as.numeric) #set variables as numeric

Catagorical_data <- c(2, 20, 21, 22, 24)
GMBApoly[Catagorical_data] <- lapply(GMBApoly[Catagorical_data], as.factor) #set variables as factor




# 3.0 Basic statistics and overview---------------------------------------------------
## 3.1 Set up a multi-panel plot for orginal data--------------------------------------------------------

# Define colors for histograms
hist_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#17becf")

# Create a list of raster datasets
raster_datasets <- list(
  Ama_dx = "dx",
  Ama_dxx = "dxx",
  Ama_dy = "dy",
  Ama_dyy = "dyy",
  Ama_slope = "slope",
  Ama_vrm = "vrm",
  Ham_Landfrom = "Landforms",
  Saph_litho = "LithoClass",
  DEM = "DEM"
)

# Set up the plotting layout
par(mfrow=c(5, 2), mar=c(3,3,1,1))

# Generate histograms for each dataset using ggplot2
for (i in seq_along(raster_datasets)) {
  dataset <- raster_datasets[[i]]
  raster_data <- get(names(raster_datasets)[i])
  
  hist(raster_data, col = hist_colors[i], main = paste("Histogram of", dataset), xlab = "", ylab = "Frequency")
}


## 3.2  Plot histograms of GMBApoly data--------------------------------------------------------

print(names(GMBApoly)[variablesHis])

Vari_Num <- c(4, 6, 20:33, 38, 39) # Define the column indices of interest Max, mean, uplift, area
Title_var_Num <- names(GMBApoly)[Vari_Num] # Extract the column names of interest



# Create histograms for each variable
par(mfrow = c(5, 4))  # Set up 4x4 layout for multiple plots

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
  abline(v = mean(GMBApoly[[Vari_Num[i]]]), col = "red", lwd = 2)  # Add vertical line for mean
  abline(v = median(GMBApoly[[Vari_Num[i]]]), col = "blue", lwd = 2, lty = 2)  # Add vertical line for median
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
  ggtitle("Violin Plots of GNBApoly All Numeric Variables with Boxplot") +
  scale_fill_manual(values = rep(rainbow(length(Title_var_Num)), each = 2)) +
  theme_minimal() +
  theme(legend.direction = "vertical", legend.box = "vertical", legend.key.size = unit(0.5, "cm"))




## 3.4 summarizing results PCA


PCAData <- GMBApoly[,c(4, 6, 33:41)] #taking only the numeric values 

res.pca <- prcomp(PCAData, scale = TRUE) #rescaled variables, which makes that they're scaled according to eachother and limiting effect of large ranges compared with smaller.

if (anyNA(PCAData)) {
  PCAData <- na.omit(PCAData)  # Remove rows with missing values, as the PCA doesn't work with that
}

# Visualize PCA results
pca_plot <- fviz_pca(res.pca, geom = "point", ggtheme = theme_minimal())







