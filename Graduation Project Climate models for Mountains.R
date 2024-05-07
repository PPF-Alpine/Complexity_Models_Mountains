# 0 setup ------------------------------------------------------------------------------------------------------------------------------------------

rm(list = ls()) #full reset

rm(dataset)

rm(list = "^Temp_|^Prec_tester", ls(), value = TRUE)

source("C:/Users/repap5991/OneDrive - University of Bergen/Data/R/Rdocs/Tester.R") #load libraries 


wd <- ("C:/Users/repap5991/OneDrive - University of Bergen/Data/R/Rdocs/")
setwd(wd) #setting the working directory

# 1.0 loading data -------------------------------------------------------------------------------------------------------------------------------------

## 1.1 GMBA ArcGIS polygon load  -----------------------------------------------------------------------------------------------------------------
GMBApoly <- st_read(paste0(wd, "DataStorage/ArcGIs Data/GMBA_Level4_ExportTable_TableToExcel.csv"))  ## arcgis dataset which contains all the variables according to the GMBA_level4.

## 1.2 orginal geomorp variables   -----------------------------------------------------------------------------------------------------------------

Ama_dx <- raster(paste0(wd, "DataStorage/Amatulli/dtm_dx_merit.dem_m_250m_s0..0cm_2018_v1.0.tif")) #Original datasets
Ama_dxx <- raster(paste0(wd, "DataStorage/Amatulli/dtm_dxx_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_dy <- raster(paste0(wd, "DataStorage/Amatulli/dtm_dy_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_dyy <- raster(paste0(wd, "DataStorage/Amatulli/dtm_dyy_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_slope <- raster(paste0(wd, "DataStorage/Amatulli/dtm_slope_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
Ama_vrm <- raster(paste0(wd, "DataStorage/Amatulli/dtm_vrm_merit.dem_m_250m_s0..0cm_2018_v1.0.tif"))
DEM <- raster(paste0(wd, "DataStorage/DEM/Chelsa V2/dem_latlong.sdat"))
Ham_Landfrom <- raster(paste0(wd, "DataStorage/NewHammond/World Ecological Facets L.tif"))
Saph_litho <- raster(paste0(wd, "DataStorage/Saphre/dtm_lithology_usgs.ecotapestry_c_250m_s0..0cm_2014_v1.0.tif"))

## 1.3 chelsa climatic variables   -----------------------------------------------------------------------------------------------------------------

Prec_CCSM_Past <- raster(paste0(wd, "DataStorage/Chelsa/Prec Bio Past/CHELSA_PMIP_CCSM4_BIO_12.tif")) #BIO 12 is the annual precipitation amount kg m-2, scale 0.1
Prec_CNRM_Past <- raster(paste0(wd, "DataStorage/Chelsa/Prec Bio Past/CHELSA_PMIP_CNRM-CM5_BIO_12.tif")) #Accumulated precipitation amount over 1 year
Prec_FGOA_Past <- raster(paste0(wd, "DataStorage/Chelsa/Prec Bio Past/CHELSA_PMIP_FGOALS-g2_BIO_12.tif")) # in the past LGM 21000 BP (before 1950)
Prec_IPSL_Past <- raster(paste0(wd, "DataStorage/Chelsa/Prec Bio Past/CHELSA_PMIP_IPSL-CM5A-LR_BIO_12.tif")) 
Prec_MIRO_Past <- raster(paste0(wd, "DataStorage/Chelsa/Prec Bio Past/CHELSA_PMIP_MIROC-ESM_BIO_12.tif")) 
Prec_MPIE_Past <- raster(paste0(wd, "DataStorage/Chelsa/Prec Bio Past/CHELSA_PMIP_MPI-ESM-P_BIO_12.tif")) 
Prec_MRIC_Past <- raster(paste0(wd, "DataStorage/Chelsa/Prec Bio Past/CHELSA_PMIP_MRI-CGCM3_BIO_12.tif")) 

Prec_Pres <- raster(paste0(wd, "DataStorage/Chelsa/Prec Bio Present/CHELSA_bio12_1981-2010_V.2.1.tif")) #BIO 12 is the annual precipitation amount kg m-2, scale 0.1 Accumulated precipitation amount over 1 year in the Present 1981 - 2010


Temp_CCSM_Past <- raster(paste0(wd, "DataStorage/Chelsa/Temp Bio Past/CHELSA_PMIP_CCSM4_BIO_01.tif")) #BIO 01 is the mean annual temperature B0C, scale 0.1, -273.15 offset
Temp_CNRM_Past <- raster(paste0(wd, "DataStorage/Chelsa/Temp Bio Past/CHELSA_PMIP_CNRM-CM5_BIO_01.tif")) #Accumulated precipitation amount over 1 year
Temp_FGOA_Past <- raster(paste0(wd, "DataStorage/Chelsa/Temp Bio Past/CHELSA_PMIP_FGOALS-g2_BIO_01.tif")) # in the past LGM 21000 YR
Temp_IPSL_Past <- raster(paste0(wd, "DataStorage/Chelsa/Temp Bio Past/CHELSA_PMIP_IPSL-CM5A-LR_BIO_01.tif"))
Temp_MIRO_Past <- raster(paste0(wd, "DataStorage/Chelsa/Temp Bio Past/CHELSA_PMIP_MIROC-ESM_BIO_01.tif"))
Temp_MPIE_Past <- raster(paste0(wd, "DataStorage/Chelsa/Temp Bio Past/CHELSA_PMIP_MPI-ESM-P_BIO_01.tif"))
Temp_MRIC_Past <- raster(paste0(wd, "DataStorage/Chelsa/Temp Bio Past/CHELSA_PMIP_MRI-CGCM3_BIO_01.tif"))

Temp_Pres <- raster(paste0(wd, "DataStorage/Chelsa/Temp Bio Present/CHELSA_bio1_1981-2010_V.2.1.tif")) #BIO 01 is the mean annual temperature B0C, scale 0.1, offset, -273.15 Accumulated Temperature amount over 1 year in the Present 1981 - 2010.

# 2.0 data type fixing -----------------------------------------------------------------------------------------------------

## 2.1 GMBApoly  ----------------------------------------------------------------------------------------------------------------------------------

desired_order <- c(1, 5, 2, 3, 4, 6, 7) # Define the desired order for the first columns
GMBApoly <- GMBApoly[, c(desired_order, setdiff(order(names(GMBApoly)), desired_order))] # Reorder the columns in the GMBApoly dataframe

str(GMBApoly) ## overview of data

numeric_Data <- c(4:19, 23, 25:41) 
GMBApoly[numeric_Data] <- lapply(GMBApoly[numeric_Data], as.numeric) #set variables as numeric

Catagorical_data <- c(2, 20, 21, 22, 24)
GMBApoly[Catagorical_data] <- lapply(GMBApoly[Catagorical_data], as.factor) #set variables as factor


## 2.2                ----------------------------------------------------------------------------------------------------------------------------------

## 2.3 CHELSA data  -------------------------------------------------------------------------------------------------------

# Define a list of temperature and precipitation raster objects
temperature_rasters <- list(
  Temp_CCSM_Past, Temp_CNRM_Past, Temp_FGOA_Past, Temp_IPSL_Past,
  Temp_MIRO_Past, Temp_MPIE_Past, Temp_MRIC_Past
)

precipitation_rasters <- list(
  Prec_CCSM_Past, Prec_CNRM_Past, Prec_FGOA_Past, Prec_IPSL_Past,
  Prec_MIRO_Past, Prec_MPIE_Past, Prec_MRIC_Past
)

# Loop through the temperature and precipitation raster lists
for (raster_list in list(temperature_rasters, precipitation_rasters)) {
  for (raster in raster_list) {
    raster[raster < -10000 | raster >= 32767] <- NA
  }
}

# 3.0 Basic statistics and overview-----------------------------------------------------------------------------------------------------------------------

## 3.1 Set up a multi-panel plot for original data--------------------------------------------------------

# Define colors for histograms
hist_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                 "#8c564b", "#e377c2", "#7f7f7f", "#17becf")
# Create a list of original raster data sets with variable descriptions
raster_datasets <- list(
  Ama_dx = "East-West Slope", Ama_dxx = "East-West Curvature Rate", Ama_dy = "North-South Slope", Ama_dyy = "North-South Curvature Rate", Ama_slope = "Terrain Slope",
  Ama_vrm = "Terrain Ruggedness Measure", Ham_Landform = "Landform Classification", Saph_litho = "Lithology Classification", DEM = "Elevation Data (m)"
)

# Set up the plotting layout
par(mfrow=c(5, 2), mar=c(3,3,1,1))

# Generate histograms for each dataset using ggplot2
for (i in seq_along(raster_datasets)) {
  dataset <- raster_datasets[[i]]
  raster_data <- get(names(raster_datasets)[i])
  
  # Extract variable description for x-axis
  x_axis_label <- switch(dataset,
                         "East-West Slope" = "East-West Slope",
                         "East-West Curvature Rate" = "East-West Curvature Rate",
                         "North-South Slope" = "North-South Slope",
                         "North-South Curvature Rate" = "North-South Curvature Rate",
                         "Terrain Slope" = "Terrain Slope",
                         "Terrain Ruggedness Measure" = "Terrain Ruggedness Measure",
                         "Landform Classification" = "Landform Classification",
                         "Lithology Classification" = "Lithology Classification",
                         "Elevation Data (m)" = "Elevation Data (m)"
  )
  
  hist(raster_data, col = hist_colors[i], main = paste(dataset), xlab = x_axis_label, ylab = "Frequency")
}




## 3.2 Set up a multi-panel plot for CHELSA data--------------------------------------------------------


# Define colors for histograms
hist_colors_chel<- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2")

# Create a list of precipitation raster data sets
precipitation_datasets <- list(
  CCSM_Past = "Prec_CCSM_Past",
  CNRM_Past = "Prec_CNRM_Past",
  FGOA_Past = "Prec_FGOA_Past",
  IPSL_Past = "Prec_IPSL_Past",
  MIRO_Past = "Prec_MIRO_Past",
  MPIE_Past = "Prec_MPIE_Past",
  MRIC_Past = "Prec_MRIC_Past",
  Prec_Pres = "Prec_Pres"
)

# Create a list of temperature raster data sets
temperature_datasets <- list(
  CCSM_Past = "Temp_CCSM_Past",
  CNRM_Past = "Temp_CNRM_Past",
  FGOA_Past = "Temp_FGOA_Past",
  IPSL_Past = "Temp_IPSL_Past",
  MIRO_Past = "Temp_MIRO_Past",
  MPIE_Past = "Temp_MPIE_Past",
  MRIC_Past = "Temp_MRIC_Past",
  Temp_Pres = "Temp_Pres"
)

# Function to plot histograms using base R hist function
plot_histograms_base <- function(datasets, title) {
  par(mfrow=c(4, 2), mar=c(3,3,1,1))
  
  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]
    raster_data <- get(dataset)  # Get raster object
    
    # Extract raster values
    raster_values <- getValues(raster_data)
    
    # Filter out NA and non-finite values
    raster_values <- raster_values[!is.na(raster_values) & is.finite(raster_values)]
    
    # Plot histogram using base R hist function
    hist(raster_values, col = hist_colors_chel[i], main = paste(dataset), xlab = "Value", ylab = "Frequency")
  }
}

# Plot precipitation histograms using base R hist function
plot_histograms_base(precipitation_datasets)

# Plot temperature histograms using base R hist function
plot_histograms_base(temperature_datasets)





for (dataset_name in names(precipitation_datasets)) {
  dataset <- get(precipitation_datasets[[dataset_name]])
  print(summary(dataset))
}


# Plot precipitation rasters
par(mfrow=c(4, 2))  # Set up a 2x4 plotting layout
for (i in 1:length(precipitation_datasets)) {
  plot(precipitation_datasets[[i]], main = names(precipitation_datasets)[i])
}

# Plot temperature rasters
par(mfrow=c(4, 2))  # Set up a 2x4 plotting layout
for (i in 1:length(temperature_datasets)) {
  plot(temperature_datasets[[i]], main = names(temperature_datasets)[i])
}

par(mfrow=c(4, 2))  # Set up a 2x4 plotting layout

plot(Prec_CCSM_Past)
plot(Prec_CNRM_Past)
plot(Prec_FGOA_Past)
plot(Prec_IPSL_Past)
plot(Prec_MIRO_Past)
plot(Prec_MPIE_Past)
plot(Prec_MRIC_Past)
plot(Prec_Pres)


par(mfrow=c(4, 2))  # Set up a 2x4 plotting layout


plot(Temp_CCSM_Past)
plot(Temp_CNRM_Past)
plot(Temp_FGOA_Past)
plot(Temp_IPSL_Past)
plot(Temp_MIRO_Past)
plot(Temp_MPIE_Past)
plot(Temp_MRIC_Past)
plot(Temp_Pres)

## 3.Prec_CCSM_Past## 3.2  Plot histograms of GMBApoly data--------------------------------------------------------------------------------------------------------------------------


print(names(GMBApoly)[variablesHis])

Vari_Num <- c(4, 6, 26:41) # Define the column indices of interest Max, mean, uplift, area
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



## 3,3 Violin plots of GMBApoly data ----------------------------------------------------------------------------------------------------------------------------
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




## 3.4 summarizing results PCA of the M---------------------------------------------------------------------------------------------


PCAData <- GMBApoly[,c(4, 6, 33:41)] #taking only the numeric values 

res.pca <- prcomp(PCAData, scale = TRUE) #rescaled variables, which makes that they're scaled according to eachother and limiting effect of large ranges compared with smaller.

if (anyNA(PCAData)) {
  PCAData <- na.omit(PCAData)  # Remove rows with missing values, as the PCA doesn't work with that
}

# Visualize PCA results
pca_plot <- fviz_pca(res.pca, geom = "point", ggtheme = theme_minimal())

pca_plot








