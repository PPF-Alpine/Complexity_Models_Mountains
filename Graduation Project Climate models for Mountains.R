# 0 setup --------------

rm(list = ls()) #full reset

install.packages("dplyr") # installing packages
install.packages( c( "dplyr", "ggplot2", "plyr" ) )

## general
library("dplyr")    #activating the package for running code for it.
library("ggplot2")
library("plyr")

library("stringr")
library("ReacTran")
library("fields")
library("animation")

## project specific

library("terra")
library("sf")
library("raster")

wd <- ("C:/Users/repap5991/Downloads/Data/R/Rdocs/")
setwd(wd) #setting the working directory

# 1.0 loading data ---------

DEM <- raster(paste0(wd, "DataStorage/DEM/dem_geo.sdat"))  ## adding dem_geo chelsa dem

GMBApoly <- st_read(paste0(wd, "DataStorage/ArcGIs Data/GMBA_Level4_ExportTable_TableToExcel.csv"))  ## arcgis dataset which contains all the variables according to the GMBA_level4.


# 2.0 data type fixing ------

str(GMBApoly) ## overview of data

numeric_dataconversion <- c("Area", 
                            "MAX_VRM",
                            "MEAN_VRM",
                            "MAX_Slope",
                            "MEAN_Slope",
                            "MAX_Dy",
                            "MEAN_Dy",
                            "Uplift_Sta",
                            "MAX_dx",
                            "MEAN_dx",
                            "MAX_DEM",
                            "MAX_Dxx",
                            "MAX_Dyy",
                            "MEAN_DEM",
                            "MEAN_Dxx",
                            "MEAN_Dyy")

Catagorical_dataconversion <- c("MAJ_LaFor",
                                "Level_04",
                                "MAJ_LaForC",
                                "MAJ_LaForP",
                                "VarC_LaFor",
                                "VAR_LaFor")


GMBApoly[numeric_dataconversion] <- lapply(GMBApoly[numeric_dataconversion], as.numeric)

GMBApoly[Catagorical_dataconversion] <- lapply(GMBApoly[Catagorical_dataconversion], as.factor)


str(GMBApoly) ## overview of data

# 3.0 Basic statistics -----------------------------

ggplot(GMBApoly, aes(x = GMBApoly$VAR_LaFor)) +
  geom_bar(stat = "count")






