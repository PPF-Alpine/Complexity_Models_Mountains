# 0 setup --------------

rm(list = ls()) #full reset

install.packages("dplyr") # installing packages

## general
library(dplyr)    #activating the package for running code for it.
library(ggplot2)
library(stringr)
library(ReacTran)
library(fields)
library(animation)

## project specific

library("terra")
library("sf")
library("raster")

wd <- ("C:/Users/repap5991/Downloads/Data/R/Rdocs/")
setwd(wd) #setting the working directory

# 1.0 loading data ---------

DEM <- raster(paste0(wd, "DataStorage/DEM/dem_geo.sdat"))  ## adding dem_geo chelsa dem

GMBApoly <- st_read(paste0(wd, "DataStorage/ArcGIs Data/GMBA_Level4_ExportTable_TableToExcel.csv"))  ## arcgis dataset which contains all the variables according to the GMBA_level4.





