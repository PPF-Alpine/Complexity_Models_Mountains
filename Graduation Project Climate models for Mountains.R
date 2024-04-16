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

wd <- ("C:/Users/repap5991/Downloads/Data/R/Rdocs")
setwd(wd) #setting the working directory

# 1.0 loading data ---------

