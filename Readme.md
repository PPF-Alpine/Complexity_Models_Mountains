---
title: "Readme.md"
output: html_document
date: "2024-03-18"
---

# Readme

This repository contains the R script and data necessary for the completion of the graduation thesis. The provided resources are essential for conducting analysis, generating results, and drawing conclusions as part of the academic project.

### **Contents**

-   **R Script:** The repository includes an R script file (**`thesis_analysis.R`**) that contains the code for data processing, statistical analysis, visualization, and any other computational tasks required for the thesis.

-   **Data:** The necessary datasets for the thesis are provided in the **`data/storage`** directory. These datasets serve as the foundation for the analysis and findings presented in the thesis.

- Dyy_250m_MERIT: contains the second derivative direction slope (N-S)

- dxx_250m_MERIT: contains the second derivative direction slope (E-W)

- VRM_250m_Merit: contains summarized ruggedness not sure what the unit is of this variable
	- Vector Ruggedness 250m Measure: quantifies terrain ruggedness by measuring the dispersion of vectors orthogonal to the terrain surface. Provides a measure of local variation in slope, providing information about the variability of terrain ruggedness within the landscape. It is dimensionless because of sine-cosine derivation, and values range from 0 to 1 in flat to rugged regions, respectively [[Vector ruggedness explanation|drawing]] [[amatulli2020]].
  - This has a resolution of 7,5 arc-seconds (~250m), under WGS84, found [250m amatulli database](https://doi.pangaea.de/10.1594/PANGAEA.899135) , 

- Slope_250m_MERIT: rate of change of elevation expressed in degrees or percentages
	- Terrain slope (slope) or rate of change of elevation is the measure of steepness in the direction of the water flow line. It is considered one of the most important terrain parameters and is often calculated first. It can be expressed in degrees or percentages, where for example, 5% means 5-m of vertical displacement over 100-m. It is especially importan for the quantification of soil erosion, water flow velocity, or agricultural suitability [18].
 - This has a resolution of 7,5 arc-seconds (~250m), under WGS84, found [250m amatulli database](https://doi.pangaea.de/10.1594/PANGAEA.899135) , 

- dx_250M_MERIT: contains directional slope (E-W) not sure what the unit is of this variable
	- First order partial derivative E-W: rate of change of the **slope** in certain aspect (direction E-W).
 - This has a resolution of 7,5 arc-seconds (~250m), under WGS84, found [250m amatulli database](https://doi.pangaea.de/10.1594/PANGAEA.899135) ,

- dy_250M_MERIT: contains directional slope (N-S) not sure what the unit is of this variable
	- First order partial derivative N-S: rate of change of the **slope** in certain aspect (direction N-S) y-axis.
 - This has a resolution of 7,5 arc-seconds (~250m), under WGS84, found [250m amatulli database](https://doi.pangaea.de/10.1594/PANGAEA.899135) ,

- dtm_lithology_usgs (15 classes): contains lithology classes
	- from [Say et al. 2014](https://www.aag.org/wp-content/uploads/2021/12/AAG_Global_Ecosyst_bklt72.pdf), the same classes as rock types above with a variation on the names and Water Bodies (WB),
 - the data can be collected from [database zenodo](https://zenodo.org/records/1464846#.Xn3P40p7lPY) USGS. used [[sayre2014]]

- LandForm_NewHammond2017 (16 classes): contains landform types classes
	- This variable contains 16 landform classes provided by [[karagulle2017]]. It uses GMTED2010 [[DEM]] and a WGS 1984 geographic coordinate system to World Equidistant Cylindrical, with a resolution of 250m. 
	- The dataset is found in [ArcGIS layer for Hammond2017_Landforms](https://www.arcgis.com/home/item.html?id=cd817a746aa7437cbd72a6d39cdb4559). 


- Mountain uplift data: contains mountain uplift starting history (MYR)
	- from [[quintero2018]] in MYR, which is added to various mountains where data was found on the uplift time in history.
 - Dataset found at [uplift sheet](https://static-content.springer.com/esm/art%3A10.1038%2Fnature25794/MediaObjects/41586_2018_BFnature25794_MOESM5_ESM.xlsx).

- CHELSA DEM: gives the elevation not sure what the unit is of this variable
	- "dem_latlong.sdat" [[CHELSA variables]], not 100% sure if these are based on GMTED2010 but in the data explanation they do refer to it a lot in terms of input data for calculations of various variables for they're bioclimate variables for instance.'
 - dataset can be found [DEMdataChelsa](https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F)
 - 
