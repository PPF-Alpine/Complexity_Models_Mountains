---
title: "Readme.md"
output: html_document
date: "2024-03-18"
---

# Readme

This repository contains the R script and data necessary for the completion of the graduation thesis. The provided resources are essential for conducting analysis, generating results, and drawing conclusions as part of the academic project.

### **Contents**

-   **R Script:** The repository includes an R script file (**`thesis_analysis.R`**) that contains the code for data processing, statistical analysis, visualization, and any other computational tasks required for the thesis.

**Data:** The necessary datasets for the thesis are provided in the **`data/storage`** directory. These datasets serve as the foundation for the analysis and findings presented in the thesis.

## Datasources

- ==Dyy_250m_MERIT==: contains the second derivative direction slope (N-S)
	- The dyy, or North-South second order partial derivative, measures the rate of change of slope in a North-South direction. It represents the curvature or concavity of a function when applied to a Digital Elevation Model (DEM). Positive values indicate convex surfaces, while negative values 	indicate concave surfaces. The unit of curvature is radians per metre.
	- This has a resolution of 7,5 arc-seconds (~250m), under WGS84,
	- Found on [250m amatulli database](https://doi.pangaea.de/10.1594/PANGAEA.899135) , 
- ==dxx_250m_MERIT==: contains the second derivative direction slope (E-W)
	- the dxx, or East-West second order partial derivative, measures the derivative of slope in an East-West direction. It also represents the curvature or concavity of a function when applied to a DEM. Positive values denote convex surfaces, while negative values denote concave surfaces. The 	unit of curvature remains radians per metre.
	- This has a resolution of 7,5 arc-seconds (~250m), under WGS84,
	- Found on [250m amatulli database](https://doi.pangaea.de/10.1594/PANGAEA.899135) , 
- ==VRM_250m_Merit==: contains summarized ruggedness not sure what the unit is of this variable
	- Vector Ruggedness 250m Measure: quantifies terrain ruggedness by measuring the dispersion of vectors orthogonal to the terrain surface. Provides a measure of local variation in slope, providing information about the variability of terrain ruggedness within the landscape. It is dimensionless 	because of sine-cosine derivation, and values range from 0 to 1 in flat to rugged regions, respectively [[Vector ruggedness explanation|drawing]] [[amatulli2020]].
	- This has a resolution of 7,5 arc-seconds (~250m), under WGS84,
 	- found [250m amatulli database](https://doi.pangaea.de/10.1594/PANGAEA.899135) , 

- ==Slope_250m_MERIT==: rate of change of elevation expressed in degrees or percentages
	- Terrain slope (slope) or rate of change of elevation is the measure of steepness in the direction of the water flow line. It is considered one of the most important terrain parameters and is often calculated first. It can be expressed in degrees or percentages, where for example, 5% means 5-m of vertical displacement over 100-m. It is especially importan for the quantification of soil erosion, water flow velocity, or agricultural suitability. 
	- This has a resolution of 7,5 arc-seconds (~250m), under WGS84,
	- Found on [250m amatulli database](https://doi.pangaea.de/10.1594/PANGAEA.899135) , 

- ==dx_250M_MERIT==: contains directional slope (E-W) not sure what the unit is of this variable
	- the dx, or East-West first order partial derivative, measures the **slope** in an **East-West** direction. It plays a role in estimating overland water flow and sediment flow using models like SIMWE. Additionally, it aids in identifying artifacts like voids, pits, sinks, and sensor stripes in DEMs due to its sensitivity to systematic noise and artifacts.
	- This has a resolution of 7,5 arc-seconds (~250m), under WGS84, 
	- found [250m amatulli database](https://doi.pangaea.de/10.1594/PANGAEA.899135) ,

- ==dy_250M_MERIT==: contains directional slope (N-S) not sure what the unit is of this variable
	- The dy, or North-South first order partial derivative, measures the **slope** in a **North-South** direction. It is useful for estimating overland water flow and sediment flow using models like SIMWE. Furthermore, it helps detect artifacts such as voids, pits, sinks, and sensor stripes in DEMs, as it is sensitive to systematic noise and artifacts.
	- This has a resolution of 7,5 arc-seconds (~250m), under WGS84, found [250m amatulli database](https://doi.pangaea.de/10.1594/PANGAEA.899135) ,

- ==dtm_lithology_usgs== (15 classes): contains lithology classes
	- from [Say et al. 2014](https://www.aag.org/wp-content/uploads/2021/12/AAG_Global_Ecosyst_bklt72.pdf), the same classes as rock types above with a variation on the names and Water Bodies (WB),
 	- the data can be collected from [database zenodo](https://zenodo.org/records/1464846#.Xn3P40p7lPY) USGS. used [[sayre2014]]

- ==LandForm_NewHammond2017== (16 classes): contains landform types classes
	- This variable contains 16 landform classes provided by [[karagulle2017]]. 
	- These landforms are classified, [[landformHammond_classes]]
	- The dataset is found in [ArcGIS layer for Hammond2017_Landforms](https://www.arcgis.com/home/item.html?id=cd817a746aa7437cbd72a6d39cdb4559). 
	- It uses GMTED2010 [[DEM]] and a WGS 1984 geographic coordinate system to World Equidistant Cylindrical, 
	- with a resolution of 250m.


- ==Mountain uplift data:== contains mountain uplift starting history (MYR)
	- from [[quintero2018]] in MYR, which is added to various mountains where data was found on the uplift time in history. 
	- Dataset found at [uplift sheet](https://static-content.springer.com/esm/art%3A10.1038%2Fnature25794/MediaObjects/41586_2018_BFnature25794_MOESM5_ESM.xlsx).

- ==GMBA.V2 level 4 Area==: in the GMBA layer they also calculated the area size of the different polygons. Calculated planimetric area of the mountain polygon (in km2) (calculated in Mollweide projection) [[snethlage2022]]. 
  
