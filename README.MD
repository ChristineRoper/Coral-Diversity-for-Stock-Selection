# Stock selection for coral restoration on the Great Barrier Reef: Does phenotypic diversity align to population genetic structure?

Assessment of coral photophysiological phenotypic performance  via a rapid, high-thoroughput phenotyping device, alongside analysis of genetic structure to determine if phenotypic diversity and population structure align.

Version 1

# Files & Folders
## `Phenotype Analysis`
Data and code for determining and visualising the QYmax at baseline temperature and ED50 of rETRmax under heat stress.  
See `Phenotype Analysis/README` for overview of files and instructions

## `Genetic Analysis`
Data and code for SNP quality control, filtering, visualisation and analysis to infer population structure.  
See `Genetic Analysis/README` for overview of files and instructions

## `Temperature`
Visualise and calculate temperature profiles from AQUA MODIS satellite data for each reef and combine with in-situ logger data. 

* Download AQUA MODIS Sea Surface Temperatur,e 2021-01-01 to 2021-01-08, Level 3, 8 Day into `AQUA_MODIS.L3m.8D.SST.x_sst`   
https://podaac.jpl.nasa.gov/dataset/MODIS_AQUA_L3_SST_MID-IR_DAILY_4KM_NIGHTTIME_V2019.0
* Run `Temperature.r` to generate plots and statistics

