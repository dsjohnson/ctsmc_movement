
library(sp)
library(sf)
library(crawl)
library(lubridate)
library(tidyverse)
library(spdep)
library(raster)
library(ncdf4)
library(ncdf4.helpers)
library(mgcv)

### Necessary package for fir seal analysis, can be omitted for others
### comment the line out if not needed, else
### install with-- 

# devtools::install_github('jmlondon/npacmaps')
library(nPacMaps) 

### At the time this code development, the development version of ggplot2 was necessary for
### the geom_sf() function 
### install with--

# devtools::install_github('tidyverse/ggplot2')
library(ggplot2)
