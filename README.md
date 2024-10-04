# Drought Index (SPI, SPEI and PDSI)

## Description

Evaluation of the `SPI`, `SPEI` and `PDSI` , using `C++` and `R` languages. 

Each section is described below:

1. The first section shows the SPI and SPEI index

2. The second section explain the PDSI index

3. The third section include animations and videos

4. The fourth and final section share the data and files

## Prerequisites and libraries

### PDSI in C++

```C++
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include <stdio.h>
#include <ctype.h>
#include <direct.h> 

```
### SPI and SPEI in R

```R
#install.packages("SPEI")
#install.packages("readr")
#install.packages("gganimate")
#install.packages("tidyverse")
library(SPEI)
library(tidyverse)
library(gganimate)
library(readr)
library(dplyr)

```

## Data acquisition and download

### PDSI 

The data used to calculate the PDSI are from stations and time series collected during the following research: [Link](https://repository.udistrital.edu.co/items/de3ecda1-01ec-4203-a938-1969240f6d24)


### SPEI and SPI 

The data used to calculate SPI and SPEI (precipitation and potential evapotranspiration) were downloaded from the Nasa Earth repository, using the Giovanni program (Geospatial Interactive Online Visualization ANd aNalysis Infrastructure), available at: [Giovanni Website](https://giovanni.gsfc.nasa.gov/giovanni/)

## Credits and more information

The SPI and SPEI index, were used with the CRAN package `SPEI Version: 1.8.1` in R language, available at: [Oficial Website](https://cran.r-project.org/web/packages/SPEI/index.html) and [GitHub Repository](https://github.com/sbegueria/SPEI)

The PDSI index, was modified following the original code, available at: [GitHub Repository](https://github.com/cszang/pdsi/blob/master/exec/scpdsi.cpp) 
## Conflict of Interest.

The author declare that there is no conflict of interest in the publication of this data annd have approved it for publication.

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate. 
