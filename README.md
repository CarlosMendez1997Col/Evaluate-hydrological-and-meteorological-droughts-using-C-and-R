# Evaluation and analysis of drought using the SPI, SPEI and PDSI index

## Use and install this repository

HTTPS
```CSS
https://github.com/CarlosMendez1997Col/Evaluate-hydrological-and-meteorological-droughts-using-C-and-R.git
```

GitHub CLI
```CSS
gh repo clone CarlosMendez1997Col/Evaluate-hydrological-and-meteorological-droughts-using-C-and-R
```

## Description

Repository builded in `C++` and `R` languages programs. 

Each section is described below:

1. The first section, calculate the `SPI` and `SPEI` index.
```HTML
https://github.com/CarlosMendez1997Col/Evaluate-hydrological-and-meteorological-droughts-using-C-and-R/tree/main/1.%20The%20SPEI%20and%20SPI%20index
```
2. The second section, calculate the `PDSI` index.
```HTML
https://github.com/CarlosMendez1997Col/Evaluate-hydrological-and-meteorological-droughts-using-C-and-R/tree/main/2.%20Palmer%20Drought%20Severity%20Index%20(PDSI)
```
3. The third section, show some animations and videos of drought.
```HTML
https://github.com/CarlosMendez1997Col/Evaluate-hydrological-and-meteorological-droughts-using-C-and-R/tree/main/3.%20Animations%20and%20videos
```
4. The fourth section, share the data and files used.
```HTML
https://github.com/CarlosMendez1997Col/Evaluate-hydrological-and-meteorological-droughts-using-C-and-R/tree/main/4.%20Data%20and%20files
```
5. The five section, plot maps and graphs of drought.
```HTML
https://github.com/CarlosMendez1997Col/Evaluate-hydrological-and-meteorological-droughts-using-C-and-R/tree/main/5.%20Maps%20and%20graphs
```
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
## Versions and releases

Version `1.0`

```HTML
https://github.com/CarlosMendez1997Col/Evaluate-hydrological-and-meteorological-droughts-using-C-and-R/releases
```

## Data acquisition and download

### PDSI 

The data used to calculate the PDSI are from stations and time series collected during the following research: [Link](https://repository.udistrital.edu.co/items/de3ecda1-01ec-4203-a938-1969240f6d24)

### SPEI and SPI 

The data used to calculate SPI and SPEI (precipitation and potential evapotranspiration) were downloaded from the Nasa Earth repository, using the Giovanni program (Geospatial Interactive Online Visualization ANd aNalysis Infrastructure), available at: [Giovanni Website](https://giovanni.gsfc.nasa.gov/giovanni/)

## Credits and acknowledgments

The SPI and SPEI index, were used with the CRAN package `SPEI Version: 1.8.1` in R language, available at: [Oficial Website](https://cran.r-project.org/web/packages/SPEI/index.html) and [GitHub Repository](https://github.com/sbegueria/SPEI)

The PDSI index, was modified following the original code, available at: [GitHub Repository](https://github.com/cszang/pdsi/blob/master/exec/scpdsi.cpp) 

## Conflict of Interest.

The author declare that there is no conflict of interest in the publication of this data and have approved it for publication.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate. 

## MIT License

Copyright (c) 2025 Carlos Mendez

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.




