# Fractional Richness 
An index for large camera trap networks

Fractional Richness is an index designed to provide intuitive and ecologically informative measures of biodiversity when working with large camera trap networks.
Camera trap networks have a number of inherent biases which make typical diversity indexes like Shannon diversity less informative. E.g. not all species are equally detectable, detections follow a poisson distribution rather than a normal distribution, and species of different body sizes or trophic levels have different maximum population densities. 
Fractional Richness accounts for these idiosycrasies of camera trap data.

**Language**: R



## Index Behavior
The Fractional Richness index behaves differently than other diversity indexes (e.g. Shannon diversity) in several key ways:
1. **All species are valued equally**. With Shannon Diversity, low-density or rare species have less impact on the overall diversity score. With Fractional Richness, all species contribute equally to the diversity score.
2. **Diversity increases linearly with species richness**. Shannon Diversity has a plateauing curvilinear distribution, so that the resolution is reduced at higher species richness. Fractional Richness has a linear relationship with richness, and performs equally well at high and low species richness. When all species are at their average population densities, Fractional Richness is equal to ½ species richness. When all species are at their maximum population densities, Fractional Richness is equal to species richness. 
3. **An increase in richness or abundance always results in an increase in Fractional Richness**. With Shannon diversity, if the abundances of species are very uneven, then a decrease in population or a complete loss of the most common species can counterintuitively result in a higher measured diversity. This is never the case with Berman diversity.
4. **Population decline results in a loss of diversity**. If all species become scarce but evenness remains the same, Fractional Richness decreases while Shannon diversity remains constant. 
5. **Species detectability does not bias the index**. Practically speaking, the measured abundance of a species, or its detection rate, is never exactly the same as the true abundance. Some species are more easily detectable than others. Because Fractional Richness values all species equally, and dsite,i and dmax,i are affected by the same set of biases, species detectability bias is effectively canceled out and does not affect the index, so long as a species that is present is detected. 
6. **Multiple sites are required**. The main drawback of Fractional Richness when compared to Shannon diversity is that it requires more information. Shannon diversity can be calculated for a single site, while Fractional Richness requires either multiple sites or multiple time points in order to calculate dmax and mean fractional population density.



# Data requirements
To calculate Fractional Richness, detection rates for multiple species at multiple sites are required. 
Each row should contain data for one site and each column should contain detection rates for one species.

For example:

ID    | Species 1  | Species 2  | Species 3
----- | ---------- | ---------- | ---------
Site1 |  0.2       |  8         |  0.001
Site2 |  0.3       |  7         |  0.007
Site3 |  0.1       |  9         |  0.004




# Equation

$FRichness_{site} = \sum_{i=1}^{n}\left ( \frac{d_{site,i}}{d_{max,i}} \right )^{\frac{log0.5}{log\overline{\left ( \frac{d_{site,i}}{d_{max,i}} \right )}}}$

Creating the following 4 functions in R allows you to calculate Fractional Richness:

```{r}
FPD <- function(species_col) {species_col / max(species_col)}
meanFPD <- function(species_col) {mean(FPD(species_col)[FPD(species_col)>0])}
NFPD <- function(species_col) {(FPD(species_col))^(log(0.5)/log(meanFPD(species_col)))}

FRichness <- function(data, species_cols, ID_col){
  data.frame(ID = data[[ID_col]],
             FRich = rowSums(sapply(data[species_cols], NFPD), na.rm=TRUE))
}
```

**FPD = Fractional Population Density**

The detection rate of a species at a site as a fraction of the maximum detection rate of that species across all sites. 



**meanFPD = mean Fractional Population Density**

A single value for each species. A measure of the skewness of the data distribution. 
The mean fractional population density of all sites where the species is present, excluding sites where the species has not been detected.



**NFPD = Normalized Fractional Population Density**

Fractional population density normalized so that the distribution of detection rates centers around 0.5. 
Normalization is species specific. For each species, sites with an average detection rate will have a NFPD of 0.5. 



**FRichness = Fractional Richness**

The sum of NFPD for each species. At sites where all species are at their typical population densities, FRichness = 1/2 species richness.



**Arguments**
 + data:             a dataframe containing species detection rates
 + ID_col:           column containing site ID
 + species_cols:     columns which contain species detection rates, one column per species
 + species_col:      a single column containing species detection rates



# Usage

Calculate Fractional Richness for a selected community of species and visualize using mapview:

```{r}
library(sf)
library(mapview)

#load dataset
DBC <- read.csv("DBC_clean.csv")

#Calculate Fractional Richness for only sensitive species
HSS <- FRichness(data=DBC,
                        species_cols = c(5:7, 15, 16, 17, 26, 33), 
                        ID_col = "CAMERA_ID")

#bind diversity values to original dataframe
HSS <- cbind(DBC, HSS)

#visualize
HSS <- st_as_sf(HSS, coords = c("LON", "LAT"), crs=4326)
mapview(HSS, zcol="FRich")

```

Intermediate steps can also be calculated independently:

```{r}
x <- FPD(DBC$BearsPerDay)

x <- meanFPD(DBC$BearsPerDay)

x <- NFPD(DBC$BearsPerDay)
```


# Provided dataset

An example dataset is available in this github which can be used to calculate Fractional Richness: Snapshot_DBC.csv

Snapshot_DBC.csv is data from the Snapshot Wisconsin camera trap network (https://dnr.wisconsin.gov/topic/research/projects/snapshot). 
This dataset includes detection rates of 33 species across 2218 camera locations in Wisconsin.
Detection rates are measured as the number of detection events / the number of days the camera location was active.
LAT/LON coordinates are coarsened for the privacy of community scientists in this publicly available dataset.

**Variables**
  + CAMERA_ID: Unique ID number of each camera location
  +	LAT: Latitude in decimal degrees, coarsened to 2 decimal places to maintain privacy of community scientists
  +	LON: Longitude in decimal degrees, coarsened to 1 decimal place to maintain privacy of community scientists
  + *Species*PerDay: One column per species. The number of times that species was detected / the number of days that camera site was active.

**Species**
 + Bear - *Ursus americanus*
 + Beaver - *Castor canadensis*
 + Bobcat - *Lynx rufus*
 + Cottontail - *Sylvilagus floridanus*
 + Coyote - *Canis latrans*
 + Deer - *Odocoileus virginianus*
 + Fisher - *Pekania pennanti*
 + Grey fox - *Urocyon cinereoargenteus*
 + Grouse - *Bonasa umbellus*
 + Mink - *Neovison vison*
 + Opossum - *Didelphis virginiana*
 + Porcupine - *Erethizon dorsatum*
 + Raccoon - *Procyon lotor*
 + Red fox - *Vulpes vulpes*
 + Sandhill crane - *Antigone canadensis*
 + Squirrel - *Sciuridae*
 + Turkey - *Meleagris gallopavo*
 + Wolf - *Canis lupus*
 + Woodchuck - *Marmota monax*


# Citation

For more details, or if you use Fractional Richness in your work, please cite:

Berman, L. M., Schneider, F. D., Pavlick, R. P., Stenglein, J., Bemowski, R., Dean, M., & Townsend, P. A. (2024). Fractional Richness: An index for camera trap networks. Ecological Indicators, 166, 112266.



