## Overview

This package generates quantitatively realistic atmospheric profiles for a given time and season. In particular, it generates the water vapor mass density, 
the air temperature, and a two-dimensional wind vector at 100m vertical increments up to 20 km.

## Methodology

weathergen randomizes the first 32 principal components of the atmospheric profile variations, fitted to each region. Eigenmodes are fitted with to the ERA5 reanalysis dataset from 1980 to 2021. The model has a diurnal resolution of an hour 
and a seasonal resolution of around two weeks. 

## Supported regions

The current version supports three regions: the Chajnantor plateau in the Atacama Desert ('chajnantor'), the South Pole ('south_pole'), and Ngari Prefecture in Tibet ('tibet'). 
