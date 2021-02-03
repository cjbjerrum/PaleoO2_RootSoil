# PaleoO2_RootSoil
PaleoO2_RootSoil model calculates root depth of (fossil) land plant associated with atmospheric oxygen levels, for soil types of non-wetland areas.  
The model calculates the oxygen diffusion into the soil, and the O2 consumption by plant root respiration as well as microbial respiration in the water film around roots and in the soil air space. 
As in Bartholomeus et al. (2008), we split the calculation into plant root O2 requirement at the macro scale, and at the micro root scale down through the soil seeking solutions so that O2 is greater than zero in the root center permitting optimal growth conditions or only survival.

The model output:
PaleoO2_RootSoil calculates soil-gas phase oxygen profile, the plant root requirement for optimal growth (eta=5), as well as the maximum root depth for optimal growth:

The right solid curve is the soil- gas phase oxygen profile, a function of: 1) macroscale root respiration and 2) microbial respiration dispersed in the soil, both with a Q10 dependence on T; and 3) the vertical diffusion of oxygen that in turn depends on soil porosity and its interconnectedness, and soil moisture.

The left dashed curve is the oxygen requirement for optimal root respiration at the micro-scale. This represents the oxygen concentration required at the interface of the gas phase of the soil and the capillary water-film around root hairs. This required oxygen level is where there is sufficient to supply oxygen to cells in the root so that there is non-zero oxygen in the root center, as well as ensuring microbial respiration in the water-film (thus implicitly including generalized mycorrhizae fungi respiration). The respiration oxygen requirement curve is a function of microscale root respiration (in turn with a Q10 dependence on T), root water film thickness where in microbial respiration is a function of  organic carbon and T, and O2 diffusivity in water a function of T. The root water-film thickness increases down through the soil as soil moisture increase. This causes the respiration requirement curve to increase with depth (cf. Bartholomeus et al., 2008; Campbell, 1985; Cook, 1995; Cook et al., 2013; de Willigen and van Noordwijk, 1987). 

Where the solid and dashed calculated O2 curves intersect is the maximum root depth that permits optimal tree growth uninhibited by O2 at a given atmospheric O2 level.

As final diagram displays the maximum root depth that is permissible for optimal growth as function of atmospheric oxygen.

Please cite model as: Christian J. Bjerrum, 2021. cjbjerrum/PaleoO2_RootSoil: First stable release of PalO2_RootSoil, doi: 10.5281/zenodo.4495296.
