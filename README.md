# HSMs_KBAs_Canada
## OVERVIEW
Screening potential areas for <a href="http://www.kbacanada.org/" target="_blank">Key Biodiversity Areas (KBAs) in Canada</a>, using Habitat Suitability Models (HSMs) for species that might trigger A1 (threatened species) and B1 (geographic restricted species) KBA criteria  (see  and <a href="https://portals.iucn.org/library/node/46259" target="_blank">KBA standard</a>).

To do that, we first evaluated habitat suitability models relating species presence-only data sets and predictor variables, using the  <a href="https://www.sciencedirect.com/science/article/pii/S030438000500267X" target="_blank">Maxent algorithm</a>  implemented in the  ENMeval v.2.0.0 package (<a href="https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13628?campaign=woletoc" target="_blank">Kass et al 2021</a>) in R. Then, we used cluster analysis to identify discrete polygons containing most of the high suitable areas.

We ran the analysis for a group of 96 species as a proof of concept aimed to provide recommendations for a full implementation (~ 1,200 species). 
