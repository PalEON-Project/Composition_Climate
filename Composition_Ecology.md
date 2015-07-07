---
author: 'Simon Goring *et al*.'
biblio-files: 'goringetal_references'
bibliography: 'goringetal_references.bib'
csl: 'ecology.csl'
date: '07 July, 2015'
output:
  md_document:
    variant: markdown
  pdf_document:
    pandoc_args: '-V geometry:vmargin=1in -V geometry:hmargin=1in'
pandoc_args: '-H margins.sty'
title: 'Composition-climate paper'
...

Land use patterns and historical contingency change climate optima for key forest genera in the Northeastern United States.
===========================================================================================================================

<!--


```
## Warning: package 'pander' was built under R version 3.2.1
```
-->
*Historic land use and forest change over the centuries since
Euroamerican settlement have significantly changed the structure and
composition of forests in the upper Midwestern United States. The extent
to which this change has affected the realized climate niche for key
taxa is less well known.*

*Predictions for future species range shifts are predicated largely on
correlative models that relate modern species distributions to climate
parameters in the modern era. Given the extent of regional forest change
and recent losses to key forest taxa including Hemlock, Elm and Chestnut
it is possible that certain taxa have shiften with respect to climate
change over the last century. The additional pressure of climatically
biased land use conversion for agriculture along the western and
southern border of the upper Midwest means that the shift in climate
space may be most pressing in a region of climate space that is
critically important for understanding future species responses to
climate change.*

      OGR data source with driver: ESRI Shapefile 
      Source: "data/input/shapes/usa/us.shp", layer: "us"
      with 51 features
      It has 51 fields
      OGR data source with driver: ESRI Shapefile 
      Source: "data/input/shapes/canada/PROVINCE.SHP", layer: "PROVINCE"
      with 13 features
      It has 2 fields

      OGR data source with driver: ESRI Shapefile 
      Source: "data/input/NaturalEarth/Lakes/ne_10m_lakes.shp", layer: "ne_10m_lakes"
      with 1353 features
      It has 10 fields
      OGR data source with driver: ESRI Shapefile 
      Source: "data/input/NaturalEarth/Rivers/ne_10m_rivers_lake_centerlines.shp", layer: "ne_10m_rivers_lake_centerlines"
      with 1454 features
      It has 7 fields
      OGR data source with driver: ESRI Shapefile 
      Source: "data/input/NaturalEarth/Coasts/ne_10m_coastline.shp", layer: "ne_10m_coastline"
      with 4132 features
      It has 2 fields

Introduction
============

Over the past 200 years centennial scale climate variability and land
use change have interacted to transform the forests of the northeastern
United States. There are clear indications that modern forests, as
represented by the Forest and Inventory Analysis show greater
homogeneity and significantly different species associations than
pre-settlement forests [@goring2015preset; @schulte2007homogenization].
This shift has occured over a period during which climate has begun to
show the effects of anthropogenic warming, and a period in which we have
increasingly used the relationship between modern tree species
distributions and climate to understand future responses of tree species
and forest assemblages to climate change
[@iverson2013tree; @iverson1998predicting]. One of the most significant
changes to forest tree distributions in the upper midwest is the loss of
extensive cover through land use change [@rhemtulla2009legacies],
particularly the .

``` {.r}
natural[[2]] + geom_tile(data = land_use, aes(x = x, y = y, fill = class)) + 
    scale_fill_manual(values = c("#ff7f00", "#e41a1c", "#33a02c", "#377eb8"))
```

![](Composition_Ecology_files/figure-markdown/unnamed-chunk-2-1.png)
**Figure 1.** *Patterns of land use change in the upper Midwestern
United States. Data from the NLCD [@jin2013comprehensive].*

If modern species distributions have been affected by land use change in
a random fashion then we should expect little effect on the resulting
predictions since, by and large, me might expect the climate envelope
encompased by the species to remain largely unchanged. In the midwestern
United States we see a strong gradient of temperature and precipitation
that results in two major ecotones, one from prairie/savanna to closed
forest, in an approixmately west to east direction, and one from conifer
dominated, sub-boreal forests to decidious forest in a north to south
direction. Goring et al [@goring2015preset] have shown that land use
change has significantly changed the strength and structure of this
ecotone, but did not examine how these shifts might affect the
relationships between climate and distribution for individual taxa.
Given that much of the land use change within the region occurs in the
south, where agricultural conversion has largely eliminated open
forests, we might expect to see that species with more southerly
distributions would show greater impacts of land use change on their
distributions in climate space, assuming no shift in regional climate.
Conversely, species in the north should show little change if we expect
that pre-settlement trees show similar regeneration patterns (with
respect to simple presence/absence) following widespread logging.

We use gridded climate data products and estimates of pre-settlement
vegetation to develop climate-vegetation relationships for 15 major tree
taxa in the upper Midwestern United States to understand how climate and
land use change have interacted since Euro-American settlement to affect
species distributions in the region, and to examine the possible
implications of these shifts to future estimates of species resillience
to climate change across the region.

Methods
=======

We pair pre-settlement vegetation data from the Public Land Survey
System aggregated to a 64km^2^ grid [@goring2015preset] with annually
resolved climate data for the region. PLSS vegetation data contains no
uncertainty, however it is spatially discontinuous, particularly in the
south where digitization of PLSS records has occured at a slower rate.
Additionally, data in the northeastern United States [@thompson2013four]
from town Proprietor Surveys is recorded arially, with no variance
within township polygons. To accomodate these different data sources a
conditional autoregressive model was applied to the data, providing the
ability to obtain uncertainty estimates across the domain
[@paciorek2015comp].

The analysis follows a blocked sampling method. We generate probability
density functions for actual PLSS presence, by taxon, at the settlement
era and for the same PLS taxa (and point locations) in the modern era.
We also generate climate data for the FIA presence, and we generate
climate space for modern land use. In this way we can test how much of
the change in shape is the result of changing land use vs changing
climate since the early part of the century.

We also develop estimates of the modern climate space that taxa would
occupy if they remained in place, with only climate changing, and vice
versa (need better explanation here).

Results
-------

![](figures/clim_plot.tiff) **Figure 1**. *Climate change in the Upper
Midwest over the last two centuries using PRISM data and North American
Drought Atlas PDSI reconstructions. Low temperatures in the 1970s result
in lower T\_max\_ during the modern era, but climate change has resulted
in higher T\_min\_ values than the early-century normals, along with
increasing A\_ppt\_. PDSI shows strong coherence with A\_ppt\_ in this
region.*

The choice of a normal for the domain is potentially problematic. The
time transgressive nature of the PLS Survey means that some of the PLS
data comes from records sampled in the early 1800s, while other come
from the early 1900s. Regardless, throughout this period climate
variability exists (Figure 1d), and trees that were on the landscape
from the 1830s to the 1880s are not expected to have been completely
extirpated by annual scale climate variability in that same time period.
As such the use of a 'pre-settlement' climate normal from 1895-1925 is
likely a reasonable choice [ADD JULY temperatures here].

Mean temperatures in the region show an average increase of 0.59^o^C
since the 1895-1925 normal. Annual precipitation has increased slightly
(55mm), and maximum annual temeprature has declined slightly
(-0.79^o^C)). The most significant change has been in annual minimum
temperatures, which have increased by 3.3 since the 1895-1925 normal.
This broad scale increase in winter temperatures is widely accepted, and
results in an overall shift in the winter 'climate space' for the
region.

The Public Lands data is spatially extensive, sampling occured across
the region in a regular pattern, while the FIA is limited to regions
with forest cover. As such, the extent of points is not overlapping.
Given the extensive use of FIA data in generating and estimating climate
responses of tree taxa and forest types to changing climate, the use of
the FIA data as an estimate for tree species distributions is acceptable
fo estimating shifts in climate space, and attributing the loss or gain
of climate niche space due to land use and climate variability and
change during the 20^th^ century.

![](figures/loss_plots.tiff)

**Figure 2**. *Gain loss and continuous presence for PLSS and FIA tree
taxa. Continuous presence is indicated by dark green. For most taxa it
is the central parts of their distribution that show continuous
presence. Few taxa show novel presence (purple), although taxa such as
poplar, ash and basswood do show gains of over 10%. Most taxa show
significant losses in distribution (orange).*

**Table 1.** *Gain and loss of the various forest tree types of
interest.*

``` {.r}
pandoc.table(loss_table, justify = "left", style = "simple")
```

      
      
      &nbsp;         Gain     Loss   Presence  
      -------------- -------- ------ ----------
      **Tamarack**   0.004191 0.8086 0.1914    
      **Pine**       0.04076  0.5403 0.4597    
      **Spruce**     0.03769  0.5424 0.4576    
      **Fir**        0.04643  0.5243 0.4757    
      **Hemlock**    0.006298 0.6779 0.3221    
      **Cedar**      0.04157  0.5512 0.4488    
      **Poplar**     0.1972   0.4371 0.5629    
      **Maple**      0.08711  0.3341 0.6659    
      **Birch**      0.03449  0.4742 0.5258    
      **Beech**      0.007812 0.7375 0.2625    
      **Ironwood**   0.03808  0.9398 0.06018   
      **Basswood**   0.1294   0.6221 0.3779    
      **Ash**        0.1457   0.5601 0.4399    
      **Elm**        0.07377  0.785  0.215     
      **Oak**        0.09938  0.506  0.494

![](Composition_Ecology_files/figure-markdown/densities-1.png)

**Figure 3** *Changes in the overall climate space in the upper
midwestern United States. The climate densities largely overlap,
although it is possible to see the shifts in precipitation, mean and
minium temperatures, where FIA data (pink) has higher values), while the
decline in maximum temperatures is also visible, largely at the very
highest values.*

While these changes in climate are visible across the whole range, we
would expect that species with significant shifts in distribution as a
result of land use change in the post-settlement era should show changes
that superceed the changes resulting from climate alone. In particular,
if the land use pressure is biased spatially or climatically then these
changes should be even more dramatic. Hellinger distance is used to
understand the difference in the shape of kernel densities. Using the
distributions in climate space for each of the key taxa we can shift
climate, or vegetation to understand the relative effects of land use
change and climatic change across the region.

![](Composition_Ecology_files/figure-markdown/climate_shifts-1.png)

**Figure X**. *For each of the four variables examined we see...*

The Hellinger distances reveal that while taxa may show the imprint of
land use change on shifting climatic niches, this shift is not uniform
across taxa not climatic variables. Maximum temperature shows the
greatest change attributed to land use change. Tamarack, elm, poplar,
spruce beech, cedar, and fir all show a greater influence for land use
change than regional climate change.

The precipitation panel (Figure X) indicates that precipitation across
the region has driven change in the species climate niche for this
variable. Only spruce and poplar show a marginally greater effect of
land use change than climate change, wheras taxa such as beech and
hemlock show a strong climate signal.

Interstingly, since the pre-settlement era maximum temperatures have
declined, while minimum temperatures have increased strongly.

Discussion
----------

Figures
-------
