# Composition_Climate

**Submitted**: 2016/09/23

## Contributors

* [Simon Goring](http://github.com/SimonGoring)
* [Jack Williams](http://github.com/IceAgeEcologist)

## About

This is the `git` repository for the paper *Effect of historic land-use and climate change on tree-climate relationships in the northern United States*.  The paper details shifts in the climate-space occupied by key taxa in the north-central United States as a result of climate change and land use change over the last 200 years.

### Key Findings

Historic land-use in the region has reinforced climatic shifts in tree species distributions toward cooler or wetter distributions than those tree species occupied in the early 1800s.

## What's in the Repository

The key file is `Composition_Ecology.Rmd`, it is an RMarkdown file that contains both text & embedded code.  The file will render to a version that is numerically equivalent to the submitted paper.  Using RStudio you simply need to load the `Rmd` file and select `Knit Word`, or `CTRL+Shift+K` or `File > Knit Document`.

From the command line you can use the command

```
> Rscript -e "rmarkdown::render('Composition_Ecology.Rmd', output_format = 'docx')
```

## What's Missing

Because some of the data, particularly in the mapping functions, uses data distributed by third parties, we have chosen to indicate where this data can be obtained, but have not chosen to make the data available directly.  In particular, the following files must be obtained from [NaturalEarthData](http://naturalearthdata.com):

* [`ne_10m_coastline.zip`](http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip)
* [`ne_10m_lakes`](http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_lakes.zip)
* [`NE1_HR_LC_SR_W_DR`](http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/raster/NE2_HR_LC_SR_W_DR.zip)
* [`ne_10m_rivers_lake_centerlines`](http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_rivers_lake_centerlines.zip)

If you use this repository and notice any other files missing, or if you notice that the paper will not compile properly, please raise an issue and we will address it quickly.

## Contributions

This paper has been submitted to *Ecology Letters*.  We welcome future contributions, and any corrections or omissions that may be required to improve this paper, however, by the standards of scientific authorship we cannot add authors to the paper at this time.
