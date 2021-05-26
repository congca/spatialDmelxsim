# spatialDmelxsim

Allelic counts from Combs & Fraser (2018), compiled into a
*SummarizedExperiment* object.

Here, `a1` assay is the *D simulans* allelic count, `a2` assay is the
*D melanogaster* allelic count.

Data downloaded from `GSE102233` as described in the publication:

> Combs PA, Fraser HB (2018) 
> Spatially varying cis-regulatory divergence in Drosophila embryos 
> elucidates cis-regulatory logic. 
> *PLOS Genetics* 14(11): e1007631. 
> https://doi.org/10.1371/journal.pgen.1007631

## R scripts

* Script to build the *SummarizedExperiment*: `inst/scripts/make-data.R`
* Demo of loading the data and calculating allelic ratio: 
  `vignettes/spatialDmelxsim.Rmd`

## Example figures of allelic ratio (simulans / total):

![slam](man/figures/slam.png)

![uif](man/figures/uif.png)
