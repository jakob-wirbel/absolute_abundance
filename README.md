# absolute_abundance

Data and code to predict absolute abundance for metagenomic samples

## Description

This project contains the analyses and results for the publication:
> [Accurate prediction of absolute prokaryotic abundance from DNA concentration](https://doi.org/10.1016/j.crmeth.2025.101030)

Please refer to the publication for more detail.

In short, we observed a strong correlation between DNA concentration and
total number of 16S rRNA copies as determined by digital droplet PCR for
stool samples from individuals undergoing allogeneic hematopoietic cell
transplantation. Based on this correlation, we trained a machine learning
model to predict 16S copies (as a proxy for absolute microbial abundance) 
from DNA concentration and high-level taxonomic information available from 
metagenomic sequencing. This model was validated on an external set of
data from stool samples of individuals with Parkinson's disease.

This repository contains:
- `data` the original data
- `src` the analysis code
- `models` a place for the models as `Rdata` files
- `figures` the raw figures for the publication

Please note that the model files are too large to be hosted on
Github, so you will have to download them from 
[this Zenodo repository](https://doi.org/10.5281/zenodo.14026970)

## Applying the model to your own data

To apply the model to your own data, you will need a data frame containing the
following columns for each sample:

- `DNA_concentration`: DNA concentration in ng/µL
- `bac_frac`: Percentage of reads retained after host-removal (range 
from 0 to 1); you can use for example
[this Nextflow pipeline](https://github.com/bhattlab/bhattlab_workflows_nf)
- `alpha`: Alpha diversity of your metagenomic sample, calculated as Shannon 
index with the `diversity` function in the vegan R package
- `k__Eukaryota`, `k__Archaea`, `k__Bacteria`, and `unclassified`: Percentage of 
relative abundance explained by various high-level taxonomic groups, determined
by [MetaPhlAn4](https://github.com/biobakery/MetaPhlAn) (range from 0 to 100)
- `Sample_type`: Information about the sample processing; can be either 
`same_day` or `next_day` processing of the sample [this feature is not
very important in the model]

There is also a dedicated function to apply the model to your data contained
in this repository, that you can use in the following way (your own data
should be in a data.frame `df`):

```{r}

library("mlr3")

# load the model
load('./models/full_model.RData`)

# load the function
source('./src/apply_model.R`)

# apply to your own data
df.predicted <- .f_apply_model(df, rr_full)

```

## Version problems

Unfortunately, some users have reported that they could not load the model
properly. The underlying problem seems to be that there was a version mismatch
between the version of `mlr3`, `paradox`, and `mlr3tuning` of my computer (where
the model was trained) and the computer where it was loaded, preventing the
proper application of the model.  
To fix this issue, the model was retrained with these package versions (not
quite sure which package is actually the culprit; these are my best guesses):

| Model version | Zenodo link | R package versions |
| ------------- | ----------- | ------------------ |
| Version 1.0   | https://zenodo.org/records/14026970 | mlr3_0.16.0, mlr3tuning_0.18.0, paradox_0.11.1, Rcpp_1.0.10 |
| Version 2.0   | https://zenodo.org/records/16969374 | mlr3_1.1.0, mlr3tuning_1.4.0, paradox_1.0.1, Rcpp_1.1.0 |

Going forward, if you run into any of these errors:
```{R}
df.predicted <- .f_apply_model(df, rr_full)
Error in private$.store_values(xs) : attempt to apply non-function
```
or 
```{R}
df.predicted <- .f_apply_model(df, rr_full)
Error in .__ResampleResult__learners(self = self, private = private, super = super,  : 
  could not find function ".__ResampleResult__learners"
```
you might want to re-train the model on your own computer, by executing the 
first 130-odd lines of the script 
[`train_models.R`](https://github.com/jakob-wirbel/absolute_abundance/blob/main/src/train_models.R) script.
On my computer, it took about 12 minutes to train the full model (see this 
[issue](https://github.com/jakob-wirbel/absolute_abundance/issues/1#issuecomment-3229301728)).


## Contact

If you have questions about the code or the analyses presented in this 
project, please feel free to contact me, 
[Jakob Wirbel](maito:wirbel@stanford.edu) 
or [open an issue in this repository](https://github.com/jakob-wirbel/absolute_abundance/issues/new)
