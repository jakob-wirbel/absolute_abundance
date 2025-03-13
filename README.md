# absolute_abundance

Data and code to predict absolute abundance for metagenomic samples

## Description

This project contains the analyses and results for the publication:
> [Accurate prediction of absolute prokaryotic abundance from DNA concentration](https://www.cell.com/cell-reports-methods/current)

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
- `models` the models as `Rdata` files
- `figures` the raw figures for the publication

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

# load the model
load('./models/full_model.RData`)

# load the function
source('./src/apply_model.R`)

# apply to your own data
df.predicted <- .f_apply_model(df, rr_full)

```

## Contact

If you have questions about the code or the analyses presented in this 
project, please feel free to contact me, 
[Jakob Wirbel](maito:wirbel@stanford.edu) 
or [open an issue in this repository](https://github.com/jakob-wirbel/absolute_abundance/issues/new)
