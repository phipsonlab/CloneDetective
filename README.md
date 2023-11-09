# CloneDetective

CloneDetective is an R package to interrogate clonal abundance data generated
using lineage tracing protocol.

It contains functions to analyse clone barcode distributions from NGS data which 
can inform the design of future scRNAseq experiments. 
For a given set of experimental parameters (e.g. the number of cells captured), 
it can predict clone barcode representation in the resulting scRNAseq data by estimating clonal abundance from NGS data.
NGS data refers to dedicated DNA barcoding data which exclusively sequences 
the synthetic lineage tracing clone barcode reads using Next Generation Sequencing.

For scRNAseq data barcoded with lineage tracing protocol, CloneDetective contains
functions to perform an initial assessment of clone barcode uptake by cells and 
to assign clones to cells.

This R package works greatly hand in hand with the 
[NextClone pipeline](https://github.com/phipsonlab/NextClone).

## Citation

If you use CloneDetective in your study, please kindly cite our preprint on bioRxiv.

## Installation

The package can be installed using `devtools`:

```
# Install devtools
install.packages("devtools")

# Install SuperCellCyto from this repository
devtools::install_github("phipsonlab/CloneDetective")
```

## Contributing

We welcome any contributions! 

Please submit a Github issue if you run into issues.

Want to suggest improvements? Kindly open a Github issue and submit a pull request referencing the issue.


