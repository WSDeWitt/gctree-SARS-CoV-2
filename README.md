# gctree-SARS-CoV-2
gctree for SARS-CoV-2 tree inference.
Uses data from:
> J Bloom, Recovery of deleted deep sequencing data sheds more light on the early Wuhan SARS-CoV-2 epidemic, bioRxiv 2021.06.18.449051; doi: https://doi.org/10.1101/2021.06.18.449051
>
>https://github.com/jbloom/SARS-CoV-2_PRJNA612766


## Environment

Set up conda environment:
```bash
conda create -f environment.yml
```
Activate it:
```bash
conda activate gctree-SARS-CoV-2
```

## Pipeline

Run snakemake pipeline:
```bash
snakemake --cores all
```
This will create a directory `build/` with subdirectories for each of three rootings shown in the paper. In each subdirectory, a `gctree.log` file will show a ranking of maximum parsimony trees (generated with PHYLIP) according to a branching process likelihood that accounts for abundance. There are also SVG renderings of each of the parsimony trees, named like `gctree.out.inference.[1,2,...].svg`, where the number indicates rank as in the log file. In these renderings, colors red, blue, and green (and counts) in pie chart nodes indicate abundance partitions into categories "Wuhan", "other China", and "outside China", as provided in the data.
