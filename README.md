# MR Corge: Sensitivity analysis of Mendelian randomization based on the core gene hypothesis for polygenic exposures


<p align="center">
  <img src="man/figures/MRCorge.png" alt="example image" width=100 height=100>
  <br>
</p>

## Overview

Mendelian randomization using polygenic exposures can be significantly biased by horizontal pleiotropy. Built on the core gene hypothesis, we developed MR Corge to identify putative core instruments that are more likely to affect genes with a direct biological role in exposures to reduce the risk of horizontal pleiotropy.

## Installation

To install the lastest version of MR Corge:

```
devtools::install_github("zhwm/MRCorge")
```

## Inputs

The input of the `mrcorge` function is a harmonized dataframe that works with the `TwoSampleMR` package. 

## Example usage

For detailed examples and usage, see the [vignette](https://zhwm.github.io/MRCorge/articles/HDL_CAD.html).


## Outputs
The `mrcorge` function outputs the MR estimates based on core instruments as well as sensitivity analysis results based on other instrument groups. Q-statistics are calculated to assess the heterogeneity of MR estimates across different groups. The `plot_mrcorge` function can be used to visualize these results. 


## Versions

* Version 0.1.0: Initial release

## Citation

If you find MR Corge useful, please cite:

[MR Corge: Sensitivity analysis of Mendelian randomization based on the core gene hypothesis for polygenic exposures.](https://doi.org/10.1101/2024.07.18.604191)

