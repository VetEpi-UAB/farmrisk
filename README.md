
<!-- README.md is generated from README.Rmd. Please edit that file -->

# farmrisk

<!-- badges: start -->
<!-- badges: end -->

This package provides tools for performing quantitative risk analysis to
prioritize farm external biosecurity measures.

## Installation

You can install the development version of farmrisk from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("VetEpi-UAB/farmrisk")
```

## Example

Below is a basic example of how to use the farmrisk package with the
farmrisk_cattle() function:

``` r
library(farmrisk)

# Run the model for example beef farm
results <- farmrisk_cattle(farm_id="eg_beef")

# Print the final result
results$node_list$total_agg$summary

# Print model expression
results$model_expression

# Look at pasture movements
results$node_list$pasture_all_mov$summary
# Look at aggregated result for all pasture movements
results$node_list$pasture_inf_agg$summary

# Examine pasture data
results$data$pasture

# Examine the references for the used pathogen parameters
unique(results$data$pasture$pathogen_ref)
```
