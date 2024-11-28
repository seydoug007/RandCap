---
title: "RandCap"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{RandCap}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


# Introduction

`RandCap` is an R package designed to streamline randomization in clinical trials. It provides robust and customizable randomization strategies, supporting blocked randomization, stratification by variables, and variable allocation ratios. The package is tailored for researchers and trial designers who require balanced allocation schemes and reproducible randomization processes.

This vignette demonstrates the key features of `RandCap` with practical examples.

# Installation

To install `RandCap`:

# install the development version from GitHub
devtools::install_github("seydoug007/RandCap")

# Load the library
```{r}
library(RandCap)
```




# Generate a Basic Randomization List

```{r}
rand_dev_list <- RandCapGen(
  n = 100,
  block_sizes = c(10, 10),
  arms = c("Control", "Treatment"),
  ratio = c(1, 1),
  seed = 1234,
  project_acronym = "RCT"
)
```

# View the full dataset
```{r}
head(rand_dev_list$tables$full_dataset)

```


# Stratified Randomization
```{r}
rand_list_strat <- RandCapGen(
  n = 100,
  block_sizes = c(20, 20),
  arms = c("Control", "Treatment"),
  strat_vars = list(Sex = c("Male", "Female"), AgeGroup = c("Young", "Old")),
  strat_vars_prefix = list(Sex = "SEX_", AgeGroup = "AGE_"),
  seed = 5678,
  project_acronym = "RCT"
)
```



# View the full dataset

```{r}
head(rand_list_strat$tables$full_dataset)
```


# Finalize Production Randomization

The production

```{r}
prod_list <- RandCapProd(randomization_object = rand_list_strat, seed = 9876)

```


<!-- # View the production-ready dataset -->
<!-- ```{r} -->
<!-- head(prod_list$tables$full_dataset) -->
<!-- ``` -->

<!-- # Assess Randomization Balance -->

<!-- ```{r} -->
<!-- RandCapBalance( -->
<!--   randomization_object = rand_list_strat, -->
<!--   output_path = "Balance_Summary.pdf" -->
<!-- ) -->
<!-- ``` -->

<!-- # Export Randomization Settings -->
<!-- ```{r} -->

<!-- RandCapSettings(rand_obj = rand_list, output_path = "Settings_Summary.pdf") -->
<!-- ``` -->

<!-- # Export Tables for REDCap -->
<!-- ```{r} -->
<!-- RandCapTable( -->
<!--   randomization_object = rand_list, -->
<!--   save_for_REDCap = TRUE, -->
<!--   save_random_table = TRUE -->
<!-- ) -->
<!-- ``` -->

<!-- # Custom Ratios -->

<!-- ```{r} -->
<!-- rand_custom <- RandCapGen( -->
<!--   n = 90, -->
<!--   block_sizes = c(18, 18), -->
<!--   arms = c("Treatment", "Control"), -->
<!--   ratio = c(2, 1), -->
<!--   seed = 3456, -->
<!--   project_acronym = "RCT" -->
<!-- ) -->
<!-- ``` -->


<!-- # View the dataset -->
<!-- ```{r} -->
<!-- head(rand_custom$tables$full_dataset) -->
<!-- ``` -->

<!-- # Visualization -->
<!-- ```{r} -->
<!-- library(ggplot2) -->

<!-- ggplot(rand_list$tables$full_dataset, aes(x = treatment_arm)) + -->
<!--   geom_bar(fill = "skyblue") + -->
<!--   labs(title = "Treatment Arm Distribution", x = "Treatment Arm", y = "Count") + -->
<!--   theme_minimal() -->
<!-- ``` -->


