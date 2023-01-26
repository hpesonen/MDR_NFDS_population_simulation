# Introduction

This repository contains source code for running the simulation study of the paper *Beta-lactam antibiotic use modulates multi-drug resistant clone success in Escherichia coli populations: a longitudinal cross-country genomic cohort study*, Pöntinen, Gladstone, Pesonen et al. 

Main simulation file : `pop_dynamic_run.R`

data-folder contains BSAC and NORM datasets

Results are in results-folder as csv-files :

- BSAC - results from BSAC
- NORM - results from NORM
- BSAC_NORM - results from BSAC (Scenario 2)
- ST131 refers to runs with ST131 variants
- sigma_f_XXX - indicates fixed value of sigma_f = XXX
- CTX_prop : CTXM prevalence every 25 time steps
- genotype_prop : prevalences of 4 most common strains every 25 time steps
