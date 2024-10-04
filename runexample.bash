#!/bin/bash

# Set the working directory (current directory)
wkd=$(pwd)

# 1. Binary Phenotype Analysis
# This command runs LTPI using binary phenotype input. 
# The input files include:
# - binput.txt: the binary phenotype data
# - genetic_covariance_bin.txt: the genetic covariance matrix for the binary traits
# - prevalence.txt: the disease prevalence information
# The results are stored in the 'test_bin' working directory.
python ${wkd}/LTPI.py --bin ${wkd}/example/binput.txt --gencov ${wkd}/example/genetic_covariance_bin.txt --prevalence ${wkd}/example/prevalence.txt --out ${wkd}/test_bin

# 2. Continuous Phenotype Analysis
# This command runs LTPI with continuous phenotype input. 
# It requires:
# - qinput.txt: the quantitative phenotype data
# - genetic_covariance_con.txt: the genetic covariance matrix for continuous traits
# - test_bin: output from the previous binary analysis (used for reference)
# Results are saved in 'test_con'.
python ${wkd}/LTPI.py --con ${wkd}/example/qinput.txt --gencov ${wkd}/example/genetic_covariance_con.txt --bout ${wkd}/test_bin --out ${wkd}/test_con

# 3. Trait Selection Optimization
# This command optimizes trait selection based on the r2_o criterion using the following:
# - trait_C: target column name for trait selection
# - genetic_covariance_bin.txt: genetic covariance matrix
# Output results are saved in 'test_pick'.
python ${wkd}/LTPI.py --pick --pi trait_C --gencov ${wkd}/example/genetic_covariance_bin.txt --out ${wkd}/test_pick
