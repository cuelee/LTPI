#!/bin/bash
#$ -l mem=8G,time=12:: -S /bin/bash -N LTPI -j y -cwd 
wkd='/ifs/scratch/msph/eigen/hl3565/01_MTGB/codes/source/LTPI_package'
python ${wkd}/LTPI.py --bin ${wkd}/example/binput.txt --gencov ${wkd}/example/genetic_covariance_bin.txt --prevalence ${wkd}/example/prevalence.txt --out ${wkd}/test_imp
#python ${wkd}/LTPI.py --con ${wkd}/example/qinput.txt --bout ${wkd}/test_imp --gencov ${wkd}/example/genetic_covariance_con.txt --out ${wkd}/test_mle 
