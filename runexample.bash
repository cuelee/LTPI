#!/bin/bash
wkd=$(pwd)
python ${wkd}/LTPI.py --bin ${wkd}/example/binput.txt --gencov ${wkd}/example/genetic_covariance_bin.txt --prevalence ${wkd}/example/prevalence.txt --out ${wkd}/test_bin
python ${wkd}/LTPI.py --con ${wkd}/example/qinput.txt --gencov ${wkd}/example/genetic_covariance_con.txt --bout ${wkd}/test_bin --out ${wkd}/test_con
