#!/usr/bin/env python3

"""
LT-PI: Liability Threshold model-based method that predicts the posterior mean genetic liability of Phenotype of Interest
Copyright(C) 2021 Cue Hyunkyu Lee 

LT-PI is a framework that estimates the probability of the phenotype of a target trait 
by leveraging the phenotypic information of multiple related traits. Our method relies 
on a liability threshold model and provides the posterior mean genetic liability scores 
to estimate the phenotypic probability. LT-PI uses the Monte Carlo integration method 
to estimate the posterior mean genetic liability. This process requires phenotypic 
information (usually Phecode; coded as 1/0) obtained from the electronic health (EMR) 
record and some additional information, including narrow-sense heritability and 
prevalence of traits, collected from literatures.


- LTPI uses a parallel operation method when performing Monte Carlo Integration that generates samples from nulls.
- LTPI only estimates configurations existed in test samples. This increases computational efficacy.
"""

import numpy as np
import pandas as pd
import os
import sys
import traceback
import argparse
import time
import warnings

# from framework.importance_sampling import LTPI_IMP
from framework.est_maxlikelihood import LTPI_MLE
from framework.ghk_algorithm import LTPI_GHK, iterated_covariance_optimization
from framework.r2_selection import run_ATSA

codename = 'LTPI'
__version__ = '0.0'
MASTHEAD = """
************************************************************
* LTPI({c})
* Version {V}
* (C)2021 Cue H. Lee
* Columbia University
* Unlicensed software
************************************************************
""".format(c=codename, V=__version__)

def sec_to_str(t):
    """Convert seconds to days:hours:minutes:seconds"""
    intervals = (('d', 86400), ('h', 3600), ('m', 60), ('s', 1))
    f = ''
    for n, c in intervals:
        v = t // c
        if v:
            t -= v * c
            if c != 1:
                f += '{}{} '.format(round(v), n)
    return f


class Logger:
    """Lightweight logging."""
    
    def __init__(self, fh):
        self.log_fh = open(fh, 'w')

    def log(self, msg):
        """Print to log file and stdout with a single command."""
        print(msg, file=self.log_fh)
        print(msg)

    def mlog(self, msg):
        """[mutelog] Print to log file without stdout with a single command."""
        print(msg, file=self.log_fh)

parser = argparse.ArgumentParser()

# Output arguments
parser.add_argument('--out', default='LTPI', type=str,
                    help='Prefix of the output files. If --out is not set, PLEIO will use out as the '
                         'default output directory (LTPI).')
parser.add_argument('--bout', default=None, type=str,
                    help='Required for --con: Prefix of the ltpi binary test output file.')

# Test mode arguments
parser.add_argument('--bin', default=None, type=str,
                    help='File prefix of the binary test input.')
parser.add_argument('--prevalence', default=None, type=str,
                    help='Required for --bin: Path of genetic disease prevalence file.')
parser.add_argument('--con', default=None, type=str,
                    help='File prefix of the quantitative test input.')
parser.add_argument('--pick', default=False, action='store_true',
                    help='Selects important non-target traits based on R2S algorithm.')
parser.add_argument('--rint', default=False, action='store_true',
                    help='Run rank-based inverse normal transformation on LTPI binary step.')
parser.add_argument('--pi', default=None, type=str,
                    help='Required for --pick: The columns name for the target')
parser.add_argument('--Q', default=5, type=int,
                    help='The number of non-targets to be selected [pick]. (Default value is 5)')

# Covariance matrix arguments
parser.add_argument('--gencov', default=None, type=str,
                    help='Required for either --bin and --con: Path of genetic covariance matrix file.')
parser.add_argument('--envcov', default=None, type=str,
                    help='Required for either --bin and --con: Path of environmental covariance matrix file.')
parser.add_argument('--shrink', default=None, type=str,
                    help='[Optional] Covariance Shrinkage on GHK input - G: GENCOV, E: ENVCOV, B: Both.')
parser.add_argument('--shrink_target', default = 1.5, type=float,
                    help='[Optional] Covariance Shrinkage target (default: 1.5).')

# Other arguments
parser.add_argument('--nsample_main', default=50000, type=int,
                    help='Number of samples to be generated from the MVN distribution. (Should be greater than 50K)')
parser.add_argument('--nsample_helper', default=300000, type=int,
                    help='The number of samples used to generate a liability score for a configuration that is '
                         'difficult to sample. You can use the default parameters for this option. (default is 300K)')
parser.add_argument('--ncore', default=1, type=int,
                    help='Number of CPU cores for parallel computing. (Default value is 1)')
parser.add_argument('--thres', default=0.1, type=float,
                    help='Threshold used in importance sampling method')


if __name__ == '__main__':
    def is_pos_def(x, name = '', thres=1e-9):
        """Check if the matrix x is positive definite"""
        c = x.columns
        w, v = np.linalg.eig(x)
        if not np.all(w > 0):
            args.log.log('The %s covariance matrix is not positive definite'%name)
            w[w < 0] = thres
        return pd.DataFrame(v.dot(np.diag(w)).dot(np.transpose(v)), index=c, columns=c)
    
    def envcov_QC(envcov, gencov):
            SIDEC = np.diag(np.sqrt(1 / np.diag(envcov)))
            OMSIDGC = np.diag(np.sqrt((1 - np.diag(gencov))))
            return(pd.DataFrame(OMSIDGC @ SIDEC @ envcov.values @ SIDEC @ OMSIDGC, 
            index = envcov.index, columns = envcov.columns)) 
            
    def cov_shrink(args, cov, keys):
        if args.shrink in keys:
            t = args.shrink_target
            res, s = iterated_covariance_optimization(cov, t)
            args.log.log('Covariance Chrinkage: Mode-%s, Sfactor-%s, Cond-%s'%(keys[0], s, np.linalg.cond(res)))
        else:
            t = None
            res = cov
        return(res)
            
    def read_prev(f):
        """Read the genetic disease prevalence file"""
        df = pd.read_csv(f, sep='\s+', header='infer', usecols=['TID', 'prev'], dtype={'TID': str, 'prev': float})
        df.set_index(['TID'], inplace=True)
        b = ~np.isnan(df.prev)
        d = {n: float(df.loc[n, 'prev']) for n in df.index[b]}
        bn = df.index[b].to_numpy(dtype='U100')
        args.log.log('Prevalence file contains {} traits'.format(len(bn)))
        return d, bn

    def read_bout(p):
        """Read the binary test output files"""
        conf = pd.read_csv(str(p) + '.conf', sep='\t', header='infer', index_col=None,
                           usecols=['CONF', 'EG', 'G_SE'], dtype={'CONF': str, 'EG': float, 'G_SE': float})
        conf.set_index(['CONF'], inplace=True)
        samp = pd.read_csv(str(p) + '.txt.gz', sep='\s+', header='infer', index_col=None,
                           usecols=['IID', 'LTPI', 'CONF'], dtype={'IID': str, 'LTPI': float, 'CONF': str})
        samp.set_index(['IID'], inplace=True)
        return conf, samp

    def read_binary_ltpiin(f, TID=None, d=str):
        """Read the binary test input file"""
        col_names = pd.read_csv(f, sep='\t', nrows=0).columns
        df = pd.read_csv(f, sep='\t', header='infer', dtype={c: d if c != 'IID' else str for c in col_names})

        # Set index, fill NaN values with 'X'
        df.set_index('IID', inplace=True)
        df.fillna('X', inplace=True)

        args.log.log('Input contains {} samples for {} traits'.format(df.shape[0], df.shape[1]))

        if TID is None:
            return df
        return df.loc[:, TID]

    def read_continuous_ltpiin(f, TID=None, d=str):
        """Read the quantitative test input file"""
        col_names = pd.read_csv(f, sep='\t', nrows=0).columns
        df = pd.read_csv(f, sep='\t', header='infer', dtype={c: d if c != 'IID' else str for c in col_names})

        # Set index and drop rows with all NaN values
        df.set_index('IID', inplace=True)
        df.dropna(axis=0, how='all', subset=None, inplace=True)

        args.log.log('Input contains {} samples for {} traits'.format(df.shape[0], df.shape[1]))

        if TID is None:
            return df
        return df.loc[:, TID]
    
    def read_cov(f, c=None):
        """Read the covariance matrix file"""
        cov_raw = pd.read_csv(f, sep='\t', header='infer', index_col=None)
        cov_raw.index = cov_raw.columns

        if c is None:
            c = cov_raw.columns

        cov_sub = cov_raw.loc[c, c]

        return cov_sub
        
    def transform_dataframe(df):
        # Get the diagonal of the DataFrame
        diag_values = np.diag(df)

        # Subtract it from 1
        transformed_values = 1 - diag_values

        # Create a diagonal matrix
        diag_matrix = np.diag(transformed_values)

        # Convert it back to DataFrame with proper index and columns
        new_df = pd.DataFrame(diag_matrix, index=df.index, columns=df.columns)

        return new_df        

    def write_r2(args):
        if len(args.selected_traits) > 1:
            selected_traits_df = pd.Series(args.selected_traits)
            selected_traits_df.to_csv(args.out + '.r2_traits', sep='\t'
            , index=False, header=False, na_rep='NA')
            best_shrinkage_df= pd.Series(args.best_shrinkage)
            best_shrinkage_df.to_csv(args.out + '.r2_best', sep='\t'
            , index=False, header=False, na_rep='NA')
            args.best_S.to_csv(args.out + '.r2_summary', sep='\t'
            , index=True, header=True, na_rep='NA')

            args.log.log('--pick summary:\nTraits selcted:\n{}\nShrinkage factor:\n{}\nSummary: \n{}'.format(args.selected_traits, args.best_shrinkage, args.best_S))

            return None
        
    def update_r2(args):
        if len(args.selected_traits) < 1:
            raise ValueError('r2_traits has the size of %s'%len(args.selected_traits))
        elif args.bin is not None:
            args.binary_traits = args.selected_traits
            args.prev = {k:args.prev[k] for k in args.binary_traits}
            args.GENCOV = args.GENCOV.loc[args.binary_traits,args.binary_traits].copy()
            args.ENVCOV = args.ENVCOV.loc[args.binary_traits,args.binary_traits].copy()
        elif args.con is not None:
            args.quantitative_traits = args.selected_traits[1:]
            args.mle_traits = args.selected_traits
            args.ltpiin = args.ltpiin.loc[:,args.quantitative_traits].copy()
            args.GENCOV = args.GENCOV.loc[args.mle_traits,args.mle_traits].copy()
            args.ENVCOV = args.ENVCOV.loc[args.mle_traits,args.mle_traits].copy()
        return args        

    args = parser.parse_args()    
    if args.out is None:
        raise ValueError('No command line arguments were provided.')
    
    log = Logger(args.out + '.log')

    try:
        defaults = vars(parser.parse_args(''))
        opts = vars(args)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        header = MASTHEAD
        header += 'Call: \n'
        header += './pleio.py \\\n'
        options = ['--' + x.replace('_', '-') + ' ' + str(opts[x]) + ' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True', '').replace('False', '')
        header = header[0:-1] + '\n'
        args.log = log
        args.log.log(header)
        args.log.log('Beginning analysis at {T}'.format(T=time.ctime()))
        start_time = time.process_time()
            
        if args.bin is not None:
            if args.prevalence is not None:
                args.prev, args.binary_traits = read_prev(args.prevalence)
            else:
                args.log.log('ValueError: Can\'t read --prev flag.')
                raise ValueError()

            if args.gencov is not None:
                args.GENCOV = read_cov(args.gencov, args.binary_traits)
                args.pi = args.GENCOV.index[0]
            else:
                args.log.log('ValueError: Can\'t read --gencov flag.')
                raise ValueError() 

            if args.envcov is not None:
                args.ENVCOV = read_cov(args.envcov, args.binary_traits)
            else:
                warnings.warn('No environmental covariance matrix input. LTPI will use 1-diag(GENCOV) as a substitute.',
                              category=UserWarning)
                args.ENVCOV = transform_dataframe(args.GENCOV)

            if args.pick:
                args.selected_traits, args.best_S, args.best_shrinkage, args.summary_data = run_ATSA(args)
                write_r2(args)
                args = update_r2(args)
                
            args.GENCOV = is_pos_def(cov_shrink(args, args.GENCOV.copy(), keys = ['G','B']),'GEN')
            args.ENVCOV = is_pos_def(cov_shrink(args, args.ENVCOV.copy(), keys = ['E','B']),'ENV')
            
            args.log.log('Condition Number: \nGEN-%s ENV-%s GEN+ENV-%s'%(np.linalg.cond(args.GENCOV),np.linalg.cond(args.ENVCOV),np.linalg.cond(args.GENCOV+args.ENVCOV)))
            
            args.ltpiin = read_binary_ltpiin(f=args.bin, TID=args.binary_traits, d=str)
            args.conf, args.samp_bin, args.time = LTPI_GHK(args)
            
            args.conf.to_csv(args.out + '.conf', sep='\t', index=True, header=True, na_rep='NA')
            args.samp_bin.to_csv(args.out + '.txt.gz', sep='\t', index=True, header=True, compression='gzip',
                                 na_rep='NA')

        elif args.con is not None:
            args.ltpiin = read_continuous_ltpiin(f=args.con, d=float)
            args.quantitative_traits = args.ltpiin.columns.to_numpy(dtype='U100')
            if args.gencov is not None:
                gencov = read_cov(args.gencov)
                args.pi = gencov.index[0]
                args.mle_traits = np.insert(args.quantitative_traits, 0, args.pi)
                args.GENCOV = gencov.loc[args.mle_traits, args.mle_traits]
            else:
                args.log.log('ValueError: Can\'t read --gencov flag.')
                raise ValueError() 
                
            if args.envcov is not None:
                envcov = read_cov(args.envcov,args.mle_traits)
                args.ENVCOV = envcov.loc[args.mle_traits, args.mle_traits]
            else:
                warnings.warn('No environmental covariance matrix input. LTPI will use 1-diag(GENCOV) as a substitute.',
                              category=UserWarning)
                args.ENVCOV = transform_dataframe(args.GENCOV)
                
            if args.pick:
                args.selected_traits, args.best_S, args.best_shrinkage, args.summary_data = run_ATSA(args)
                write_r2(args)
                args = update_r2(args)
            
            if args.bout is not None:
                args.conf, args.samp_bin = read_bout(args.bout)
            else:
                args.log.log('ValueError: Can\'t read --bout flag.')
                raise ValueError()

            args.GENCOV = is_pos_def(args.GENCOV.copy(),'GEN')
            args.ENVCOV = is_pos_def(args.ENVCOV.copy(),'ENV')
            
            args.samp_mle, args.time = LTPI_MLE(args)
            args.samp_mle.to_csv(args.out + '.txt.gz', sep='\t', index=True, header=True, compression='gzip',
                                 na_rep='NA')
        
        elif args.pick:  
            if args.pi is None:
                args.log.log('ValueError: --pick requires pi column name (--pi [PI]).')
                raise ValueError() 

            if args.gencov is not None:
                    args.GENCOV = read_cov(args.gencov)
            else:
                args.log.log('ValueError: Can\'t read --gencov flag.')
                raise ValueError() 
            
            args.log.log('Start R2 selection')
            args.selected_traits, args.best_S, args.best_shrinkage, args.summary_data = run_ATSA(args)
            args.log.log('Start Writing Result')
            write_r2(args)

        else:
            args.log.log('ValueError: Can\'t read test flag.')
            raise ValueError()            

    except Exception:
        args.log.mlog(traceback.format_exc())
        raise

    finally:
        args.log.log('Analysis finished at {T}'.format(T=time.ctime()))
        time_elapsed = round(time.process_time() - start_time, 2)
        args.log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))
