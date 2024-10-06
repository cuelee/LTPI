#!/usr/bin/env python3

"""
LT-PI: Liability Threshold Model for Disease Risk Prediction
Copyright (C) 2024 Cue Hyunkyu Lee

LT-PI is a statistical framework that estimates the conditional expected genetic liability of a target disease by leveraging phenotypic data from multiple related traits. Built on a liability threshold model, LT-PI calculates disease probabilities using phenotypic information from electronic health records (EHR), combined with heritability and prevalence data from the literature. The model integrates binary and continuous traits using the GHK algorithm and maximum likelihood estimation.
"""

import numpy as np
import pandas as pd
import argparse
import time
import warnings
import traceback

from framework.est_maxlikelihood import LTPI_MLE
from framework.ghk_algorithm import LTPI_GHK, iterated_covariance_optimization
from framework.r2_selection import run_ATSA
from framework.ltpi_utils import (
    is_pos_def, envcov_QC, cov_shrink, read_prev, read_bout,
    read_binary_ltpiin, read_continuous_ltpiin, read_cov,
    transform_dataframe, write_r2, update_r2
)

codename = 'LTPI'
__version__ = '1.0'
MASTHEAD = f"""
************************************************************
* LTPI ({codename})
* Version {__version__}
* (C) 2024 Cue H. Lee, Columbia University
* MIT license
************************************************************
"""

def sec_to_str(t):
    """Convert seconds to days:hours:minutes:seconds."""
    intervals = [('d', 86400), ('h', 3600), ('m', 60), ('s', 1)]
    f = ''
    for label, sec in intervals:
        v = t // sec
        if v:
            t -= v * sec
            f += f'{round(v)}{label} '
    return f

class Logger:
    """Lightweight logging utility."""
    
    def __init__(self, filepath):
        self.log_fh = open(filepath, 'w')

    def log(self, msg):
        """Log message to both file and stdout."""
        print(msg, flush=True)
        print(msg, file=self.log_fh)

    def mlog(self, msg):
        """Log message only to file."""
        print(msg, file=self.log_fh)


parser = argparse.ArgumentParser()

# Output arguments
parser.add_argument('--out', default='LTPI', type=str,
                    help='Output file prefix. Default is "LTPI".')

# Test mode arguments
## Binary test input
parser.add_argument('--bin', default=None, type=str,
                    help='Path to the binary phenotype input.')
parser.add_argument('--prevalence', default=None, type=str,
                    help='Path to disease prevalence information (required for --bin).')

## Continuous test input
parser.add_argument('--con', default=None, type=str,
                    help='Path to the quantitative phenotype input.')
parser.add_argument('--bout', default=None, type=str,
                    help='Same as the --out argument used in the --bin test. Required for --con.')

# Trait selection (Optimization with --pick)
parser.add_argument('--pick', default=False, action='store_true',
                    help='Optimize non-target trait selection based on the r2_o criterion.')
parser.add_argument('--pi', default=None, type=str,
                    help='Target column name for the trait selection (required for --pick).')
parser.add_argument('--Q', default=30, type=int,
                    help='Number of non-target traits to select (required for --pick; default: 30).')

## Rank-based inverse normal transformation
parser.add_argument('--rint', default=False, action='store_true',
                    help='Apply rank-based inverse normal transformation on LTPI scores (optional for --bin, --con.')

# Covariance matrix arguments
parser.add_argument('--gencov', default=None, type=str,
                    help='Path to genetic covariance matrix (required for --bin or --con).')
parser.add_argument('--envcov', default=None, type=str,
                    help='Path to environmental covariance matrix (required for --bin or --con).')
parser.add_argument('--shrink', default=None, type=str,
                    help='Apply covariance shrinkage: G for GENCOV, E for ENVCOV, B for both.')
parser.add_argument('--shrink_target', default=1.5, type=float,
                    help='Target condition number for covariance shrinkage (default: 1.5).')

# Parameters specific to GHK algorithm
parser.add_argument('--nsample', default=50000, type=int,
                    help='Number of samples for the GHK algorithm (default: 50K).')

# Parameters specific to R2 selection (ATSA)
parser.add_argument('--r2', default=0.0, type=float,
                    help='R2_o threshold for ATSA (default: 0).')


if __name__ == '__main__':

    args = parser.parse_args()    
    if not args.out:
        raise ValueError('No output file prefix (--out) provided.')
   
    log = Logger(f'{args.out}.log')

    try:
        # Generate command call header
        defaults = vars(parser.parse_args(''))
        opts = vars(args)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        header = MASTHEAD
        header += 'Call:\n./pleio.py \\\n'
        options = ['--' + x.replace('_', '-') + ' ' + str(opts[x]) + ' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True', '').replace('False', '') + '\n'
        
        # Log the start of analysis
        args.log = log
        args.log.log(header)
        args.log.log(f'Beginning analysis at {time.ctime()}')
        start_time = time.process_time()

        if args.bin:
            if args.prevalence:
                args.prev, args.binary_traits = read_prev(args, args.prevalence)
            else:
                args.log.log('ValueError: Can\'t read --prev flag.')
                raise ValueError()

            if args.gencov:
                args.GENCOV = read_cov(args.gencov, args.binary_traits)
                args.pi = args.GENCOV.index[0]
            else:
                args.log.log('ValueError: Can\'t read --gencov flag.')
                raise ValueError()

            args.ENVCOV = read_cov(args.envcov, args.binary_traits) if args.envcov else transform_dataframe(args.GENCOV)

            if args.pick:
                args.selected_traits, args.best_S, args.best_shrinkage, args.summary_data = run_ATSA(args)
                write_r2(args)
                args = update_r2(args)
                
            args.GENCOV = is_pos_def(args, cov_shrink(args, args.GENCOV.copy(), keys=['G', 'B']), 'GEN')
            args.ENVCOV = is_pos_def(args, cov_shrink(args, args.ENVCOV.copy(), keys=['E', 'B']), 'ENV')
            
            args.log.log(f'Condition Number:\nGEN-{np.linalg.cond(args.GENCOV)} ENV-{np.linalg.cond(args.ENVCOV)} '
                         f'GEN+ENV-{np.linalg.cond(args.GENCOV + args.ENVCOV)}')
            
            args.ltpiin = read_binary_ltpiin(args, args.bin, args.binary_traits, str)
            args.conf, args.samp_bin, args.time = LTPI_GHK(args)
            
            args.conf.to_csv(f'{args.out}.config_summary', sep='\t', index=True, header=True, na_rep='NA')
            args.samp_bin.to_csv(f'{args.out}.ltpi_scores.gz', sep='\t', index=True, header=True, compression='gzip', na_rep='NA')

        elif args.con:
            args.ltpiin = read_continuous_ltpiin(args, args.con, None, float)
            args.quantitative_traits = args.ltpiin.columns.to_numpy(dtype='U100')
            
            if args.gencov:
                gencov = read_cov(args.gencov)
                args.pi = gencov.index[0]
                args.mle_traits = np.insert(args.quantitative_traits, 0, args.pi)
                args.GENCOV = gencov.loc[args.mle_traits, args.mle_traits]
            else:
                args.log.log('ValueError: Can\'t read --gencov flag.')
                raise ValueError()
                
            args.ENVCOV = read_cov(args.envcov, args.mle_traits) if args.envcov else transform_dataframe(args.GENCOV)
                
            if args.pick:
                args.selected_traits, args.best_S, args.best_shrinkage, args.summary_data = run_ATSA(args)
                write_r2(args)
                args = update_r2(args)
            
            if args.bout:
                args.conf, args.samp_bin = read_bout(args.bout)
            else:
                args.log.log('ValueError: Can\'t read --bout flag.')
                raise ValueError()

            args.GENCOV = is_pos_def(args, args.GENCOV.copy(), 'GEN')
            args.ENVCOV = is_pos_def(args, args.ENVCOV.copy(), 'ENV')

            args.samp_mle, args.time = LTPI_MLE(args)
            args.samp_mle.to_csv(f'{args.out}.ltpi_scores.gz', sep='\t', index=True, header=True, compression='gzip', na_rep='NA')

        # Handle --pick option without --bin or --con
        elif args.pick:
            if not args.pi:
                args.log.log('ValueError: --pick requires pi column name (--pi [PI]).')
                raise ValueError()

            if args.gencov:
                args.GENCOV = read_cov(args.gencov)
            else:
                args.log.log('ValueError: Can\'t read --gencov flag.')
                raise ValueError()

            # Sanity check: Ensure PI is in the column names of the covariance matrix
            if args.pi not in args.GENCOV.columns:
                args.log.log(f"Error: PI '{args.pi}' is not found in the covariance matrix columns.")
                raise ValueError(f"PI '{args.pi}' is not a valid column name in the covariance matrix.")
            
            args.log.log('Start R2 selection')
            args.selected_traits, args.best_S, args.best_shrinkage, args.summary_data = run_ATSA(args)
            write_r2(args)
        
        else:
            args.log.log('ValueError: No test flag provided.')
            raise ValueError()      

    except Exception:
        args.log.mlog(traceback.format_exc())
        raise
        
    finally:
        args.log.log(f'Analysis finished at {time.ctime()}')
        time_elapsed = round(time.process_time() - start_time, 2)
        args.log.log(f'Total time elapsed: {sec_to_str(time_elapsed)}')