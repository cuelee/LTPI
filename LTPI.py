#!/usr/bin/env python3

"""
LT-PI: Liability Threshold Model for Disease Risk Prediction
Copyright (C) 2024 Cue Hyunkyu Lee

LT-PI is a statistical framework that estimates the conditional expected genetic liability of a target disease by leveraging phenotypic data from multiple related traits. Built on a liability threshold model, LT-PI calculates disease probabilities using phenotypic information from electronic health records (EHR), combined with heritability and prevalence data from the literature. The model integrates binary and continuous traits using the GHK algorithm and maximum likelihood estimation.
"""

import numpy as np
import pandas as pd
import argparse
import time, sys
import warnings
import traceback

from framework.est_maxlikelihood import LTPI_MLE
from framework.ghk_algorithm import LTPI_GHK, iterated_covariance_optimization
from framework.r2_selection import run_ATSA
from framework.ltpi_utils import (
    is_pos_def, envcov_QC, cov_shrink, read_prev, read_bout,
    read_binary_ltpiin, read_continuous_ltpiin, read_cov,
    transform_dataframe, write_r2, update_r2, put_pi_first
)

codename = 'LTPI'
__version__ = '1.2'
MASTHEAD = f"""
************************************************************
* LTPI ({codename})
* Version {__version__}
* (C) 2024-2026 Cue H. Lee, Columbia University
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

# Output
parser.add_argument('--out', default='LTPI', type=str,
    help='Output file prefix (default: LTPI).')

# Modes
parser.add_argument('--bin', type=str,
    help='Path to binary phenotype input.')
parser.add_argument('--con', type=str,
    help='Path to continuous phenotype input.')

# Inputs
parser.add_argument('--prevalence', type=str,
    help='Path to disease prevalence file (required for --bin).')
parser.add_argument('--bout', type=str,
    help='Prefix of binary LTPI output (required for --con).')

# Target trait
parser.add_argument('--pi', type=str,
    help='Target phenotype (PI). Required for all modes. Will be placed first and used as prediction target.')

# Trait selection
parser.add_argument('--pick', action='store_true',
    help='Enable trait selection using ATSA.')
parser.add_argument('--Q', default=30, type=int,
    help='Number of non-target traits to select (default: 30).')
parser.add_argument('--r2', default=0.0, type=float,
    help='R2 threshold for ATSA (default: 0).')

# Covariance
parser.add_argument('--gencov', type=str,
    help='Path to genetic covariance matrix (required).')
parser.add_argument('--envcov', type=str,
    help='Path to environmental covariance matrix (optional).')
parser.add_argument('--shrink', type=str,
    help='Covariance shrinkage mode: G, E, or B.')
parser.add_argument('--shrink_target', default=1.5, type=float,
    help='Target condition number for shrinkage (default: 1.5).')

# GHK
parser.add_argument('--nsample', default=50000, type=int,
    help='Number of GHK samples (default: 50K).')

# Optional
parser.add_argument('--rint', action='store_true',
    help='Apply rank-based inverse normal transform (if supported).')

def require_pi(args):
    if not args.pi:
        args.log.log('ValueError: --pi is required.')
        raise ValueError('--pi is required.')
    return args.pi

def validate_pi_in_columns(df, pi_name, df_name):
    if pi_name not in df.columns:
        raise ValueError(f"PI '{pi_name}' not found in {df_name} columns")

def validate_same_traits(df1, df2, name1, name2):
    cols1 = list(df1.columns)
    cols2 = list(df2.columns)
    if set(cols1) != set(cols2):
        only1 = sorted(set(cols1) - set(cols2))
        only2 = sorted(set(cols2) - set(cols1))
        msg = (
            f"Trait mismatch between {name1} and {name2}. "
            f"Only in {name1}: {only1}. Only in {name2}: {only2}."
        )
        raise ValueError(msg)
        
def setup_header(args):
    defaults = vars(parser.parse_args([]))
    opts = vars(args)
    non_defaults = [x for x in defaults.keys() if opts.get(x) != defaults[x]]

    header = MASTHEAD
    header += f"Call:\n{' '.join(sys.argv)}\n"
    options = ['--' + x.replace('_', '-') + ' ' + str(opts[x]) + ' \\' for x in non_defaults]
    header += '\n'.join(options).replace('True', '').replace('False', '') + '\n'
    return header

def require_arg(args, arg_name, message=None):
    value = getattr(args, arg_name)
    if value is None:
        if message is None:
            message = f"Missing required argument --{arg_name}"
        args.log.log(f'ValueError: {message}')
        raise ValueError(message)
    return value

def load_gencov_and_set_pi(args, traits=None):
    require_arg(args, 'gencov', "Can't read --gencov flag.")
    require_pi(args)

    cov = read_cov(args.gencov, traits)
    validate_pi_in_columns(cov, args.pi, 'GENCOV')

    cov, order = put_pi_first(cov, args.pi)
    return cov, order

def load_envcov(args, reference_cov, traits=None):
    if args.envcov:
        envcov = read_cov(args.envcov, traits)
        validate_same_traits(reference_cov, envcov, 'GENCOV', 'ENVCOV')
        validate_pi_in_columns(envcov, args.pi, 'ENVCOV')
    else:
        envcov = transform_dataframe(reference_cov)

    envcov, _ = put_pi_first(envcov, args.pi)
    return envcov

def maybe_run_pick(args):
    if args.pick:
        args.selected_traits, args.best_S, args.best_shrinkage, args.summary_data = run_ATSA(args)
        write_r2(args)
        args = update_r2(args)
    return args

def run_binary_mode(args):
    require_pi(args)
    require_arg(args, 'prevalence', "Can't read --prevalence flag.")

    args.prev, args.binary_traits = read_prev(args, args.prevalence)

    args.GENCOV, _ = load_gencov_and_set_pi(args, args.binary_traits)
    args.ENVCOV = load_envcov(args, args.GENCOV, args.binary_traits)

    args.ltpiin_bin = read_binary_ltpiin(args, args.bin, args.binary_traits, str)
    validate_pi_in_columns(args.ltpiin_bin, args.pi, 'binary phenotype input')
    args.ltpiin_bin, _ = put_pi_first(args.ltpiin_bin, args.pi)

    args = maybe_run_pick(args)

    args.GENCOV = is_pos_def(args, cov_shrink(args, args.GENCOV.copy(), keys=['G', 'B']), 'GEN')
    args.ENVCOV = is_pos_def(args, cov_shrink(args, args.ENVCOV.copy(), keys=['E', 'B']), 'ENV')

    args.log.log(
        f'Condition Number:\n'
        f'GEN-{np.linalg.cond(args.GENCOV)} '
        f'ENV-{np.linalg.cond(args.ENVCOV)} '
        f'GEN+ENV-{np.linalg.cond(args.GENCOV + args.ENVCOV)}'
    )

    args.conf, args.samp_bin, args.time = LTPI_GHK(args)

    args.conf.to_csv(
        f'{args.out}.config_summary',
        sep='\t',
        index=True,
        header=True,
        na_rep='NA'
    )
    args.samp_bin.to_csv(
        f'{args.out}.ltpi_scores.gz',
        sep='\t',
        index=True,
        header=True,
        compression='gzip',
        na_rep='NA'
    )

def run_continuous_mode(args):
    require_pi(args)

    args.ltpiin_con = read_continuous_ltpiin(args, args.con, None, float)

    gencov, _ = load_gencov_and_set_pi(args)

    # --- enforce PI exists in covariance ---
    validate_pi_in_columns(gencov, args.pi, 'GENCOV')

    # --- intersect ONLY continuous traits ---
    continuous_traits = [
        t for t in args.ltpiin_con.columns
        if t in gencov.columns and t != args.pi
    ]

    ordered_traits = [args.pi] + continuous_traits
    args.GENCOV = gencov.loc[ordered_traits, ordered_traits]
    args.ltpiin_con = args.ltpiin_con.loc[:, continuous_traits]

    args.mle_traits = np.array(ordered_traits)
    args.quantitative_traits = np.array(continuous_traits)

    args.ENVCOV = load_envcov(args, args.GENCOV, ordered_traits)

    args = maybe_run_pick(args)

    require_arg(args, 'bout', "Can't read --bout flag.")
    args.conf, args.samp_bin = read_bout(args.bout)

    args.GENCOV = is_pos_def(args, args.GENCOV.copy(), 'GEN')
    args.ENVCOV = is_pos_def(args, args.ENVCOV.copy(), 'ENV')

    args.samp_mle, args.time = LTPI_MLE(args)

    args.samp_mle.to_csv(
        f'{args.out}.ltpi_scores.gz',
        sep='\t',
        index=True,
        header=True,
        compression='gzip',
        na_rep='NA'
    )

def run_pick_only_mode(args):
    require_pi(args)
    require_arg(args, 'gencov', "Can't read --gencov flag.")

    args.GENCOV = read_cov(args.gencov)
    validate_pi_in_columns(args.GENCOV, args.pi, 'GENCOV')
    args.GENCOV, _ = put_pi_first(args.GENCOV, args.pi)

    args.log.log('Start R2 selection')
    args.selected_traits, args.best_S, args.best_shrinkage, args.summary_data = run_ATSA(args)
    write_r2(args)

if __name__ == '__main__':

    args = parser.parse_args()
    if not args.out:
        raise ValueError('No output file prefix (--out) provided.')

    log = None
    start_time = None

    try:
        log = Logger(f'{args.out}.log')
        args.log = log

        header = setup_header(args)
        args.log.log(header)
        args.log.log(f'Beginning analysis at {time.ctime()}')
        start_time = time.process_time()

        if args.bin:
            run_binary_mode(args)

        elif args.con:
            run_continuous_mode(args)

        elif args.pick:
            run_pick_only_mode(args)

        else:
            args.log.log('ValueError: No test flag provided.')
            raise ValueError('No test flag provided.')

    except Exception:
        if log is not None:
            log.mlog(traceback.format_exc())
        raise

    finally:
        if log is not None:
            log.log(f'Analysis finished at {time.ctime()}')
            if start_time is not None:
                time_elapsed = round(time.process_time() - start_time, 2)
                log.log(f'Total time elapsed: {sec_to_str(time_elapsed)}')