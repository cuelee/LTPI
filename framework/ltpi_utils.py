import numpy as np
import pandas as pd
import warnings

def is_pos_def(args, x, name='', threshold=1e-9):
    """Check if the matrix is positive definite."""
    w, v = np.linalg.eig(x)
    if not np.all(w > 0):
        args.log.log(f'The {name} covariance matrix is not positive definite')
        w[w < 0] = threshold
    return pd.DataFrame(v @ np.diag(w) @ v.T, index=x.columns, columns=x.columns)

def envcov_QC(envcov, gencov):
    """Perform QC for the environmental covariance matrix using the genetic covariance matrix."""
    SIDEC = np.diag(np.sqrt(1 / np.diag(envcov)))
    OMSIDGC = np.diag(np.sqrt((1 - np.diag(gencov))))
    adjusted_envcov = OMSIDGC @ SIDEC @ envcov.values @ SIDEC @ OMSIDGC
    return pd.DataFrame(adjusted_envcov, index=envcov.index, columns=envcov.columns)

def cov_shrink(args, cov, keys):
    """Apply covariance shrinkage if specified."""
    if args.shrink in keys:
        target = args.shrink_target
        res, s = iterated_covariance_optimization(cov, target)
        args.log.log(f'Covariance Shrinkage: Mode-{keys[0]}, Sfactor-{s}, Cond-{np.linalg.cond(res)}')
    else:
        res = cov
    return res

def read_prev(args, f):
    """Read genetic disease prevalence file."""
    df = pd.read_csv(f, sep='\s+', usecols=['TID', 'prev'], dtype={'TID': str, 'prev': float})
    df.set_index('TID', inplace=True)
    valid = df['prev'].notna()
    prevalence = df.loc[valid, 'prev'].to_dict()
    args.log.log(f'Prevalence file contains {len(prevalence)} traits.')
    return prevalence, df.index[valid].to_numpy(dtype='U100')

def read_bout(p):
    """Read the binary test output files."""
    conf = pd.read_csv(f'{p}.config_summary', sep='\t', usecols=['CONF', 'EG', 'G_SE'], dtype={'CONF': str, 'EG': float, 'G_SE': float}, index_col=None)
    conf.set_index('CONF', inplace=True)
    samp = pd.read_csv(f'{p}.ltpi_scores.gz', sep='\s+', usecols=['IID', 'LTPI', 'CONF'], dtype={'IID': str, 'LTPI': float, 'CONF': str}, index_col=None)
    samp.set_index('IID', inplace=True)
    return conf, samp

def read_binary_ltpiin(args, f, TID=None, d=str):
    """Read the binary test input file."""
    col_names = pd.read_csv(f, sep='\t', nrows=0).columns
    df = pd.read_csv(f, sep='\t', dtype={col: d if col != 'IID' else str for col in col_names})
    df.set_index('IID', inplace=True)
    df.fillna('X', inplace=True)
    args.log.log(f'Input contains {df.shape[0]} samples for {df.shape[1]} traits')
    if TID is None:
        return df
    return df.loc[:, TID]

def read_continuous_ltpiin(args, f, TID=None, d=str):
    """Read the quantitative test input file."""
    col_names = pd.read_csv(f, sep='\t', nrows=0).columns
    df = pd.read_csv(f, sep='\t', dtype={col: d if col != 'IID' else str for col in col_names})
    df.set_index('IID', inplace=True)
    df.dropna(axis=0, how='all', inplace=True)
    args.log.log(f'Input contains {df.shape[0]} samples for {df.shape[1]} traits')
    if TID is None:
        return df
    return df.loc[:, TID]

def read_cov(f, traits=None):
    """Read covariance matrix file."""
    cov_raw = pd.read_csv(f, sep='\t', index_col=None)
    cov_raw.index = cov_raw.columns
    if traits is None:
        traits = cov_raw.columns
    return cov_raw.loc[traits, traits]

def transform_dataframe(df):
    """Transform the diagonal elements of a DataFrame by subtracting them from 1."""
    warnings.warn('No environmental covariance matrix input. LTPI will use 1-diag(GENCOV) as a substitute.', category=UserWarning)
    
    diag_values = np.diag(df.values)
    transformed_values = 1 - diag_values
    diag_matrix = np.diag(transformed_values)
    return pd.DataFrame(diag_matrix, index=df.index, columns=df.columns)

def write_r2(args):
    """Write selected traits, shrinkage factors, and summary statistics to output files."""
    if len(args.selected_traits) > 1:
        selected_traits_df = pd.Series(args.selected_traits)
        selected_traits_df.to_csv(f'{args.out}.selected_traits', sep='\t', index=False, header=False, na_rep='NA')
        best_shrinkage_df = pd.Series(args.best_shrinkage)
        best_shrinkage_df.to_csv(f'{args.out}.shrinkage_value', sep='\t', index=False, header=False, na_rep='NA')
        args.best_S.index.name = 'TID'
        args.best_S.to_csv(f'{args.out}.r2_summary', sep='\t', index=True, header=True, na_rep='NA')
        args.log.log(f"--pick summary:\nTraits selected: {args.selected_traits}\n"
                     f"Shrinkage factor: {args.best_shrinkage}\nSummary: {args.best_S}")
    return None

def update_r2(args):
    """Update R2-related attributes in 'args' based on selected traits."""
    if len(args.selected_traits) < 1:
        raise ValueError(f'r2_traits has the size of {len(args.selected_traits)}')
    if args.bin is not None:
        args.binary_traits = args.selected_traits
        args.prev = {k: args.prev[k] for k in args.binary_traits}
        args.GENCOV = args.GENCOV.loc[args.binary_traits, args.binary_traits].copy()
        args.ENVCOV = args.ENVCOV.loc[args.binary_traits, args.binary_traits].copy()
    elif args.con is not None:
        args.quantitative_traits = args.selected_traits[1:]
        args.mle_traits = args.selected_traits
        args.ltpiin = args.ltpiin.loc[:, args.quantitative_traits].copy()
        args.GENCOV = args.GENCOV.loc[args.mle_traits, args.mle_traits].copy()
        args.ENVCOV = args.ENVCOV.loc[args.mle_traits, args.mle_traits].copy()
    return args
