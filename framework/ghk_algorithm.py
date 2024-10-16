#!/bin/python3

# dependencies
import numpy as np
import pandas as pd
import sys
import os
import time
import re
import random
import multiprocessing as mp
import warnings
import scipy.stats as stats
import scipy.linalg as linalg
from scipy import optimize 
from scipy.stats import multivariate_normal as MVN
from sklearn.covariance import shrunk_covariance
from decimal import *
from itertools import product

# Define helper function
def make_pos_def(x, name = '', thres=1e-9):
    w, v = np.linalg.eig(x)
    if not np.all(w > 0):
        w[w < 0] = thres
    return np.dot(v, np.dot(np.diag(w), v.T))   

def optimize_covariance(cov, x0, t=10):
    def objective_function(x, cov_matrix, target):
        shrinkage_value = x if isinstance(x, float) else x[0]
        reg_cov = make_pos_def(shrunk_covariance(cov_matrix, shrinkage=shrinkage_value))
        cond = np.linalg.cond(reg_cov)
        return (np.abs(cond) - target) ** 2

    res = optimize.minimize(fun=objective_function, x0=x0, args=(cov, t), 
                            method='L-BFGS-B', bounds=[(0.0001, 0.9999)])
    
    # Ensure shrinkage is a float
    s = res.x if isinstance(res.x, float) else res.x[0]

    # Regularize the covariance matrix using the best regularization parameter
    res_arr = shrunk_covariance(cov, shrinkage=s)
    
    return res_arr, s

def iterated_covariance_optimization(cov, t=10, iterations=10):
    # Save original type, index, and columns if input is a DataFrame
    original_type = type(cov)
    original_index = cov.index if isinstance(cov, pd.DataFrame) else None
    original_columns = cov.columns if isinstance(cov, pd.DataFrame) else None

    # Convert DataFrame to numpy array if necessary
    cov = cov.values if isinstance(cov, pd.DataFrame) else cov

    x0_values = np.linspace(0.1, 0.9, num=iterations)  # compute different initial guesses
    
    best_diff = float('inf')
    best_matrix = None
    best_s = None

    for x0 in x0_values:
        res_arr, s = optimize_covariance(cov, x0, t)

        diff = np.abs(np.linalg.cond(res_arr) - t)
        if diff < best_diff:
            best_diff = diff
            best_matrix = res_arr
            best_s = s

    # Convert back to DataFrame if original input was a DataFrame
    if original_type == pd.DataFrame:
        best_matrix = pd.DataFrame(best_matrix, index=original_index, columns=original_columns)

    return best_matrix, best_s

def ldlt_decomposition(A):
    def is_pos(x, thres=1e-6):
        w, v = np.linalg.eig(x)
        if not np.all(w > 0):
            w[w < 0] = thres
        return v.dot(np.diag(w)).dot(np.transpose(v)) 

    l, d, perm = linalg.ldl(is_pos(A), lower=True, hermitian=True, overwrite_a=False, check_finite=True)
    d[d < 0] = 1e-18
    L = l.dot(np.sqrt(d))
    T_ind = [i for i in range(A.shape[0]) if perm[i] == 0][0]
    return np.real(L), perm, T_ind

def estimate_ab(x, p):
    def get_thresholds(loc, scale):
        return stats.norm.ppf(loc, loc=0, scale=1), stats.norm.ppf(scale, loc=0, scale=1)
    
    mapping = {
        '1': (lambda: get_thresholds(1-p, 1)),
        '0': (lambda: get_thresholds(0, 1-p)),
        'X': (lambda: get_thresholds(0, 1)),
    }
    
    if x in mapping:
        a, b = mapping[x]()
    else:
        raise ValueError('Critical error: Cannot read Disease Prevalence')

    return {'a': a, 'b': b}

def estim_su(S, u, i):
    # This is faster than below
    su = S[i, :i] @ u[:, :i].T
    return su

# Function to estimate U
def sample_U(S, omega_S, N_samp, ns, key, prev, mu=None):
    # Initialize array
    u = np.zeros((N_samp, ns))

    # Calculation for the first element
    res_0 = estimate_ab(x=key[0], p=prev[0])
    low_p = stats.norm.cdf(x=(res_0['a'] - mu[0]) / S[0, 0], loc=0, scale=1)
    high_p = stats.norm.cdf(x=(res_0['b'] - mu[0]) / S[0, 0], loc=0, scale=1)
    p_0 = stats.uniform.rvs(loc=low_p, scale=high_p - low_p, size=N_samp)
    u[:, 0] = stats.norm.ppf(p_0, loc=0, scale=1)

    # Calculation for the rest of the elements
    for i in range(1, ns):
        res_i = estimate_ab(x=key[i], p=prev[i])
        su = estim_su(S, u, i)
        low_p_i = stats.norm.cdf(x=((res_i['a'] - mu[i] - su) / S[i, i]), loc=0, scale=1)
        high_p_i = stats.norm.cdf(x=((res_i['b'] - mu[i] - su) / S[i, i]), loc=0, scale=1)
        p_i = stats.uniform.rvs(loc=low_p_i, scale=(high_p_i - low_p_i), size=N_samp)
        u[:, i] = stats.norm.ppf(p_i, loc=0, scale=1)

    # Compute liabilities
    liab = np.dot(S, u.T).T
    genliab = np.dot(omega_S, u.T).T
        
    # Return dictionary of results
    return {
        'liab': liab, 'genliab': genliab, 'omega_S': omega_S, 'S': S, 'u': u, 
        'prev': prev, 'N_samp': N_samp, 'key': key, 'ns': ns
    }

# Function to estimate w
def estim_w(data, mu=None):
    # Extract data
    S, u, key, prev, ns = data['S'], data['u'], data['key'], data['prev'], data['ns']

    # Calculation for the first element
    res_0 = estimate_ab(x=key[0], p=prev[0])
    w = np.array(
        stats.norm.cdf(x=((res_0['b'] - mu[0]) / S[0, 0]), loc=0, scale=1) -
        stats.norm.cdf(x=((res_0['a'] - mu[0]) / S[0, 0]), loc=0, scale=1)
    )

    # Calculation for the rest of the elements
    for i in range(1, ns):
        res_i = estimate_ab(x=key[i], p=prev[i])
        su = estim_su(S, u, i)
        w = w * (
            stats.norm.cdf(x=((res_i['b'] - mu[i] - su) / S[i, i]), loc=0, scale=1) -
            stats.norm.cdf(x=((res_i['a'] - mu[i] - su) / S[i, i]), loc=0, scale=1)
        )
    
    return w

# Monte Carlo simulation function
def monte_carlo(key2test, prev, eg_df, U, E):
    Nsam = 1000000
    tval = stats.norm.ppf(1-prev[0], loc=0, scale=1)

    for k in key2test:
        # Generate liabilities
        genliab = stats.norm.rvs(loc=0, scale=np.sqrt(U), size=Nsam)
        envliab = stats.norm.rvs(loc=0, scale=np.sqrt(E), size=Nsam)
        liab = genliab + envliab

        # Condition based on key
        if k == '1':
            ind = liab >= tval
        elif k == '0':
            ind = liab < tval
        elif k == 'X':
            ind = np.array([True] * Nsam)
        else:
            raise ValueError(f'GHK: Invalid configuration code found - {k}')

        # Compute and store results
        eg_df.loc[k, 'EG'] = np.mean(genliab[ind])
        eg_df.loc[k, 'G_SE'] = np.sqrt(np.var(genliab[ind]))

    return eg_df

# Function to estimate L from a covariance matrix
def estimate_L(COV):
    max_retries = 5  # Maximum number of retries before failing
    for attempt in range(max_retries):
        try:
            L = np.linalg.cholesky(COV)
            return L, 'NP:Cholesky'
        except:
            # If the error occurs, slightly adjust random.seed(i)
            if attempt < max_retries - 1:  # Be slient until it reaches the final attempt
                COV = make_pos_def(COV)
            else:  # This block will be executed if the loop completed all iterations without breaking
                raise ValueError("Error: GHK failed to estimate L %s"%k)
#                 try: ## I will later implement LDLT decomposition
#                     L = ldlt_decomposition(COV)
#                     return L, 'Scipy:LDLT'
#                 except:
#                     raise ValueError('Fatal Error: matrix decomposition of the genetic covariance matrix failed')
    return 'ERROR', 'ERROR'

# Function to estim posterior mean genetic liability and its standard error
def estim_eg(ns, S, omega_S, N_samp, k, prev):
    mu = np.array([0.0] * ns)
    eg_data = sample_U(S=S, omega_S=omega_S, N_samp=N_samp, ns=ns, key=k, prev=prev, mu=mu)
    eg_w = estim_w(eg_data, mu=mu)
    eg_T = 1 / N_samp * np.sum(eg_w)

    mu = np.array([1 / eg_T / N_samp * np.sum(eg_data['genliab'][:, i] * eg_w) for i in range(ns)])
    mu_liab = np.array([1 / eg_T / N_samp * np.sum(eg_data['liab'][:, i] * eg_w) for i in range(ns)])
    EG = mu[0]

    gse_data = sample_U(S=S, omega_S=omega_S, N_samp=N_samp, ns=ns, key=k, prev=prev, mu=mu_liab)
    gse_w = estim_w(gse_data, mu=mu_liab)
    gse_T = 1 / N_samp * np.sum(gse_w)
    G_SE = np.sqrt(1 / gse_T / N_samp * np.sum((gse_data['genliab'][:, 0]) ** 2 * gse_w))
    return(EG, G_SE)

def gen_Lsub(U, E, omega_S, max_retries = 5):
    E_cond = np.linalg.cond(E)
    U_cond = np.linalg.cond(U)
    if E_cond > 9:
        ts=[1,2,4]
        Ssub = {i:estimate_L(U+iterated_covariance_optimization(E,t=ts[i],iterations=10)[0])[0] for i in range(max_retries - 2)}
    else:
        E_indp = np.diag(1-np.diag(U))
        Ssub = {i:estimate_L(U + E_indp)[0] for i in range(max_retries-2)}
    if U_cond > 10:
        omega_Ssub,_ = estimate_L(iterated_covariance_optimization(U,t=10,iterations=10)[0])
    else: 
        omega_Ssub = omega_S
    return Ssub, omega_Ssub
    
def ghk_algorithm(
    log, U, E, prev, key2test, 
    N_samp=100000):
    """
    U: genetic covariance
    tval: inverse cdf of prevalence
    """

    log.log('\nRunning GHK algorithm...\n')

    ns = len(prev)
    
    if np.sum(prev < 0.01) > 0:
        warnings.warn(
            'The input contains traits with prevalence less than 0.01.'
            ' This can cause problems in the estimation process.'
        )

    # Perform MonteCarlo if Qb is 1
    if ns == 1:
        posterior_mean_genetic_liability = pd.DataFrame(index=key2test, columns=['EG', 'G_SE'])
        eg_df = monte_carlo(key2test, prev, posterior_mean_genetic_liability, U, E)
        return eg_df

    # U: gencov, E: envcov
    cov = U + E
    S, M_S = estimate_L(cov)
    omega_S, M_O = estimate_L(U)
    
    Ssub, omega_Ssub = gen_Lsub(U, E, omega_S)
    
    log.log(f'Decomposition method: S-{M_S}, U-{M_O}')

    posterior_mean_genetic_liability = pd.DataFrame(index=key2test, columns=['EG', 'G_SE'])
    eg_df = posterior_mean_genetic_liability

    for k in key2test:
        _S = S
        _omega_S = omega_S 
        
        max_retries = 5  # Maximum number of retries before failing
        for attempt in range(max_retries):
            try:
                eg_df.loc[k, 'EG'], eg_df.loc[k, 'G_SE'] = estim_eg(ns, _S, _omega_S, N_samp, k, prev)
                break  # If this line is reached, no error was raised, so exit the loop
            except ValueError:
                # If the error occurs, slightly adjust random.seed(i)
                if attempt < max_retries - 2:  # Be slient until it reaches the final attempt
                    if attempt < 3:
                        _S = Ssub[attempt]
                    elif attempt == 3:
                        _omega_S = omega_Ssub
                else:  # This block will be executed if the loop completed all iterations without breaking
                    raise ValueError("Error: GHK failed to process variant %s"%k)
        
        log.log(f'Estimated posterior score of conf {k}: MEAN:{eg_df.loc[k,"EG"]}, SE:{eg_df.loc[k,"G_SE"]}')

    return eg_df

def RINT(df, c=3.0/8):
    def rank_to_normal(rank, c, n):
        """Convert rank to normal distribution using specified constant."""
        return stats.norm.ppf((rank - c) / (n - 2*c + 1))
    """
    Processes the LTPIb output using RINT_D and returns the processed DataFrame 
    with additional columns 'R' and 'LTPI' and a summary DataFrame with mean 'EG' 
    and unique 'LTPI' values for each 'CONF'.
    
    :param df: Input DataFrame with at least columns 'CONF' and 'EG'.
    :return: A tuple containing the processed DataFrame with new columns 'R' and 'LTPI', 
             and a summary DataFrame.
    """
    df['R'] = df['EG'].rank(method='average', ascending=True)
    df['LTPI'] = df['R'].apply(rank_to_normal, c=c, n=df.shape[0])
      
    summary = df.groupby('CONF')[['EG', 'LTPI']].first()
    
    return df.reindex(['EG', 'LTPI', 'CONF'], axis=1), summary

def LTPI_GHK(args):
    """
    Estimates the posterior mean genetic liability using the Geweke-Hajivassiliou-Keane (GHK) Algorithm.

    Parameters: 
    - args: dictionary of arguments necessary for the GHK estimation.

    Returns: 
    - configuration_info: DataFrame of unique configurations from GHK algorithm.
    - sample_info: DataFrame with LTPI scores and CONF per sample.
    - runtime: execution time of the process.
    """
    def sanity_check(df):
        """Check if the input DataFrame has valid elements."""
        if not isinstance(df, pd.DataFrame):
            raise ValueError('Input must be a Pandas DataFrame.')

        valid_elements = set(['0', '1', 'X'])
        Z = df.values.astype(str)
        invalid_elements = set(Z.flatten()) - valid_elements
        if invalid_elements:
            raise ValueError('Input has an element other than 0, 1, or X')

        keys = np.apply_along_axis(lambda l: ''.join(l.astype(str)), axis=1, arr=Z)
        return {'keys': keys}

    # Perform sanity check
    binary_traits = args.binary_traits
    binary_phenotype = sanity_check(args.ltpiin.loc[:, binary_traits])
    
    # Prepare parameters for GHK
    log = args.log
    prev = np.array([args.prev[t] for t in binary_traits], dtype=float)
    gencov = np.asarray(args.GENCOV)
    envcov = np.asarray(args.ENVCOV)
    sample_size = args.nsample
    run_rint = args.rint
    key2t = np.unique(binary_phenotype['keys'])
    
    # Start timing
    start_time = time.time()
        
    # Generate samples using GHK algorithm
    samples_df = ghk_algorithm(log=log, U=gencov, E=envcov, prev=prev, key2test=key2t, N_samp=sample_size)
    samples_df.index.names = ['CONF']
    
    # Prepare sample information
    sample_info = pd.DataFrame(index=args.ltpiin.index)
    sample_info.index.names = ['IID']
    sample_info['EG'] = np.array([samples_df.EG[k] for k in binary_phenotype['keys']], dtype=float)
    sample_info['CONF'] = np.array(binary_phenotype['keys'], dtype=str)
    
    if run_rint:
        sample_info, summary = RINT(sample_info)
        samples_df.loc[samples_df.index,'LTPI'] = summary.loc[samples_df.index,'LTPI']
        samples_df.rename(columns={'EG':'EG_obs','LTPI':'EG'}, inplace = True)
    else:
        sample_info.rename(columns={'EG':'LTPI'}, inplace = True)
    
    # Calculate runtime
    runtime = time.time() - start_time

    return samples_df, sample_info, runtime