#!/bin/python3

import multiprocessing as mp
import pandas as pd
import numpy as np
import numba
import scipy.stats as stats
from scipy.optimize import minimize
from scipy import stats
import sys, os, time, random

def get_meanTN(p):
    t = stats.norm.ppf(1-p, loc=0, scale=1)    
    meanTN = {
        '0':stats.truncnorm.mean(-np.Inf, t, loc=0, scale=1),
        '1':stats.truncnorm.mean(t, np.Inf, loc=0, scale=1)
    }
    return(meanTN)

@numba.jit(nopython=True)
def LL_fun(x, q, l, c, m):
    '''
    Objective: minimize - LL 
    '''
    X = np.empty(q.shape[0]+1, dtype=q.dtype)
    X[0] = x[0] - m
    X[1:] = q
    n = len(X)
    etas = np.transpose(c).dot(X)
    etasqs = etas ** 2
    LL = -0.5 * (n * np.log(2 * np.pi) + np.sum(np.log(l)) + np.sum(etasqs / l))
    return -LL 

def estim_liab_from_quantitative_traits(q, l, c, m):
    '''
    q: n-1 values of quantitative phenotypes vector
    x: phenotype value estimate of PI
    K: U(gencov)+E(envcov). Its diagonal elements should be ones [1,1,1,...,1]
    '''
    try:
        x0 = np.random.uniform(low=0.0, high=1.0, size=1)
        res = minimize(LL_fun, x0, method='nelder-mead', options={'xatol': 1e-8, 'disp': False}, args=(q, l, c, m))
        return res.x[0]
    except ValueError:
        print("LTPIq: could not find local minimum.")
    except BaseException as err:
        print(f"Unexpected {err=}, {type(err)=}")
        raise    

def eigen_decomposition(K):
    K = (K + K.T)/2  # Ensure K is symmetric
    l, c = np.linalg.eigh(K)
    ind = np.argsort(l)[::-1][:len(l)]
    return l[ind], c[:,ind]

def estimate_NA_GSE(conf):
    mean_value = np.mean(conf.G_SE)
    conf.G_SE.fillna(value = mean_value, inplace=True)
    return(conf)

def argmax_mle(param):
    random.seed(220)
    gencov = param['gencov']
    envcov = param['envcov']
    envcov[0, :] = 0
    envcov[:, 0] = 0
    sample_keys = param['sample_keys']
    q_pheno = param['phenotypes']
    configuration_info = param['configuration_info']

    Nsam = len(sample_keys)
    key_df = pd.Series(sample_keys, index = q_pheno.index)

    # Preallocate sample_estimate DataFrame
    sample_estimate = pd.DataFrame({'LTPI':np.empty([Nsam]), 'CONF':sample_keys}, index = q_pheno.index)
    PI_ind = 0
    h2_T = gencov[PI_ind,PI_ind]
    
    configuration_info_set = configuration_info.loc[np.unique(sample_keys)]
    for idx in q_pheno.index:
        s = q_pheno.loc[idx,:]
        non_null = ~pd.isnull(s)
        cov_inc = np.insert(non_null.to_numpy(bool), 0, True)
        q = s.loc[non_null].to_numpy()
        key = key_df[idx]
        eg = configuration_info_set.loc[key, 'EG']
        d = [1.0]*np.sum(cov_inc)
        d[0] = configuration_info_set.loc[key, 'G_SE']/np.sqrt(h2_T)
        L = np.diag(d)
        gcov = gencov[:,cov_inc][cov_inc,:]
        ecov = envcov[:,cov_inc][cov_inc,:]
        K = L.dot(gcov).dot(L) + ecov
        res = eg + K[1:,0] @ np.linalg.inv(K[1:,1:]) @ q
#         l, c = eigen_decomposition(K)
#         res = np.nan if np.isnan(eg) else estim_liab_from_quantitative_traits(q=q, l=l, c=c, m=eg)
        sample_estimate.loc[idx,'LTPI'] = res

    return sample_estimate

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
    df['obs_LTPI'] = df['LTPI'].copy()
    df['R'] = df['obs_LTPI'].rank(method='average', ascending=True)
    df['LTPI'] = df['R'].apply(rank_to_normal, c=c, n=df.shape[0])
    
#     # Min-Max scaling the 'LTPI' column to match 'EG' min and max
#     df['LTPI'] = (df['LTPI'] - df['LTPI'].min()) * (df['EG'].max() - df['EG'].min()) / (df['LTPI'].max() - df['LTPI'].min()) + df['EG'].min()
    
    return df.reindex(['LTPI', 'CONF', 'obs_LTPI'], axis=1)

def LTPI_MLE(args):
    "LTPI_MLE - dependencies: argmax_mle. "
    
    # Read data
    mle_traits = args.mle_traits
    PI_name = args.pi
    quantitative_traits = args.quantitative_traits
    phenotype_matrix = args.ltpiin
    gencov_mle = np.asarray(args.GENCOV)
    envcov_mle = np.asarray(args.ENVCOV)
    configuration_info = estimate_NA_GSE(args.conf)
    sample_info = args.samp_bin
    run_rint = args.rint
    ncore = args.ncore
    
    # starting time
    time_var = time.time()
    
    IID = phenotype_matrix.index
    sample_key = sample_info.loc[IID,'CONF']

    # maximum likelihood estimation
    param_mle = {
        'PI_name':PI_name,
        'mle_traits':mle_traits,
        'gencov':np.asarray(gencov_mle),
        'envcov':np.asarray(envcov_mle),
        'sample_keys':sample_key,
        'phenotypes':phenotype_matrix,
        'configuration_info':configuration_info
    }
    sample_info = argmax_mle(param_mle)
    
    if run_rint:
        sample_info = RINT(sample_info)

    # end time    
    runtime = time.time() - time_var

    return(sample_info, runtime)