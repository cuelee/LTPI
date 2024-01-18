# #!/bin/python3
# # dependency: numpy,scipy.stats, scipy.optimize,import numpy as np
# from scipy.stats import multivariate_normal as MVN
# import scipy.stats as stats
# import numpy as np
# import pandas as pd
# import sys, os, time, re, random
# import multiprocessing as mp
# from itertools import product
# from decimal import *

# try:
#     from framework.utilities import *
# except ImportError:
#     print('failed to import framework.utilities.py')
#     pass
   
# ### generate multi sampling distributions 
# def generate_P(n, factors, gencov):
#     class Pj(object):
#         '''
#         Pj is a sampling density function of the deterministic importance sampling procedure 
#         The class define the covariance matrix and means of Pj
#         '''
#         def __init__(self, means, cov):
#             self.means = means;
#             self.cov = cov;
#             self.pdf = None;

#     P = [Pj([0]*n, np.diag( [f] * n ).dot( gencov ).dot( np.diag([f] * n) )) for f in factors ];
#     return(P) 

# ### sampling from mixture densities
# def mixture_sampling (N, ns, alpha, P):
#     K = len(alpha); choices = [i for i in range(K)]; count = [0] * K;
#     comp = np.random.choice(choices, N, replace = True, p = alpha); samp = np.empty(shape = [0,ns])
#     for i in range(len(comp)): count[comp[i]] += 1;
#     for j in range(K):
#         Pj_mean = P[j].means; Pj_cov = P[j].cov; Pj_N = count[j];
#         samp_pj = np.random.multivariate_normal(mean = Pj_mean, cov = Pj_cov, size = Pj_N);
#         samp = np.concatenate((samp,samp_pj),axis = 0)
#     return(samp);

# def P_density_estimation (P, input_df):
#     nP = len(P); pdf_P = [];
#     for i in range(nP):
#         Pj_mean = P[i].means;
#         Pj_cov = P[i].cov;
#         Pj_pdf = MVN.pdf( x = input_df, mean = Pj_mean, cov = Pj_cov );
#         pdf_P.append( np.array( Pj_pdf ) );
#     return( pdf_P );

# ### These definitions are necessary for estimating I 
# def estim_cov_tm(pdf_Pj, m):
#     l = len(pdf_Pj);
#     tm_vec = [0]*l;
#     for i in range(l):
#         tm_vec[i] = np.cov(m, pdf_Pj[i])[0][1];
#     return(np.array(tm_vec));

# def estim_cov_t(pdf_Pj, Palpha):
#     l = len(pdf_Pj);
#     t_mat = [];
#     for i in range(l):
#         array = np.array([pdf_Pj[i][j]/Palpha[j] for j in range(len(pdf_Pj[i]))])
#         t_mat.append(array)
#     return(np.cov(np.array(t_mat)));

# def svd_inv(cov_t):
#     u,s,v = np.linalg.svd(cov_t);
#     ds = np.diag([1/s[j] for j in range(len(s)-1)]);
#     us = np.matrix(u)[:,:-1];
#     vs = np.matrix(np.transpose(v))[:,:-1];
#     inv_cov_t = vs.dot(ds).dot(np.transpose(us));
#     return(inv_cov_t) 

# def const_mul(array, pdf_Pj):
#     alist = [];
#     for i in range(len(array)):
#         amul = [array[i] * pdf_Pj[i][j] for j in range(len(pdf_Pj[i]))];
#         alist.append(np.array(amul));
#     return(alist);

# def vector_sum(alist):
#     sumvec = [0]*len(alist[0]);
#     for i in range(len(alist)):
#         sumvec = [sumvec[j] + alist[i][j] for j in range(len(alist[i]))];
#     return(sumvec);

# def key2eg(key, tval, sample_liab, gen_liab_pi):
#     N = len(gen_liab_pi)
#     ind = np.array([True] * N)
#     for i in range(len(key)):
#         test = np.array([True if ind[j] == True else False for j in range(N)])
#         x = sample_liab.iloc[:,i].to_numpy(dtype = float)[test]
#         if key[i] == '1':
#             ind[test] = np.array([x[j] >= tval[i] for j in range(len(x))])
#         else:
#             ind[test] = np.array([x[j] < tval[i] for j in range(len(x))])
#     return(np.array([gen_liab_pi[i] if ind[i] == True else 0.0 for i in range(N)]))

# def key2f(key, tval, sample_liab):
#     N = sample_liab.shape[0]
#     ind = np.array([True] * N)
#     for i in range(len(key)):
#         test = np.array([True if ind[j] == True else False for j in range(N)])
#         x = sample_liab.iloc[:,i].to_numpy(dtype = float)[test]
#         if key[i] == '1':
#             ind[test] = np.array([x[j] >= tval[i] for j in range(len(x))])
#         else:
#             ind[test] = np.array([x[j] < tval[i] for j in range(len(x))])
#     return(np.array([1 if ind[i] == True else 0 for i in range(N)]))    

# def gen_index(keys, samp_keys):
#     dat = dict()
#     n = len(samp_keys)
#     for k in keys:
#         p = k.replace('X','\d')
#         dat[k] = np.array([True if re.match(p,samp_keys[i]) else False for i in range(n)])
#     return(dat)

# def _boostrap(h,s,N_fold = 5):
#     h = np.array(h)
#     s = np.array(s)
#     ind_nz = h != 0
#     ind_z = h == 0
#     set_A = s[ind_nz]
#     set_B = s[ind_z]
#     set_A_ind = np.random.choice(range(N_fold), size = len(set_A), replace = True)
#     set_B_ind = np.random.choice(range(N_fold), size = len(set_B), replace = True)
#     mean_I = np.array([0.0]*N_fold)
#     for i in range(N_fold):
#         set_I = np.concatenate([set_A[set_A_ind==i],set_B[set_B_ind==i]])
#         N_I = len(set_I)
#         mean_I[i] = sum(set_I) / N_I;
#     B_mean = 1 / N_fold * np.sum(mean_I)
#     B_var = 1 / N_fold * np.sum(mean_I[i] - B_mean)**2
#     return(B_var)

# ## Importance sampling with control covariate
# def iscc(h, d_Q, d_P, alpha, Palpha, nPj, N):
#     cov_t = estim_cov_t(d_P, Palpha)
#     inv_cov_t = svd_inv(cov_t);
#     denominator = vector_sum(const_mul(alpha, d_P));
#     m = [h[i] * d_Q[i] / Palpha[i] for i in range(N)]
#     cov_tm = estim_cov_tm(d_P, m)
#     betas = [inv_cov_t.dot(cov_tm)[0,i] for i in range( nPj )];
#     control_variate = vector_sum(const_mul(betas, d_P));
#     nominator = np.array([ h[i] * d_Q[i] - control_variate[i] for i in range(N) ]);
#     beta_sum = np.sum( betas )
#     samp = [nominator[i] / denominator[i] + beta_sum for i in range(N)]
#     I_hat = sum(samp)/N;
#     ## estimate the variance using Boostrap
#     I_var = _boostrap(h,samp)
#     return(I_hat, I_var, np.array(samp, dtype = float))

# ## GenCor and RECor: np.matrix, N: int, outfn: str
# def importance_sampling(
#     U, tval, key2test, 
#     N_importance_sampling = 100000, 
#     QC_thres = 1.0,
#     ncores = 1):
#     '''
#     U: genetic covariance
#     tval: inverse cdf of prevalence
#     QC_thres: Check if standard error of EG is greater than abs(QC_thres*EG)
#     '''
#     print('\nRun importance sampling\n')

#     ### we set random seed 
#     np.random.seed(1)

#     N = N_importance_sampling
#     ns = len(tval)

#     ### set multi processing options 
#     if(ncores == 0):
#         cores = mp.cpu_count() - 1; partitions = cores;
#     else:
#         cores = ncores; partitions = cores;

#     # U: gencov, E: envcov
#     E = np.diag(np.diag(1-U))
    
#     h2_pi = np.diag(U)[0]
    
#     ## LTPI's importance sampling method reqiures probability densities to generate samples. They have means of [0] * n and the covariance matrix of c_Pj * Ce
#     c_Pj = [1, 1.1, 1.2, 1.3, 1.4, 1.5];
#     nPj = len(c_Pj); mean_P = [0] * nPj; alpha = [1 / nPj] * nPj;

#     COV = U + E
#     P = generate_P(ns, c_Pj, COV)

#     gen_liab = mixture_sampling(N, ns, alpha, generate_P(ns, c_Pj, U))
#     env_liab = mixture_sampling(N, ns, alpha, generate_P(ns, c_Pj, E))

#     choices = [i for i in range(nPj)]; count = [0] * nPj;
#     comp = np.random.choice(choices, N, replace = True, p = alpha); 
#     for i in range(len(comp)): count[comp[i]] += 1;
#     k=0; inverse_scaling = np.empty(shape = [N]); 
#     for i in range(nPj): inverse_scaling[k:(count[i]+k)] = (c_Pj[i] * (1-h2_pi) + h2_pi) ; k+=count[i];
    
#     sample_liab = gen_liab + env_liab
#     sample_key = [''.join(['1' if sample_liab[i,j] >= tval[j] else '0' for j in range(ns)]) for i in range(N)]

#     gen_liab_pi = gen_liab[:,0] 
#     index_dict = gen_index(key2test,sample_key)

#     d_Q = MVN.pdf( sample_liab, [0] * ns, COV );
#     d_P = P_density_estimation( P, sample_liab );

#     Palpha = vector_sum(const_mul(alpha,d_P));
#     ones = np.array([1]*N)
    
#     posterior_mean_genetic_liability = pd.DataFrame(index = key2test, columns = ['EG','EG_SE','CDF','CDF_SE','N_valid'])
#     eg_df = posterior_mean_genetic_liability
#     eg_samples = dict()
#     for k in key2test:
#         index = index_dict[k]
#         N_valid = sum(index)
#         eg_df.loc[k,'N_valid'] = N_valid
#         if N_valid <= 5:
#             eg_df.loc[k,'QC'] = False
#             continue

#         cdf_estim, cdf_var, _ = iscc(np.array([1 if i else 0 for i in index]), d_Q, d_P, alpha, Palpha, nPj, N)
#         eg_df.loc[k,'CDF'] = cdf_estim
#         eg_df.loc[k,'CDF_SE'] = np.sqrt(cdf_var) 

#         eg_estim, _ , samp = iscc(np.array([gen_liab_pi[i] if index[i] else 0 for i in range(N)]), d_Q, d_P, alpha, Palpha, nPj, N) 
#         eg_var = sum([(samp[i]/np.sqrt(inverse_scaling[i])-eg_estim)**2 for i in range(N)])/N**2
#         eg_df.loc[k,'EG'] = eg_estim / np.abs(cdf_estim)
#         eg_df.loc[k,'EG_SE'] = np.sqrt( eg_var / np.abs(cdf_estim)**2)
    
#         # _ , _ , samp = iscc(np.array([(gen_liab_pi[i]-eg_df.loc[k,'EG'])**2 / inverse_scaling[i] if index[i] else 0 for i in range(N)]), d_Q, d_P, alpha, Palpha, nPj, N) 
#         _ , _ , samp = iscc(np.array([(gen_liab_pi[i]-eg_df.loc[k,'EG'])**2 / inverse_scaling[i] if index[i] else 0 for i in range(N)]), d_Q, d_P, alpha, Palpha, nPj, N) 
#         eg_df.loc[k,'G_SE'] = np.sqrt(np.mean([samp[i] for i in range(N)]))/np.sqrt(np.abs(cdf_estim))
#         eg_df.loc[k,'QC'] =  (QC_thres > (eg_df.loc[k,'EG_SE']/np.abs(eg_df.loc[k,'EG']))) or (eg_df.loc[k,'EG_SE'] < 0.05)
#         print('estimated posterior score of conf {}: MEAN:{},SE:{},SEM:{},Pass:{}'.format(k, eg_df.loc[k,'EG'], eg_df.loc[k,'G_SE'], eg_df.loc[k,'EG_SE'],eg_df.loc[k,'QC'] ))

#     return(eg_df, eg_samples)

# def mean_eg_imp(param):
#     ### we set random seed 
#     random.seed(10)
    
#     def count_similarity(a,b,c=0):
#         n = len(a)
#         if a[0] != b[0]:
#             return(0)
#         for i in range(n):
#             if a[i] == b[i]:
#                 c+=1
#         return(c)

#     ts = time.time()
#     prev=param['prev']
#     gencov = np.array(param['gencov'])
#     sample_size=param['size_main']
#     sample_size_helper=param['size_helper']
#     ltpiin_keys = param['k2t']
#     iter_max = param['iter_max']
#     ncore=param['ncore']
#     QC_thres=param['QC_thres']
    
#     # first try a raw sample! -- only simulate from trucaated if sem is too big
#     tval = stats.norm.ppf(1-prev, loc=0, scale=1)
#     n = len(tval)
    
#     eg_df, eg_samples = importance_sampling(U=gencov,tval=tval,key2test=ltpiin_keys,N_importance_sampling=sample_size,QC_thres=QC_thres,ncores=ncore)
       
#     helper_keys = eg_df.index[~(eg_df.QC.to_numpy(dtype = bool))]
    
#     mean_eg_helper_On = len(helper_keys) > 0
#     if not mean_eg_helper_On:
#         return(eg_df)
    
#     print('run helper function')
#     factor_keys = eg_df.index[eg_df.QC.to_numpy(dtype = bool)]

#     # estimate factor
#     factor_param = {'gencov':gencov,'tval':tval,'keys':factor_keys,'iter':1,'size':sample_size_helper,'ncore':ncore,'correction_factor':{'ca':1.0,'co':1.0,'X':1.0}}
#     factor_sample = mean_eg_helper(factor_param)
#     h_means = get_sample_mean(factor_sample, factor_keys)
#     imp_means = {k:eg_df.loc[k,'EG'] for k in factor_keys}

#     correction_factor = get_factor(imp_means, h_means, factor_keys)
#     print('genetic liability estimates from TN approx will be adjusted: case:{} control:{}'.format(correction_factor['ca'],correction_factor['co']))
#     for k in helper_keys:
#         eg_samples[k]  = np.array([], dtype = float)
#     sems = {k:0.0 for k in helper_keys}
    
#     i=0;iter_max=100;
#     HELPER_param = {'gencov':gencov,'tval':tval,'keys':helper_keys,'iter':i,'size':sample_size_helper,'ncore':ncore,'correction_factor':correction_factor}
#     while(mean_eg_helper_On):
#         i+=1
#         HELPER_param['iter'] = i
#         add_sample = mean_eg_helper(HELPER_param)
#         for key, values in add_sample.items():
#             if key[0] == '1':
#                 eg_samples[key] = np.append(eg_samples[key],values*correction_factor['ca'])
#             else:
#                 eg_samples[key] = np.append(eg_samples[key],values*correction_factor['co'])
            
#             if len(eg_samples[key]) > 1:
#                 sems[key] = stats.sem(eg_samples[key])
#             else:
#                 sems[key] = 1
#         mean_eg_helper_On, HELPER_param['keys'], _ = test_sem(sems)
#         if(i >= iter_max):
#             mean_eg_helper_On = False
            
#     eg_mean = get_sample_mean(eg_samples, eg_df.index[~(eg_df.QC.to_numpy(dtype = bool))])
#     for k in eg_df.index[~(eg_df.QC.to_numpy(dtype = bool))]:
#         eg_df.loc[k,'EG'] = eg_mean[k]
#         eg_df.loc[k,'EG_SE'] = sems[k]
#         print('{}:{}'.format(eg_df.loc[k,'EG'],eg_df.loc[k,'EG_SE']))
    
#     return(eg_df)

# def LTPI_IMP(args):
#     "LTPI_IMP - dependencies: mean_eg, mean_eg_helper, integer2binary, generate_phecode"
    
#     def sanity_check(ltpiin):
#         "This parse the input matrix in a pandas format into numpy binaryVector matrix and key per individual"

#         def Vector_to_key(l):
#             return(''.join(map(str,l)))

#         def check_Z(Z):
#             if any(~((Z.flatten()=='1') | (Z.flatten()=='0') | (Z.flatten()=='X'))):
#                 raise ValueError('LTPI: LTPI input has an element other than 1 or 0 (or NA)')
#             return(Z)   

#         Z = check_Z(ltpiin.to_numpy(dtype = str, copy = False))
#         keys = np.apply_along_axis(func1d = Vector_to_key, axis=1, arr = Z) # key is a vector of N
#         return({'keys':keys})
    
#     # importance sampling
#     binary_traits = args.binary_traits
#     binary_phenotype = sanity_check(args.ltpiin.loc[:,binary_traits])
#     prev = np.array([args.prev[t] for t in binary_traits], dtype = float)
#     gencov_imp = np.asarray(args.cov)
#     sample_size = args.nsample_main
#     sample_size_helper = args.nsample_helper 
#     QC_thres = args.thres
#     ncore = args.ncore
    
#     key2t = np.unique(binary_phenotype['keys'])

#     # starting time
#     time_var = time.time()
    
#     IID = args.ltpiin.index 
#     param_imp = {
#         'prev':prev,
#         'gencov':gencov_imp,
#         'size_main':sample_size,
#         'size_helper':sample_size_helper,
#         'k2t':key2t,
#         'ncore':ncore,
#         'iter_max':100,
#         'QC_thres':QC_thres
#     }
    
#     configuration_info = mean_eg_imp(param_imp)
#     configuration_info.index.names = ['CONF']
#     sample_info = pd.DataFrame(index = IID)
#     sample_info.index.names = ['IID']
#     sample_info['LTPI'] = np.array([configuration_info.EG[k] for k in binary_phenotype['keys']], dtype = float)
#     sample_info['CONF'] = np.array(binary_phenotype['keys'], dtype = str)

#     # end time    
#     runtime = time.time() - time_var

#     return(configuration_info,sample_info,runtime)