import multiprocessing as mp
import pandas as pd
import numpy as np
import scipy.stats as stats
import sys, os, time
import random
try:
    from framework.utilities import *
except ImportError:
    print('failed to import framework.utilities.py')
    pass

def mean_eg_MVN(param):
    ts = time.time()
    prev=param['prev']
    gencov = np.array(param['gencov'])
    envcov = np.diag(1 - np.diag(gencov))
    sample_size=param['size_main']
    sample_size_helper=param['size_helper']
    ltpiin_keys = param['k2t']
    iter_max = param['iter_max']
    ncore=param['ncore']
    random.seed(10)
    
    # first try a raw sample! -- only simulate from trucaated if sem is too big
    h2_vec = np.abs(np.diag(gencov))
    h2_pi = h2_vec[0]
    tval = stats.norm.ppf(1-prev, loc=0, scale=1)
    n = len(tval)
    #######
    
    mp_nsample = [int(sample_size/ncore)]*ncore
    sample_size = sum(mp_nsample)

    global genliab_pi, sample_key

    genliab_pi = np.array([],dtype = float)
    sample_key = np.array([],dtype = str)
    
    if(ncore == 1):
        gen_liab = stats.multivariate_normal.rvs(mean=[0]*n, cov = gencov, size= sample_size)
        env_liab = stats.multivariate_normal.rvs(mean=[0]*n, cov = envcov, size= sample_size)
        liab = gen_liab + env_liab
        genliab_pi = gen_liab[:,0]
        del gen_liab, env_liab
        sample_key = [''.join(['1' if liab[i,j] >= tval[j] else '0' for j in range(n)]) for i in range(sample_size)]
        del liab
    else:    
        print('CPU counts:{}'.format(ncore))
        pool = mp.Pool(mp.cpu_count())
        for i in range(ncore):
            args=(n, gencov + envcov, h2_pi, tval, mp_nsample[i])
            pool.apply_async(sampling_from_MVN, args = args, callback=get_result)
        pool.close()
        pool.join()
    
    print('Time in MC-MVN step:', time.time() - ts)   
    ts = time.time()
    if (len(genliab_pi) < 1 | len(sample_key) < 1):
        raise ValueError('MonteCarlo simulation: Number of samples < 1')
    keys_unique = np.unique(np.concatenate((ltpiin_keys, sample_key))    )
    eg_samples = {key: [] for key in set(keys_unique)}
    for i in range(len(sample_key)):
        eg_samples[sample_key[i]] += [genliab_pi[i]]
    sems = dict()
    for key in keys_unique:
        if len(eg_samples[key]) > 1:
            sems[key] = stats.sem(eg_samples[key])
        else:
            sems[key] = 1
    del genliab_pi, sample_key
    print('Time in parsing step:', time.time() - ts)    
    
    mean_eg_helper_On, mean_eg_helper_keys, factor_keys = test_sem(sems)
    correction_factor = None
    ts = time.time()
    i=0
    HELPER_param = {'gencov':gencov,'tval':tval,'keys':mean_eg_helper_keys,'iter':i,'size':sample_size_helper,'ncore':ncore}
    while(mean_eg_helper_On):
        print('LTPI Detected \'SEM > 0.01\' from MC step')

        # estimate factor
        factor_param = {'gencov':gencov,'tval':tval,'keys':factor_keys,'iter':-9,'size':sample_size_helper,'ncore':ncore}
        factor_sample = mean_eg_helper(factor_param)

        MC_means = get_sample_mean(eg_samples, factor_keys)
        h_means = get_sample_mean(factor_sample, factor_keys)
        correction_factor = get_factor(MC_means, h_means, factor_keys)
        print('genetic liability estimates from TN approx will be adjusted: case:{} control:{}'.format(correction_factor['ca'],correction_factor['co']))
        
        i+=1
        add_sample = mean_eg_helper(HELPER_param)
        for key, values in add_sample.items():
            if key[0] == '1':
                eg_samples[key] = np.append(eg_samples[key],values*correction_factor['ca'])
            else:
                eg_samples[key] = np.append(eg_samples[key],values*correction_factor['co'])
            
            if len(eg_samples[key]) > 1:
                sems[key] = stats.sem(eg_samples[key])
            else:
                sems[key] = 1
        mean_eg_helper_On, HELPER_param['keys'], _ = test_sem(sems)
        if(i >= iter_max):
            mean_eg_helper_On = False        
    
    print('Time in help step:', time.time() - ts)    
    return(sems, eg_samples, correction_factor)

def LTPI_MVN(args):
    "LTPI main - dependencies: mean_eg, mean_eg_helper, integer2binary, generate_phecode"
    trait_names = args.prev.keys()
    
    prev = np.array([args.prev[t] for t in trait_names])
    gencov = np.asarray(args.cov.to_numpy(dtype = float))
    res = generate_ltpiin(args.ltpiin.loc[:,trait_names])
    key2t = np.unique(res['keys'])
    sample_size = args.nsample_main
    sample_size_helper = args.nsample_helper 
    ncore = args.ncore
    
    # starting time
    time_var = time.time()
    
    IID = args.ltpiin.index 
    param = {'prev':prev,'gencov':gencov,'size_main':sample_size,'size_helper':sample_size_helper,'k2t':key2t,'ncore':ncore, 'iter_max':100}
    res['sems'], res['eg_samples'], res['correction_factor'] = mean_eg_MVN(param)
    res['eg_means'] = get_sample_mean(res['eg_samples'])
    res['liab'] = [res['eg_means'][k] for k in res['keys']]
    
    # end time
    res['run_time'] = time.time() - time_var

    # generate output
    configuration_info = pd.DataFrame(index = key2t)
    configuration_info.index.names = ['CONF']
    configuration_info['pmliab'] = np.array([res['eg_means'][k] for k in key2t], dtype = float)
    configuration_info['sems'] = np.array([res['sems'][k] for k in key2t], dtype = float)
    sample_info = pd.DataFrame(index = IID.to_numpy())
    sample_info.index.names = ['IID']
    sample_info['liab'] = res['liab']
    sample_info['conf'] = res['keys']
    runtime = res['run_time']
    correction_factor = res['correction_factor']

    return(configuration_info,sample_info,runtime, res['eg_samples'])

if(__name__ == "__main__"):   
    print(__name__) 