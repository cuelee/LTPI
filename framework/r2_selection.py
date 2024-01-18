#!/bin/python3

import numpy as np
import pandas as pd
import warnings
from scipy import stats,optimize
from sklearn.model_selection import GridSearchCV
from sklearn.covariance import ShrunkCovariance, shrunk_covariance

def contribution_based_greedy_algorithm(gencov = None, 
                                        PI = None, 
                                        N = 30, 
                                        shrinkages = np.round(np.linspace(0.01, 0.99, num = 12),2)):
    
    def gen_R2(gencov, PI, O):
        Omega_O = gencov.loc[O, O].to_numpy()
        Omega_OT = gencov.loc[O, [PI]].to_numpy()
        Omega_OTP = gencov.loc[[PI], O].to_numpy()
        h2_T = float(gencov.loc[PI, PI])

        cond_O = np.log10(np.linalg.cond(Omega_O))
        Omega_O_inv = np.linalg.pinv(Omega_O)

        tau_r = h2_T - Omega_OTP.dot(Omega_O_inv).dot(Omega_OT)
        tau = tau_r/h2_T
        R2 = 1-tau
        return(float(R2), cond_O)    
    
    def gencor_filter(gencov, thres = 0.3):
        h2 = np.diag(gencov)
        inv_sqrt_h2 = np.sqrt(1/h2)
        gencor = pd.DataFrame(inv_sqrt_h2.dot(gencov).dot(inv_sqrt_h2), index = gencov.index, columns = gencov.columns)
        ind = gencor.loc[:,PI] >= 0.3 
        return(gencov.loc[ind,ind])
    
    def gencov_max(gencov):
        count = pd.Series(0, index = gencov.index )
        match = pd.Series('', index = gencov.index )
        for c in gencov.columns:
            for i in gencov.index:
                if c != i:
                    c_sqrt_var = np.sqrt(gencov.loc[c,c])
                    i_sqrt_var = np.sqrt(gencov.loc[i,i])
                    ic_max_gencov = c_sqrt_var * i_sqrt_var
                    if np.abs(ic_max_gencov) < np.abs(gencov.loc[i,c]):
                        count[i] += 1
                        match[i] = c
                        gencov.loc[i,c] = ic_max_gencov
        sorted_count = count.sort_values()
        remove = []
        for i in sorted_count.index:
            if (sorted_count.loc[i] > 0) and i != PI:
                remove += ['%s'%i]
                match_cur = match.loc[match.index == i]
                for j in range(match_cur.shape[0]):
                    count[match_cur.iloc[j]] -= 1
        s = [i for i in gencov.index if (i not in remove) or (i == PI)]
            
        return(gencov.loc[s,s])
    
    def define_max_r2_set(S,PI):
        s = [PI]
        for x in S.index:
            s+=[x]
            if x == S.R2.idxmax(axis=0):
                break
        return(s)
    
    def greedy_algorithm(gencov, PI):
        S = pd.DataFrame(np.array([]), index = [], columns = ['R2'])
        S.index.name = 'TID'
        
        for iter in range(N):
            l = [nt for nt in gencov.columns if (nt not in S.index) and (nt not in PI)]
            if len(l) < 1:
                break
            val = pd.Series([0] * len(l), index = l)
            cond = pd.Series([0] * len(l), index = l)
            for temp in l:
                O = np.append(S.index,temp)
                R2, cond_O = gen_R2(gencov, PI, O)
                val.loc[temp] = R2
                cond.loc[temp] = cond_O
            idx = val.idxmax()
            R2 = val.loc[idx]
            condition_number = cond.loc[idx] 
            S.loc[idx, 'R2'] = R2
            S.loc[idx, 'CN'] = condition_number
        trait_list = np.insert(S.index.to_numpy(), 0, PI)
        return(trait_list, S)

    def test_multiple_shrinkages(shrinkages, gencov, PI, N):
        summary_R2 = pd.DataFrame(index = range(N), columns = shrinkages) 
        summary_TID = pd.DataFrame(index = range(N), columns = shrinkages) 
        summary_COND = pd.DataFrame(index = range(N), columns = shrinkages) 

        for shrinkage in shrinkages:
            reg_gencov = pd.DataFrame(shrunk_covariance(gencov, shrinkage = shrinkage), index = gencov.index, columns = gencov.columns)
            trait_list, S = greedy_algorithm(reg_gencov,  PI)
            S['TID'] = S.index
            S.index = range(N)
            summary_R2[shrinkage] = S.R2
            summary_TID[shrinkage] = S.TID
            summary_COND[shrinkage] = S.CN
        return(summary_R2, summary_TID, summary_COND)
    
    def select_best_shrinkage(summary_R2, summary_TID, summary_COND):
        def test_non_decreasing(l):
            l = np.array(l)
            if len(l) < 2:
                return(True)
            
            ind = np.array(range(len(l)-1))
            if np.all(l[ind] < l[ind+1]):
                return(True)
            else:
                return(False)

        def test_R2_sanity(l):
            l = np.array(l)
            if np.all((l >= 0) & (l <= 1)):
                return(True)
            else:
                return(False)
        
        test_ND = summary_R2.apply(test_non_decreasing, axis=0, raw=True)
        test_sanity = summary_R2.apply(test_R2_sanity, axis=0, raw=True)    

        if not any(test_ND & test_sanity):
            raise ValueError('Matrix Regularization failed\nCause:\n\t1.R2 is not ranged between (0,1).\n\tSelected R2 is not non-decreasing vector.\n\n We suggest the followings:\n\t1. check the sanity of covariance matrix\n\tor 2. Use option shrinkages with values higher than 0.9')

        best_shrinkage = summary_R2.columns[test_ND & test_sanity][0]

        best_S = pd.DataFrame(
            {'R2':summary_R2.loc[:,best_shrinkage].to_numpy(float),
             'COND':summary_COND.loc[:,best_shrinkage].to_numpy(float)}
            , index = summary_TID.loc[:,best_shrinkage].to_numpy(str))
        return(best_S, best_shrinkage)

   
    if not (isinstance(shrinkages, (np.ndarray,list) ) & (np.array(shrinkages).dtype == 'float')):
        raise ValueError('Shrinkages should be a float vector')
    else:
        shrinkages = np.array(shrinkages)
        
    if gencov.shape[0] == 1 and gencov.index[0] == PI:
        warnings.warn("The GeneticCovariance matrix contains no non-target traits to be tested", category=Warning)
        return(np.array([PI]), None, None, None)
    
    if (gencov.shape[0] <= 1) or (PI not in gencov.index) or (PI not in gencov.columns):
        raise ValueError('Genetic Covariance Matrix should have the size > N and must contain target trait:%s'%PI)

    if gencov.shape[0] < N:
        raise ValueError('Genetic Covariance matrix is smaller than N'%gencov.shape[0])
    
    if N <= 0:
        return(np.array([PI]), None, None, None)
        
    summary_R2, summary_TID, summary_COND = test_multiple_shrinkages(shrinkages, gencov, PI, N)
    best_S, best_shrinkage = select_best_shrinkage(summary_R2, summary_TID, summary_COND)
    selected_traits = np.insert(best_S.index,0,PI)
    summary_data = {'R2':summary_R2, 'TID':summary_TID, 'COND':summary_COND, 'shrinkages':shrinkages}
    
    return(selected_traits, best_S, best_shrinkage, summary_data)

def ATSA(Gennetic_covariance, PI, N = 5, shrinkage = np.round(np.linspace(0.01, 0.99, num = 12),2)):
    """
    Automatic Trait Selection Algorithm(ATSA)
    The algorithm employs genetic covariance and iterative selection to identify the most relevant non-target traits that significantly influence the genetic variance of the target disease. By iteratively selecting traits, the algorithm aims to build a trait subset that maximizes the impact of the target trait in prediction models. The algorithm is designed to handle datasets with a large number of traits efficiently.
    
    Parameters: 
    - gencov: Genetic Covariance Matrix. 
    - PI: Name of Phenotype of Interest. 
    - N[Default = 5]: Number of traits to be selected(#Q-1).
    - shrinkages[Default = 12 points from 0.01 to 0.99]: parameter space of srinkage factors Default is 12 points (0.01, 0.99)
    
    Returns: 
    
    - selected_traits: The selected traits (TIDs in set Q)
    - best_S: A dataframe object which gives R and COND for best_shrinkage
    - best_shrinkage: The selected shrinkage factor
    - summary_data: A dictionary obejcts, containing the following variables {R2-R2[Per Shrinkage Factor], TID-selected_trait_ID[Per Shrinkage Factor], COND-condition_number[Per Shrinkage Factor] of matrix, shrinkages-shrinkage_factors[Per Shrinkage Factor]}
    """

    def check_n_by_n_matrix(df):

        if not isinstance(df, pd.DataFrame):
            raise ValueError("Input is not a Pandas DataFrame.")

        if df.shape[0] != df.shape[1] or not df.index.equals(df.columns):
            raise ValueError("Input is not an N-by-X matrix.")
            
    def check_PI(variable, df):

        if variable not in df.index and variable not in df.columns:
            raise ValueError("Variable is not found in the DataFrame's index or columns.")
    
    def check_numeric_value(value):
        if not isinstance(value, (int, float)):
            raise ValueError("Input must be a numeric value.")
            
    def validate_numeric_vector(vector):
        if not isinstance(vector, np.ndarray):
            raise ValueError("Input must be a numpy array.")

        if len(vector) >= 100:
            raise ValueError("Vector size must be less than 100.")

        for value in vector:
            if not (0 <= value <= 1):
                raise ValueError("Values in the vector must be between 0 and 1.")

        return vector
      
    
    check_n_by_n_matrix(Gennetic_covariance)
    check_PI(PI, Gennetic_covariance)
    check_numeric_value(N)
    validate_numeric_vector(shrinkage)
    
    selected_traits, best_S, best_shrinkage, summary_data = contribution_based_greedy_algorithm(Gennetic_covariance,PI,N,shrinkage)
    
    return(selected_traits, best_S, best_shrinkage, summary_data)

def run_ATSA(args):
    args.shrinkage = np.round(np.linspace(0.01, 0.99, num = 12),2)
    Gennetic_covariance=args.GENCOV
    print(Gennetic_covariance)
    PI=args.pi
    N = args.Q - 1
    
    # If N+1 is not less than the first dimension of Gennetic_covariance
    if N+1 > Gennetic_covariance.shape[0]:
        # Update N to be the maximum possible value
        N = Gennetic_covariance.shape[0] - 1
        args.log.log(f"Adjusted N to be {N} due to the size of Gennetic_covariance, which is {Gennetic_covariance.shape[0]} x {Gennetic_covariance.shape[1]}.")
    
    shrinkage = args.shrinkage
    selected_traits, best_S, best_shrinkage, summary_data = ATSA(Gennetic_covariance, PI, N, shrinkage)
    
    return selected_traits, best_S, best_shrinkage, summary_data