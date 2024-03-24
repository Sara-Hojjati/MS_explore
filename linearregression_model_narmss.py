#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import scipy.stats as sts
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import r2_score
import statsmodels.api as sm
from scipy.signal import argrelextrema
from scipy.stats import t

def read_csf_data(protein_set, data_file):
    csf_data = pd.read_csv(data_file, sep='\t', index_col='Unnamed: 0')
    
    narmss_scores = csf_data.loc['nARMSS', :].dropna()
            
    csf_earlyms = csf_data[narmss_scores.index]
    
    csf_earlyms.loc['sex', :] = csf_earlyms.loc['sex', :].astype(int)
    csf_earlyms.loc['age', :] = csf_earlyms.loc['age', :].astype(int)
    
    # Training data
    csf_dis = csf_earlyms.loc[:, csf_earlyms.columns.str.contains('Dis')]
    narmss_dis = narmss_scores[narmss_scores.index.str.contains('Dis')].astype('float')
    csf_dis = csf_dis.loc[protein_set, narmss_dis.index]
    
    # Test data
    csf_rep = csf_earlyms.loc[:,csf_earlyms.columns.str.contains('Rep')]
    narmss_rep = narmss_scores[narmss_scores.index.str.contains('Rep')].astype('float')
    csf_rep = csf_rep.loc[protein_set, narmss_rep.index]
    
    return csf_dis, narmss_dis.values, csf_rep, narmss_rep.values

def loo_regression(X, y, model=LinearRegression()):
    loo = LeaveOneOut()
    
    y_pred = []
    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]    
        reg = model.fit(X_train, y_train)
        y_pred += [reg.predict(X_test)]
    
    r2 = r2_score(y, y_pred)
    return np.array(y_pred), r2

def forward_selection(X_train, y_train, features, model=LinearRegression()):  
    selected_features = []
    r2_features = []
    
    for i in range(len(features)):
        
        r2_loo_scores = []
        for prot in features:
            
            pred, r2 = loo_regression(X=X_train.loc[selected_features + [prot],:].T.values, 
                                      y=y_train, model=model)
            r2_loo_scores += [r2]
        
        selected_features += [features[np.argmax(r2_loo_scores)]]
        features.pop(np.argmax(r2_loo_scores))
        r2_features += [r2_loo_scores[np.argmax(r2_loo_scores)]]
    
    return selected_features, r2_features

def coefficient_significance(X_train, y_train, features):
    X2 = sm.add_constant(X_train.loc[features,:].T)
    est = sm.OLS(y_train, X2.astype(float))
    est2 = est.fit()
    return est2

def backward_selection_coefficient_significance(X_train, y_train, features):
    significant_features = features.copy()    
    for i in range(len(features)):
        est2 = coefficient_significance(X_train, y_train, features=significant_features)
        if np.max(est2.pvalues[1:]) >= 0.05:
            significant_features = np.delete(np.array(significant_features), np.argmax(est2.pvalues[1:]))
        else:
            break
    return significant_features

def concordance_correlation_coefficient(y_true, y_pred):
    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)

    centered_true = y_true - mean_true
    centered_pred = y_pred - mean_pred

    covariance = np.mean(centered_true * centered_pred)

    true_var = np.mean(centered_true ** 2)
    pred_var = np.mean(centered_pred ** 2)

    numerator = 2 * covariance
    denominator = true_var + pred_var + (mean_true - mean_pred) ** 2

    ccc = numerator / denominator
    
    n = len(y_true)
    df = n - 1
    t_value = ccc * np.sqrt(df / (1 - ccc**2))
    p_value = 2 * (1 - t.cdf(np.abs(t_value), df))
    return ccc, p_value

def calculate_performance(X_train, X_test, y_train, y_test, features):
    est2 = coefficient_significance(X_train, y_train, features=features)

    reg = LinearRegression().fit(X_train.loc[features,:].T.values, y_train)
    
    pred_train = reg.predict(X_train.loc[features,:].T.values)
    
    performance_train = [r2_score(y_train, pred_train), sts.spearmanr(y_train, pred_train), concordance_correlation_coefficient(y_train, pred_train)]
    
    pred_test = reg.predict(X_test.loc[features,:].T.values)
    performance_test = [r2_score(y_test, pred_test), sts.spearmanr(y_test, pred_test), concordance_correlation_coefficient(y_test, pred_test)]
    return est2, pred_train, pred_test, performance_train, performance_test

def select_plasma_features(csf_dis, narmss_dis, plasma_dis, significant_features):
    corr_prot = []
    for prot in significant_features:
        pcc, p = sts.pearsonr(csf_dis.loc[prot, :].astype(float), plasma_dis.loc[prot, :].astype(float))
        if p < 0.05:
            corr_prot += [prot]
            
    plasma_features = backward_selection_coefficient_significance(csf_dis, narmss_dis, features=np.array(corr_prot))
    return plasma_features

'''
Train nARMSS model for CSF samples
'''
# Read file with differentially expressed proteins in discovery cohort
csf_dep_disc = pd.read_csv("results/Discovery_DEP.txt", header=None).values.flatten()

# Read CSF data
csf_dis, narmss_dis, csf_rep, narmss_rep = read_csf_data(protein_set=np.concatenate((csf_dep_disc, ['sex', 'age'])),
                                                         data_file='data/CSF_corrected_data.tsv')

#Feature selection: forward selection and backward selection
selected_features, r2_features = forward_selection(X_train=csf_dis,
                                                   y_train=narmss_dis, 
                                                   features=list(csf_dis.index),
                                                   model=LinearRegression())

features = selected_features[:argrelextrema(np.array(r2_features), np.greater, order=5)[0][0]+1]

significant_features = backward_selection_coefficient_significance(csf_dis, narmss_dis, features=features)

print("CSF model features:", significant_features)

# Evaluate performance
est2, pred_train, pred_test, performance_train, performance_test = calculate_performance(X_train=csf_dis, X_test=csf_rep,
                                                                                         y_train=narmss_dis, y_test=narmss_rep,
                                                                                         features=significant_features)

print("CSF model performance:", 
      "Discovery SCC =", performance_train[1][0], 
      "Replication SCC =", performance_test[1][0],
      "Discovery CCC =", performance_train[2][0],
      "Replication CCC =", performance_test[2][0])

# Test effect of treatment
csf_dis_treat, narmss_dis_treat, csf_rep_treat, narmss_rep_treat = read_csf_data(protein_set=np.concatenate((csf_dep_disc, ['sex', 'age', 'treatment_duration_index'])), 
                                                                                 data_file='data/CSF_corrected_data.tsv')

features_treatment = np.concatenate((significant_features, ['treatment_duration_index']))

results_treat = calculate_performance(X_train=csf_dis_treat, X_test=csf_rep_treat,
                                      y_train=narmss_dis_treat, y_test=narmss_rep_treat,
                                      features=features_treatment)
print("Treatment effect on CSF model:", 
      "Discovery SCC =", results_treat[3][1][0], 
      "Replication SCC =", results_treat[4][1][0],
      "Discovery CCC =", results_treat[3][2][0],
      "Replication CCC =", results_treat[4][2][0])

'''
Train nARMSS model for plasma samples
'''
# Read Plasma data
plasma_dis, narmss_pla_dis, plasma_rep, narmss_pla_rep = read_csf_data(np.concatenate((csf_dep_disc, ['sex', 'age'])),
                                                                       data_file='data/Plasma_corrected_data.tsv')

#Feature selection: correlation CSF vs Plasma and backward selection
plasma_features = select_plasma_features(csf_dis, narmss_dis, plasma_dis, significant_features)

# Evaluate performance
est2, pred_train, pred_test, performance_train, performance_dis = calculate_performance(X_train=csf_dis, X_test=plasma_dis,
                                                                                         y_train=narmss_dis, y_test=narmss_pla_dis,
                                                                                         features=plasma_features)

est2, pred_train, pred_test, performance_train, performance_rep = calculate_performance(X_train=csf_dis, X_test=plasma_rep,
                                                                                         y_train=narmss_dis, y_test=narmss_pla_rep,
                                                                                         features=plasma_features)

print("Plasma model features:", plasma_features)
print("Plasma model performance:", 
      "Discovery SCC =", performance_dis[1][0], 
      "Replication SCC =", performance_rep[1][0],
      "Discovery CCC =", performance_dis[2][0],
      "Replication CCC =", performance_rep[2][0])
