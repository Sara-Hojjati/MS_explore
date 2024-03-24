#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

def format_data(data, below_lod, clinical_data, dataset):
    data = data.pivot(index='OlinkID', columns='SampleID', values='NPX')
    data = data.loc[:,  clinical_data.index[clinical_data.index.isin(data.columns)]]
    
    data = data.loc[~below_lod.loc[data.index, f'{dataset}_belowLOD']]
    data.index = below_lod.loc[data.index, 'Symbol']
    data.index.name = None
    data = data.groupby(data.index).mean()
    return data

def control_correct_per_assay(data):
    dis_hc_mean = data.loc[:, data.columns.str.contains("Dis_HC_")].mean(1)
    dis_hc_std = data.loc[:, data.columns.str.contains("Dis_HC_")].std(1)
    rep_hc_mean = data.loc[:, data.columns.str.contains("Rep_HC_")].mean(1)
    rep_hc_std = data.loc[:, data.columns.str.contains("Rep_HC_")].std(1)
    
    dis = data.loc[:, data.columns.str.contains("Dis")]
    rep = data.loc[:, data.columns.str.contains("Rep")]
    
    dis_adjust = dis.subtract(dis_hc_mean, axis=0).divide(dis_hc_std, axis=0)
    rep_adjust = rep.subtract(rep_hc_mean, axis=0).divide(rep_hc_std, axis=0)
    
    data_adjust = pd.merge(dis_adjust, rep_adjust, left_index=True, right_index=True)
    return data_adjust

# Read input data
below_lod = pd.read_csv("data/proteins_below_LOD.csv", sep='\t', index_col='OlinkID')
clinical_data = pd.read_csv("data/clinical_data.csv", sep='\t', index_col = 'sample_ID')
csf = pd.read_csv('data/CSF_proteindata.csv')
plasma = pd.read_csv('data/plasma_proteindata.csv')

# Format data from long to wide format
csf_mat = format_data(csf, below_lod, clinical_data, dataset='CSF')
plasma_mat = format_data(plasma, below_lod, clinical_data, dataset='Plasma')

# Save NPX data
csf_mat.to_csv("data/CSF_data.tsv", sep='\t')
plasma_mat.to_csv("data/Plasma_data.tsv", sep='\t')

# Control correction
csf_data_cadj = control_correct_per_assay(csf_mat)
plasma_data_cadj = control_correct_per_assay(plasma_mat)

# Save corrected data
csf_data_cadj.to_csv("data/CSF_control_corrected.tsv", sep='\t')
plasma_data_cadj.to_csv("data/Plasma_control_corrected.tsv", sep='\t')
