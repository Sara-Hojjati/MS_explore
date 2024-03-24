#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

def add_clinical_info_to_data(data, clinical_data):
    
    clinical_data = clinical_data.loc[:,['sex', 'age', 'NEDA_EDA_2years', 'nARMSS', 'treatment_duration_index']]
    clinical_data['sex'][clinical_data['sex'] == "f"] = 1
    clinical_data['sex'][clinical_data['sex'] == "m"] = 0

    clinical_data = clinical_data.loc[data.columns, :]

    data = pd.concat([data, clinical_data.T])
    
    data.loc['ms_control'] = 0
    data.loc['ms_control'][data.columns.str.contains('MS')] = 1

    return data

# Read data
clinical_data = pd.read_csv("data/clinical_data.csv", sep='\t', index_col = 'sample_ID')
csf_data = pd.read_csv("data/CSF_batch_corrected.tsv", sep='\t')
plasma_data = pd.read_csv("data/Plasma_batch_corrected.tsv", sep='\t')

# Add clinical information to protein data sets
csf_data = add_clinical_info_to_data(csf_data, clinical_data)
plasma_data = add_clinical_info_to_data(plasma_data, clinical_data)

# Save data
csf_data.to_csv('data/CSF_corrected_data.tsv', sep='\t')
plasma_data.to_csv('data/Plasma_corrected_data.tsv', sep='\t')
