#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import scipy.stats as sts
import pandas as pd
import numpy as np

def fisher_test(geneset1, geneset2, background):
    geneset1 = [gene for gene in geneset1 if gene in background]
    geneset2 = [gene for gene in geneset2 if gene in background]
    overlap = np.intersect1d(geneset1, geneset2)
    
    table = [[len(overlap), 
              len(geneset1)-len(overlap)],
             [len(geneset2)-len(overlap), 
              len(background)-len(geneset2)-len(geneset1)+len(overlap)]]
    oddsratio, pvalue = sts.fisher_exact(table, alternative='two-sided')
    return oddsratio, pvalue, table

# Background proteins
background = pd.read_csv('data/CSF_data.tsv', sep='\t', index_col='Unnamed: 0').index.values

# Differentially expressed proteins in CSF in discovery cohort
csf_dep_disc = pd.read_csv("results/Discovery_DEP.txt", header=None).values.flatten()

# Known MS biomarkers
ms_markers = pd.read_csv("data/known_ms_biomarkers.txt", header=None).values

# Run fishers exact test
oddsratio, pvalue, table = fisher_test(csf_dep_disc, ms_markers, background)
print("OR =", oddsratio, "p =", pvalue)
