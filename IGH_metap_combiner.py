#!/usr/bin/python3
#########################################################################
# Author : Porter B. Howland                                            #
# Date   : April 23, 2025                                               #
# Email  : pbh2@illinois.edu                                            #
#                                                                       #
# Description :                                                         #
#                                                                       #
# Parse the raw output file of Seurat's FindConservedMarkers() adding a #
# patient label to each row.  A new dataframe is then made with columns #
#                                                                       #
# 1_avg_log2FC ... 4_avg_log2FC  1_p_val ... 4_p_val                    #
#                                                                       #
# in addition to the additional columns                                 #
#                                                                       #
# mean_log2FC  {method}_p_val  adj_{method}_p_val                       #
#                                                                       #
# which are computed.  Gene lists are saved to files according to the   #
# thresholds set at the beginning of the file.                          #
#########################################################################

import pandas as pd
import numpy as np
from scipy.stats import combine_pvalues
from scipy.stats import false_discovery_control
import re


# Thresholds and method for p-value combination
pval_thresh = 1e-5
log2FC_thresh = 0.2
method = "stouffer"
label = method + "_padj_5_log2FC_02"


def slicer(header_lines, header_matches):
    """
    Returns a function which takes patient integer and returns a slice
    of the indices (0-based, left/right inclusive) for use with the
    .loc[] method.
    """

    patient_splits = {}
    header_lines = [header_lines[n]-n for n in range(len(header_lines))]

    for i, m in enumerate(header_matches[:-1]):
        patient_number = int(m.split('_')[0])
        patient_splits[patient_number] = (header_lines[i], header_lines[i+1] - 1)

    return lambda p: slice(*patient_splits[p])


# find the line number locations of each experimental group and parse
# group id

input_file = 'path_to_de_file_from_Seurat.txt'
header_lines = []
header_matches = []
with open(input_file, 'r') as f:
    for n, line in enumerate(f.readlines()):
        match = re.match(r'.*p_val', line)
        if match:
            header_lines.append(n)
            header_matches.append(match.group())

# create the slicer (by patient)
patient = slicer(header_lines, header_matches)

# read in the dataframe skipping header sections
df = pd.read_csv(input_file,
                 delimiter=r'\s+',
                 skiprows=header_lines,
                 header=None,
                 names=["gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"])

# add a patient column and drop the combined analysis rows
for i in range(1, 5):
    df.loc[patient(i), "patient"] = i
df = df.dropna()
df['patient'] = df['patient'].astype(int)


# combine p-values using the Fisher method and adjust using
# Benjamini-Hochberg. Raw pvals are combined then adjusted; use mean
# log2FC instead of 'mostly' up / down conventions

combined_pvals = pd.DataFrame()

for k in range(1, 5):
    combined_pvals[str(k) + '_avg_log2FC'] = \
        df.query('patient == @k')[['gene', 'avg_log2FC']].set_index('gene')

for k in range(1, 5):
    combined_pvals[str(k) + '_p_val'] = \
        df.query('patient == @k')[['gene', 'p_val']].set_index('gene')

combined_pvals['mean_log2FC'] = \
    combined_pvals.loc[:,'1_avg_log2FC':'4_avg_log2FC'].apply(np.mean, axis=1)

combined_pvals[method + '_p_val'] = \
    combined_pvals.loc[:, '1_p_val':'4_p_val']\
                  .apply(
                      lambda x: combine_pvalues(x, method=method)[1],
                      axis=1
                  )

combined_pvals['adj_' + method + '_p_val'] = \
    false_discovery_control(combined_pvals[method + '_p_val'])


with open(f"genes_igh_combined_adjusted_up_{label}.txt", 'w') as f:
    for gene in combined_pvals.query(f"adj_{method}_p_val < @pval_thresh and mean_log2FC > @log2FC_thresh").index:
        print(gene, file=f)

with open(f"genes_igh_combined_adjusted_down_{label}.txt", 'w') as f:
    for gene in combined_pvals.query(f"adj_{method}_p_val < @pval_thresh and mean_log2FC < -@log2FC_thresh").index:
        print(gene, file=f)

combined_pvals.to_csv(f"{method}_combined_pvals.txt", index=True)
