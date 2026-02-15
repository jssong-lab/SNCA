
############################################################################
## Author  : Porter B. Howland                                            ##
## Email   : pbh2@illinois.edu                                            ##
## Date    :                                                              ##
## Updated :                                                              ##
##                                                                        ##
## Description :                                                          ##
##                                                                        ##
## Creates scatter plots of labeled GO categories (enrichment score) vs.  ##
## -log_10(adjusted p-value) for labeled GO categories, and plots a       ##
## background of other GO terms shown in grey.                            ##
##                                                                        ##
############################################################################


import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


# Set default plot fonts and fontsizes
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['font.size'] = 14
mpl.rcParams['text.latex.preamble'] = r'\usepackage{helvet}'

# Data extracted from DAVID GO analysis was saved in an excel file
# from which specific GO terms were selected for vieweing (line below)
# and background (next read below).
#
# The format of the csv file should be:
#
# GO category name, enrichment score, adjusted p-value
#
df = pd.read_csv("path_to/go_scores.csv")

# Remove small filled point from marker list and gray from color list
# since these will be used for the background points
markers = Line2D.filled_markers[1:]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors.remove("#7f7f7f")
num_markers = len(markers)
num_colors = len(colors)

# Load in points (adjusted p-values, enrichment scores) for GO
# categories of the background
mat = np.loadtxt("path_to/background_case_vs_control.txt")

mat[:, 0] = -np.log10(mat[:, 0])

# specify figure and inset axes
fig, ax = plt.subplots(1, 1, figsize=(8, 6.5))

# plot annotated points on both main plot and inset
for i, row in df.iterrows():
    ax.scatter(-np.log10(row.iloc[3]),
               row.iloc[2],
               label=row.iloc[1],
               marker=markers[i % num_markers],
               color=colors[i % num_colors])

# Scatter plot background points for both main and inset
# BUG : put this before the annotations results in the first
# annotation point having the same plot style as the background
# otherwise labels get miscolored
ax.scatter(mat[:, 0], mat[:, 1], alpha=0.1, marker='.', color='grey')

ax.legend(list(df["GO Tag"]), bbox_to_anchor=(1.05, 1.0))
ax.set_title("Gene Ontologies: case vs. control")
ax.set_xlabel(r'$- \log_{10}(\text{adj. }p)$')
ax.set_ylabel("Enrichment score")
fig.savefig('go-scatter-case-vs-control.png', bbox_inches='tight')
fig.savefig('go-scatter-case-vs-control.pdf', bbox_inches='tight')
