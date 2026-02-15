import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from functools import reduce
import os
import glob
from matplotlib.colors import LinearSegmentedColormap


# Set matplotlib defaults to use Arial
mpl.rcParams['font.family'] = 'sans-serif'          # use sans-serif family
mpl.rcParams['font.sans-serif'] = ['Arial']         # prefer Arial
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Arial'
mpl.rcParams['mathtext.it'] = 'Arial:italic'
mpl.rcParams['mathtext.bf'] = 'Arial:bold'

# You can also set default sizes if desired:
mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['legend.fontsize'] = 11
mpl.rcParams['figure.titlesize'] = 14

# Tell seaborn to use matplotlib's rc (this is usually automatic)
sns.set_theme()  # respects mpl.rcParams

# Load in the combined patients csv generated from the
# 'metap_igh_melanocyte.py' script, which combines DE
# data from all patients (output from Seurat pipeline).
data = pd.read_csv("path_to_combined_patient_data.tsv", sep="\t", index_col=0)


def reader(filename):
    """Returns a set containing the gene names in the GO category"""
    return set(np.loadtxt(filename, delimiter='\t', skiprows=1, usecols=0, dtype=str))


# Set to directory where GO genes are located
# Each file in this directory should contain gene
# names associated with the category. All GO categories
# which we want represented should be included in
# 'go_directory'. Each file therein should begin with
# 'GO' preferrably followed by some kind of unique identifier
# (e.g., the GO accession number).
go_directory = 'Neurodegenerative_GO_categories_IGH'
os.chdir(go_directory)
genes = list(reduce(lambda x, y: x | y, map(reader, glob.glob("GO*"))))
os.chdir("..")

patients = [1,2,3,4]
avg_log2FC_cols = [f"{i}_avg_log2FC" for i in patients]
pvalue_cols = [f"{i}_p_val" for i in patients]

# Store sliced dataframe
df = data.loc[genes][avg_log2FC_cols]
df['mean_log2FC'] = data.loc[genes][avg_log2FC_cols].apply(lambda x: np.mean(x[avg_log2FC_cols]), axis=1)

# Build colormap
cmap = LinearSegmentedColormap.from_list("pure_bwr", ["#0000FF", "#FFFFFF", "#FF0000"])
maxabs = 1.7
vmin, vmax = -maxabs, maxabs

sns.heatmap(df.sort_values(by=['mean_log2FC'], ascending=False)[avg_log2FC_cols],
            cmap=cmap, center=0, vmin=vmin, vmax=vmax, annot=True, fmt=".2f",
            xticklabels=[f"Patient {i}" for i in range(1, 5)],
            cbar_kws={"label": r"avg. $\log_2$-FC"})

# plt.savefig("IGH_heatmap_neurodegeneration.png", dpi=600, bbox_inches='tight')
plt.show()
