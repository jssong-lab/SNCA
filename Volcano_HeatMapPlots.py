import math
import re
import os,sys
import numpy as np
import scipy as spy
import scipy.stats as stats
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import pandas as pd
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from collections import Counter
import itertools

import seaborn as sns


import warnings


warnings.filterwarnings("ignore")


seed=0;
np.random.seed(seed);

plt.style.use('seaborn-v0_8-whitegrid');
plt.ion();



ColorScheme=[
(0.368417, 0.506779, 0.709798), 
(0.880722, 0.611041, 0.142051), 
(0.560181, 0.691569, 0.194885), 
(0.922526, 0.385626, 0.209179), 
(0.528488, 0.470624, 0.701351), 
(0.772079, 0.431554, 0.102387), 
(0.363898, 0.618501, 0.782349), 
(1, 0.75, 0), 
(0.647624, 0.37816, 0.614037), 
(0.571589, 0.586483, 0.), 
(0.915, 0.3325, 0.2125), 
(0.40082222609352647, 0.5220066643438841, 0.85), 
(0.9728288904374106, 0.621644452187053, 0.07336199581899142), 
(0.736782672705901, 0.358, 0.5030266573755369), 
(0.28026441037696703, 0.715, 0.4292089322474965),
];


pd.set_option('display.max_rows', 100)


Sample_List_siTFE=['4A','4B','C1','C2-1','C2-2','C2-3','5A','5B','M','3-1','3-2','3-3','B1-1','B1-2','B1-3',];
SampleToSapleTag_dict_siTFE={'4A':'siControl1_1','4B':'siControl1_2','C1':'siControl1_3','C2-1':'siControl2_1','C2-2':'siControl2_2','C2-3':'siControl2_3','5A':'siMITF_1','5B':'siMITF_2','M':'siMITF_3','3-1':'siTFE3_1','3-2':'siTFE3_2','3-3':'siTFE3_3','B1-1':'siTFEB_1','B1-2':'siTFEB_2','B1-3':'siTFEB_3',};

genes_fpkm_table_df_dict_siTFE={Sample:pd.read_table("Exp/"+"genome_hg38"+"/"+Sample+"_gene_abund.tab",sep='\t',header=0) for Sample in Sample_List_siTFE};

for Sample in Sample_List_siTFE: genes_fpkm_table_df_dict_siTFE[Sample]=genes_fpkm_table_df_dict_siTFE[Sample][[Sample+':Gene ID',Sample+':'+ExpressionType]].groupby(Sample+':Gene ID').max()

    
ExpressionTable_DESeq_df_siTFE=pd.concat([genes_fpkm_table_df_dict_siTFE[Sample] for Sample in Sample_List_siTFE],axis=1)
ExpressionTable_DESeq_df_siTFE=ExpressionTable_DESeq_df_siTFE.rename(columns={Sample+':'+ExpressionType:SampleToSapleTag_dict_siTFE[Sample] for Sample in Sample_List_siTFE})


transcript_fpkm_table_df_dict_MITF={Sample:pd.read_table(Dir_siTFE+"/../Exp/"+"genome_hg38"+"/"+Sample+"_t_data.ctab",sep='\t',header=0) for Sample in Sample_List_siTFE};

for Sample in Sample_List_siTFE: transcript_fpkm_table_df_dict_MITF[Sample]=transcript_fpkm_table_df_dict_MITF[Sample][[Sample+':t_name',Sample+':'+'FPKM']].groupby(Sample+':t_name').max()

TranscriptTable_DESeq_df_siTFE=pd.concat([transcript_fpkm_table_df_dict_MITF[Sample] for Sample in Sample_List_siTFE],axis=1)
TranscriptTable_DESeq_df_siTFE=TranscriptTable_DESeq_df_siTFE.rename(columns={Sample+':'+'FPKM':SampleToSapleTag_dict_siTFE[Sample] for Sample in Sample_List_siTFE})



MIT_Table_df=pd.concat([
ExpressionTable_DESeq_df_siTFE.loc[['MITF','TFE3','TFEB']][['siControl1_1','siControl1_2','siControl1_3']].mean(axis=1).rename('Control1'),
ExpressionTable_DESeq_df_siTFE.loc[['MITF','TFE3','TFEB']][['siControl2_1','siControl2_2','siControl2_3']].mean(axis=1).rename('Control2'),
ExpressionTable_DESeq_df_siTFE.loc[['MITF','TFE3','TFEB']][['siMITF_1','siMITF_2','siMITF_3']].mean(axis=1).rename('siMITF'),
ExpressionTable_DESeq_df_siTFE.loc[['MITF','TFE3','TFEB']][['siTFE3_1','siTFE3_2','siTFE3_3']].mean(axis=1).rename('siTFE3'),
ExpressionTable_DESeq_df_siTFE.loc[['MITF','TFE3','TFEB']][['siTFEB_1','siTFEB_2','siTFEB_3']].mean(axis=1).rename('siTFEB'),
ExpressionTable_DESeq_df_siMITF.loc[['MITF','TFE3','TFEB']][['siControl_0','siControl_1']].mean(axis=1).rename('Control Akinori'),
ExpressionTable_DESeq_df_siMITF.loc[['MITF','TFE3','TFEB']][['siMITF_0','siMITF_1']].mean(axis=1).rename('siMITF Akinori'),],axis=1)
 
MIT_Table_df

MIT_Family_List=['MITF','TFE3','TFEB']

fig=plt.figure(figsize=(8, 6));
ax=fig.add_axes([0.2, 0.2, 0.8, 0.8]);
LegendHandles=[];

W=0.1;
Ind=np.array([i*W*3.0 for i in range(len(MIT_Family_List))]);
Width_Data=[W]*len(MIT_Family_List);


n=-0.5;
Bar_Data=np.array([MIT_Table_df.loc['MITF']['Control1'],MIT_Table_df.loc['TFE3']['Control2'],MIT_Table_df.loc['TFEB']['Control2']]);
Color_Data=[ColorScheme[0]]*len(MIT_Family_List);
ax.barh(Ind+W*n, Bar_Data, Width_Data, color=Color_Data,align='center');
LegendHandles.append(mpatches.Patch(color=Color_Data[0], label='Control'));

n=+0.5;
Bar_Data=np.array([MIT_Table_df.loc[g]['si'+g] for g in MIT_Family_List]);
Color_Data=[ColorScheme[3]]*len(MIT_Family_List);
ax.barh(Ind+W*n, Bar_Data, Width_Data, color=Color_Data,align='center');
LegendHandles.append(mpatches.Patch(color=Color_Data[0], label='Knock down'));


plt.gca().invert_yaxis()

ax.set_ylabel('Gene',fontsize=22);
ax.set_xlabel('Expression',fontsize=22);


ax.set_yticks(np.array(Ind));
ax.set_yticklabels(MIT_Family_List,fontsize=22,rotation=0)


for tick in ax.xaxis.get_major_ticks(): tick.label1.set_fontsize(18) ;
ax.legend(handles=LegendHandles,fontsize=20,loc=[0.45,0.25]);

Axis=ax;
Plot=Axis.tick_params(axis='both', labelsize=24, width=2, length=4);
Plot=Axis.tick_params(axis='x', pad=6);
Plot=Axis.tick_params(axis='y', pad=8);
Plot=Axis.xaxis.set_major_locator(plt.MaxNLocator(5));
Plot=Axis.yaxis.set_major_locator(plt.MaxNLocator(6));
Plot=Axis.set_frame_on(True);
Axis.spines['top'   ].set_linewidth(1.75);
Axis.spines['bottom'].set_linewidth(1.75);
Axis.spines['left'  ].set_linewidth(1.75);
Axis.spines['right' ].set_linewidth(1.75);


Axis.spines['right'].set_visible(False)
Axis.spines['top'].set_visible(False)


plt.savefig("MIT_Genes_KD.pdf",format='pdf',bbox_inches='tight');



MITF_peaks_Dir='MITF_peaks'
GenesSet_Dir='Genes'

BoundByMITF_list=ReadFile(MITF_peaks_Dir+'/'+'genes_2.txt'); #genes bound in at least 2 replicates


[ActivatedByMITF_list,RepressedByMITF_list]=[{},{}]
    
for Type in ['Autophagy','Lysosome']: ActivatedByMITF_list[Type]=ReadFile(GenesSet_Dir+'/'+Type+'_MITF_up'+'.txt');
for Type in ['Autophagy','Lysosome']: RepressedByMITF_list[Type]=ReadFile(GenesSet_Dir+'/'+Type+'_MITF_down'+'.txt');
ActivatedByMITF_list['All']=ReadFile(GenesSet_Dir+'/'+Type+'_MITF_up'+'.txt')    
RepressedByMITF_list['All']=ReadFile(GenesSet_Dir+'/'+Type+'_MITF_down'+'.txt');

    

for Type in ['Autophagy','Lysosome','All']: ActivatedByMITF_list[Type]=['C9orf72' if item=='C9ORF72' else item for item in ActivatedByMITF_list[Type]]    #fix for upper case 
for Type in ['Autophagy','Lysosome', 'All']: RepressedByMITF_list[Type]=['C9orf72' if item=='C9ORF72' else item for item in RepressedByMITF_list[Type]]    #fix for upper case 



DiTag='siControl1_siMITF';
#DiTag='siControl2_siTFE3';

Type='Autophagy'

Expression_df=gene_exp_diff_df_DESeq_dict[DiTag].dropna()

fig=plt.figure(figsize=(5, 5));
Axis=fig.add_axes([0.2, 0.2, 0.8, 0.8]);

Plot_data_x=Expression_df['log2FoldChange'].values
Plot_data_y=-np.log10(Expression_df['padj'].values)
Plot_data_padj=Expression_df['padj'].values
Plot_data_name=Expression_df.index


Plot_data_size=5

Plot_data_color=[(0.25,0.25,0.25,0.05) if Plot_data_padj[i]<0.05 else (0.0,0.0,0.0,1.0) for i in range(len(Plot_data_padj))]

Axis.scatter(Plot_data_x, Plot_data_y,s=Plot_data_size,c=Plot_data_color);

Axis.set_xlabel(r'$\rm {log_{2}} $(fold-change)',size=20),Axis.set_ylabel(r'$\rm {-log_{10}}$ (adj. p-value)',size=20,rotation=90);
Axis.set_title(Type+' Genes',size=20,pad=10);

Plot=Axis.tick_params(axis='both', labelsize=20, width=2, length=4);
Plot=Axis.tick_params(axis='x', pad=6);
Plot=Axis.tick_params(axis='y', pad=8);
Plot=Axis.xaxis.set_major_locator(plt.MaxNLocator(5));
if(DiTag=='siControl2_siTFE3'): Plot=Axis.xaxis.set_major_locator(plt.MaxNLocator(6));
Plot=Axis.yaxis.set_major_locator(plt.MaxNLocator(5));
Plot=Axis.set_frame_on(True);
Axis.spines['top'   ].set_linewidth(1.75);
Axis.spines['bottom'].set_linewidth(1.75);
Axis.spines['left'  ].set_linewidth(1.75);
Axis.spines['right' ].set_linewidth(1.75);


[Range_x_min,Range_x_max]=[-10,10];
[Range_y_min,Range_y_max]=[-5,165];

if(DiTag=='siControl1_siMITF'):
    [Range_x_min,Range_x_max]=[-10,10];
    [Range_y_min,Range_y_max]=[-2,80];
    
if(DiTag=='siControl2_siTFE3'):
    [Range_x_min,Range_x_max]=[-3,3];
    [Range_y_min,Range_y_max]=[-0.5,5];

    
Line_data_x=[0.0,0.0];
Line_data_y=[Range_y_min,Range_y_max];
Axis.plot(Line_data_x, Line_data_y, color=(0.5, 0.5, 0.5),linewidth=1,linestyle='-');
Line_data_x=[1.0,1.0];
Line_data_y=[Range_y_min,Range_y_max];
Axis.plot(Line_data_x, Line_data_y, color=(0.5, 0.5, 0.5),linewidth=1,linestyle='dashed');
Line_data_x=[-1.0,-1.0];
Line_data_y=[Range_y_min,Range_y_max];
Axis.plot(Line_data_x, Line_data_y, color=(0.5, 0.5, 0.5),linewidth=1,linestyle='dashed');


Line_data_x=[Range_x_min,Range_x_max];
Line_data_y=[0.0,0.0];
Axis.plot(Line_data_x, Line_data_y, color=(0.5, 0.5, 0.5),linewidth=1,linestyle='-');

if(DiTag=='siControl2_siTFE3' or True):
    Line_data_x=[Range_x_min,Range_x_max];
    Line_data_y=[-np.log10(0.05),-np.log10(0.05)];
    Axis.plot(Line_data_x, Line_data_y, color=(0.5, 0.5, 0.5),linewidth=1,linestyle='dashed');


Axis.set_xlim([Range_x_min,Range_x_max]);
Axis.set_ylim([Range_y_min,Range_y_max]);


for i in range(len(Plot_data_name)):
    if(Plot_data_name[i] in (ActivatedByMITF_list[Type]+RepressedByMITF_list[Type])):
        Axis.scatter(Plot_data_x[i], Plot_data_y[i],s=Plot_data_size*2,c=ColorScheme[3]);



fig.savefig("Volcano_"+DiTag+'_'+Type+".pdf", bbox_inches='tight');



DiTag='siControl1_siMITF';
#DiTag='siControl2_siTFE3';

Type='Lysosome'

Expression_df=gene_exp_diff_df_DESeq_dict[DiTag].dropna()

fig=plt.figure(figsize=(5, 5));
Axis=fig.add_axes([0.2, 0.2, 0.8, 0.8]);

Plot_data_x=Expression_df['log2FoldChange'].values
Plot_data_y=-np.log10(Expression_df['padj'].values)
Plot_data_padj=Expression_df['padj'].values
Plot_data_name=Expression_df.index


Plot_data_size=5

Plot_data_color=[(0.25,0.25,0.25,0.05) if Plot_data_padj[i]<0.05 else (0.0,0.0,0.0,1.0) for i in range(len(Plot_data_padj))]

Axis.scatter(Plot_data_x, Plot_data_y,s=Plot_data_size,c=Plot_data_color);

Axis.set_xlabel(r'$\rm {log_{2}} $(fold-change)',size=20),Axis.set_ylabel(r'$\rm {-log_{10}}$ (adj. p-value)',size=20,rotation=90);
Axis.set_title(Type+' Genes',size=20,pad=10);

Plot=Axis.tick_params(axis='both', labelsize=20, width=2, length=4);
Plot=Axis.tick_params(axis='x', pad=6);
Plot=Axis.tick_params(axis='y', pad=8);
Plot=Axis.xaxis.set_major_locator(plt.MaxNLocator(5));
if(DiTag=='siControl2_siTFE3'): Plot=Axis.xaxis.set_major_locator(plt.MaxNLocator(6));
Plot=Axis.yaxis.set_major_locator(plt.MaxNLocator(5));
Plot=Axis.set_frame_on(True);
Axis.spines['top'   ].set_linewidth(1.75);
Axis.spines['bottom'].set_linewidth(1.75);
Axis.spines['left'  ].set_linewidth(1.75);
Axis.spines['right' ].set_linewidth(1.75);


[Range_x_min,Range_x_max]=[-10,10];
[Range_y_min,Range_y_max]=[-5,165];

if(DiTag=='siControl1_siMITF'):
    [Range_x_min,Range_x_max]=[-10,10];
    [Range_y_min,Range_y_max]=[-2,80];
    
if(DiTag=='siControl2_siTFE3'):
    [Range_x_min,Range_x_max]=[-3,3];
    [Range_y_min,Range_y_max]=[-0.5,5];

    
Line_data_x=[0.0,0.0];
Line_data_y=[Range_y_min,Range_y_max];
Axis.plot(Line_data_x, Line_data_y, color=(0.5, 0.5, 0.5),linewidth=1,linestyle='-');
Line_data_x=[1.0,1.0];
Line_data_y=[Range_y_min,Range_y_max];
Axis.plot(Line_data_x, Line_data_y, color=(0.5, 0.5, 0.5),linewidth=1,linestyle='dashed');
Line_data_x=[-1.0,-1.0];
Line_data_y=[Range_y_min,Range_y_max];
Axis.plot(Line_data_x, Line_data_y, color=(0.5, 0.5, 0.5),linewidth=1,linestyle='dashed');


Line_data_x=[Range_x_min,Range_x_max];
Line_data_y=[0.0,0.0];
Axis.plot(Line_data_x, Line_data_y, color=(0.5, 0.5, 0.5),linewidth=1,linestyle='-');

if(DiTag=='siControl2_siTFE3' or True):
    Line_data_x=[Range_x_min,Range_x_max];
    Line_data_y=[-np.log10(0.05),-np.log10(0.05)];
    Axis.plot(Line_data_x, Line_data_y, color=(0.5, 0.5, 0.5),linewidth=1,linestyle='dashed');


Axis.set_xlim([Range_x_min,Range_x_max]);
Axis.set_ylim([Range_y_min,Range_y_max]);


for i in range(len(Plot_data_name)):
    if(Plot_data_name[i] in (ActivatedByMITF_list[Type]+RepressedByMITF_list[Type])):
        Axis.scatter(Plot_data_x[i], Plot_data_y[i],s=Plot_data_size*2,c=ColorScheme[3]);



fig.savefig("Volcano_"+DiTag+'_'+Type+".pdf", bbox_inches='tight');



Type='Autophagy'

UsedSamples=['siControl1_1','siControl1_2','siControl1_3','siMITF_1','siMITF_2','siMITF_3']
if (DiTag=='siControl2_siTFE3'): UsedSamples=['siControl2_1','siControl2_2','siControl2_3','siTFE3_1','siTFE3_2','siTFE3_3']
if (DiTag=='siControl2_siTFEB'): UsedSamples=['siControl2_1','siControl2_2','siControl2_3','siTFEB_1','siTFEB_2','siTFEB_3']

N_Genes_Cutoff=40;

RegulatedGenes_list=[g for g in ActivatedByMITF_list[Type] if g in BoundByMITF_list]

SelectedGenes=list(gene_exp_diff_df_DESeq_dict[DiTag].loc[RegulatedGenes_list].sort_values('padj')[0:N_Genes_Cutoff].index)


ClusterData_Raw=gene_exp_DESEq_df[UsedSamples].loc[SelectedGenes];

ClusterData=((ClusterData_Raw.T-ClusterData_Raw.T.mean())/ClusterData_Raw.T.std()).T

ax=sns.clustermap(ClusterData,method='average', metric='cosine' ,vmin=-1.0, vmax=1.0, cmap="coolwarm", figsize=(3,8), linewidths=0.0, yticklabels=True, cbar_kws={"ticks":[-1.0,-0.5,0.0,0.5,1]});
ax.ax_heatmap.set_title(Type+" Genes",pad=110)
ax.ax_cbar.set_position((1.0, 0.58, 0.08, 0.2))
ax.ax_heatmap.set_xticklabels(ax.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
ax.ax_heatmap.set_yticklabels(ax.ax_heatmap.get_ymajorticklabels(), fontsize = 8)

ax.ax_cbar.set_xlabel("z-score",size=12,labelpad=-135);
plt.savefig("ClusterMap_"+Type+'_'+DiTag+'_'+str(N_Genes_Cutoff)+".pdf", bbox_inches='tight');


Type='Lysosome'

UsedSamples=['siControl1_1','siControl1_2','siControl1_3','siMITF_1','siMITF_2','siMITF_3']
if (DiTag=='siControl2_siTFE3'): UsedSamples=['siControl2_1','siControl2_2','siControl2_3','siTFE3_1','siTFE3_2','siTFE3_3']
if (DiTag=='siControl2_siTFEB'): UsedSamples=['siControl2_1','siControl2_2','siControl2_3','siTFEB_1','siTFEB_2','siTFEB_3']

N_Genes_Cutoff=40;

RegulatedGenes_list=[g for g in ActivatedByMITF_list[Type] if g in BoundByMITF_list]

SelectedGenes=list(gene_exp_diff_df_DESeq_dict[DiTag].loc[RegulatedGenes_list].sort_values('padj')[0:N_Genes_Cutoff].index)


ClusterData_Raw=gene_exp_DESEq_df[UsedSamples].loc[SelectedGenes];

ClusterData=((ClusterData_Raw.T-ClusterData_Raw.T.mean())/ClusterData_Raw.T.std()).T

ax=sns.clustermap(ClusterData,method='average', metric='cosine' ,vmin=-1.0, vmax=1.0, cmap="coolwarm", figsize=(3,8), linewidths=0.0, yticklabels=True, cbar_kws={"ticks":[-1.0,-0.5,0.0,0.5,1]});
ax.ax_heatmap.set_title(Type+" Genes",pad=110)
ax.ax_cbar.set_position((1.0, 0.58, 0.08, 0.2))
ax.ax_heatmap.set_xticklabels(ax.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
ax.ax_heatmap.set_yticklabels(ax.ax_heatmap.get_ymajorticklabels(), fontsize = 8)

ax.ax_cbar.set_xlabel("z-score",size=12,labelpad=-135);
plt.savefig("ClusterMap_"+Type+'_'+DiTag+'_'+str(N_Genes_Cutoff)+".pdf", bbox_inches='tight');


cg=sns.clustermap(CutAndRun_SignalData_df_all,method='average', metric='euclidean',row_cluster=True,col_cluster=False,vmin=0.0, vmax=200.0, cmap="coolwarm", figsize=(8,8), linewidths=0.0, xticklabels=False, yticklabels=True, cbar_kws={"ticks":[0,100,200]});

pos_x=[0.28+0.23*i for i in range(len(TagList))]
pos_y=0.80
for i in range(len(TagList)): cg.fig.text(x=pos_x[i], y=pos_y, s=TagList[i]);
cg.ax_row_dendrogram.set_visible(False)


BoarderColor=(0.0,0.75,0.0)
BoarderColor=ColorScheme[7]
BoarderColor=(0.90,0.90,0.0)

[Line_data_x,Line_data_y]=[[0+10,0+10],[0.0,40]];
cg.ax_heatmap.plot(Line_data_x, Line_data_y, color=BoarderColor,linewidth=2,linestyle='-');
[Line_data_x,Line_data_y]=[[2*Slop,2*Slop],[0.0,40]];
cg.ax_heatmap.plot(Line_data_x, Line_data_y, color=BoarderColor,linewidth=2,linestyle='-');
[Line_data_x,Line_data_y]=[[4*Slop,4*Slop],[0.0,40]];
cg.ax_heatmap.plot(Line_data_x, Line_data_y, color=BoarderColor,linewidth=2,linestyle='-');
[Line_data_x,Line_data_y]=[[6*Slop-40,6*Slop-40],[0.0,40]];
cg.ax_heatmap.plot(Line_data_x, Line_data_y, color=BoarderColor,linewidth=2,linestyle='-');

plt.savefig("CutAndRun_Tracks"+".pdf", bbox_inches='tight');

