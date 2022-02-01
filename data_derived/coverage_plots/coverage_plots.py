#Importing modules

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd
from pathlib import Path
from sklearn import preprocessing
import os
import seaborn as sns; sns.set()
from matplotlib import pyplot

###################################################
# Using ANhui files as an example, but code applies to all of the strains:

folder='TSV_files/'
#Anhui rep1
for file in Path(folder).glob('*/rep1/*/Anhui*_coverage.tsv'):
    df_anhui=pd.read_csv(file, header=0, sep='\t')
    
    #Adding a new column ("GENE") containing gene segment name only, obtained by deleting characters from #ID column
    df_anhui['Gene'] = df_anhui['#ID'].str[15:]
    df_anhui['Sample'] = os.path.basename(file)
    df_anhui['Read_depth'] = df_anhui[['Plus_reads', 'Minus_reads']].mean(axis=1)
    df_anhui.append(df_anhui)
    #We added this step to clean up the sample names and keep them consistent across all timepoints:
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret1_':'_ferret01_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret2_':'_ferret02_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret3_':'_ferret03_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret4_':'_ferret04_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret5_':'_ferret05_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret6_':'_ferret06_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret7_':'_ferret07_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret8_':'_ferret08_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret9_':'_ferret09_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day1_':'_day01_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day2_':'_day02_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day3_':'_day03_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day4_':'_day04_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day5_':'_day05_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day6_':'_day06_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day7_':'_day07_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day8_':'_day08_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day9_':'_day09_'}, regex=True)
    df_anhui['Sample'] = df_anhui['Sample'].str[:25]
    # drop columns that are not needed (i.e. "Length")
    df_anhui = df_anhui.drop(['#ID', 'Length', 'Ref_GC', 'Covered_bases', 'Plus_reads', 'Minus_reads', 'Read_GC', 'Median_fold'], axis=1)
    #reordering folumns for easier access
    df_anhui = df_anhui[['Gene','Sample', 'Covered_percent',"Read_depth", 'Avg_fold', 'Std_Dev']]
 
    #saving for later reference/use
    df_anhui.to_csv(file.with_suffix('.csv'), index = False)

#Anhui rep2
for file in Path(folder).glob('*/rep2/*/Anhui*_coverage.tsv'):
    df_anhui=pd.read_csv(file, header=0, sep='\t')
     #Adding a new column ("GENE") containing gene segment name only, obtained by deleting characters from #ID column
    df_anhui['Gene'] = df_anhui['#ID'].str[15:]
    df_anhui['Sample'] = os.path.basename(file)
    df_anhui['Read_depth'] = df_anhui[['Plus_reads', 'Minus_reads']].mean(axis=1)
    df_anhui.append(df_anhui)
    #We added this step to clean up the sample names and keep them consistent across all timepoints:
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret1_':'_ferret01_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret2_':'_ferret02_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret3_':'_ferret03_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret4_':'_ferret04_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret5_':'_ferret05_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret6_':'_ferret06_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret7_':'_ferret07_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret8_':'_ferret08_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_ferret9_':'_ferret09_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day1_':'_day01_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day2_':'_day02_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day3_':'_day03_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day4_':'_day04_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day5_':'_day05_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day6_':'_day06_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day7_':'_day07_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day8_':'_day08_'}, regex=True)
    df_anhui['Sample'] = df_anhui.Sample.replace({'_day9_':'_day09_'}, regex=True)
    df_anhui['Sample'] = df_anhui['Sample'].str[:25]
    # drop columns that are not needed (i.e. "Length")
    df_anhui = df_anhui.drop(['#ID', 'Length', 'Ref_GC', 'Covered_bases', 'Plus_reads', 'Minus_reads', 'Read_GC', 'Median_fold'], axis=1)
    #reordering
    df_anhui = df_anhui[['Gene','Sample', 'Covered_percent',"Read_depth", 'Avg_fold', 'Std_Dev']]
   
    #saving for later reference/use
    df_anhui.to_csv(file.with_suffix('.csv'), index = False)
  
# We used the clean csv files containing Anhui data:

dfs = []
for file in Path(folder).glob('*/*/*/Anhui*.csv'):
    dfs.append(pd.read_csv(file))
all_anhui = pd.concat(dfs, ignore_index=True)
    
#create matrix using pivot
anhui_matrix = all_anhui.pivot("Sample", "Gene", "Read_depth")
    
    
#print(df_anhui)
#print(anhui_matrix)
fig = plt.figure(figsize=(5,12))
color=sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=True)
r = sns.heatmap(anhui_matrix, cmap=color, linewidths=0.5,  cbar_kws={'label': 'Read Depth'}, vmin=0, vmax=220000)
r.set_title("Read depth across gene segments per sample \n", fontdict= { 'fontsize': 20})
r.figure.savefig("Anhui_coverage_plot.pdf", bbox_inches="tight", dpi=300)
r.figure.savefig("Anhui_coverage_plot.png", bbox_inches="tight", dpi=300)

#end
