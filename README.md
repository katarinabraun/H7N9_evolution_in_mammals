# H7N9 avian influenza virus spillover infections are evolutionarily constrained by stochastic processes

### Authors: 

Katarina M. Braun 1#, Luis A. Haddock III 1#, Chelsea M. Crooks 1, Gabrielle L. Barry 1, Joe Lalli 1, Gabrielle Neumann 1,2, Tokiko Watanabe 2,3,4, Masaki Imai 2,5, Seiya Yamayoshi 2,5, Mutsumi Ito 2, Yoshihiro Kawaoka1,2,5, Thomas C. Friedrich1*

1 Department of Pathobiological Sciences, University of Wisconsin-Madison, Madison, WI, United States of America  
2 Division of Virology, Institute of Medical Science, University of Tokyo, Japan   
3 Department of Molecular Virology, Research Institute for Microbial Diseases, Osaka University  
4 Center for Infectious Disease Education and Research (CiDER), Osaka University   
5 The Research Center for Global Viral Diseases, National Center for Global Health and Medicine Research Institute, Tokyo, Japan.  

\# These authors contributed equally.   
*Corresponding author; tfriedri@wisc.edu 

\*These authors contributed equally

\#These authors contributed equally 

## Abstract
H7N9 avian influenza viruses (AIV) have caused over 1,500 dead end infections in humans since emerging in 2013. Although wild-type H7N9 AIV can transmit by respiratory droplets in ferrets, they have yet to adapt to and cause widespread outbreaks in humans. Previous studies have revealed molecular determinants of H7N9 AIV virus host-switching, but little is known about potential evolutionary constraints on this process. Here we compare patterns of sequence evolution for H7N9 AIV and mammalian H1N1 viruses transmitting in ferrets. We show that three main factors – purifying selection, stochasticity, and very narrow bottlenecks – combine to severely constrain the ability of H7N9 AIV to effectively adapt to mammalian hosts. We only find rare evidence of natural selection favoring new or mammalian-adapting mutations within ferrets and no evidence of natural selection acting during transmission. We conclude that human-adapted H7N9 viruses are unlikely to emerge during typical acute spillover infections. Our findings are instead consistent with a model in which the emergence of a human-transmissible virus would be a rare, though highly consequential, jackpot event. Strategies to limit the total number of spillover infections will also limit opportunities for the virus to win this evolutionary lottery.  

## Repository 

Primary data generated and analyzed in this study have been deposited in the Sequence Read Archive under Bioproject ID: ​**​PRJNA758865**. 

This repository contains the `data_derived` and `code` to replicate the analyses and figures presented in the manuscript.

The scripts are generally in the form of Jupyter Notebooks in the `code` directory. `snpgenie.pl` and the `Bottleneck_size_estimation_exact.r` scripts can also be found in this directory. `Sniffles` and its dependencies can be found in the `Sniffles_pipeline` directory. 

The jupyter notebook should generally be run in the following order: 
1. `data_cleaning_and_intersection_plots.ipynb`
2. `iSNVs_counts.ipynb`
3. `iSNV_over_time.ipynb`
    - requires the user to run for each ferret and transmission pair, however this has already been done and the output files can be found in the `data_derived/iSNVs-over-time` folder. 
4. `iSNVs_over_time_plots.ipynb`
5. `iSNVs_frequencies_and_location.ipynb`
6. `iSNVs_frequencies_and_location_stock.ipynb`
7. `iSNV_frequency_spectrums.ipynb`
8. `pairwise_diversity_plots.ipynb`
    - shows the command line code used to run the `Bottleneck_size_estimation_exact.r`. 
9. `transmission_bottleneck_input_files_and_TV_plots.ipynb`
10. `surveillance_analysis_with_JSON_files_from_nextstrain.ipynb`

## Dependencies
1. pandas
2. numpy
3. sklearn
4. matplotlib
5. seaborn
6. scipy 

The jupyter notebooks generally write data files to the `data_derived/*` directory and figures to the `figures` direcotry. 

The figures formatted for publication that have been edited in Adobe Illustrator can be found in `ms_figures`. The supplementary figures formatted publication that have been edited in Adobe Illustrator can be found in `supplementary_figures`.
