#!/bin/bash

#Running the full joined_families pipeline
python shared_nbr_candidates.py ../data/gene_families/modified_ALL_GENE_file 6 ../results/pipeline_results/candidatePairs.txt 

python silhouette.py all_families ../data/GenFamClust_files/translatedGenFam.syc ../data/gene_families/modified_ALL_GENE_file ../results/pipeline_results/nc_silhouette.txt

python silhouette.py joined ../data/GenFamClust_files/translatedGenFam.syc ../data/gene_families/modified_ALL_GENE_file ../results/pipeline_results/nc_silhouette.txt ../results/pipeline_results/candidatePairs.txt ../results/pipeline_results/joined_all_candidates.txt

python silhouette.py shuffled ../data/GenFamClust_files/translatedGenFam.syc ../data/gene_families/modified_ALL_GENE_file ../results/pipeline_results/shuffled_weights.txt ../results/pipeline_results/shuffled_nc_silhouette.txt ../results/pipeline_results/shuffled_joined_families.txt 1

python FDR_analysis.py ../results/pipeline_results/joined_all_candidates.txt ../results/pipeline_results/shuffled_joined_families.txt 0.01


