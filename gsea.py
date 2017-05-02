import csv
import numpy as np
import random
from scipy import stats

import gsea_functions as gf

'''Program to calculate the normalized enrichment scores and pvalues according to the publication: 
http://software.broadinstitute.org/gsea/doc/subramanian_tamayo_gsea_pnas.pdf
The program accepts 2 files, which user can specify in the respective input step after the program 
is launched. First file should contain the gene expression data in the format of the example file 
"leukemia.txt". The second file should contain the gene sets data, in the format of the example 
file "pathways.txt". Further input includes 2 phenotype lables (case sensitive, should be the same 
as in the gene expression data). The advantage is, if the data contains multiple phenotypes, they 
can be examined by pairs. Number of permutations defines how many permutations will be run to 
calculate the pvalues and the normalized enrichment scores. The file "gsea_functions.py" should be 
in the same folder as this file. It contains functions required in the calculattion. 
Returns 2 files with results, each for one of the phenotypes. The results files contain gene sets 
names, normalized enrichment scores and pvalues based on the user inputs. The results are sorted 
by the normalized enrichment score. '''

#asking for user input
print('Thank you for starting the GSEA program.')
gene_list = input('Enter the path to the file containing gene expression data\n>')
pathways = input('Enter the path to the file containing gene sets\n>')
phenotype_1 = input('Enter the name of the first phenotype to be examined\n>')
phenotype_2 = input('Enter the name of the second phenotype to be examined\n>')
number_permutations = input('Enter the number of permutations for the calculation\n>')
number_permutations = int(number_permutations)
phenotype_1_columns = []
phenotype_2_columns = []
unranked_data = []
gene_sets = []

#reading in the list of genes and their expression profiles
with open(gene_list, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
#reading the file header and columns for different phenotypes   
#the columns indicating the data for specific phenotypes will be needed for the ranking function 
    header = next(reader)
    for index, column_header in enumerate(header):
        if column_header == phenotype_1:
            phenotype_1_columns.append(index)
        elif column_header == phenotype_2:
            phenotype_2_columns.append(index)
#transforming the values from str to float
    for line in reader:
        new_line = []
        new_line.append(line[0])        
        for element in line[1:]:
            new_line.append(float(element))
        unranked_data.append(new_line)
        
#reading in the list of gene sets
with open(pathways, newline='') as f:
    reader = csv.reader(f, delimiter='\t')  
    for line in reader:
        gene_sets.append(line)

#ranking the gene expression profiles by the signal to noise metric
ranked_list = gf.build_ranked_list(unranked_data, phenotype_1_columns, phenotype_2_columns)

#following lists will contain the calculation results for both phenotypes in question
results_1 = []
results_2 = []

#this loop contains the calculation of the original enrichment scores, enrichment scores for each 
#permutation as well as the respective p values and the normalized enrichment scores
for gene_set in gene_sets:
    print('Analyzing the gene set: ' + gene_set[0])
#calculating the original enrichment score for a given gene set
    enrichment_score = gf.calculate_es(ranked_list, gene_set)
#preparing a list to store the results of the permutations    
    null_hypotheses_scores = [[],[]]
    for permutation in range(number_permutations):
#randomly permuting the phenotype labels, however the number of the samples for each phenotype is 
#preserved
        new_columns_1 = []
        new_columns_2 = []
        new_columns_1 = random.sample(range(1,len(header)), len(phenotype_1_columns))
        for index in range(1,len(header)):
            if index not in new_columns_1:
                new_columns_2.append(index)
#ranking the gene expression profiles according to the new metric                
        permuted_ranked_list = gf.build_ranked_list(unranked_data, new_columns_1, new_columns_2)    
#calculating the enrichment score for a given permutation
        permuted_enrichment_score = gf.calculate_es(permuted_ranked_list, gene_set)
        null_hypotheses_scores[0].append(max(permuted_enrichment_score[0]))
        null_hypotheses_scores[1].append(max(permuted_enrichment_score[1]))
#calculating the pvalues for each phenotype based on the results of the permutations
    pvalue_1 = stats.ttest_1samp(null_hypotheses_scores[0], max(enrichment_score[0]))
    pvalue_2 = stats.ttest_1samp(null_hypotheses_scores[1], max(enrichment_score[1]))
#calculating the normalized enrichment scores based on the results of the permutations
    normalized_es_1 = max(enrichment_score[0])/np.mean(null_hypotheses_scores[0])
    normalized_es_2 = max(enrichment_score[1])/np.mean(null_hypotheses_scores[1])
#storing the data for later output    
    results_1.append([gene_set[0], normalized_es_1, pvalue_1[1]])
    results_2.append([gene_set[0], normalized_es_2, pvalue_2[1]])
#sorting the data by the normalized enrichment score
sorted_results_1 = sorted(results_1, key=lambda x: x[1], reverse=True)
sorted_results_2 = sorted(results_2, key=lambda x: x[1], reverse=True)
#creating the output in form of 2 files, one for each phenotype
with open(phenotype_1 + '_output.txt', 'w') as output:
    output.write('name of the gene set\tnormalized ES\tP-Value\n')
    for row in sorted_results_1:
        output.write(str(row[0]) + '\t' + str(row[1]) + '\t' + str(row[2]) + '\n')
with open(phenotype_2 + '_output.txt', 'w') as output:
    output.write('name of the gene set\tnormalized ES\tP-Value\n')
    for row in sorted_results_2:
        output.write(str(row[0]) + '\t' + str(row[1]) + '\t' + str(row[2]) + '\n')
print('The calculation is completed.\nThank you for using the GSEA program!')
