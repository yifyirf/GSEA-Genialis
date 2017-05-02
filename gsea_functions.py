import numpy as np

def build_ranked_list(data, columns_1, columns_2):
    '''Function to sort the list of gene expression profiles by the signal to noise ratio 
    between 2 phenotypes. Uses the unsorted data and 2 lists with the column numbers 
    associated with the phenotypes. Returns the list of lists containing:
    gene name    expression mean value phenotype1   standard deviation phenotype1   
    expression mean value phenotype2   standard deviation phenotype2    signal to noise ratio'''
    processed_data = []    
    for row in data:
        processed_line = []
        processed_line.append(row[0])
        values_of_interest_1 = [row[index] for index in range(len(row)) if index in columns_1]
        values_of_interest_2 = [row[index] for index in range(len(row)) if index in columns_2]
        processed_line.append(np.mean(values_of_interest_1))
        processed_line.append(np.std(values_of_interest_1))
        processed_line.append(np.mean(values_of_interest_2))
        processed_line.append(np.std(values_of_interest_2))
        processed_line.append(np.abs(
            (np.mean(values_of_interest_1)-np.mean(values_of_interest_2))/
            (np.std(values_of_interest_1)+np.std(values_of_interest_2))))
        processed_data.append(processed_line)
    processed_data.sort(key=lambda x: x[5], reverse=True)
    return processed_data
    
    
def calculate_es(ranked_gene_list, gene_set, exponent=1):
    '''Function to calculate the enrichment score of a given gene set using the ranked gene list. 
    Accepts the ranked gene list, containing the gene name and the mean expression values for 
    both phenotypes in each row; and a gene set, containing the name of the set, 
    a comment about the gene set, and the genes in the set starting at position [2]. The exponent 
    to weight the running sum steps is set to 1 by default. Returns an list with 2 lists, each 
    for one of the phenotypes, containing the values of the running sum.'''
    running_sum = [[0],[0]]
#the value being subtracted, in case of a gene not being a part of the gene set in question, is 
#constant for a given gene set        
    miss_value = 1/(len(ranked_gene_list)-len(gene_set[2:]))
#calculation of N_R is required prior to the calculation of the running sum
#as the running sum is calculated for both phenotypes in question, so does the N_R 
    weight_phenotype_1 = []
    weight_phenotype_2 = []  
    for gene in gene_set[2:]:
        for row in ranked_gene_list:
            if row[0] == gene:
                weight_phenotype_1.append(np.abs(row[1])**exponent)
                weight_phenotype_2.append(np.abs(row[3])**exponent)
    N_R_1 = np.sum(weight_phenotype_1)
    N_R_2 = np.sum(weight_phenotype_2)
#calculation of the running sum
    for row in ranked_gene_list:
        if row[0] not in gene_set[2:]:
            running_sum[0].append(running_sum[0][-1]-miss_value)
            running_sum[1].append(running_sum[1][-1]-miss_value)
        else: 
            running_sum[0].append(running_sum[0][-1]+np.abs(row[1])**exponent/N_R_1)
            running_sum[1].append(running_sum[1][-1]+np.abs(row[3])**exponent/N_R_2)
    return running_sum
