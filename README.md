# GSEA-Genialis
This repository contains the excersize from the Genialis team. The goal is to implement the GSEA method, described here: 
http://software.broadinstitute.org/gsea/doc/subramanian_tamayo_gsea_pnas.pdf

The program accepts 2 files, which user can specify in the respective input step after the program 
is launched. First file should contain the gene expression data in the format of the example file 
"leukemia.txt". The second file should contain the gene sets data, in the format of the example 
file "pathways.txt". Further input includes 2 phenotype lables (case sensitive, should be the same 
as in the gene expression data). The advantage is, if the data contains multiple phenotypes, they 
can be examined by pairs. Number of permutations defines how many permutations will be run to 
calculate the pvalues and the normalized enrichment scores. 
The file "gsea_functions.py" contains functions required in the calculattion. The file "gsea.py" is the main. Both should be in the 
same folder to run the script. 
Returns 2 files with results, each for one of the phenotypes. The results files contain gene sets 
names, normalized enrichment scores and pvalues based on the user inputs. The results are sorted 
by the normalized enrichment score.

The repository also contains 2 output files, "ALL.txt" and "AML.txt", created during a test run using the first 100 gene sets 
from the "pathways.txt" file. 
During the test run 100 permutations were used for each gene set. The test run took around 3 hours to complete.
