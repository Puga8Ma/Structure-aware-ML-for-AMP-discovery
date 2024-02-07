# Structure-aware-machine-learning-strategies-for-antimicrobial-peptide-discovery

This repository contains the scripts that were used to build a binary and ternary model, which can classify between action mechanisms MDPs, MPPs, and PAPs, as stipulated in the article. 
This repository is divided into data, descriptors, results, and scripts used for the article's methodology. 

The data carpet includes all the fasta files with the sequences collected from the databases (DBAASP, APD3, PDBe, CPPsite 2.0). 

The results folder includes the raw data and clean datasets with the physicochemical properties derived from the peptide sequences, which were used to train the binary and ternary models, as well as the external validation data. the subsets into which the original dataset was divided along with their physicochemical properties. 
The prediction results for each binary, ternary, and subsets model. 

The scripts folder is divided into scripts used to develop the project in the R and Python programming languages. 
 
The first stage contains 3 scripts written in R 3.3.1, ordered in the development of the research project. 

Calculation of physical-chemical properties based on the peptide sequence.R 

Statistical analysis of the comparison of distributions of physical-chemical properties between the MDPs and MPPs groups.R 

Multicollinearity and reduction of redundant information.R 

The second stage contains 3 scripts written in Python 3.6+ 

Calculator of Descriptors (modlamp) from peptide primary sequences.ipynb   

Classification binary models on datasets with 49 descriptors.ipynb   

Classification multilabel models on datasets with 56 descriptors.ipynb  
