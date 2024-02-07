# Structure-aware-machine-learning-strategies-for-antimicrobial-peptide-discovery

This repository contains the scripts that were used to build binary and ternary models capable of classifying action mechanisms MDPs, MPPs, and PAPs as described in the article. The repository is divided into data, descriptors, results, and scripts used for the article's methodology.

The data folder includes all the fasta files with the sequences collected from various databases such as DBAASP, APD3, PDBe, and CPPsite 2.0.

The results folder includes the raw and clean datasets with physicochemical features derived from peptide sequences used to train the binary and ternary models, along with external validation data. It also contains subsets into which the original dataset was divided, along with their physicochemical properties. Finally, the prediction results for each binary, ternary, and subsets model are also included.

Scripts folder is divided into scripts used to develop the project in R 3.4 and Python 3.6+ programming languagesorder as follows: 

1. Calculation of physicochemical properties based on peptide sequence.R
2. Statistical analysis of the comparison of distributions of physical-chemical properties between the MDPs and MPPs groups.R
3. Multicollinearity and reduction of redundant information.R
4. Calculator of Descriptors (modlamp) from peptide primary sequences.ipynb
5. Classification binary models on datasets with 49 descriptors.ipynb
6. Classification multilabel models on datasets with 56 descriptors.ipynb  
