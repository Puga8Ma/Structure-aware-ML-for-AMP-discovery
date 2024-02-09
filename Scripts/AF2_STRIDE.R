# Read .pdb files and predict their secondary structures with STRIDE
#https://search.r-project.org/CRAN/refmans/bio3d/html/dssp.html

#Working directory
setwd("/Users/(yourname)/Desktop/")

# Load packages
install.packages("bio3d")
install.packages("DSSP")
library(bio3d)
library(DSSP)

#Load 2-3 examples 
pdb <- read.pdb('./AF2_predictions/AF2_results_MDPs/MDPs_1.pdb')
pdb2 <- read.pdb('./AF2_predictions/AF2_results_MDPs/MDPs_10.pdb')
pdb3 <- read.pdb('./AF2_predictions/AF2_results_MDPs/MDPs_115.pdb')

#Run STRIDE to obtain secondary structures
sse <- stride(pdb)
[1] "C" "H" "H" "H" "H" "H" "H" "H" "H" "H" "H"
[12] "H" "H" "H" "H" "H" "H" "H" "H" "C"

sse <- stride(pdb2)
[1] "H" "H" "H" "H" "H" "H" "H" "H" "H" "H" "H"
[12] "H" "H" "H" "H" "H" "H" "H" "H" "H" "H" "H"
[23] "H" "H" "H" "H" "H" "C"

sse <- stride(pdb3)
[1] "C" "E" "E" "E" "C" "C" "H" "H" "H" "H" "H"
[12] "H" "H" "H" "H" "H" "H" "T" "T" "T" "T" "C"
[23] "C" "H" "H" "H" "H" "H" "H" "H" "H" "C" "E"
[34] "E" "E" "C" "C" "C" "C" "C" "C" "T" "T" "T"
[45] "T" "T" "C"

# Percentages
sum(sse$helix$length)/sum(pdb$calpha) * 100
[1] 40.42553
sum(sse$sheet$length)/sum(pdb$calpha) * 100
[1] 12.76596
sum(sse$turn$length)/sum(pdb$calpha) * 100
[1] 19.14894
#Note: turn (T) and coil (C) are equivalent in 3-state secondary structures system


# Define the R function "calculate_sse_percentages"
calculate_sse_percentages <- function(path_to_pdbs){
  
  # Initialize empty vectors to store results
  pdb_names <- character()
  sse_sequences <- character()
  helix_percent <- numeric()
  sheet_percent <- numeric()
  coil_percent <- numeric()
  
  # List all PDB files in the specified path
  pdb_files <- list.files(path_to_pdbs, pattern = "\\.pdb$", full.names = TRUE)
  
  # Loop over each PDB file
  for (pdb_file in pdb_files) {
    
    # Read the PDB structure
    pdb <- read.pdb(pdb_file)
    
    # Calculate secondary structure using stride function from Bio3D package
    stride_result <- stride(pdb)
    print(pdb_file)
    
    # Extract the secondary structure elements (SSE) from the stride result as a string
    sse_seq <- paste0(stride_result$sse, collapse = "")
    
    # Calculate percentages of helix, sheet, and turn
    helix_percentage <- sum(stride_result$helix$length)/sum(pdb$calpha) * 100
    sheet_percentage <- sum(stride_result$sheet$length)/sum(pdb$calpha) * 100
    coil_percentage <- 100 - (helix_percentage + sheet_percentage)
    
    # Append results to vectors
    pdb_names <- c(pdb_names, basename(pdb_file))
    sse_sequences <- c(sse_sequences, sse_seq)
    helix_percent <- c(helix_percent, helix_percentage)
    sheet_percent <- c(sheet_percent, sheet_percentage)
    coil_percent <- c(coil_percent, coil_percentage)
  }
  
  # Create a dataframe to store the results
  result_df <- data.frame(
    PDB_Name = pdb_names,
    SSE_Sequence = sse_sequences,
    '%_Helix' = helix_percent,
    '%_Sheet' = sheet_percent,
    '%_Coil' = coil_percent
    )
  
  return(result_df)
}

# Warnings
options(warn=-1)

# Example usage:
path_to_pdbs <- "../AF2_predictions/test/"
result_dataframe <- calculate_sse_percentages(path_to_pdbs)
print(result_dataframe)
write.csv(result_dataframe, "./STRIDE_results/test_stride_results.csv", row.names = FALSE)

# Note : the connection to measure SSE for each PDB structure is lost after 260 runs. We divided our datasets so these characteristics are measured for a list up to 260 PDB structures every time.

#MDPs
path_to_pdbs <- "../AF2_predictions/AF2_results_MDPs/"
MDPs_dataframe <- calculate_sse_percentages(path_to_pdbs)
print(MDPs_dataframe)
write.csv(MDPs_dataframe, "./STRIDE_results/MDPs_stride_results.csv", row.names = FALSE)

#MPPs
path_to_pdbs <- "../AF2_predictions/AF2_results_MPPs"
MPPs_dataframe <- calculate_sse_percentages(path_to_pdbs)
print(MPPs_dataframe)
write.csv(MPPs_dataframe, "./STRIDE_results/MPPs_stride_results.csv", row.names = FALSE)

#PAPs
path_to_pdbs <- "../AF2_predictions/AF2_results_PAPs"
PAPs_dataframe <- calculate_sse_percentages(path_to_pdbs)
print(PAPs_dataframe)
write.csv(PAPs_dataframe, "./STRIDE_results/PAPs_stride_results.csv", row.names = FALSE)

#External_validation
path_to_pdbs <- "../AF2_predictions/AF2_results_externalvalidation"
EV_dataframe <- calculate_sse_percentages(path_to_pdbs)
print(EV_dataframe)
write.csv(EV_dataframe, "./STRIDE_results/EV_stride_results.csv", row.names = FALSE)

#Alpha_helices
path_to_pdbs <- "./AF2_predictions/AF2_results_alphahelices/"
EAHs_dataframe <- calculate_sse_percentages(path_to_pdbs)
print(EAHs_dataframe)
write.csv(EAHs_dataframe, "./AF2_predictions/STRIDE_results/EAHs_stride_results.csv", row.names = FALSE)

# Some AF2-generated PDB structures (.pdb) produced some errors ("No hydrogen found in temporary files") using STRIDE. 
# Coincidentally, many of these structures presented poor pLDDT scores (<50).

# PDB structure references
path_to_pdbs <- "../AF2_predictions/PDB_examples/"
result_dataframe <- calculate_sse_percentages(path_to_pdbs)
print(result_dataframe)
write.csv(result_dataframe, "./STRIDE_results/pdb_stride_results.csv", row.names = FALSE)

