# PATOPAI

### 1.Input Data Sets:

“Input”: All the matrices represent the structure of phylogenetic trees output by PhyloWGS for 31 patients of TCGA pancreatic Adenocarcinoma dataset. The matrices for each patient are under separated files. For each matrix, each column and row represents a specific SNV in the patient, and each element of the matrix suggests whether the corresponding mutation in the row is a downstream mutation of the corresponding mutation in the column (0) or not (1) in the corresponding phylogenetic tree.

“TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf”: All the  mutations in TCGA pancreatic Adenocarcinoma dataset. Which contains the sample names, mutated genes and Polyphen-2 scores information.

"pathway_database14": Contains the mapping of pathways to genes (genes entrez id). The name of the list are the name of pathways, the elements of each list are the genes included in the pathways.

"mapping": Contains the mapping of Chromosome ID to mutations. The first column is the Chromosome ID and start position, the second column is the Hugo symbol of the mutation and the protein sequence of the variant in HGVS recommended format(HGVSp).


### 2.Example files:

The file "example.R” is an example of R codes which read the all the files mentioned above together and analyze the Cell cycle and MAPK signaling pathway for TCGA ancreatic Adenocarcinoma dataset. The file returns “pairprob.txt” which contains the probabilities of A>B, A<B and A=B for all pairs of pathways.


### 3.Function library:

The file "function_library.r" contains R codes used in estimating the order of mutations and simulations. It calls C code "loglik-call2_2.c”. One needs to compile "loglik-call2_2.c” using the following command before calls it:

R CMD SHLIB loglik-call2_2.c
