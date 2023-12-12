########Data Processing########
mapping="mapping"
filewd="input"
tcgafile="TCGA.PAAD.mutect.fea333b5-78e0-43c8-bf76-4c78dd3fac92.DR-10.0.somatic.maf"
target="PAAD"
nmc=5000
strlen=12

source("function_library_dataprocessing.R")
patient_list=read.table("patient_list_PAAD")
patient_list=as.vector(as.matrix(patient_list))
datalist_mapping=get_mappingfile_fun(mapping,strlen,patient_list)  # Get mapping files which map from Chromosome IDs (IDs used in the mutation ordering matrix) to mutations for all patients.
result=convert_fun(datalist_mapping,patient_list,strlen,target,tcgafile,filewd) #Get a mutation list with Entrez Gene Ids,functional impact scores,Hugo Symbols for all patients.
convert_allmatrix(patient_list,strlen,nmc,filewd,result$gene_index,result$gene_names,target) #Convert all the matrixes' with Chromosome IDs as row and column names to mutation IDs.

######Calculation of Pivor Probabilities Matrix######
source("function_library_calculation.R")
dyn.load("loglik-call2_2.so")
library(gtools)
library(alabama)
library(doMC) # for parallel computing
registerDoMC()
require(methods)
parallel <- F
N <- 10 # N is the number of times that the optimization with different inital values is repeated. If the variable "parallel" is TRUE, then the repeated optimization is done using parallel computing R package doMC. If it is not possible, let "parallel" be FALSE, then the optimization is repeated using for loop. We recommend N to be at least 10.
MAX=3 # assume the same distribution after MAX events. In this case, assume P_{k,i}=P_{4,i} for k>4
p1=1 #index of pathways to analyze  
p2=2 #index of pathways to analyze

load("pathway_database14")
load("pathwayindex")
pathway.database=pathway.database[pathwayindex[[target]]] #get all pathway mapping information for all core pathways for a cancer type.
load("input/mutated.genes.matrix")
mutated.genes=unique(mutated.genes.matrix2[,1])

result1=reduce_index.fun(patient_list,strlen,target,filewd,pathway.database[c(p1,p2)],mutated.genes.matrix2) #Get the index for all mutations within selected pathways.
index.list=result1$index.list
patient_list=result1$patient_list
no.pathway=result1$no.pathway
pathway=result1$pathway

result2=reduce_matrix.fun(index.list,patient_list,no.pathway,pathway,target,filewd,nmc) #reduce matrix to selected pathways, and only keep matrix with top probabilities with a cumulative probability reaches 95%.

matrix.list2=result2$matrix.list2
prob=result1$prob
nmatrix=result2$nmatrix
prob_matrix=result2$prob_matrix

max.mut=max(sapply(1:length(matrix.list2),function(i) nrow(matrix.list2[[i]])))  #calculate the number of maximal mutations among all patients
prob_nondammut=probnon(pathway,no.pathway,prob,nmatrix) #calculate the probabilities of a random non-functional mutations occurs in each pathways
mle<-order_estimate(matrix.list2,N,parallel,prob,pathway,no.pathway,nmatrix,prob_matrix,prob_nondammut,0.1,10) # mle of P_{k,i}, the probability of kth mutation occurring in pathway/geneset i

genesets=result1$genesets
round.mle=round(mle[[1]],3)
row.names(round.mle)=names(genesets)
driver.geneset=names(pathway.database)
pair.order.prob=pair.order.prob.fun(driver.geneset[c(p1,p2)],round.mle) #Calculate the probabilities of A<B, A>B and A=B
write.table(pair.order.prob,"pairprob.txt",row.names=F,col.names=T)


