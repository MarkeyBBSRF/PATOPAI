#' Get mapping files which map from IDs(IDs used in the mutation ordering matrix) to gene_mutations for all patients.
#'
#' @param path_mapping The directory of the mapping files. Each mapping files has two column, the first column is the  IDs used in the mutation ordering matrix, The second column is the HugoSymbol_HGVSp.
#' @param strlen A number indicate the length of the patients ID we want to keep.
#' @param patient_list A vector of all patient names.
#' @return All mapping files for all patients as a list.
get_mappingfile_fun=function(path_mapping,strlen,patient_list){
      patient_ids=substr(patient_list,1,strlen)
      print(patient_ids)
      filelist_mapping <- paste(patient_ids,"filtered.mapping.txt",sep = ".")
      filelist_mapping <- paste(path_mapping,filelist_mapping,sep = "/")
      datalist_mapping <- lapply(filelist_mapping,FUN=read.table,header=TRUE,
                           stringsAsFactors=FALSE)
      return(datalist_mapping)
}


#' Get a mutation list with Entrez Gene Ids,functional impact scores,Hugo Symbols for all patients.
#'
#' @param datalist_mapping All mapping files for all patients as a list.
#' @param patient_list A vector of all patient names.
#' @param strlen A number indicate the length of the patients ID we want to keep.
#' @param target TCGA project ID.
#' @param maffile TCGA maf file.
#' @param filewd The directory for mutation ordering matrices.
#' @return A list with two part. i)gene_names: the gene names ordered as in the mutation ordering matrix; ii)mutated.genes.matrix: a list with four columns: 1) Entrez Gene Id 2)functional impact score of each mutations; 3)HugoSymbol_HGVSc(HGVS coding sequence name); 4)HugoSymbol

convert_fun=function(datalist_mapping,patient_list,strlen,target,maffile,filewd){
  patient_ids=substr(patient_list,1,strlen)
  data=read.delim(maffile,header=T,stringsAsFactors = F,comment.char = "#")
  data$PolyPhen_score = gsub(".+\\((.+)\\)",'\\1',data$PolyPhen)
  data$Tumor_Sample_Barcode=substr(data$Tumor_Sample_Barcode,1,strlen)
  gene_index=list()
  gene_names=list()
  for(n in 1:length(patient_list)){
    filename=paste("ssm_data_ordering_matrix_0.csv",sep="")
    filepath=paste(filewd,patient_ids[n],filename,sep="/")
    tmatrix=read.csv(filepath,row.names = 1)
    mutated.genes=row.names(tmatrix)
    for(mpf in 1:nrow(datalist_mapping[[n]])){

	 mutated.genes[mpf] <- datalist_mapping[[n]]$Mutation[mpf]
    }
    mutated.genes.list=sapply(1:length(mutated.genes),function(i) strsplit(mutated.genes[i],"_"))
    mutated.genes.matrix=as.data.frame(matrix(NA,nrow=length(mutated.genes),ncol=3))
    for(i in 1:length(mutated.genes)){
      if(length(mutated.genes.list[[i]])==3){mutated.genes.list[[i]][2]=paste(mutated.genes.list[[i]][2],mutated.genes.list[[i]][3],sep="_")}
      index=which(data$Tumor_Sample_Barcode==patient_ids[n] & data$Hugo_Symbol==mutated.genes.list[[i]][1] & data$HGVSp_Short == mutated.genes.list[[i]][2] &(data$Variant_Type!="INS")&(data$Variant_Type!="DEL")&(data$Variant_Classification!="RNA")&(data$Variant_Classification!="Silent")& !(data$Variant_Classification=="Missense_Mutation" & data$PolyPhen_score =="") )
           if(length(index)>0){
        if(length(index)>1){
        mutated.genes.matrix[i,1]=as.character(data$Entrez_Gene_Id[index[1]])
        mutated.genes.matrix[i,2]=as.numeric(data$PolyPhen_score[index[1]])
        mutated.genes.matrix[i,3]=paste(mutated.genes.list[[i]][1],data$HGVSc[index[1]],sep="_")
        mutated.genes.matrix[i,4]=mutated.genes.list[[i]][1]
      }else{
        mutated.genes.matrix[i,1]=as.character(data$Entrez_Gene_Id[index])
        mutated.genes.matrix[i,2]=as.numeric(data$PolyPhen_score[index])
        mutated.genes.matrix[i,3]=paste(mutated.genes.list[[i]][1],data$HGVSc[index],sep="_")
        mutated.genes.matrix[i,4]=mutated.genes.list[[i]][1]
      }
      }
    }

#print(mutated.genes.matrix)
      index.noscore=which(is.na(mutated.genes.matrix[,1]))
      gene_index[[n]]=1:length(mutated.genes)
      if(length(index.noscore)>0){
        mutated.genes.matrix=mutated.genes.matrix[-index.noscore,]
        gene_index[[n]]=gene_index[[n]][-index.noscore]
      }
      gene_names[[n]]=mutated.genes.matrix[,3]
      index.score1=which(is.na(mutated.genes.matrix[,2]))


      if(length(index.score1)>0){
        mutated.genes.matrix[index.score1,2]=1
      }
    if(n==1){
      mutated.genes.matrix2=mutated.genes.matrix
      }else{
      mutated.genes.matrix2=rbind(mutated.genes.matrix2,mutated.genes.matrix)
      }
  }
  filename2=paste(filewd,"mutated.genes.matrix",sep="/")
  save(mutated.genes.matrix2,file=filename2)
  result=list( gene_names= gene_names,gene_index=gene_index,mutated.genes.matrix=mutated.genes.matrix2)
  return(result)
}

#' Convert all the matrixes' with Chromosome IDs as row and column names to mutations.
#'
#' @param patient_list A vector of all patient names.
#' @param strlen A number indicate the length of the patients ID we want to keep.
#' @param num  number of ordering matrices.
#' @param filewd The directory for mutation ordering matrices.
#' @param gene_names i)gene_names: the gene names ordered as in the mutation ordering matrix.
#' @param target TCGA project ID.
#' @return All matrixes with row and column names converted.
convert_allmatrix=function(patient_list,strlen,num,filewd,gene_index,gene_names,target){
  patient_ids=substr(patient_list,1,strlen)
  for(n in 1:length(patient_list)){
    for(i in 0:(num-1)){
      filename=paste("ssm_data_ordering_matrix_",i,".csv",sep="")
      filepath=paste(filewd,patient_ids[n],filename,sep="/")
      newfilename=paste(patient_ids[n],"_matrix_",i,sep="")
      newfilepath=paste(filewd,patient_ids[n],newfilename,sep="/")

      tmatrix=read.csv(filepath,row.names = 1)
      tmatrix=as.matrix(tmatrix)
      newmatrix=tmatrix[gene_index[[n]],gene_index[[n]]]
      row.names(newmatrix)=gene_names[[n]]
      colnames(newmatrix)=gene_names[[n]]
      save(newmatrix,file=newfilepath)
    }
  }
}


#' Delete highly mutated samples.
#' 
#' @param patient_list A vector of all patient names.
#' @param strlen A number indicate the length of the patients ID we want to keep.
#' @param target TCGA project ID.
#' @param filewd The directory for mutation ordering matrices.
#' @param percent percentage of highly mutated sample to be deleted.
#' @return All matrixes with row and column names converted.
delete_fun=function(patient_list,strlen,target,filewd,percent){
  nsample=length(patient_list)
  nmut=rep(NA,nsample)
  patient_ids=substr(patient_list,1,strlen)
  for(n in 1:length(patient_list)){
    filename=paste("ssm_data_ordering_matrix_0.csv",sep="")
    filepath=paste(filewd,patient_ids[n],filename,sep="/")
    tmatrix=read.csv(filepath,row.names = 1)
    nmut[n]=nrow(tmatrix)
  }
  sort_nmut=sort(nmut,decreasing = F,index.return=T)
  nsample=ceiling(nsample*percent)
  patient_list_del=patient_list[sort_nmut$ix[1:nsample]]
  nor.sample=substr(patient_list_del,1,strlen)
}



#' Reduce matrix to selected pathways--get the index for all mutations within selected pathways
#'
#' @param patient_list A vector of all patient names.
#' @param strlen A number indicate the length of the patients ID we want to keep.
#' @param target TCGA project ID.
#' @param filewd The directory for mutation ordering matrices with IDs converted.
#' @param pathway.list A list, the name of each elements is the name of pathways/genesets, the vectors are genes included in the pathways/genesets.
#' @param mutated.genes.matrix A list with four columns: 1) Entrez Gene Id for each mutations 2)functional impact score of each mutations; 3)HugoSymbol_HGVSc(HGVS coding sequence name); 4)HugoSymbol
#' @return A list with six parts, 1) index.list: order index for all mutations in the selected pathways; 2)patient_list: list for all patients who has mutaions in selected pathways; 3)no.pathway: number of redefined gene sets which has no overlap with each other. 4)pathway: pathway list for each mutations. 5)prob: functional impact score ordered as the index for corresponding mutations. 6)genesets: A list which contains the redefined gene sets which has no overlap with each other. The name of each elements is the name of pathways/genesets, the vectors are genes included in the pathways/genesets.

reduce_index.fun=function(patient_list,strlen,target,filewd,pathway.list,mutated.genes.matrix){
  patient_ids=substr(patient_list,1,strlen)
  index.list=list()
  prob=list()
  pathway=list()
  index_row=0
  for(n in 1:length(patient_list)){
    filename=paste(patient_ids[n],"_matrix_",0,sep="")
    filepath=paste(filewd,patient_ids[n],filename,sep="/")
    load(filepath)
    mutated.genes=row.names(newmatrix)
    
    
    mutated.index=sapply(1:length(mutated.genes),function(i) which(mutated.genes.matrix[,3]==mutated.genes[i])[1])
    
    mutated.genes_id=mutated.genes.matrix[mutated.index,1]
    
    mutated.genes_score=mutated.genes.matrix[mutated.index,2]
    unique.mutated.genes_id=unique(mutated.genes.matrix[,1])
    genesets=redefine.geneset.fun2(pathway.list,unique.mutated.genes_id)
    index=which(mutated.genes_id %in% unlist(genesets))
    
    no.pathway=length(genesets)
    
    
    if(length(index)>0){
      index_row = index_row+1
      index.list[[index_row]]=index
      prob[[index_row]]=mutated.genes_score[index]
      
      pathway[[index_row]]=rep(NA,length(index))
      for(j in 1:length(genesets)){
        pathway[[index_row]][which(mutated.genes_id[index] %in% genesets[[j]])]=j
      }
      
      
    }else{
      patient_list[n]=NA
    }
  }
  
  patient_list=patient_ids[-which(is.na(patient_list))]
  result=list(index.list=index.list,patient_list=patient_list,no.pathway=no.pathway,pathway=pathway,prob=prob,genesets=genesets)
  return(result)
}

#' Reduce matrix to selected pathways--keep matrix with top probabilities until accumulated probabilities reached 95%
#'
#' @param index.list order index for all mutations in the selected pathways.
#' @param patient_list A vector of all patient names.
#' @param no.pathway number of redefined gene sets which has no overlap with each other.
#' @param pathway pathway list for each mutations.
#' @param target TCGA project ID.
#' @param filewd The directory for mutation ordering matrices with IDs converted.
#' @param nmc number of ordering matrices.
#' @return A list with three parts, 1)A list with all matrices kept for all samples; 2)nmatrix: A vector indicate the number of reduced matrices for each patients; 3)prob_matrix: the corresponding posterior probabilities for each matrices kept.  

reduce_matrix.fun=function(index.list,patient_list,no.pathway,pathway,target,filewd,nmc){
  matrix.list2=list()
  prob_matrix=list()
  nmatrix_cumu=0
  nmatrix=rep(NA,length(patient_list))
  for(n in 1:length(patient_list)){
    combine_result=combine.fun(patient_list[n],target,filewd,nmc,index.list[[n]])
    prob_matrix[[n]]=combine_result$matrix.prob
    matrix.list1=combine_result$matrix.list1
    summation=sum(prob_matrix[[n]])
    sortmatrix=sort(prob_matrix[[n]],decreasing = T,index.return=T)
    nmatrix2=length(matrix.list1)
    cum=0
    for(i in 1:nmatrix2){
      cum=cum+sortmatrix$x[i]
      if(n>1){
        matrix.list2[[i+nmatrix_cumu]]=as.matrix(matrix.list1[[sortmatrix$ix[i]]])
      }else{
        matrix.list2[[i]]=as.matrix(matrix.list1[[sortmatrix$ix[i]]])
      }
      if((cum/summation)>0.95){
        break
      }
    }
    prob_matrix[[n]]=prob_matrix[[n]][sortmatrix$ix[1:i]]
    nmatrix[n]=i
    
    nmatrix_cumu=nmatrix_cumu+i
  }
  result=list(matrix.list2=matrix.list2,nmatrix=nmatrix,prob_matrix=prob_matrix)
  return(result)
}


#' Reduce matrix to selected pathways--combine the same matrices within selected pathways.
#'
#' @param patient_ids Id for one patient.
#' @param target TCGA project ID.
#' @param filewd The directory for mutation ordering matrices with IDs converted.
#' @param nmc number of ordering matrices.
#' @param index order index for a patient.
#' @return A list with two parts, 1)matrix.list1: A list with all reduced matrices for all samples; 2)matrix.prob: the corresponding posterior probabilities for each reduced matrices.  
combine.fun=function(patient_ids,target,filewd,nmc,index){
  matrix.list1=list()
  filename=paste(patient_ids,"_matrix_",0,sep="")
  filepath=paste(filewd,patient_ids,filename,sep="/")
  load(filepath)
  
  tmatrix=newmatrix[index,index]
  tmatrix[which(tmatrix==0)]=1
  tmatrix[which(tmatrix==-1)]=0
  matrix.list1[[1]]=tmatrix
  matrix.prob=rep(0,nmc)
  for(i in 1:(nmc-1)){
    ncheck=0
    filename=paste(patient_ids,"_matrix_",i,sep="")
    filepath=paste(filewd,patient_ids,filename,sep="/")
    load(filepath)
    tmatrix=newmatrix[index,index]
    tmatrix[which(tmatrix==0)]=1
    tmatrix[which(tmatrix==-1)]=0
    for(j in 1:length(matrix.list1)){
      ncheck=ncheck+1
      if(sum(tmatrix==matrix.list1[[j]])==length(index)*length(index)){
        matrix.prob[j]=matrix.prob[j]+1
        break
      }
      if(ncheck==length(matrix.list1)){
        matrix.list1[[length(matrix.list1)+1]]=tmatrix
        matrix.prob[length(matrix.list1)+1]=matrix.prob[length(matrix.list1)+1]+1
      }
    }
  }
  matrix.prob=matrix.prob[1:length(matrix.list1)]
  matrix.prob=matrix.prob/nmc
  result=list(matrix.list1=matrix.list1,matrix.prob=matrix.prob)
  return(result)
}


#####generate pathway.database##########
pathwaylist.fun=function(groupgene){
  pathway=unique(groupgene[,1])
  pathwaylist=list()
  for(i in 1:length(pathway)){
    pathwaylist[[i]]=groupgene[which(groupgene[,1]==pathway[i]),2]
  }
  names(pathwaylist)=pathway
  return(pathwaylist)
}

###########redefine geneset############
redefine.geneset.fun2=function(driver.database,mutated.genes){
  driver.geneset=names(driver.database)
  no.dri=length(driver.geneset)
  mutated.genes=as.integer(mutated.genes)
  indicate=sapply(1:no.dri,function(j) mutated.genes %in% driver.database[[j]] )
  indicate2=t(sapply(1:length(mutated.genes),function(i) indicate[i,]*(1:no.dri)))
  gene.score=sapply(1:length(mutated.genes),function(i) sum(2^indicate2[i,which(indicate2[i,]!=0)]))
  unique.score=sort(unique(gene.score))
  
  unique.score=unique.score[which(unique.score!=0)]
  unique.score.index=sapply(1:length(unique.score),function(i) which(gene.score==unique.score[i]))
  genes.per.genesets=sapply(1:length(unique.score),function(i) mutated.genes[unique.score.index[[i]]])
  pathway.per.genesets.index=t(sapply(1:length(unique.score),function(i) indicate[unique.score.index[[i]][1],]))
  pathway.per.genesets=sapply(1:length(unique.score),function(i) driver.geneset[which(pathway.per.genesets.index[i,]==1)])
  names(genes.per.genesets)=sapply(1:length(pathway.per.genesets), function(i) paste(pathway.per.genesets[[i]],collapse = ","))
  return(genes.per.genesets)
}

