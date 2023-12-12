#' Function for estimating the pivot mutation probabilities matrix--with multiple initial values.
#'
#' @param matrix A list of (reduced and combined) mutation ordering matrices for all patients.
#' @param N The number of times that the optimization with different inital values is repeated.
#' @param parallel A Boolean value indicate whether parallel computing is to be used.
#' @param score.list A list contains the functional impact scores for each mutations in each sample.
#' @param pathway Pathway list for each mutations.
#' @param no.pathway Number of redefined gene sets which has no overlap with each other.
#' @param nmatrix A vector indicates the number of reduced matrices for each patients.
#' @param prob_matrix The corresponding posterior probabilities for each matrices.  
#' @param prob_nondammut The probabilities a non-functional mutation occurs at each pathways.
#' @param cutoff A number indicates the cutoff value of functional impact score for treating an mutations as determinated functional/non-functional. For example, a cutoff value of 0.05 means that all impact scores which are less than 0.05 will be set as 0 and all scores greater than 0.95 will be set as 1.The default value for this parameter is 0.
#' @param max_unfix A number indicates the maximum value of undeterminated mutations.
#' @return A matrix which indicate the probabilities of each pathways/genesets mutated as the k-th mutational event.

order_estimate <- function(matrix,N,parallel,score.list,pathway,no.pathway,nmatrix,prob_matrix,prob_nondammut,cutoff=0,max_unfix=0)
{
  
  weight=weightcalculation(matrix,score.list,pathway,no.pathway,nmatrix,prob_matrix,prob_nondammut ,cutoff,max_unfix,p1,p2)
  save(weight,file=paste(p1,p2,"weight",sep=""))
  ## Run optimization with different inital values N times using parallel computing if the variable "parallel" is TRUE :
  if(parallel)
    tmp <- foreach (kk = 1:N) %dopar%  main.function(weight$weight,weight$al,MAX,weight$max.dmut,no.pathway)  
  else ## Otherwise, run for loop
  {
    tmp <- vector("list", N)
    for(i in 1:N) tmp[[i]] <- main.function(weight$weight,weight$al,MAX,weight$max.dmut,no.pathway)  
  }
  
  minusloglik=rep(Inf,N)
  for(l in 1:N)
    if(is.list(tmp[[l]])) minusloglik[l]=tmp[[l]][[2]][length(tmp[[l]][[2]])]
  
  result <- tmp[[which(minusloglik==min(minusloglik))[1]]]#[[1]] #find the one giving the maximum likelihood
  
  return(result)
  
}


#' Function for estimating the pivot mutation probabilities matrix--main matrix for one initial value.
#'
#' @param weight The accumulated weight for each possible pathway/geneset mutational orders with all patients. 
#' @param al all possible mutational orders with MAX predefined. 
#' @param MAX Assume the same distribution after MAX events.
#' @param max.dmut  maximal mumber of functional mutations.
#' @param no.pathway Number of redefined gene sets which has no overlap with each other.
#' @return A list with three parts, 1)AllPtrue: The estimated pivot mutation probabilities matrix; 2) loglikelihood: loglikelihood for the estimated matrix; 3)AllPseq: the updated matrix in each iteration step.

main.function <- function(weight,al,MAX,max.dmut,no.pathway)  
{
  #    nomut <- max(rowSums(Y))
  nparam <- min(MAX+1,max.dmut)
  division <- vector("list",nparam)
  for(i in 1:(nparam-1)) division[[i]] <- i
  division[[nparam]] <- nparam:max.dmut
  
  
  ###############################################################
  
  AllPtrue <- matrix(0,nr=no.pathway,nc=max.dmut)
  
  
  for(j in 1:nparam)
  {
    init=runif(no.pathway) ;init=init/sum(init)
    AllPtrue[,division[[j]]] <-  init
    #     AllPtrue[,division[[j]]] <-  c(0.5,0.5,0)
  }
  
  prev0=10^30
  prev=10^30
  n=0
  decrease=1
  
  no.repeat=50
  loglikelihood=rep(NA,60)
  AllPseq=list()
  decrease.limit=10^(-6)
  
  
  while((decrease>decrease.limit || decrease<0) & n<=no.repeat )
  {
    n=n+1
    print(n)        
    prev0=prev              # For each k, the length of \vec{P_{k}} optimized is N-1 where N is the number of driver genes. This is because \sum_{i=1}^{N} P_{k,i} = 1 and the value of P_{k,N} is determeined by the other N-1 P_{k,i}, thus we optimize only for P_{k,i} for i=1,..N-1. we let P_{k,N}=1-\sum_{i=1}^{N-1}{P_{k,i}}
    for(sel in 1:nparam)    # find \vec{P_{k}} maximizing the log likelihood  for each k in turn with the constraint that 0<P_{k,i} for i=1,..N-1 and 0< 1-\sum_{i=1}^{N-1}{P_{k,i}}
    {
      #          repeat{
      no.fail=0
      error=T
      init=runif(no.pathway-1,-5,5);  # initial values for \vec{u_{k}}
      
      if(sel==1 & n==1){
        x <- try(optim(par=init,fn=loglik,control=list(maxit=10000),AllP=AllPtrue,division=division,sel=sel,weight=weight,al=al),silent=T)
      }else{
        
        hin <- function(x,AllPtrue,division,sel,weight,al){     # a vector function specifying inequality constraints such that hin[j] > 0 for all j used in the constrained optimization function auglag
          P=c(exp(x)/(sum(exp(x))+1),1/(sum(exp(x))+1))
          if(sel>1){
            sign=-1
          }else{
            sign=1
          }
          return((entropy+sum(P*log(P)))*sign)
        }
        
        x <- try(auglag(par=init,fn=loglik,hin=hin,control.optim=list(maxit=10000),control.outer=list(trace=F),AllP=AllPtrue,division=division,sel=sel,weight=weight,al=al),silent=T)
      }
      error=is(x,"try-error")
      if(is(x,"try-error")){
        print(x)

        no.fail=no.fail+1
        
      }
      if(!is(x,"try-error")){
        fval=x$value
        decrease=prev-fval
        prev=fval
        AllPtrue[-no.pathway,division[[sel]]] <- exp(x$par)/(sum(exp(x$par))+1)
        AllPtrue[no.pathway,division[[sel]]] <- 1/(sum(exp(x$par))+1)
        entropy=-sum(AllPtrue[,sel]*log(AllPtrue[,sel]))
      }
    }
    decrease=prev0-prev
    loglikelihood[n]=prev   
    AllPseq[[n]]=AllPtrue
  }
  
  loglikelihood=loglikelihood[1:n]    
  return(list(AllPtrue,loglikelihood,AllPseq))
}

#' Function for calculating the "accumulated weight" for each possible pathway/geneset mutational orders.
#'
#' @param matrix A list of (reduced and combined) mutation ordering matrices for all patients.
#' @param prob A list contains the functional impact scores for each mutations in each sample.
#' @param pathway Pathway list for each mutations.
#' @param no.pathway Number of redefined gene sets which has no overlap with each other.
#' @param nmatrix A vector indicates the number of reduced matrices for each patients.
#' @param prob_matrix The corresponding posterior probabilities for each matrices.  
#' @param prob_nondammut The probabilities a non-functional mutation occurs at each pathways.
#' @param cutoff A number indicates the cutoff value of functional impact score for treating an mutations as determinated functional/non-functional. For example, a cutoff value of 0.05 means that all impact scores which are less than 0.05 will be set as 0 and all scores greater than 0.95 will be set as 1.The default value for this parameter is 0.
#' @param max_unfix A number indicates the maximum value of undeterminated mutations.
#' @param p1 The name/index of pathway1(only for file names).
#' @param p2 The name/index of pathway2(only for file names).
#' @return A list with three parts, 1)weight: The accumulated weight for each possible pathway/geneset mutational orders with all patients. 2)al: all possible mutational orders with MAX(assume the same distribution after MAX events) predefined. 3)max.dmut: maximal mumber of functional mutations.

weightcalculation <- function(matrix,prob,pathway,no.pathway,nmatrix,prob_matrix,prob_nondammut ,cutoff,max_unfix,p1,p2)  
{
  #    nomut <- max(rowSums(Y))
  if(cutoff>0){
    for(i in 1:length(prob)){
      for(k in 1:length(prob[[i]])){
        if(prob[[i]][k]<=cutoff || prob[[i]][k]>=1-cutoff ) {prob[[i]][k]=round(prob[[i]][k])}
      }
    }
  }
  
  
  if(max_unfix!=0){
    for(i in 1:length(prob)){
      if(sum(prob[[i]]!=0 & prob[[i]]!=1)>max_unfix){
        sort.index=sort(abs(prob[[i]]-0.5),decreasing = T,index.return=T)[[2]]
        round.index=sort.index[1:(length(prob[[i]])-max_unfix)]
        prob[[i]][round.index]=round(prob[[i]][round.index]) 
      }
    }
  }
  
  nonc_index=sapply(1:length(prob),function(i) which(prob[[i]]==0))
  nnonc=sapply(1:length(prob),function(i) length(nonc_index[[i]]))
  #  print(prob)
  #  print(nmatrix)
  #  print(nonc_index)
  #  print(nnonc)
  cnmatrix=0
  for(i in 1:length(prob)){
    if(nnonc[i]>0){ # & (length(prob[[i]])>nnonc[i])){
      prob[[i]]=prob[[i]][-nonc_index[[i]]]
      pathway[[i]]=pathway[[i]][-nonc_index[[i]]]
      #print("i")
      #print(i)
      for(j in 1:nmatrix[i]){
        #print(j)
        matrix[[j+cnmatrix]]=as.matrix(matrix[[j+cnmatrix]][-nonc_index[[i]],-nonc_index[[i]]])
      }
    }
    cnmatrix=cnmatrix+nmatrix[i]
  }
  #print(matrix)
  
  max.dmut=max(sapply(1:length(prob), function(i) sum(prob[[i]]!=0 )))
  #  print(max.dmut)
  #  print(sum(sapply(1:length(prob), function(i) 2^sum(prob[[i]]!=0 & prob[[i]]!=1))))
  
  #  print(max.dmut)
  base=max.dmut+1
  a=generatea2(MAX,max.dmut,no.pathway)
  cl=unlist(a$c)
  ncl=length(cl)
  al=list()
  indexa=1
  for(i in 1:length(a$a)){
    for(j in 1:nrow(a$a[[i]])){
      al[[indexa]]=a$a[[i]][j,]
      indexa=indexa+1
    }
  }
  weight=weight_cal(cl, matrix, prob,prob_nondammut,nmatrix,prob_matrix,pathway,nnonc,base,p1,p2)
  weightmatrix=matrix(weight,ncol=ncl+1,byrow=F)
  result=list("weight"=weightmatrix,"al"=al,"max.dmut"=max.dmut)
}

#' Function to caculate the probabilities of a random non-functional mutations occurs in each pathways/datasets.
#'
#' @param pathway pathway list for each mutations.
#' @param npathway Number of redefined gene sets which has no overlap with each other.
#' @param prob The list of functional impact score for each mutations.
#' @param nmatrix A vector indicates the number of reduced matrices for each patients.
#' @return A vector with probabilities of a random non-functional mutations occurs in each pathways/datasets.
probnon <- function(pathway, npathway,prob,nmatrix){
  for (i in 1:length(matrix))
    for (i in 1:length(prob))
      storage.mode(prob[[i]]) <- "double"
  storage.mode(nmatrix) <- "integer"
  for (i in 1:length(pathway))
    storage.mode(pathway[[i]]) <- "integer"
  storage.mode(npathway) <- "integer"
  .Call("probnon",pathway, npathway,prob,nmatrix)
}

#' Function for calculating the "accumulated weight" for each possible pathway/geneset mutational orders.
weight_cal<-function(cl, matrix, prob,prob_nondammut,nmatrix,prob_matrix,pathway,nnonc,base,p1,p2){
  storage.mode(cl) <- "integer"
  storage.mode(p1) <- "integer"
  storage.mode(p2) <- "integer"
  storage.mode(MAX) <- "integer"
  for (i in 1:length(matrix))
    storage.mode(matrix[[i]]) <- "integer"
  for (i in 1:length(prob))
    storage.mode(prob[[i]]) <- "double"
  for (i in 1:length(prob_nondammut))
    storage.mode(prob_nondammut[[i]]) <- "double"
  storage.mode(nmatrix) <- "integer"
  for (i in 1:length(prob_matrix))
    storage.mode(prob_matrix[[i]]) <- "double"
  for (i in 1:length(pathway))
    storage.mode(pathway[[i]]) <- "integer"
  storage.mode(nnonc) <- "integer"
  storage.mode(base) <- "integer"
  
  .Call("weight_cal", cl, matrix, MAX,prob,prob_nondammut,nmatrix,prob_matrix,pathway,nnonc,base,p1,p2)
}


#' Function for generating all permutations with MAX.mut mutations and calculate the 

generatea2<-function(MAX,MAX.mut,no.pathway){
  base=MAX.mut+1
  a <- vector("list",MAX.mut)
  b <- vector("list",MAX.mut)
  c <- vector("list",MAX.mut)
  for(N in 1:min(MAX,MAX.mut)){
    a[[N]]=permutations(n=no.pathway,r=N,repeats.allowed = T)
    b[[N]]=a[[N]]
    c[[N]]=sapply(1:nrow(b[[N]]),function(i) order2int(b[[N]][i,],MAX.mut+1))
  }
  if(MAX.mut>MAX){
    for(N in (MAX+1):MAX.mut){
      y=combinations(n=no.pathway,r=N-MAX,repeats.allowed = T)
      #      a[[N]]=matrix(NULL,nrow=nrow(y)*nrow(a[[MAX]]),ncol=N)
      y2=NULL
      y3=NULL
      y4=NULL#matrix(NULL,nrow=nrow(y)*nrow(a[[MAX]]),ncol=N)
      for(i in 1:nrow(y)){
        y2=rbind(y2,a[[MAX]])
        y3=rbind(y3,matrix(rep(y[i,],nrow(a[[MAX]])),ncol=ncol(y)))
        y4t=rep(NA,no.pathway)
        for(j in 1:no.pathway){
          y4t[j]=sum(y[i,]==j)
        }
        #print("y4t")
        #print(y4t)
        y4=rbind(y4,t(matrix(rep(y4t,nrow(a[[MAX]])),nrow=length(y4t))))
        #print("y4")
        #print(y4)
      }
      a[[N]]=cbind(y2,y3)
      #      b[[N]]=matrix(NULL,nrow=nrow(y)*nrow(a[[MAX]]),ncol=N)
      b[[N]]=cbind(y2,y4)
      c[[N]]=sapply(1:nrow(b[[N]]),function(i) order2int(b[[N]][i,],MAX.mut+1))
    }
  }
  result=list("a"=a,"b"=b,"c"=c)
  return(result)
}



loglik <- function( xx, AllPtrue, division, sel,weight,al)
{
    storage.mode(xx) <- "double"
    storage.mode(AllPtrue) <- "double"
    for (i in 1:length(division))
        storage.mode(division[[i]]) <- "integer"
    storage.mode(sel) <- "integer"
    storage.mode(MAX) <- "integer"
    divsel=division[[sel]]
    storage.mode(divsel) <- "integer"
    for (i in 1:length(al))
       storage.mode(al[[i]]) <- "integer"
    storage.mode(weight) <- "double"
    
    .Call("loglik", xx, AllPtrue, divsel, weight,al,MAX)
}


order2int<-function(order,base){
  n=length(order)
  result=0
  for(i in 0:(n-1)){
    result=result+order[i+1]*base^i
  }
  return(result)
}



