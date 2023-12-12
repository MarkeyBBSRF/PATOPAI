#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>


#define DIMS(x)							\
  INTEGER(coerceVector(getAttrib((x), R_DimSymbol), INTSXP))   \

#define e 2.718281828



int transorder2int(int *order,int nmut,int base){
    int tresult=0;
    for(int i=0;i<nmut;i++){
        tresult += order[i]*pow(base,i);
    }
    return(tresult);
}

int transorder2b(int *order,int nmut,int npathway,int MAX,int *b){
    for(int j=0;j<MAX;j++){
            b[j] = order[j];
    }
    for(int i=0;i<npathway;i++){
        b[i+MAX]=0;
        for(int j=MAX;j<nmut;j++){
            if(order[j]==(i+1)){
                b[i+MAX] += 1;
            }
        }
    }
}

/*void transtoperm(int *result,int tresult, int nmut,int basen){
    int index=0;{
    int dividend=tresult;
    for(int i=nmut-1;i>0;i--){
        int divider=pow(basen,i);
        result[index]= dividend%divider;
       dividend=dividend/divider;
    }
}
*/



double sum(double *AllP, int *resultcal, int *nAll, int nmut){
    double tmp;
    tmp=1.0;
    for(int i=0;i<nmut;i++){
        tmp *= AllP[resultcal[i]-1+i*nAll[0]];
    }
    return tmp;
}




/*double sum(double *AllP, int *resultcal, int *nAll, int *dim_result){
    double tmp,sum;
    sum=0.0;
    for(int i=0;i<dim_result[0];i++){
        tmp=1.0;
        for( int j=0;j<dim_result[1];j++){
            //printf("AllP\t%d\t",resultcal[i*dim_result[1]+j]-1+j*nAll[0]);
            tmp *= AllP[resultcal[i*dim_result[1]+j]-1+j*nAll[0]];
        }
        sum += tmp;
    }
    return sum;
}*/

int irowsum(int *x, int* nx, int row)
{
    int i, sum = 0;
    for (i = 0; i < nx[1]; i++) {
        sum += x[row + i*nx[0]];
    }
    return sum;
}

int  rowsum(int *x,int nrow,int ncol,int *sumN1,int *N2){
    for(int i=0;i<nrow;i++){
        sumN1[i]=0;
        for(int j=0;j<ncol;j++){
            sumN1[i] += x[j+i*ncol];
        }
        N2[i]=ncol-sumN1[i];
    }
}

int  partrowsum(int *x, int nx,int start){     //sum of a vector x from a start place , nx is the length of the part of x we want to calculate.
    int i,sum=0;
    for (i = start; i < start+nx; i++) {
        sum += x[i];
    }
    return sum;
}

double partsum(double *P, int *a, int *nP, int *na, int beginrow, int endrow)  //caculate part sum of multiplication of probability values correspoing to A matrix from beginrow to endrow.
{
    int i, j;
    double tmp, sum;
    sum = 0.0;
    for (i = beginrow; i < endrow; i++) {
        tmp = 1.0;
        for (j = 0; j < na[1]; j++)
            tmp *= P[a[i + j*na[0]] - 1 + j*nP[0]];
        sum += tmp;
    }
    return sum;
}

//void indexpathway(int *order,int *pathway,int nmut){
    //*order:vector of a possible mutation order; *pathway: vector indicating the corresponding pathways each mutation belongs to; nmut: # of mutations.
//    for (int i =0;i<nmut;i++){
//        order[i]=pathway[order[i]-1];
//    }
//}

void del_cal( int *order, int *del_order, int N, int *rcausal, int c, int N1, int caus, int *pathway){
    //*order:vector of a possible mutation order; *del_order: output for all possible functional/non-functional orders; N: # of mut; N1: # of functional  mut; rcausal: possible causal matrix for a tree; cuas:0 for non-functional,1 for functional.
    int cN=c*N;
    if(N1==N){//no nondam mut or dammut
        for(int i=0;i<N;i++){
            del_order[i]=pathway[order[i]-1];
        }
    }else{ //there both nondam mut and dam mut
        int k=0;
        for(int i=0;i<N;i++){
            if (rcausal[order[i]-1+cN]==caus){
                del_order[k++]=pathway[order[i]-1];
            }
        }
    }
}


/*void *cal_lik(int *order,int nmut, int *rcausal, int nunsign, double *pweight, int *pathway,int npathway, double *AllP, int *nAll,double *lik_order,int *N1,int *N2,int num_nonc){
    //   *order:vector of a possible mutation order; rcausal: possible causal matrix for a tree.; nc: # of non-causal mut; nunsign: # of mut which has functional impact score between 0 and 1.; pweight: the probabilities for each case; *pathway: vector indicating the corresponding pathways each mutation belongs to; npathway: # of pathways; AllP:the i,j-th element indicate the probabilities that the j-th mutations occurs in pathway i; nAll: dimension of AllP; lik_order: updated result of the total likelihood for a patient.
      int dim_result[2];
//        int *cal_order=(int *)malloc(N1*sizeof(int));
      double lik_temp,summation;
      for(int c=0;c<pow(2,nunsign);c++){
//        N1=partrowsum(rcausal,nmut,nmut*c);  //# of functional mut for c-th combination;
// N2=nmut-N1;

        lik_temp=pweight[c];
        if(N1[c]!=0){
        int *cal_order=(int *)malloc(N1[c]*sizeof(int));
        del_cal(order,cal_order,nmut,rcausal,c,N1[c],1,pathway);
        summation=sum(AllP, cal_order, nAll,N1[c]);
        lik_temp*=summation;
        free(cal_order);
        cal_order=NULL;
        }
        int sN2=N2[c]+num_nonc;
        if(sN2!=0){
            lik_temp *= pow(1.0/npathway,(float)sN2);
        }
        *lik_order += lik_temp;
     }
}*/

void *count_order(int *order,int nmut, int *rcausal, int nunsign, double *pweight,double prob_matrix_ij, int *pathway,int npathway,int *N1,int *N2,int num_nonc,int base,int *rcl, int ncl, double *weight,int MAX){
    //*order:vector of a possible mutation order; rcausal: possible causal matrix for a tree.; nc: # of non-causal mut; nunsign: # of mut which has functional impact score between 0 and 1.; pweight: the probabilities for each case; *pathway: vector indicating the corresponding pathways each mutation belongs to; npathway: # of pathways; AllP:the i,j-th element indicate the probabilities that the j-th mutations occurs in pathway i; nAll: dimension of AllP; lik_order: updated result of the total likelihood for a patient.
    //        int *cal_order=(int *)malloc(N1*sizeof(int));
    int tresult;
    double count_weight;
    int *tempweight=(int *)malloc(ncl*sizeof(int));
//FILE *pfile;
//pfile=fopen("output2.txt","a");
    for(int c=0;c<pow(2,nunsign);c++){
        count_weight=pweight[c]*prob_matrix_ij;
        int sN2=N2[c]+num_nonc;
        if(sN2!=0){
            count_weight *= pow(1.0/npathway,(float)sN2);
        }
        if(N1[c]!=0){
            int *cal_order=(int *)malloc(N1[c]*sizeof(int));
            del_cal(order,cal_order,nmut,rcausal,c,N1[c],1,pathway);
            
//            fprintf(pfile,"cal_order\t");
//            for(int j=0;j<N1[c];j++){
//                fprintf(pfile,"%d\t",cal_order[j]);
//            }
//            fprintf(pfile,"\n");

            
            if(N1[c]<=MAX){
                tresult=transorder2int(cal_order,N1[c],base);
            }else{
                int *b=(int *)malloc((MAX+npathway)*sizeof(int));
                transorder2b(cal_order,N1[c],npathway,MAX,b);
                tresult=transorder2int(b,MAX+npathway,base);
                
            
//            for(int j=0;j<MAX+npathway;j++){
//                fprintf(pfile,"%d\t",b[j]);
//            }
//            fprintf(pfile,"\n");
                
                free(b);
                b=NULL;
 
            }
            
            free(cal_order);
            cal_order=NULL;
            for(int i=0;i<ncl;i++){
                if(tresult==rcl[i]){
//                    weight[i] += count_weight;
                      tempweight[i]=1;
                }
            }
        }else{//no dammut
            weight[ncl] += count_weight;
        }
    }
    for(int i=0;i<ncl;i++){
      if(tempweight[i]==1){
           weight[i] += prob_matrix_ij;
      }
    }
//fclose(pfile);
}




void *tree_perm(int *rtmatrix, int *A, int nmut, int m, int *rcausal,int nunsign, double *pweight, double prob_matrix_ij,int *pathway, int npathway,int *index,int *N1,int *N2,int num_nonc,int *index2, int *rcl, int ncl, double *weight,int base,int MAX){
    //tmatrix:matrix indicate the tree structure; A: 1234...nmut  nmut: No. of damaging mutation,equals to dimension of tmatrix.
    int temp;
time_t t;
struct tm *timeinfo; 
//FILE *pfile2;
    if(m==nmut){
//            for(int i=0;i<nmut;i++){
//            result[index[0]*nmut+i]=A[i];
//            }
         (*index)++;
        count_order(A,nmut,rcausal,nunsign,pweight,prob_matrix_ij,pathway,npathway,N1,N2,num_nonc,base,rcl,ncl,weight,MAX);

//         cal_lik(A,nmut,rcausal,nunsign,pweight,pathway,npathway,AllP,nAll,lik_order,N1,N2,num_nonc);
  if((*index)%100000000==0){
(*index2)++;
//pfile2=fopen("output.txt","a");
printf("index2%d\t",(*index2));
time(&t);
timeinfo = localtime(&t);
printf("time%s\n", asctime(timeinfo));
//fclose(pfile2);
   }
    }else{
        for(int i=m;i<nmut;i++){
            int flag=1;
            for(int j=m;j<nmut;j++){
                flag *= rtmatrix[(A[j]-1)*nmut+(A[i]-1)];
            }
            if(flag==1){
                temp=A[m];
                A[m]=A[i];
                A[i]=temp;
                tree_perm(rtmatrix,A,nmut,m+1,rcausal,nunsign,pweight,prob_matrix_ij,pathway,npathway,index,N1,N2,num_nonc,index2,rcl,ncl,weight,base,MAX);
                temp=A[m];
                A[m]=A[i];
                A[i]=temp;
            }
        }
    }
}


/*void *tree_perm(int *result,int *rtmatrix, int *A, int nmut,int m,int *index){
    //tmatrix:matrix indicate the tree structure; A: 1234...nmut  nmut: No. of damaging mutation,equals to dimension of tmatrix.
    int temp;
    int basen=10;
//    int static *result;
//    int *newresult;

//    index=0；
    if(m==nmut){
//        newresult = (int *)realloc(result, nmut*sizeof(int));
//        if(newresult==NULL){
//            printf("fail\n");
//            return(result);
//        }
//        printf("newresult%p\n",newresult);
//        result=newresult;
//        printf("result%p\n",result);
 

//        printf("A");
//        printf("%d\t",index[0]);
//        for (int i=0; i<nmut; i++) {
//            printf("%d\t",A[i]);
//        }
//        printf("\n");
        for(int i=0;i<nmut;i++){
            result[index[0]]=trans(&A[i],nmut,basen);
        for(int i=0;i<nmut;i++){
            result[index[0]]=trans(&A[i],nmut,basen);
//            result[index[0]*nmut+i]=A[i];
        
//            printf("%d\t",index[0]*nmut+i);
 //           printf("%d\t",result[index[0]*nmut+i]);
       }
//        printf("\n");
        (*index)++;
//        if((*index)*nmut>1e7){
//       printf("index %d\t", *index);
//        }

    }
    else{
        for(int i=m;i<nmut;i++){
            int flag=1;
            for(int j=m;j<nmut;j++){
                flag *= rtmatrix[(A[j]-1)*nmut+(A[i]-1)];
            }
//            printf("%d\t",flag);
            if(flag==1){
                temp=A[m];
                A[m]=A[i];
                A[i]=temp;
                tree_perm(result,rtmatrix,A,nmut,m+1,index);
                temp=A[m];
                A[m]=A[i];
                A[i]=temp;
            }
        }
    }
//    printf("resultf%p\n",result);
//    return(result);
}
*/
                
/*
SEXP test(SEXP tmatrix){
 int *rtmatrix,nmut;
 double *rans;
 rtmatrix=INTEGER(tmatrix);
 nmut=DIMS(tmatrix)[0];
 int A[nmut];
 for(int i=0;i<nmut;i++){
 A[i]=i+1;
 }
// printf("%d\n",nmut);
 //for(int i=0;i<nmut;i++){
 //for(int j=0;j<nmut;j++){
 //printf("%d\t",rtmatrix[i+j*nmut]);
 //}
 //printf("\n");
 //}
 int index=0;
 int result[100];
 
for( int i1 = 0;i1<nmut;i1++){
   for( int i2 = 0;i2<nmut;i2++){
  printf("%d\t",rtmatrix[i1+nmut*i2]);
  }
  printf("\n");
  }


  tree_perm(rtmatrix,A,nmut,0,&index);
     printf("%d\n",index);
  SEXP ans;
 
  PROTECT(ans = allocVector(REALSXP, 12));
  rans = REAL(ans);
     for(int i=0;i<12;i++){
         rans[i]=result[i];
     }
 
  UNPROTECT(1);
 
  return ans;
 
 
  }

*/


void causind(int *rcausal,double *prob, int nmut ,int N) /*generate causal matrix for a sample with nmut mutaitons.nunde_mut:*/
{
    int num,j;
    num=N;
    for(j=0;j<nmut;j++){
        if(prob[j]<1e-6) {rcausal[j]=0;}
        else if(prob[j]>1-(1e-6) ) {rcausal[j]=1;}
        else{
            if(num>0){
                rcausal[j]=num % 2;
                num=num-rcausal[j];
                num /= 2;
            }else{
                rcausal[j]=0;
            }
        }
    }
}


void causind2(int *rcausal,double *prob, int nmut, int nunsign) /*generate causal matrix for a sample with nmut mutaitons. The output rcausal is a matrix that each row represents a possible combination*/
{
    for(int c=0;c<pow(2,nunsign);c++){
        int num=c;
        for(int j=0;j<nmut;j++){
            if(prob[j]<1e-6) {rcausal[j+c*nmut]=0;}
            else if(prob[j]>1-(1e-6) ) {rcausal[j+c*nmut]=1;}
            else{
                if(c>0){
                    rcausal[j+c*nmut]=num % 2;
                    num=num-rcausal[j+c*nmut];
                    num /= 2;
                    
                }else{
                    rcausal[j+c*nmut]=0;
                }
            }
        }
    }
}

double conprob( int *rcausal, double *rprob_j, int nmut) /*caculate p(n_j|m_j),rcausal is the indicate of causal vec for a situation  （pnm is the output，rprob_j is the input probability）*/
{
    int j;
    double pnm;
    pnm=1.0;
    for(j = 0; j< nmut;j++){
        if(rcausal[j]==0){
            pnm *= (1-rprob_j[j]);
        }else{
            pnm *= rprob_j[j];
        }
    }
    return(pnm);
}


void conprob2( double *pnm, int *rcausal, double *rprob_j, int nmut, int nunsign) /*caculate p(n_j|m_j),rcausal is the indicate of causal vec for a situation  （pnm is the output，rprob_j is the input probability）*/
{
    for(int c=0;c<pow(2,nunsign);c++){
        pnm[c]=1.0;
        for(int j = 0; j< nmut;j++){
            if(rcausal[j+c*nmut]==0){
                pnm[c] *= (1-rprob_j[j]);
            }else{
                pnm[c] *= rprob_j[j];
            }
        }
    }
}


/*void delnoncausal(int *result,int *resultcal, int N,int index,int *rcausal,int c){
    int indexcal=0;
    for (int j =0;j<N*index;j++){
        if (rcausal[result[j]-1+c*N]==1) {
            resultcal[indexcal++]=result[j];
        }
    }
}
 */


void pathway_cumu(int *result,int *pathway,int *pathwaymutcumu, int index,int N,int npathway){
    for (int i =0;i<index;i++){
        for( int j=0;j<N;j++){
            pathwaymutcumu[pathway[result[i+j*index]-1]-1]++; //+j*npathway
        }
    }
}


double logfactorial(int n){   //compute log(n!)
    double result=0.0;
    for(int i=1;i<n+1;i++){
        result += log(i);
    }
    return result;
}

int factorial(int n){   //compute n!
    int result=1;
    for(int i=2;i<n+1;i++){
        result *= i;
    }
    return result;
}


void count(double *prob, int nmut,int* nunsign, int* nc){//, double *loglik,int j){
    for(int i=0;i<nmut;i++){
        if((prob[i] > 1e-6) & (prob[i]< 1-(1e-6))){
            (*nunsign)++;
        }else if(prob[i]> 1-(1e-6)){
            (*nc)++;
        }
    }
}



void caustable(int *rcausal, int *rY, int *newrY_j, int i, int nY0, int nY1)/*casual mut table of sample j corresponding to causal matrix(ncausal=2^n) */
/* rcausal: causal vector; rY_j: mut vector for jth sample; nY_j: length of rY_j ; newrY_j: output;
 newY_j: first row is for damage mut and second row is for non_damage mut*/
{
    int index=0;
    for(int j=0; j < nY1; j++){
        newrY_j[j]=0;
        newrY_j[j+nY1]=0;
        for(int k=0;k<rY[i+j*nY0];k++){
            newrY_j[j] += rcausal[index];
            newrY_j[j+nY1] += 1-rcausal[index];
            index += 1;
        }
    }
}


/*double sum(double *P, int *a, int *nP, int *na)
{
    int i, j;
    double tmp, sum;
    sum = 0.0;
    for (i = 0; i < na[0]; i++) {
        tmp = 1.0;
        for (j = 0; j < na[1]; j++)
//            printf("a\t%f\t",P[a[i + j*na[0]] - 1 + j*nP[0]]);
//            printf("%d\t",i);
            tmp *= P[a[i + j*na[0]] - 1 + j*nP[0]];
        sum += tmp;
    }
    return sum;
}*/



double sumnon(double *Pnon, int *anon, int nPnon, int *nanon)
{
    int i, j;
    double tmp, sum;
    sum = 0.0;
    for (i = 0; i < nanon[0]; i++) {
        tmp = 1.0;
        for (j = 0; j < nanon[1]; j++)
            tmp *= Pnon[anon[i + j*nanon[0]] - 1];
        sum += tmp;
    }
    return sum;
}







double logsum(double *P, int *a, int *nP, int *na)
{
    int i, j;
    double tmp, sum;
    sum = 0.0;
    for (i = 0; i < na[0]; i++) {
        tmp = 1.0;
        for (j = 0; j < na[1]; j++)
            tmp *= P[a[i + j*na[0]] - 1 + j*nP[0]];
        sum += tmp;
    }
    return log(sum);
}



double dvecsum(double *x, int n)
{
    int i;
    double sum = 0.0;
    for (i = 0; i < n; i++)
        sum += x[i];
    return sum;
}





void dsubmat(double *old, double *new, int *nold, int *nnew,
             int* rowold, int* colold)
{
    int i, j;
    if (rowold != NULL && colold == NULL) {
        for (i = 0; i < nnew[0]; i++)
            for (j = 0; j < nnew[1]; j++)
                new[i + j*nnew[0]] = old[rowold[i]-1 + j*nold[0]];  //edit: -1
    } else if (rowold == NULL && colold != NULL) {
        for (i = 0; i < nnew[0]; i++)
            for (j = 0; j < nnew[1]; j++)
                new[i + j*nnew[0]] = old[i + colold[j]*nold[0]];
    } else if (rowold != NULL && colold != NULL) {
        for (i = 0; i < nnew[0]; i++)
            for (j = 0; j < nnew[1]; j++)
                new[i + j*nnew[0]] = old[rowold[i] + colold[j]*nold[0]];
    } else {
        for (i = 0; i < nnew[0]; i++)
            for (j = 0; j < nnew[1]; j++)
                new[i + j*nnew[0]] = old[i + j*nold[0]];
    }
}

void isubmat(int *old, int *new, int *nold, int *nnew,
             int* rowold, int* colold)
{
    int i, j;
    if (rowold != NULL && colold == NULL) {
        for (i = 0; i < nnew[0]; i++)
            for (j = 0; j < nnew[1]; j++)
                new[i + j*nnew[0]] = old[rowold[i] + j*nold[0]];
    } else if (rowold == NULL && colold != NULL) {
        for (i = 0; i < nnew[0]; i++)
            for (j = 0; j < nnew[1]; j++)
                new[i + j*nnew[0]] = old[i + colold[j]*nold[0]];
    } else if (rowold != NULL && colold != NULL) {
        for (i = 0; i < nnew[0]; i++)
            for (j = 0; j < nnew[1]; j++)
                new[i + j*nnew[0]] = old[rowold[i] + colold[j]*nold[0]];
    } else {
        for (i = 0; i < nnew[0]; i++)
            for (j = 0; j < nnew[1]; j++)
                new[i + j*nnew[0]] = old[i + j*nold[0]];
    }
}



SEXP weight_cal(SEXP cl, SEXP matrix, SEXP MAX, SEXP prob, SEXP prob_nondammut, SEXP nmatrix, SEXP prob_matrix, SEXP pathway, SEXP nnonc,SEXP base)
//pathway:a list has the same length as the # of sample, indicate the mutation in the matrix belong to which pathway.
//nmc: # of mc trees for each sample.
//nmatrix: a vector indicate the number of matrix for each sample
//prob_matrix: a list has the same length as the # of sample, indicate the prob for each matrix.

{
    int *nP,nsample,npathway,ncl;
    int *rMAX,*rbase, *matrix_j, *pathway_i, *num_matrix, *num_nonc, *rcl;
    int N,cumu_nmatrix;
    
    double *rans,*rprob_nondammut;
    double *weight;
    double *prob_i,*prob_matrix_i;
    //  int *newY_i;
    
    
    
    
    nP = (int *)R_alloc(2, sizeof(int));
    
    num_matrix=INTEGER(nmatrix);
    rMAX = INTEGER(MAX);
    rbase=INTEGER(base);
    nsample=length(nmatrix);
    rprob_nondammut = REAL(prob_nondammut);
    npathway=length(prob_nondammut);
    num_nonc=INTEGER(nnonc);
//    loglik = (double *)R_alloc(nsample, sizeof(double));
    ncl=length(cl);
    rcl=INTEGER(cl);

    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, (ncl+1)*nsample));
    rans = REAL(ans);

    
    weight=(double *)R_alloc(ncl+1, sizeof(double));
    
    cumu_nmatrix=0;
    
    //   for (i = 0; i < nAll[0] ; i++){
    //       printf("%f\t",AllP[i + (rdivsel[0] - 1)*nAll[0]]);
    //       }
    //  printf("\n");
    
    /*    printf("Allp");
     for(int j=0; j<nAll[0] ;j++)
     {
     for (int k = 0; k < nAll[1]; k++){
     printf("%f\t",AllP[j + k*nAll[0]]);
     }
     printf("\n");
     }
     */
    
    /*  for(j=0; j<ndivsel ;j++)
     {
     for (i = 0; i < nAll[0] - 1; i++)
     AllP[i + (rdivsel[j] - 1)*nAll[0]] = rx[i];
     AllP[i + (rdivsel[j]-1)*nAll[0]] = 1 - dvecsum(rx, nx);
     }
     */
    
    //i sample;j matrix num for each sample; c index for each causal situation
    
    //printf("nsample%d\t",nsample);

FILE *pfile;
//pfile=fopen("output2.txt","w");
//fclose(pfile);
    for (int i = 0; i < nsample; i++) {
        //printf("nmatrix%d\t",num_matrix[i]);
        prob_i=REAL(VECTOR_ELT(prob,i));
        pathway_i=INTEGER(VECTOR_ELT(pathway,i));
        N=DIMS(VECTOR_ELT(matrix,cumu_nmatrix))[1]; //# of mut
        prob_matrix_i=REAL(VECTOR_ELT(prob_matrix,i));
        
        int nc=0;
        int nunsign=0;
        for(int i=0;i<(ncl+1);i++){
            weight[i]=0.0;
        }
        count(prob_i,N,&nunsign,&nc);
        int pow2nunsign=pow(2,nunsign);
        //printf("unsign%d\n",nunsign);
        int rcausal[N*pow2nunsign];
        causind2(&rcausal[0],prob_i,N,nunsign);
        //printf("pow\t%d\n",pow2nunsign);
        //for( int k=0;k<N*pow2nunsign;k++){
        //printf("%d\t",rcausal[k]);
        //}
        //printf("\n");
        double pweight[pow2nunsign];
        conprob2(pweight,&rcausal[0],prob_i,N,nunsign);
        //rcuasal:indicate matrix for causal/noncausal;pweight: probs for each situation in rcausal.
        int N1[pow2nunsign];
        int N2[pow2nunsign];
        rowsum(rcausal,pow2nunsign,N,N1,N2);
        
  
pfile=fopen("output.txt","a");
fprintf(pfile,"i%d\t",i);
fprintf(pfile,"N%d\n",N);
        int A[N];
        for(int j=0;j<N;j++){
            A[j]=j+1;
//fprintf(pfile,"%d\t",A[j]);
        }
//fprintf(pfile,"\n");
fclose(pfile);
      
        for( int j = 0; j < num_matrix[i]; j++){
pfile=fopen("output.txt","a");
fprintf(pfile,"j%d\n",j);
fclose(pfile);
            
            matrix_j=INTEGER(VECTOR_ELT(matrix,j+cumu_nmatrix));
            int index=0;
            int index2=0;
            tree_perm(matrix_j, A,  N, 0, rcausal, nunsign, pweight,prob_matrix_i[j],pathway_i,npathway,&index,N1,N2,num_nonc[i],&index2,rcl,ncl,weight,rbase[0],rMAX[0]);
/*pfile=fopen("output2.txt","a");
fprintf(pfile,"liki_j%f\n",lik_i);
fclose(pfile);*/      
         }//j

        
//        loglik[i]=log(lik_i);
//pfile=fopen("output2.txt","a");
//fprintf(pfile,"logliki%f\n",loglik[i]);
//fclose(pfile);
        cumu_nmatrix += num_matrix[i];
        for(int j=0;j<ncl+1;j++){
            rans[(ncl+1)*i+j]=weight[j];
        }
    }//i
    
    
    UNPROTECT(1);
    
    return ans;
}

SEXP probnon(SEXP pathway, SEXP npathway, SEXP prob,SEXP nmatrix){
    int *pathway_i;
    int nsample,rnpathway;
    rnpathway=INTEGER(npathway)[0];
    double *prob_i;
    
    double *rans, *prob_nonmut;
    long double pweight_c,sum_prob;
    
    nsample=length(nmatrix);

    SEXP ans;
    PROTECT(ans = allocVector(REALSXP,rnpathway));
    rans = REAL(ans);
    
    
    
    prob_nonmut =(double *)R_alloc(rnpathway, sizeof(double));
    for (int i = 0; i < rnpathway; i++){
        prob_nonmut[i]=0.0;
    }
    
    for(int i = 0; i < nsample; i++){
        prob_i=REAL(VECTOR_ELT(prob,i));
        pathway_i=INTEGER(VECTOR_ELT(pathway,i));
        int N=length(VECTOR_ELT(prob,i));
        int nc=0;
        int nunsign=0;
        count(prob_i,N,&nunsign,&nc);
        for(int c=0;c<pow(2,nunsign);c++){
            int rcausal[N];
            causind(&rcausal[0],prob_i,N,c);
            pweight_c=conprob(&rcausal[0],prob_i,N);
            for(int k=0;k<N;k++){
                prob_nonmut[pathway_i[k]] += rcausal[k]*pweight_c;
            }
        }
    }
    sum_prob=0.0;
    
    for(int k=0;k<rnpathway;k++){
        sum_prob += prob_nonmut[k];
    }
    
    
    printf("%Lf\t",sum_prob);
    
    for(int k=0;k<rnpathway;k++){
        rans[k] = prob_nonmut[k]/sum_prob;
    }
    
    
    UNPROTECT(1);
    
    return ans;

}

SEXP loglik(SEXP x, SEXP AllPtrue, SEXP divsel,
            SEXP weight, SEXP al, SEXP MAX)
//pathway:a list has the same length as the # of sample, indicate the mutation in the matrix belong to which pathway.
//nmc: # of mc trees for each sample.
//nmatrix: a vector indicate the number of matrix for each sample
//prob_matrix: a list has the same length as the # of sample, indicate the prob for each matrix.

{
    int nx, ndivsel, *nAll, *nP,*nweight,nal_j,*ral_j;
    int *rdivsel, *rMAX;
    
    double *rx, *rAllPtrue, *rans,summation, *rweight;
    double *AllP, *loglik;
    //  int *newY_i;
    
    
    
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, 1));
    rans = REAL(ans);
    nx   = length(x);
    
    ndivsel = length(divsel);
    nAll = DIMS(AllPtrue);
    nP = (int *)R_alloc(2, sizeof(int));
    
    rx = REAL(x);
    rAllPtrue = REAL(AllPtrue);
    rdivsel = INTEGER(divsel);
    rMAX = INTEGER(MAX);
    nweight=DIMS(weight);
    rweight=REAL(weight);
    loglik = (double *)R_alloc(nweight[0], sizeof(double));

//printf("%d\t",nweight[0]);
//printf("%d\t",nweight[1]);
    
    

    
    AllP = (double *)R_alloc(nAll[0]*nAll[1], sizeof(double));
    for (int i = 0; i < nAll[0]; i++)
        for (int j = 0; j < nAll[1]; j++)
            AllP[i + j*nAll[0]] = rAllPtrue[i + j*nAll[0]];
    
    double sumP=0.0;
    
    for(int j=0; j<ndivsel ;j++)
    {
        for (int i = 0; i < nAll[0] - 1; i++){
            AllP[i + (rdivsel[j] - 1)*nAll[0]] = exp(rx[i]);
            sumP += AllP[i + (rdivsel[j] - 1)*nAll[0]] ;
        }
    }
    
    sumP /= ndivsel;
    
    
    for(int j=0; j<ndivsel ;j++)
    {
        for (int i = 0; i < nAll[0] - 1; i++){
            AllP[i + (rdivsel[j] - 1)*nAll[0]] /= sumP+1;
            
        }
        AllP[nAll[0]-1 + (rdivsel[j]-1)*nAll[0]] = 1.0/(sumP+1);
    }
    
 

   
    for(int i=0;i<nweight[0];i++){//i:sample
        loglik[i]=0.0;
        summation=0.0;
        for(int j=0;j<nweight[1]-1;j++){  //j:order
            if(rweight[i*nweight[1]+j]>1e-6){
                nal_j = length(VECTOR_ELT(al,j));
                ral_j = INTEGER(VECTOR_ELT(al,j));
                summation += sum(AllP, ral_j, nAll,nal_j)*rweight[i*nweight[1]+j];

//printf("i%d\t",i);
//printf("j%d\t",j);
//printf("%f\t\t",sum(AllP, ral_j, nAll,nal_j));
            }
       }
       summation += rweight[(i+1)*nweight[1]-1];
            

//printf("sumation%f\n",summation);
       
       loglik[i] += log(summation);
    }

    *rans = -dvecsum(loglik, nweight[0]);
    
    UNPROTECT(1);
    
    return ans;
}


/*long double factorial(long double initial,int n){   //compute initial*n!
    long double result=initial;
    for(int i=1;i<n+1;i++){
        result *= n;
            }
//    if(result <0.0){
//        printf("%s\n","overflow");
//    }

    return result;
}
*/
    
/*double sum_tree(double *AllP, int *rcausal, int *result, int index, int nmut){
    double result=0.0;
    for(int i=0;i<index;i++){
        for(int j=0;j<nmut;j++){
            if(rcausal[j]==1){
                
            }
        }
    }
}*/

 /*SEXP loglik(SEXP x, SEXP AllPtrue, SEXP divsel,
	    SEXP Y, SEXP a, SEXP MAX,SEXP prob,SEXP prob_nondammut)
    
{
  int nx, ndivsel, *nAll, *nY, *nP, *na_N;
  int *rY, *rdivsel, *rMAX, *ra_N;
  int i, j, k, l, N,N1,c;

  double *rx, *rAllPtrue, *rans,*rprob_nondammut;
  double *AllP, *loglik;
  double *prob_i;
  double pweight_c;
  //  int *newY_i;
    
  
 
  SEXP ans;

  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  nx	 = length(x);
  ndivsel = length(divsel);
  nAll = DIMS(AllPtrue);
  nY	 = DIMS(Y);
  nP = (int *)R_alloc(2, sizeof(int));

  rx = REAL(x);
  rAllPtrue = REAL(AllPtrue);
  rdivsel = INTEGER(divsel); 
  rY = INTEGER(Y);
  rMAX = INTEGER(MAX);
    rprob_nondammut = REAL(prob_nondammut);

   loglik = (double *)R_alloc(nY[0], sizeof(double));

  AllP = (double *)R_alloc(nAll[0]*nAll[1], sizeof(double));
  for (i = 0; i < nAll[0]; i++)
    for (j = 0; j < nAll[1]; j++)
      AllP[i + j*nAll[0]] = rAllPtrue[i + j*nAll[0]];
    
    double sumP=0.0;
  for(j=0; j<ndivsel ;j++)
     {
         for (i = 0; i < nAll[0] - 1; i++){
              AllP[i + (rdivsel[j] - 1)*nAll[0]] = exp(rx[i]);
              sumP += AllP[i + (rdivsel[j] - 1)*nAll[0]] ;
         }
     }
    
    sumP /= ndivsel;

    
    for(j=0; j<ndivsel ;j++)
    {
        for (i = 0; i < nAll[0] - 1; i++){
            AllP[i + (rdivsel[j] - 1)*nAll[0]] /= sumP+1;
            
        }
            AllP[nAll[0]-1 + (rdivsel[j]-1)*nAll[0]] = 1.0/(sumP+1);
    }

    
 
  for (i = 0; i < nY[0]; i++) {
      loglik[i]=0.0;
      N = irowsum(rY, nY, i);
      prob_i=REAL(VECTOR_ELT(prob,i));
      int nc=0;
      int nunsign=0;
      count(prob_i,N,&nunsign,&nc);
      
      for(c=0;c<pow(2,nunsign);c++){
          int rcausal[N];
          causind(&rcausal[0],prob_i,N,c);
          pweight_c=conprob(&rcausal[0],prob_i,N);
          
          int newY_i[nY[1]*2];
          caustable(rcausal, rY , newY_i, i, nY[0],nY[1] );
          
          N1=partrowsum(&newY_i[0],nY[1],0); //# of damage mut
          int N2=partrowsum(&newY_i[0],nY[1],nY[1]); //# of non dam mut
          
          double logprob_nonmut_sum=0.0;
          if (N2!=0){
              for(j=0;j<nY[1];j++){
                  for(k = 0; k < newY_i[nY[1]+j]; k++){
                      if(rprob_nondammut[j]!=0.0){
                          logprob_nonmut_sum += log(rprob_nondammut[j]);
                      }else{

printf("jump\t");
printf("%s\t",k);
printf("%f\n",rprob_nondammut[k]);
                          goto end;
                      }
                  }
              }
              logprob_nonmut_sum = logprob_nonmut_sum+logfactorial(N2);
              
              
           }
          
          if(N1==0){
              loglik[i] += pweight_c*pow(e,logprob_nonmut_sum);
          }else{
              nP[0] = N1;
              nP[1] = N1;
              double P[nP[0]*nP[1]];
              int rowAll[nP[0]];
              
          
              l = 0;
              for (j = 0; j < nY[1]; j++)
                  for (k = 0; k < newY_i[j]; k++) {
                      rowAll[l++] = j;
                  }
              
        
          
              dsubmat(AllP, &P[0], nAll, nP, &rowAll[0], NULL);
          
              na_N = DIMS(VECTOR_ELT(a, N1-1));
              ra_N = INTEGER(VECTOR_ELT(a, N1-1));
              
              double summation=sum(P, ra_N, nP, na_N);
              
 printf("%.9f\t",logfactorial(N1-rMAX[0])+log(summation));
  printf("%.9f\t",pweight_c);
   printf("%f\t",logprob_nonmut_sum);        
          //    double temp1=summation*pweight_c*prob_nonmut_sum;
              double logprob_c=logfactorial(N1-rMAX[0])+log(summation)+log(pweight_c)+logprob_nonmut_sum;
//printf("%s\t","temp");
//printf("%.9f\n",temp);

              double prob_c=pow(e,logprob_c);
printf("%s\t","prob_c");
printf("%.9f\n",prob_c);
     //         double temp = temp1*temp2*temp3;//sum(P, ra_N, nP, na_N)*pweight_c*pow(1.0/nY[1],N-N1)*factorial(N1-rMAX[0]);
             
              loglik[i] += prob_c;
  //            loglik[i] += sum(P, ra_N, nP, na_N)*pweight_c*pow(1.0/nY[1],N-N1)*factorial(N1-rMAX[0]);
              }
               end: ;
//printf("jump!\n");
          }

  printf("%.10f\n",loglik[i]);
    
      loglik[i]=log(loglik[i]);

      }
  
    
  *rans = -dvecsum(loglik, nY[0]);

  UNPROTECT(1);

  return ans;
}
  */

void chindex(int* seq, int nseq, int a, int omit)
{
  int i, j;

  j = 0;
  for (i = 0; i < nseq; i++) {
    if (i == omit)
      j++;
    seq[i] = j + a;
    j++;
  }
}


void makeimat(int *old, int *new, int *no, int *nn,
	      int radd, int cadd, int romit, int comit)
{
  int i, j, r, c;
  for (i = 0, r = 0; i < nn[0]; i++, r++) {
    if (i == romit)
      r++;
    for (j = 0, c = 0; j < nn[1]; j++, c++) {
      if (j == comit)
	c++;
      new[i + j*nn[0]] = old[r + radd + (c + cadd)*no[0]];
    }
  }
}

void rowreduce(double *x, int nx)
{
  int i;
  for (i = 0; i < nx; i++)
    x[i] -= x[nx-1];
}

/*SEXP predict(SEXP Yj, SEXP AllPesti, SEXP a,SEXP prob,SEXP MAX)  //function to caculate P(A<B),P(A>B),P(A=B)

{
    int nYj, *na_N,*nAll;
    int *rYj,*rMAX, *ra_N;
    int N,N1;
    
    double *rAllPesti, *rans, *rprob;
    long double pweight_c;
    
    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, 3)); //Probability of P(A=B),P(A<B),P(A>B)
    rans = REAL(ans);
    
 
    nYj	 = length(Yj);
    rYj = INTEGER(Yj);
    N=rowsum(rYj,nYj);
    nAll = DIMS(AllPesti);
    rMAX=INTEGER(MAX);
    int nP[2];
    
    rAllPesti = REAL(AllPesti);
    rprob=REAL(prob);
    
    int nunsign=0;
    int nc=0;
    count(rprob,N,&nunsign,&nc);
    
    long double jointprob[3];  //P(Y_i,j intersection A<B)
    
    for (int i=0;i<3;i++ ) {
        jointprob[i] += 0.0;
    }

    
    for(int c=0;c<pow(2,nunsign);c++){
        int rcausal[N];
        causind(&rcausal[0],rprob,N,c);
        pweight_c=conprob(&rcausal[0],rprob,N);
        int newYj[nYj];
        caustable(rcausal, rYj , newYj, 0, 1,nYj);
        
        N1=rowsum(&newYj[0],nYj);   //N1: number of causal mut corresponding to c
        
        if (N1>0){
            nP[0] = N1;
            nP[1] = N1;
            double P[nP[0]*nP[1]];
            int rowAll[nP[0]];
            
            
            printf("rporb\t");
            for (int i=0; i<N; i++) {
                printf("%f\t",rprob[i]);
            }
            printf("\n");
            
            printf("N\t%d\n",N);
            
            
            printf("pweight\t%Lf\n",pweight_c);
            
            printf("rYj\t");
            for (int i=0; i<nYj; i++) {
                printf("%d\t",rYj[i]);
            }
            printf("\n");
            
            printf("rcausal\t");
            for (int i=0; i<N; i++) {
                printf("%d\t",rcausal[i]);
            }
            printf("\n");
            
            printf("newY\t");
            for (int i=0; i<nYj; i++) {
                printf("%d\t",newYj[i]);
            }
            printf("\n");
            
            int l = 0;
            for (int j = 0; j < nYj; j++)
                for (int k = 0; k < newYj[j]; k++) {
                    rowAll[l++] = j;
                }
            
            dsubmat(rAllPesti, P, nAll, nP, rowAll, NULL);
            
            printf("rowAll\t");
            for (int i=0; i<N1; i++) {
                printf("%d\t",rowAll[i]);
            }
            printf("\n");
            
            
            na_N = DIMS(VECTOR_ELT(a, N1-1));
            ra_N = INTEGER(VECTOR_ELT(a, N1-1));
            int arows[3];
            int partnum=na_N[0]/N1;
            arows[0]=partnum*newYj[0];
            printf("arrow\t%d\t",arows[0]);
            
            for( int i =1;i<3;i++){
                arows[i]=(na_N[0]/N1)*newYj[i]+arows[i-1];
                printf("%d\t",arows[i]);
            }
            
            
            
            long double temp[3];
            for(int i=0;i<3;i++){
                temp[i]=0.0;
            }
            
            
            
            temp[0]=partsum(P, ra_N, nP, na_N,0,arows[0])*pweight_c*pow(1.0/3,N-N1);
            
            temp[1]=partsum(P, ra_N, nP, na_N,arows[0],arows[1])*pweight_c*pow(1.0/3,N-N1);
            temp[2]=partsum(P, ra_N, nP, na_N,arows[1],arows[2])*pweight_c*pow(1.0/3,N-N1);
            printf("temp\t");
            for (int i=0; i<3; i++) {
                printf("%Lf\t",temp[i]);
            }
            printf("\n");
            
             temp[0]=factorial(temp[0],N1-rMAX[0]);
             temp[1]=factorial(temp[1],N1-rMAX[0]);
             temp[2]=factorial(temp[2],N1-rMAX[0]);
             
             printf("temp\t");
             for (int i=0; i<3; i++) {
             printf("%Lf\t",temp[i]);
             }
             printf("\n\n");
             
             for (int i=0;i<3;i++ ) {
             jointprob[i] += temp[i];
             }

        }
    }
    
    for( int i=0;i<3;i++){
        rans[i]=jointprob[i]*rAllPesti[i];
    }
    
    double sumrans=dvecsum(rans,3);
    for( int i=0;i<3;i++){
        rans[i]=rans[i]/sumrans;
    }
    
//    rans[0]=1.0;
//    rans[1]=1.0;
//    rans[2]=1.0;

    UNPROTECT(1);
    
    return ans;
    
} */

/*SEXP prob(SEXP Y, SEXP prob)

{
    int *nY, *nP;
    int *rY;
    int i, k, N,c;
    
    double *rans, *prob_nonmut;
    double *prob_i;
    long double pweight_c,sum_prob;
    //  int *newY_i;
    
    
    
    nY	 = DIMS(Y);
    nP = (int *)R_alloc(2, sizeof(int));
    rY = INTEGER(Y);
    
    
    SEXP ans;
    PROTECT(ans = allocVector(REALSXP,nY[1]));
    rans = REAL(ans);


    
    prob_nonmut =(double *)R_alloc(nY[1], sizeof(double));
    for (i = 0; i < nY[1]; i++){
        prob_nonmut[i]=0.0;
    }
  
    
     for (i = 0; i < nY[0]; i++) {
        N = irowsum(rY, nY, i);
        prob_i=REAL(VECTOR_ELT(prob,i));
        int nc=0;
        int nunsign=0;
        count(prob_i,N,&nunsign,&nc);
        
        for(c=0;c<pow(2,nunsign);c++){
            int rcausal[N];
            causind(&rcausal[0],prob_i,N,c);
            pweight_c=conprob(&rcausal[0],prob_i,N);
   //         printf("%Lf\t",pweight_c);

      
         int newY_i[nY[1]*2];
         caustable(rcausal, rY , newY_i, i, nY[0],nY[1] );
            
       
            for(k=0;k<nY[1];k++){
                prob_nonmut[k] += newY_i[k]*pweight_c;
            }
        }

    }
  
 
    sum_prob=0.0;
    
    for(k=0;k<nY[1];k++){
        sum_prob += prob_nonmut[k];
    }
 
 
    printf("%Lf\t",sum_prob);

    for(k=0;k<nY[1];k++){
        rans[k] = prob_nonmut[k]/sum_prob;
    }


    UNPROTECT(1);
    
    return ans;
}
*/



/*    for(int c=0;c<pow(2,nunsign);c++){
        int rcausal[N];
        causind(&rcausal[0],rprob,N,c);
    
        
        printf("rporb\t");
        for (int i=0; i<N; i++) {
            printf("%f\t",rprob[i]);
        }
        printf("\n");
        
        printf("N\t%d\n",N);

 
        printf("pweight\t%Lf\n",pweight_c);
        
        printf("rYj\t");
        for (int i=0; i<nYj; i++) {
            printf("%d\t",rYj[i]);
        }
        printf("\n");

        printf("rcausal\t");
        for (int i=0; i<N; i++) {
            printf("%d\t",rcausal[i]);
        }
        printf("\n");

        for (int i=0; i<nYj; i++) {
            printf("%d\t",newYj[i]);
        }
        
        N1=rowsum(&newYj[0],nYj);   //N1: number of causal mut corresponding to c
        nP[0] = N1;
        nP[1] = N1;
        double P[nP[0]*nP[1]];
        int rowAll[nP[0]];
        
        int l = 0;
        for (int j = 0; j < nYj; j++)
            for (int k = 0; k < newYj[j]; k++) {
                rowAll[l++] = j;
            }
        
        dsubmat(rAllPesti, P, nAll, nP, rowAll, NULL);
        
        na_N = DIMS(VECTOR_ELT(a, N1-1));
        ra_N = INTEGER(VECTOR_ELT(a, N1-1));
        int arows[3];
        for( int i =0;i<3;i++)
            arows[i]=na_N[i]/N1*newYj[i];
        
        long double temp[3];
        
        temp[0]=partsum(P, ra_N, nP, na_N,1,arows[0])*pweight_c*pow(1.0/3,N-N1);
        temp[0]=factorial(temp[0],N1-rMAX[0]);
        temp[1]=partsum(P, ra_N, nP, na_N,arows[0],arows[1])*pweight_c*pow(1.0/3,N-N1);
        temp[1]=factorial(temp[1],N1-rMAX[0]);
        temp[2]=partsum(P, ra_N, nP, na_N,arows[1],arows[2])*pweight_c*pow(1.0/3,N-N1);
        temp[2]=factorial(temp[2],N1-rMAX[0]);
            
        double jointprob[3];  //P(Y_i,j intersection A<B)
            
        for (int i=0;i<3;i++ ) {
                jointprob[i] += temp[i];
        }
        
        
        for( int i=0;i<3;i++){
            rans[i]=rans[i]*rAllPesti[i];
        }
        for( int i=0;i<3;i++){
            rans[i]=dvecsum(rans,3);
        }

        
        
        
   }


  
    
 //   rans[0]=1.0;
 //   rans[1]=1.0;
 //   rans[2]=1.0;
    UNPROTECT(1);
    
    return ans;

    
} */




/*
SEXP gradient(SEXP x, SEXP AllPtrue, SEXP divsel,
	      SEXP Y, SEXP a, SEXP MAX, SEXP prob,SEXP aa)
{
  int nx, ndivsel, nAll[2], nY[2], nP[2], nP1[2], na_N[2], na_N1[2], ngrad,N1;
  int *rY, *rdivsel, *rMAX, *ra_N, *ra_N1;
  int i, j, k, l, cnt, N;

  double *rx, *rAllPtrue, *rans;
  double *AllP, *lik, *grad,*prob_i,pweight_c;
 
  SEXP ans;

  nx = length(x);
  ndivsel=length(divsel);
  nAll[0] = DIMS(AllPtrue)[0]; 
  nAll[1] = DIMS(AllPtrue)[1];
  nY[0] = DIMS(Y)[0]; 
  nY[1] = DIMS(Y)[1];

  rx = REAL(x);
  rAllPtrue = REAL(AllPtrue);
  rdivsel = INTEGER(divsel);
  rY = INTEGER(Y);
  rMAX = INTEGER(MAX);

  lik = (double *)R_alloc(nY[0], sizeof(double));

  AllP = (double *)R_alloc(nAll[0]*nAll[1], sizeof(double));
  for (i = 0; i < nAll[0]; i++)
    for (j = 0; j < nAll[1]; j++)
      AllP[i + j*nAll[0]] = rAllPtrue[i + j*nAll[0]];
 

  for(j=0; j<ndivsel ;j++)
    {
      for (i = 0; i < nAll[0] - 1; i++)
	AllP[i + (rdivsel[j] - 1)*nAll[0]] = rx[i];
      AllP[i + (rdivsel[j]-1)*nAll[0]] = 1 - dvecsum(rx, nx);
    }

  ngrad = nY[1];
  grad = (double *)R_alloc(ngrad, sizeof(double));
  for (i = 0; i < ngrad; i++)
    grad[i] = 0.0;

  PROTECT(ans = allocVector(REALSXP, ngrad - 1));
  rans = REAL(ans);
    


  for (i = 0; i < nY[0]; i++) {
      
      lik[i]=0.0;
      N = irowsum(rY, nY, i);
      prob_i=REAL(VECTOR_ELT(prob,i));
      int nc=0;
      int nunsign=0;
      count(prob_i,N,&nunsign,&nc);
      
      double grad_i[ngrad];
    //  grad_i=(double *)R_alloc(ngrad, sizeof(double));
      for (int m = 0; m < ngrad; m++)
          grad_i[m] = 0.0;
      
      for(int c=0;c<pow(2,nunsign);c++){
          int rcausal[N];
          causind(&rcausal[0],prob_i,N,c);
          pweight_c=conprob(&rcausal[0],prob_i,N);
          
          int newY_i[nY[1]];
          caustable(rcausal, rY , newY_i, i, nY[0],nY[1] );
          
          N1=rowsum(&newY_i[0],nY[1]);
          if(N1==0){
              lik[i] += pow(1.0/nY[1],N)*pweight_c;
          }else{
              nP[0] = N1;
              nP[1] = N1;
              double P[nP[0]*nP[1]];
              int rowAll[nP[0]];
              
              l = 0;
              for (j = 0; j < nY[1]; j++)
                  for (k = 0; k < newY_i[j]; k++) {
                      rowAll[l++] = j;
                  }
              
              dsubmat(AllP, &P[0], nAll, nP, &rowAll[0], NULL);
              
              na_N[0] = DIMS(VECTOR_ELT(a, N1-1))[0];
              na_N[1] = DIMS(VECTOR_ELT(a, N1-1))[1];
              ra_N = INTEGER(VECTOR_ELT(a, N1-1));
              lik[i] += sum(P, ra_N, nP, na_N)*pweight_c*pow(1.0/nY[1],N-N1)*factorial(N1-rMAX[0]);


    if (rdivsel[0] <= N1) {

      if (N1 > 1) {


		
	nP1[0] = nP[0] - 1;
	nP1[1] = nP[1] - 1;
	//P1 = (double *)R_alloc(nP1[0]*nP1[1], sizeof(double));
	//rowP = (int *)R_alloc(nP1[0], sizeof(int));
	//colP = (int *)R_alloc(nP1[1], sizeof(int));
          double P1[nP1[0]*nP1[1]];
          int rowP[nP1[0]];
          int colP[nP1[0]];

	cnt = 0;
	for (j = 0; j < nP[0]; j++) {
	  cnt++;

	  for(l=rdivsel[0];l<=imin2(N,rdivsel[ndivsel-1]);l++)
	    {
	      chindex(rowP, nP1[0], 0, cnt - 1);
	      chindex(colP, nP1[1], 0, l - 1);
	      dsubmat(P, P1, nP, nP1, rowP, colP);

	      if (l>rMAX[0]) {

		na_N1[0] = DIMS(VECTOR_ELT(a, N1-2))[0];
		na_N1[1] = DIMS(VECTOR_ELT(a, N1-2))[1];
		ra_N1 = INTEGER(VECTOR_ELT(a, N1-2));
		grad_i[rowAll[j]] +=
		  //exp(logsum(P1, ra_N1, nP1, na_N1) - loglik[i])/(N-rMAX[0]);
              sum(P1, ra_N1, nP1, na_N1)*pweight_c*pow(1.0/nY[1],N-N1)*factorial(N1-rMAX[0]);

	      } else
		{
		  na_N1[0] = DIMS(VECTOR_ELT(aa, N1-2))[0];
		  na_N1[1] = DIMS(VECTOR_ELT(aa, N1-2))[1];
		  ra_N1 = INTEGER(VECTOR_ELT(aa, N1-2));
		  grad_i[rowAll[j]] +=
            sum(P1, ra_N1, nP1, na_N1)*pweight_c*pow(1.0/nY[1],N-N1)*factorial(N1-rMAX[0]) ;

		}  
	    }
	}
	rowreduce(grad_i, ngrad);

	  if (rY[i + j*nY[0]] == 1)
	    grad[j] += (1./P[0])*pweight_c*pow(1.0/nY[1],N-N1);

	rowreduce(grad_i, ngrad);
      }
    }
    
    }
    }
    for (j=0;j< ngrad -1;j++){
        grad_i[j] /=lik[i];
        grad[j] += grad_i[j];
    }
    
  }
    
    for (i = 0; i < ngrad - 1; i++)
    rans[i] = -grad[i];
    
  UNPROTECT(1);
  return ans;
}

*/

