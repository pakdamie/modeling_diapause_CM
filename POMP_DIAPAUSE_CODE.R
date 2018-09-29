library(here)
library(pomp)

###This loads the covariate table
load(here("Data","covar.RData"))

###This sets the number of subcompartments, or the number within each array 
n = 8                          
nE <- as.integer(n)
nL1 <-  as.integer(n)
nL2 <-  as.integer(n)
nL3 <-  as.integer(n)
nL4 <-  as.integer(n)
nL5 <-  as.integer(n)
nP <-  as.integer(n)
nAr <-  as.integer(n)
nDL4 <-  as.integer(n)
nDL5  <-  as.integer(n)

globs <- c(paste0("static int nE = ",nE,";"),
           paste0("static int nL1 = ",nL1,";"),
           paste0("static int nL2 = ",nL2,";"),
           paste0("static int nL3 = ",nL3,";"),
           paste0("static int nL4 = ",nL4,";"),
           paste0("static int nL5 = ",nL5,";"),
           paste0("static int nDL4 = ",nDL4,";"),
           paste0("static int nDL5 = ",nDL5,";"),
           paste0("static int nP = ",nP,";"),
           paste0("static int nAr = ",nAr,";"))


DIA_SKEL <- Csnippet("
                      //I need this function for the birth rate and diapause rate
                      // example- the number of eggs is dependent on all the reproductive adults 
                     // in the subcompartments 



                      const double  sum(const double arr[], int n) 
                      { 
                   double  sum = 0.0;  
                         for (int i = 0; i < n; i++) 
                        sum += arr[i]; 
  
                           return sum; 
                      } 

                  // Here are the states that will require pointers to form 
                 // an array of vectors. For all stages except for senescent adult

                    double  *E = &E1;
                    double *DE = &DE1;
                    
                    double *L1 = &L11;
                    double *DL1 = &DL11;
                    
                    double *L2 = &L21;
                    double *DL2 = &DL21;
                    
                    double *L3 = &L31;
                    double *DL3 = &DL31;
                    
                    double  *L4= &L41;
                    double *DL4 = &DL41;
                    
                    double *L5= &L51;
                    double *DL5 = &DL51;
                    
                    double *DIAL4= &DIAL41;
                    double *DDIAL4 = &DDIAL41;
                    
                    double *DIAL5= &DIAL51;
                    double *DDIAL5 = &DDIAL51;
                    
                    double *P= &P1;
                    double *DP = &DP1;
                    
                    double *Ar= &Ar1;
                    double *DAr = &DAr1;

                     int i; 
                     
                  //Birth Rate
                    long  double B;
                   B =  10.52 * exp(-0.50* powl ((TT- 24.629), 2.0) / powl (3.427 , 2.0));

                   //Development Rate
                   double alphaE, alphaL1, alphaL2, alphaL3, 
                   alphaL4, alphaL5, alphaP, alphaA, alphaDL;
                     
                     
                     alphaE  =  0.1927/ (1+ exp(-0.3039*(TT - 18.5929)));
                     alphaL1  =  0.30/(1+ exp(-0.327*(TT - 17.60  )));
                     alphaL2 =  0.457/(1+ exp(- 0.2301 * (TT - 21.23)));
                     alphaL3  =  0.338/(1+ exp(- 0.350*(TT - 18.45)));
                     alphaL4 = 0.305/(1+ exp(- 0.429*(TT - 21.54)));
                     alphaL5 =  0.33/(1+ exp(- 0.2144*(TT - 20.94)));
                     alphaP  = 0.09287/(1+ exp(-0.28966*(TT - 18.48736 )));
                     alphaA =    0.1707/(1+ exp(-0.1830*(TT - 20.2206 )));
                     alphaDL =  0.7543 / (1+ exp(0.2392 *(TT + 7.1729)));
                       
                     
                      //Mortality Rate
                    double muE, muL, muP, muA, muDL;
                   
                     muE = 0.00056 * pow(TT, 2)  - (0.038* TT) + 0.65;
                  
                     
                     muL =  0.00067 * pow(TT, 2)  - (0.017 * TT) + 0.04;
             
                     
                     muP =  -1.598e-05 * pow(TT, 2)  - (6.787e-05  * TT) + 1.606e-02;
                    
                     muA =3.328e-04 * pow(TT, 2)  - ( 1.179e-03 * TT) +  2.452e-02 ;
                     
                     muDL =3.328e-05 * pow(TT, 2)  - ( 1.179e-03 * TT) + 2.452e-02;
             
                     
                     //Diapause induction
                     
                      double DIA_INDUC_L4;
                      double DIA_INDUC_L5;
                     
                     DIA_INDUC_L4=0.2130/( 1+exp(80 *PP_DIFF));
                     DIA_INDUC_L5 = 0.5130/( 1+exp(80 *PP_DIFF));
                     
                     
                     //BALANACE THE EQUATION
                     
                     
                     //EGG
                     DE[0] = B*sum(Ar,nAr)  - (nE * alphaE* E[0]) - (muE* E[0]);
                     for ( i = 1; i  < nE;  i ++){
                     DE[i] = nE * alphaE * (E[i-1]- E[i]) - muE * E[i];}
                     
                     // LARVAE 1
                     DL1[0] = nE * alphaE * E[nE-1] - (nL1 * alphaL1  * L1[0])-  (muL*L1[0]);
                     for ( i = 1;  i < nL1; i++){
                     DL1[i] = nL1 * alphaL1 * (L1[i-1]- L1[i])- muL * L1[i];}
                     
                     //LARVAE 2
                     DL2[0] =  nL1 * alphaL1 * L1[nL1-1]  - (nL2 * alphaL2 * L2[0])- (muL*L2[0]);
                     for (i = 1;  i < nL2 ;  i++){
                     DL2[i] = nL2 * alphaL2 * (L2[i-1]- L2[i])- muL*L2[i];}
                     
                     //LARVAE 3
                     DL3[0] =  nL2 * alphaL2 * L2[nL2-1] -  (nL3 * alphaL3 * L3[0])- (muL*L3[0]);
                     for (i = 1; i < nL3; i++){
                     DL3[i]= nL3 * alphaL3 *(L3[i-1]- L3[i])- muL* L3[i];
                     }
                     
                     //LARVAE 4
                     DL4[0] =   (nL3 * alphaL3 * L3[nL3-1]) +
                     (nDL4 * alphaDL * DIAL4[nDL4-1]) - 
                     (nL4 * alphaL4 * L4[0]) -
                     (DIA_INDUC_L4 * L4[0])-
                     (muL * L4[0]);

                     for ( i = 1;  i <nL4; i++){
                     DL4[i]=  nL4 * alphaL4 * (L4[i-1] - L4[i]) - (DIA_INDUC_L4*L4[i]) - 
                     (muL* L4[i]);
                     }
                     
                     //LARVAE 5
                     DL5[0] = (nL4 * alphaL4 * L4[nL4-1]) +
                     (nDL5 * alphaDL * DIAL5[nDL5-1]) -           
                     (nL5 * alphaL5 * L5[0])  -
                      (DIA_INDUC_L5 * L5[0]) - 
                     (muL*L5[0]);
                     
                     for ( i = 1;  i < nL5; i++){
                     DL5[i]=  nL5 * alphaL5 * (L5[i-1] - L5[i]) - 
                    (DIA_INDUC_L5*L5[i]) - 
                    (muL*L5[i]);
                     }
                     
                     //DIAPAUSING LARVAL STAGE 4 
                     
                     DDIAL4[0] =  DIA_INDUC_L4 * sum(L4,nL4) -
                     (nDL4 * alphaDL* DIAL4[0])-(muDL *DIAL4[0]);
                   
                      for ( i = 1;  i <nDL4; i++){
                     DDIAL4[i] =  nDL4 * alphaDL *(DIAL4[i-1] -DIAL4[i]) - muDL * DIAL4[i];
                     }
                     
                     
                     //DIAPAUSING LARVAL STAGE 5 
                     DDIAL5[0] = DIA_INDUC_L5*sum(L5,nL5)  - 
                     (nDL5 * alphaDL*DIAL5[0])-  (muDL * DIAL5[0]);

                     for ( i = 1;  i <  nDL5; i++){
                     DDIAL5[i]=  nDL5* alphaDL*(DIAL5[i-1] - DIAL5[i]) - muDL * DIAL5[i];
                     }
                     
                     //PUPAE 
                     
                     DP[0] = nL5 * alphaL5 * L5[nL5-1] - nP * alphaP *  P[0] - muP * P[0] ;
                     for ( i = 1;  i <nP; i++){
                     DP[i] =  nP * alphaP * (P[i-1] - P[i]) - muP*P[i];
                     }

                     //REPRODUCTIVE ADULT
                     
                     DAr[0] = nP * alphaP * P[nP-1] -nAr * alphaA*  Ar[0] - muA * Ar[0] ;
                     
                      for ( i = 1;  i < nAr; i++){
                     DAr[i] =  nAr * alphaA * (Ar[i-1] - Ar[i]) - muA*Ar[i];
                      }


                     //SENESCENT ADULT
                     DAs =  nAr * alphaA * Ar[nAr-1] - muA * As;
                     
                     
                     

                     ")

###This initializes the function. I wrote this is as an R code instead of a Csnippet
init <-function(params,t0,...)
  {
  c(E = rep(0, nE),
  L1 = rep(0,nL1),
  L2 = rep(0,nL2),
  L3 = rep(0,nL3),
  L4 = rep(0,nL4),
  L5 = rep(0,nL5),
  DIAL4 =rep(2,nDL4),
  DIAL5 =rep(0,nDL5),
  P = rep(0,nP),
  Ar = rep(0,nAr),
  As = 0)}




POMP_CM<- pomp(covartable[,c(1,4)],
         times="time",
         t0=1,
         globals=globs,
         initializer=init,
         covar = covartable ,
         tcovar = "time",
         skeleton = vectorfield(DIA_SKEL),
         statenames=c(sprintf("E%1d",seq_len(nE)),
                      sprintf("L1%1d",seq_len(nL1)),
                      sprintf("L2%1d",seq_len(nL2)),
                      sprintf("L3%1d",seq_len(nL3)),
                      sprintf("L4%1d",seq_len(nL4)),
                      sprintf("L5%1d",seq_len(nL5)),
                      sprintf("DIAL4%1d",seq_len(nDL4)),
                      sprintf("DIAL5%1d",seq_len(nDL5)),
                      sprintf("P%1d",seq_len(nP)),
                      sprintf("Ar%1d",seq_len(nAr)),
                      "As"))
###The warning messages shouldn't really mateter

###Make a dummy parameter, because I put all the parameter 
###values unto the extra code for debugging easier-ness
params = c('alphaTmp'=1.000)


###This runs the desolver and we integrate
TRAJ_MODEL <- trajectory(POMP_CM,params,times=seq(1,365*5),as.data.frame=TRUE)

###To get the total number of individuals within a life-stage, 
###We have to sum up the whole subcompartments 

EGG<-rowSums(TRAJ_MODEL [,1:(nE)])

LARVAE1<- rowSums(TRAJ_MODEL [,(nE+1): (nE+nL1)])

LARVAE2<- rowSums(TRAJ_MODEL[,(nE+nL1+1): (nE+nL1+nL2)])
LARVAE3<- rowSums(TRAJ_MODEL[,(nE+nL1+nL2+1): (nE+nL1+nL2+nL3)])
LARVAE4<- rowSums(TRAJ_MODEL[,(nE+nL1+nL2+nL3+1): (nE+nL1+nL2+nL3+nL4)])
LARVAE5<- rowSums(TRAJ_MODEL[,(nE+nL1+nL2+nL3+nL4+1): (nE+nL1+nL2+nL3+nL4+nL5)])

DIAL4<-  rowSums(TRAJ_MODEL[,(nE+nL1+nL2+nL3+nL4+nL5+1): (nE+nL1+nL2+nL3+nL4+nL5+nDL4)])

DIAL5<- rowSums(TRAJ_MODEL[,(nE+nL1+nL2+nL3+nL4+nL5+nDL4+1): (nE+nL1+nL2+nL3+nL4+nL5+nDL4+nDL5)])

Pupae<- rowSums(TRAJ_MODEL[,(nE+nL1+nL2+nL3+nL4+nL5+nDL4+nDL5+1): (nE+nL1+nL2+nL3+nL4+nL5+nDL4+nDL5+nP)])

R.adult<- rowSums(TRAJ_MODEL[,(nE+nL1+nL2+nL3+nL4+nL5+nDL4+nDL5+nP+1): (nE+nL1+nL2+nL3+nL4+nL5+nDL4+nDL5+nP+nAr)])


###IF YOU RUN THE DESOLVE CODE AS WELL
###PLOT THESE

###WIth the model right now, the only stages that would be changing would be
### L4, L5, DIAL4,DIAL5.and P
plot(EGG,col='red',type='l')
plot(LARVAE1,col='red')
plot(LARVAE2,col='red')
plot(LARVAE3,col='red')
plot(LARVAE4,col='red')
plot(LARVAE5,col='red')
plot(DIAL4,col='red')
plot(DIAL5,col='red')
plot(Pupae,col='red')
plot(R.adult,type='l',col='red')


