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

globs <- c(paste0("int nE = ",nE,";"),
           paste0("int nL1 = ",nL1,";"),
           paste0("int nL2 = ",nL2,";"),
           paste0("int nL3 = ",nL3,";"),
           paste0("int nL4 = ",nL4,";"),
           paste0("int nL5 = ",nL5,";"),
           paste0("int nDL4 = ",nDL4,";"),
           paste0("int nDL5 = ",nDL5,";"),
           paste0("int nP = ",nP,";"),
           paste0("int nAr = ",nAr,";"))


DIA_SKEL <- Csnippet("
                      //I need this function for the birth rate and diapause rate
                      // example- the number of eggs is dependent on all the reproductive adults 
                     // in the subcompartments 

                      const double  sum(const double arr[], int n) 
                      { 
                      long double  sum = 0.0;  
                         for (int i = 0; i < n; i++) 
                        sum += arr[i]; 
  
                           return sum; 
                      } 

                  // Here are the states that will require pointers to form 
                 // an array of vectors 

                  const double  *E = &E1;
                  double *DE = &DE1;
                  
                  const double *L1 = &L11;
                  double *DL1 = &DL11;
                  
                  const double *L2 = &L21;
                  double *DL2 = &DL21;
                  
                  const double *L3 = &L31;
                  double *DL3 = &DL31;
                  
                  const  double  *L4= &L41;
                  double *DL4 = &DL41;
                  
                  const double *L5= &L51;
                  double *DL5 = &DL51;
                  
                  const double *DIAL4= &DIAL41;
                  double *DDIAL4 = &DDIAL41;
                  
                  const  double *DIAL5= &DIAL51;
                  double *DDIAL5 = &DDIAL51;
                  
                  const double *P= &P1;
                  double *DP = &DP1;
                  
                  const double *Ar= &Ar1;
                  double *DAr = &DAr1;
                     
                  int i; 
                     
                  //Birth Rate
                  double B;
                   B =  10.52 * exp(-0.50* pow ((TT- 24.629), 2.0) / pow (3.427 , 2.0));

                   //Development Rate
                   double alphaE, alphaL1, alphaL2, alphaL3, 
                   alphaL4, alphaL5, alphaP, alphaA, alphaDL;
                     
                     
                   alphaE  = 0;
                   alphaL1  = 0;
                   alphaL2 =  0;
                   alphaL3  =  0;
                   alphaL4 = 0.305/(1+ exp(- 0.429*(TT - 21.54)));
                   alphaL5 =  0.33/(1+ exp(- 0.2144*(TT - 20.94)));
                   alphaP  = 0;
                   alphaA =   0;
                   alphaDL =   0.7543 / (1+ exp(0.2392 *(TT + 7.1729)));
                       
                     
                      //Mortality Rate
                     double muE, muL, muP, muA, muDL;
                     
                     muE =0.0;
                     if(muE < 0.0){
                     muE = 0;
                     }else if (muE > 1){
                     muE = 1.0;
                     }else (muE = muE);
                     
                     muL = 0.0 ;
                     if(muL< 0){
                     muL = 0.0;
                     }else if (muL > 1){
                     muL = 1.0;
                     }else (muL = muL);
                     
                     muP =0.0 ;
                     if(muP < 0.0){
                     muP = 0.0;
                     }else if (muP > 1){
                     muP = 1.0;
                     }else (muP = muP);
                     
                     muA = 0.0  ;
                     if(muA < 0.0){
                     muA = 0;
                     }else if (muA > 1.0){
                     muA = 1.0;
                     }else (muA = muA);
                     
                     muDL =0.0 ;
                     if(muDL< 0.0){
                     muDL = 0.0;
                     }else if (muDL > 1.0){
                     muP = 1.0;
                     }else (muDL = muDL);
                     
                     
                     //Diapause induction
                     
                     double DIA_INDUC_L4;
                     double DIA_INDUC_L5;
                     
                    DIA_INDUC_L4=0.2130/( 1+exp(80 *PP_DIFF));
                    DIA_INDUC_L5 = 0.5130/( 1+exp(80 *PP_DIFF));
                     
                     
                     //BALANACE THE EQUATION
                     
                     
                     //EGG
                     DE[0] = B*sum(Ar,nAr)  - (nE * alphaE + muE) * E[0];
                     for ( i = 1; i  < nE;  i ++){
                     DE[i] = nE * alphaE * (E[i-1]- E[i]) - muE * E[i];}
                     
                     // LARVAE 1
                     DL1[0] = nE * alphaE * E[nE-1] - (nL1 * alphaL1 +  muL) * L1[0];
                     for ( i = 1;  i < nL1; i++){
                     DL1[i] = nL1 * alphaL1 * (L1[i-1]- L1[i])- muL * L1[i];}
                     
                     //LARVAE 2
                     DL2[0] =  nL1 * alphaL1 * L1[nL1-1]  - (nL2 * alphaL2 +  muL) * L2[0];
                     for (i = 1;  i < nL2 ;  i++){
                     DL2[i] = nL2 * alphaL2 * (L2[i-1]- L2[i])- muL*L2[i];}
                     
                     //LARVAE 3
                     DL3[0] =  nL2 * alphaL2 * L2[nL2-1] - (nL3 * alphaL3 + muL) * L3[0];
                     for (i = 1; i < nL3; i++){
                     DL3[i]= nL3 * alphaL3 *(L3[i-1]- L3[i])- muL* L3[i];
                     }
                     
                     //LARVAE 4
                     DL4[0] =   nL3 * alphaL3 * L3[nL3-1] +
                     nDL4 * alphaDL * DIAL4[nDL4-1] - 
                     nL4 * alphaL4 * L4[0] -
                     DIA_INDUC_L4 * L4[0]-
                     (muL * L4[0]);

                     for ( i = 1;  i <nL4; i++){
                     DL4[i]=  nL4 * alphaL4 * (L4[i-1] - L4[i]) - (DIA_INDUC_L4*L4[i]) - 
                     (muL* L4[i]);
                     }
                     
                     //LARVAE 5
                     DL5[0] = nL4 * alphaL4 * L4[nL4-1] +
                     nDL5 * alphaDL * DIAL5[nDL5-1] -           
                     nL5 * alphaL5 * L5[0]  -
                      DIA_INDUC_L5 * L5[0] - 
                     (muL*L5[0]);
                     
                     for ( i = 1;  i < nL5; i++){
                     DL5[i]=  nL5 * alphaL5 * (L5[i-1] - L5[i]) - 
                    (DIA_INDUC_L5*L5[i]) - 
                    (muL*L5[i]);
                     }
                     
                     //DIAPAUSING LARVAL STAGE 4 
                     
                     DDIAL4[0] =  DIA_INDUC_L4 * sum(L4,nL4) -
                     (nDL4 * alphaDL + muDL) * DIAL4[0];
                   
                      for ( i = 1;  i <nDL4; i++){
                     DDIAL4[i] =  nDL4 * alphaDL *(DIAL4[i-1] -DIAL4[i]) - muDL *DIAL4[i];
                     }
                     
                     
                     //DIAPAUSING LARVAL STAGE 5 
                     DDIAL5[0] = DIA_INDUC_L5*sum(L5,nL5)  - 
                     (nDL5 * alphaDL +muDL) * DIAL5[0];

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
  DIAL4 =rep(1,nDL4),
  DIAL5 =rep(1,nDL5),
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
TRAJ_MODEL <- trajectory(POMP_CM,params,times=seq(1,365),as.data.frame=TRUE)

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

plot(L4,type='l') ###R code simulation
lines(LARVAE4, col='red')
plot(diff(LARVAE4-L4),type='l')

plot(L5,type='l') ###R code simulation
lines(LARVAE5, col='red')
plot(diff(LARVAE5-L5),type='l')

plot(DL4,type='l') ###R code simulation
lines(DIAL4, col='red')
plot(diff(DIAL4-DL4),type='l')

plot(DL5,type='l') ###R code simulation
lines(DIAL5, col='red')
plot(diff(DIAL5-DL5),type='l')

plot(P,type='l') ###R code simulation
lines(Pupae, col='red')
plot(diff(Pupae-P),type='l')
