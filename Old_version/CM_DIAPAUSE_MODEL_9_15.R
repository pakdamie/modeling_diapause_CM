
library(deSolve)
library(geosphere)

PHOTO_PERIOD_GEOSPHERE <- geosphere::daylength(40, 1:366)
DIFF <- diff(PHOTO_PERIOD_GEOSPHERE)

MODEL_DDE_CM_DIA = function(t, y, parms) {
  
  
  ###This is for controlling the erlang distribution in each stage.
  
  nE=8#The number of subcompartments in the egg-stage
  nL1 =8 #The number of subcompartments in larval instar 1
  nL2 =8 #The number of subcompartments in larval instar 2
  nL3 = 8#The number of subcompartments in larval instar 3
  nL4 =8#The number of subcomppartmetns in larval instar 4
  nL5=8#The number of subcompartments in larval instar 5
  nDL4 =8  #The number of subcompartments in diapausing larval instar 4
  nDL5 =8 #The number of subcomparments in diapausing larval instar 5
  nP =8 # The number of subcompartments in pupae
  nA.r =8 #The number of subcompartments in reproductive adults
  
  
  ###The environmental variable- this is the daily mean temperature 
  ###from Biglerville, PA. It contains 33 years of data 
  ### (January 1st, 1984- December 31, 2016)
  TT =(SUBSET_TEMP)$avg[t]
  
  
  ###The parameters that I am interested in 
  
  ###Development
  ###The parameter associated with the "diapause termination rate"
  ###for individuals in the diapausing instar 4 and diapausing instar 5
  alphaDL_a = parms[1]; alphaDL_b = parms[2];  alphaDL_c = parms[3];
  
  ###Mortality
  ###The parameters associated with the mortality of the Adults
  deltaAa= parms[4];deltaAb= parms[5];deltaAc= parms[6];
  ###The parameters associated with the mortality of the Diapausing Larvae
  deltaDLa= parms[7]; deltaDLb= parms[8]; deltaDLc = parms[9];
  
  ###The parameter associated with diapause induction:
  ###Diapause induction rate 1 is for larvae 4
  ###Diapause induction rate 2 is for larvae 5
  DIA_INDUC1_A= parms[10];DIA_INDUC1_B = parms[11];
  DIA_INDUC2_A=parms[12];DIA_INDUC2_B = parms[13]
  
  #########
  #STAGES##
  #########
  ###The stages
  E = y[(1:nE)]
  L1 = y[(nE+1):(nE + nL1)]
  L2=  y[(nE+nL1+1):(nE+nL1 + nL2 )]
  L3 = y[(nE+nL1+nL2+1):(nE+nL1+nL2+nL3)]
  L4 =  y[(nE+nL1+nL2+nL3+ 1):(nE+nL1+nL2+nL3+nL4)]
  L5 =  y[(nE+nL1+nL2+nL3+nL4 +1):(nE+nL1+nL2+nL3+nL4+nL5)]
  DL4 = y[(nE+nL1+nL2+nL3+nL4+nL5+1):(nE+nL1+nL2+nL3+nL4+nL5+nDL4)]
  DL5 = y[(nE+nL1+nL2+nL3+nL4+nL5+nDL4+ 1):(nE+nL1+nL2+nL3+nL4+nL5+nDL4+nDL5)]
  P =  y[(nE+nL1+nL2+nL3+nL4+nL5+nDL4+nDL5+1):(nE+nL1+nL2+nL3+nL4+nL5+nDL4+nDL5+nP)]
  A.r= y[(nE+nL1+nL2+nL3+nL4+nL5+nDL4+nDL5+nP+1):(nE+nL1+nL2+nL3+nL4+nL5+nDL4+nDL5+nP+nA.r) ]
  A.s = y[(nE+nL1+nL2+nL3+nL4+nL5+nDL4+nDL5+nP+nA.r+1) ]
  
  
  #Birth Rate
  B =0* exp(-1 / 2 * (TT - 24.629) ^ 2 / ( 3.427  ^ 2))
  
  #Development rate (including the diapause termination rate)
  alphaE <- 0.1927/(1+ exp(-0.3039*(TT - 18.5929)))
  alphaL1 <- 0.30/(1+ exp(- 0.327*(TT - 17.60  )))
  alphaL2 <- 0.457/(1+ exp(- 0.2301*(TT - 21.23)))
  alphaL3 <- 0.338/(1+ exp(- 0.350*(TT - 18.45)))
  alphaL4 <- 0.305/(1+ exp(- 0.429*(TT - 21.54)))
  alphaL5 <- 0.33/(1+ exp(- 0.2144*(TT - 20.94)))
  alphaP <- 0.09287/(1+ exp(-0.28966 * (TT - 18.48736 )))
  alphaA <-   0.1707 /(1+ exp(-0.1830*(TT - 20.2206)))
  alphaDL <-   alphaDL_a /(1+ exp(alphaDL_b*(TT + alphaDL_c )))
  
  #Mortality Function (all quadratic, if it's over 1, set to 1. If it's negative,
  ### set to 0)
  deltaE <- 0
  deltaE[deltaE <0 ]<-0
  deltaE[deltaE > 1] <- 1
  
  deltaL = 0
  deltaL[deltaL <0 ]<-0
  deltaL[deltaL > 1] <- 1
  
  deltaP =0
  deltaP[deltaP <0] <- 0
  deltaP[deltaP > 1] <- 1
  
  deltaA = 0
  deltaA[deltaA <0] <- 0
  deltaA[deltaA > 1] <- 1
  
  deltaDL <-   0
  deltaDL[deltaDL <0] <- 0
  deltaDL[deltaDL > 1] <- 1
  
  ###DIAPAUSE INDUCTION RATES- I used the change in photoperiod to get
  ###the shape I want.
  
  DIA_INDUC1 <-0.5
  DIA_INDUC2 <-  0
  
  
  ######
  #Egg##
  ######
  dE= rep(0,nE)
  ###For the first substage 
  dE[1] = B * sum(A.r) - (nE*alphaE + deltaE)*E[1] 
  for (i in seq(2, nE)){
    dE[i] = nE*alphaE*(E[i-1]- E[i]) - deltaE*E[i]
  } 
  ###LARVAL INSTAR 1
  dL1= rep(0,nL1)
  dL1[1] = nE*alphaE*E[nE] - (nL1 *alphaL1 +deltaL) * L1[1]
  
  for (i in seq(2, nL1)){
    dL1[i] = nL1 * alphaL1  * (L1[i - 1] - L1[i]) - deltaL * L1[i]   
  } 
  
  ###LARVAL INSTAR 2
  dL2= rep(0,nL2)
  dL2[1] = nL1*alphaL1*L1[nL1] - (nL2 * alphaL2 + deltaL) * L2[1] 
  
  for (i in seq(2, nL2)){
    dL2[i] = nL2 *   alphaL2  * (L2[i-1]-L2[i]) - deltaL*L2[i] 
  } 
  ###LARVAE 3
  dL3= rep(0,nL3)
  dL3[1] = (nL2*alphaL2*L2[nL2]) - (nL3 * alphaL3 * L3[1]) - deltaL * L3[1]  
  
  for (i in seq(2, nL3)){
    dL3[i] = nL3 *   alphaL3  * (L3[i - 1] - L3[i]) -  
      deltaL * L3[i] #Background mortality 
  } 
  
  ###LARVAE 4
  dL4= rep(0,nL4)
  
  dL4[1] = nL3*alphaL3*L3[nL3]  + #Recruitment from the previous stage
    nDL4 * alphaDL* DL4[nDL4] - #Recruitment from the breaking diapause stage
    (nL4 * alphaL4 * L4[1])- #Recruitment into the next developmental stage 
    DIA_INDUC1*L4[1]-                       #Recruitment into the diapausing stage
    (deltaDL+((0e-04*sum(L4)))) * L4[1]                        #Background mortality 
  
  for (i in seq(2, nL4)){
    dL4[i] =  nL4 *   alphaL4  * (L4[i - 1]- L4[i])-   
      (DIA_INDUC1*L4[i])- #Recruitment into the diapausing stage 
      (deltaL +((0e-04*sum(L4)))) * L4[i] #Background mortality 
  } 
  ###LARVAE 5
  
  dL5= rep(0,nL5)
  dL5[1] = nL4*alphaL4*L4[nL4] + #Recruitment from the previous stage 
    nDL5 * alphaDL* DL5[nDL5] -           #Recruitment from the diapause stage
    (nL5 * alphaL5  * L5[1]) #Recruitment into the next developmental stage
  - DIA_INDUC2*L5[1] -                 # Recriutment into the diapausing stage
    (deltaL+((0e-04*sum(L5))))*L5[1]                     #Background mortality 
  
  for (i in seq(2, nL5)){
    dL5[i] = nL5 *   alphaL5  * (L5[i - 1]-L5[i])-  #Recruitment from the previous stage
      #Recruitment intot he next substage      
      (DIA_INDUC2* L5[i])-                       #Recruitment into the diapausing stage
      (deltaL+(0e-04*sum(L5))) * L5[i]                             #Background mortality 
  } 
  
  ###The diapausing stage - 4th instar
  dDL4 = rep(0,nDL4)
  dDL4[1] =DIA_INDUC1 * sum(L4)-  #Recruitment from the fourth instar larvae
    alphaDL *nDL4*DL4[1] - #Recruitment out into the next substage
    deltaDL*DL4[1] #Backgroudn mortality 
  
  for (i in seq(2, nDL4))
  {
    dDL4[i] = nDL4 *   alphaDL  * (DL4[i - 1]-DL4[i])-  
      deltaDL * DL4[i] #Background mortality 
  } 
  ###The diapausing stage - 5th instar
  dDL5 = rep(0,nDL5)
  dDL5[1] = DIA_INDUC2*sum(L5) - #Diapause in
    alphaDL * nDL5 * DL5[1]- #Development out 
    deltaDL*DL5[1] # Mortality
  
  for (i in seq(2, nDL5)){
    dDL5[i] =nDL5 *   alphaDL  * (DL5[i - 1]-DL5[i]) - 
      deltaDL * DL5[i]  }
  
  #######
  #Pupae#
  #######
  dP = rep(0,nP)
  dP[1] =  nL5*alphaL5* L5[nL5] - #Recruitment in from previous stage
    (nP * alphaP * P[1]) - deltaP * P[1] #
  for (i in seq(2, nP))
  {
    dP[i] = nP * alphaP * (P[i - 1] - P[i])-  #Recruitment in from prevoius substage
      deltaP* P[i] #Backgroudn mortality 
  }
  
  #Reproductive Adult 
  dA.r = rep(0,nA.r)
  dA.r[1] = nP * alphaP* P[nP] - #Recruitment in from previous stage
    (nA.r* alphaA + deltaA ) * A.r[1]  #Recruitment out
  for (i in seq(2, nA.r))
  {
    dA.r[i] = nA.r* alphaA * (A.r[i - 1] - A.r[i]) -deltaA* A.r[i]
  }
  
  
  #senscent adult
  dA.s = nA.r*alphaA * A.r[nA.r]  - deltaA*A.s
  
  res = c(dE,dL1,dL2,dL3,dL4,dL5,dDL4,dDL5,dP,dA.r,dA.s)
  list(res)
  
}
tstep = seq(1,365,1)

nE=8 #The number of subcompartments in the egg-stage
nL1 =8 #The number of subcompartments in larval instar 1
nL2 =8 #The number of subcompartments in larval instar 2
nL3 = 8 #The number of subcompartments in larval instar 3
nL4 =8#The number of subcomppartmetns in larval instar 4
nL5=8 #The number of subcompartments in larval instar 5
nDL4 =8  #The number of subcompartments in diapausing larval instar 4
nDL5 =8 #The number of subcomparments in diapausing larval instar 5
nP =8 # The number of subcompartments in pupae
nA.r =8


yinit = c(
  E = rep(0, nE),
  L1 = rep(0,nL1),
  L2 = rep(0,nL2),
  L3 = rep(0,nL3),
  L4 = rep(0,nL4),
  L5 = rep(0,nL5),
  DL4 =rep(1,nDL4),
  DL5 =rep(1,nDL5),
  P = rep(0,nP),
  A.r = rep(0,nA.r),
  A.s = 0
)



p <- c( 
        alphaDL_a =0.7543 , alphaDL_b =0.2392,  alphaDL_c =7.1729,
        deltaAa=  3.328e-05,deltaAb= 1.179e-03,deltaAc=2.452e-02 ,
        deltaDLa=  3.328e-05, deltaDLb= 1.179e-03 , deltaDLc = 2.452e-02 ,
        DIA_INDUC1_A=0.2130,DIA_INDUC1_B = 80.0004,
        DIA_INDUC2_A= 0.5130,DIA_INDUC2_B = 80.0004)

yout<- ode(
  y = yinit ,
  times = tstep,
  func = MODEL_DDE_CM_DIA,
  parms=p)

E<-rowSums(yout[,(2:(nE+1))])
plot(E,type='l',main='E')

L1<- rowSums(yout[,(nE+2):(nE+nL1+1)])
plot(L1,type='l',main='L1')


L2 <- rowSums(yout[,(nE+nL1+2):(nE+nL1+nL2 +1)])
plot(L2, type='l',main ='L2')

L3 <- rowSums(yout[,(nE+nL1+nL2+2):(nE+nL1+nL2 + nL3+1)])
plot(L3, type='l',main ='L3')

L4 <- rowSums(yout[,(nE+nL1+nL2+nL3+2):(nE+nL1+nL2 + nL3+nL4+1)])
plot(L4, type='l',main ='L4')


L5 <- rowSums(yout[,(nE+nL1+nL2+nL3+nL4+2):(nE+nL1+nL2 + nL3+nL4+nL5+1)])
plot(L5, type='l',main ='L5')

DL4 <- rowSums(yout[,(nE+nL1+nL2+nL3+nL4+nL5+2):(nE+nL1+nL2 + nL3+nL4+nL5+nDL4 + 1)])
plot(DL4,main ='DL4',type='b')

DL5<- rowSums(yout[,(nE+nL1+nL2 + nL3+nL4+nL5+nDL4 + 2):(nE+nL1+nL2 + nL3+nL4+nL5+nDL4 +nDL5+ 1)])
plot(DL5, type='l',main ='DL5')


P<- rowSums(yout[,(nE+nL1+nL2 + nL3+nL4+nL5+nDL4 + nDL5+2):(nE+nL1+nL2 + nL3+nL4+nL5+nDL4 +nDL5+nP+ 1)])
plot(P, type='l',main ='P')

A.r<- rowSums(yout[,(nE+nL1+nL2 + nL3+nL4+nL5+nDL4 + nDL5+nP+2):
                     (nE+nL1+nL2 + nL3+nL4+nL5+nDL4 +nDL5+nP+nA.r +1)])
plot(A.r,col='red')
lines(FULLCM$julian, Fu)
plot(seq(6028,12054),log(A.r+1)[6028:12054], type='l',main ='A')
lines(FullCM$julian[FullCM$julian > 6028],
      log(FullCM$Mean[FullCM$julian > 6028]+1),col='red')
