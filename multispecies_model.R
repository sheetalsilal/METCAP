# MULTIPATCH MULTI-SPECIES MODEL FOR FALCIPARUM AND VIVAX MALARIA TRANSMISSION AND CONTROL # 
# THE MODEL IS WRITTEN IN R AND HAS SUPPLEMENTARY C CODE THAT MUST BE COMPILED AND LOADED 
#
# The data time step is one month and can be adapted by altering the dtout parameter
# 
# INSTRUCTIONS
# Set initial set up parameters (lines 25-32)
# Load in data set and define sheets accordingly
# Set names for RData image and top parameter sets object.
# Set number of cores for parallelisation (default=24)
# Set number of parameter sets to generate. 

# ************************************************************************************* #
#Load packages
library(adaptivetau)
library(deSolve)
library(XLConnect)               # load XLConnect package

#setwd("")  # Set appropriate directory for specific analysis
setwd("/Volumes/SHEETAL-RES/APLMA/Manuscripts")  # Set appropriate directory for specific analysis

system("R CMD SHLIB eq0.c")     # compiles c code
dyn.load("eq0.so")              # loads compiled code to memory (.so - OS X, .dll - Windows)


N<-22   # number of patches
B<-23   # number of variables per patch
A<-261  # number of transitions per patch
V<-N*B # total number of variables
L<-N*A #total number of transitions
startyear=1995 # starting year of simulation
tyears<-20 # total years of simulation
dtout<-1/12 # output timestep
tsteps<-round(tyears/dtout) # number of time steps
time<-startyear+seq(0,tyears,dtout) # time vector

# Call C function from eq0.so in R
EQ<-function(L, N, oldeq, transit,transitionsiu1,transitionsiu2,transitionsiv1,transitionsiv2){
  len<-length(oldeq)
  .C("EQ",
     as.integer(L),
     as.integer(N),
     as.double(oldeq),
     as.double(transit),
     as.integer(transitionsiu1),
     as.integer(transitionsiu2),
     as.integer(transitionsiv1),
     as.integer(transitionsiv2),
     as.double(vector("double",length = len))
  )[[9]]
  
}

# ************************************************************************************* #
# import data
# ************************************************************************************* #

alldata = loadWorkbook("data_template.xlsx") 

# Population data
pvxy = readWorksheet(alldata, sheet="Pop")
pvxy<-pvxy[1:N,]
yrprim<-pvxy[,10] #year primaquine adopted

#ITN distribution over time (No. of LLIN+ITN)
itndat = readWorksheet(alldata, sheet="ITN")
itn_dis<-itndat[,(1:N)+2]
itn_time<-itndat[,1]+ itndat[,2]/12
cov_itn_dat<-array(0,c(length(itn_time),N))
for (n in 1:N){ 
  for (i in 1:length(itn_time)){
    cov_itn_dat[i,n]<-100*min(as.numeric(itn_dis[i,n]), pvxy[n,6])/pvxy[n,6]
  }}

#IRS coverage over time (No. of people protected by IRS)
irsdat = readWorksheet(alldata, sheet="IRS")
irs_dis<-irsdat[,(1:N)+2]
irs_time<-irsdat[,1]+irsdat[,2]/12
cov_irs_dat<-matrix(0,nrow=length(irs_time),ncol=N)
for (n in 1:N){
  for (i in 1:length(itn_time)){
    cov_irs_dat[i,n]<-100*min(pvxy[n,6],as.numeric(irs_dis[i,n]))/pvxy[n,6]
  }}

# VMW coverage over time
vmwdat = readWorksheet(alldata, sheet="VMW")
vmw_dis<-vmwdat[,(1:N)+2]
vmw_time<-vmwdat[,1] + vmwdat[,2]/12
cov_vmw_dat<-matrix(0,nrow=length(vmw_time),ncol=N)
for (n in 1:N){
  for (i in 1:length(vmw_time)){
    cov_vmw_dat[i,n]<-100*min(pvxy[n,3],as.numeric(vmw_dis[i,n]))/pvxy[n,3]
  }}

# El Nino index over time
elnino_index = readWorksheet(alldata, sheet="Climate")
elnino<-elnino_index[,5] # column 5 for average value until 2000; else column 4 for data values from start
elnino_t<-elnino_index[,1]+elnino_index[,2]/12 


# ************************************************************************************* #
# define variables
# ************************************************************************************* #
# FALCIPARUM
# 1=S: uninfected non-immune
# 2=In: infected submicro
# 3=Ia: infected asymp
# 4=Ic: infected clinical
# 5=Is: infected severe
# 6=To: treated unrecorded
# 7=Tv: treated VMW
# 8=Th: treated HIS
# 9=R: uninfected immune
# 10=H: uninfected and HRP2 positive

#VIVAX
# 11=S: uninfected non-immune
# 12=In: infected submicro
# 13=Ia: infected asymp
# 14=Ic: infected clinical
# 15=Is: infected severe
# 16=To: treated unrecorded
# 17=Tv: treated VMW without G6PDd
# 18=Th: treated HIS without G6PDd
# 19=R: uninfected immune with no hypnozoites
# 20=L: uninfected immune with hypnozoites 
# 21=Togd: treated unrecorded with G6PDd
# 22=Tvgd: treated VMW with G6PDd
# 23=Thgd: treated HIS with G6PDd

falpop<-1:10
vivpop<-11:23

# ************************************************************************************* #
# define indices
# ************************************************************************************* #
varind<-matrix(0,nrow=B,ncol=N)
traind<-matrix(0,nrow=A,ncol=N)
for (n in 1:N){
  for (b in 1:B){
    varind[b,n]<-(n-1)*B+b
  }
  for (a in 1:A){
    traind[a,n]<-(n-1)*A+a
  }
}

# ************************************************************************************* #
# define transitions
# ************************************************************************************* #
# first transition is given without index
transitions =ssa.maketrans(V,rbind(varind[4,1], +1)) # birth mos patch 1
for (n in 1:N){
  # Falciparum indep flows (1:59)
  transitions[traind[1,n]]<-ssa.maketrans(V,rbind(varind[1,n],0,varind[1,n], +1)) # birth humans
  transitions[traind[2,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[1,n],0)) # death S=1 
  transitions[traind[3,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1,varind[1,n],0)) # death In=2 
  transitions[traind[4,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1,varind[1,n],0)) # death Ia=3
  transitions[traind[5,n]]<-ssa.maketrans(V,rbind(varind[4,n], -1,varind[1,n],0)) # death Ic=4
  transitions[traind[6,n]]<-ssa.maketrans(V,rbind(varind[5,n], -1,varind[1,n],0)) # death Is=5  
  transitions[traind[7,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1,varind[1,n],0)) # death To=6 
  transitions[traind[8,n]]<-ssa.maketrans(V,rbind(varind[7,n], -1,varind[1,n],0)) # death Tv=7
  transitions[traind[9,n]]<-ssa.maketrans(V,rbind(varind[8,n], -1,varind[1,n],0)) # death Th=8
  transitions[traind[10,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1,varind[1,n],0)) # death R=9
  transitions[traind[11,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1,varind[1,n],0)) # death H=10
  transitions[traind[12,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[2,n],+1)) # incidence S to In
  transitions[traind[13,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[3,n],+1)) # incidence S to Ia
  transitions[traind[14,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[4,n],+1)) # incidence S to Ic 
  transitions[traind[15,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[6,n],+1)) # incidence S to To
  transitions[traind[16,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[7,n],+1)) # incidence S to Tv
  transitions[traind[17,n]]<-ssa.maketrans(V,rbind(varind[1,n], -1,varind[8,n],+1)) # incidence S to Th
  transitions[traind[18,n]]<-ssa.maketrans(V,rbind(varind[5,n], -1, varind[4,n], +1)) # recovery Is to Ic
  transitions[traind[19,n]]<-ssa.maketrans(V,rbind(varind[4,n], -1, varind[3,n], +1)) # recovery Ic to Ia
  transitions[traind[20,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1, varind[2,n], +1)) # recovery Ia to In
  transitions[traind[21,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[9,n], +1)) # recovery In to R
  transitions[traind[22,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1, varind[10,n], +1)) # recovery To to H
  transitions[traind[23,n]]<-ssa.maketrans(V,rbind(varind[7,n], -1, varind[10,n], +1)) # recovery Tv to H
  transitions[traind[24,n]]<-ssa.maketrans(V,rbind(varind[8,n], -1, varind[10,n], +1)) # recovery Th to H
  transitions[traind[25,n]]<-ssa.maketrans(V,rbind(varind[5,n], -1, varind[10,n], +1)) # recovery Is to H
  transitions[traind[26,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[3,n], +1)) # incidence In to Ia
  transitions[traind[27,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1, varind[4,n], +1)) # incidence Ia to Ic
  transitions[traind[28,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1, varind[6,n], +1)) # incidence Ia to To
  transitions[traind[29,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1, varind[7,n], +1)) # incidence Ia to Tv
  transitions[traind[30,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1, varind[8,n], +1)) # incidence Ia to Th
  transitions[traind[31,n]]<-ssa.maketrans(V,rbind(varind[4,n], -1, varind[5,n], +1)) # incidence Ic to Is
  transitions[traind[32,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[4,n], +1)) # incidence In to Ic
  transitions[traind[33,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[6,n], +1)) # incidence In to To
  transitions[traind[34,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[7,n], +1)) # incidence In to Tv
  transitions[traind[35,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[8,n], +1)) # incidence In to Th
  transitions[traind[36,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[2,n], +1)) # incidence R to In 
  transitions[traind[37,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[3,n], +1)) # incidence R to Ia 
  transitions[traind[38,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[4,n], +1)) # incidence R to Ic 
  transitions[traind[39,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[6,n], +1)) # incidence R to To 
  transitions[traind[40,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[7,n], +1)) # incidence R to Tv
  transitions[traind[41,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[8,n], +1)) # incidence R to Th
  transitions[traind[42,n]]<-ssa.maketrans(V,rbind(varind[9,n], -1, varind[1,n], +1)) # loss imm R to S 
  transitions[traind[43,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[2,n], +1)) # incidence H to In 
  transitions[traind[44,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[3,n], +1)) # incidence H to Ia
  transitions[traind[45,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[4,n], +1)) # incidence H to Ic
  transitions[traind[46,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[6,n], +1)) # incidence H to To
  transitions[traind[47,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[7,n], +1)) # incidence H to Tv
  transitions[traind[48,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[8,n], +1)) # incidence H to Th
  transitions[traind[49,n]]<-ssa.maketrans(V,rbind(varind[10,n], -1, varind[9,n], +1)) # loss HRP2 H to  R
  transitions[traind[50,n]]<-ssa.maketrans(V,rbind(varind[5,n], -1, varind[1,n],0)) # death Is fatal malaria
  transitions[traind[51,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1, varind[3,n],+1)) # failed trt To to Ia
  transitions[traind[52,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1, varind[4,n],+1)) # failed trt To to Ic
  transitions[traind[53,n]]<-ssa.maketrans(V,rbind(varind[6,n], -1, varind[8,n],+1)) # failed trt To to Th
  transitions[traind[54,n]]<-ssa.maketrans(V,rbind(varind[7,n], -1, varind[3,n],+1)) # failed trt Tv to Ia
  transitions[traind[55,n]]<-ssa.maketrans(V,rbind(varind[7,n], -1, varind[4,n],+1)) # failed trt Tv to Ic
  transitions[traind[56,n]]<-ssa.maketrans(V,rbind(varind[7,n], -1, varind[7,n],+1)) # failed trt Tv to Th
  transitions[traind[57,n]]<-ssa.maketrans(V,rbind(varind[8,n], -1, varind[3,n],+1)) # failed trt Th to Ia
  transitions[traind[58,n]]<-ssa.maketrans(V,rbind(varind[8,n], -1, varind[4,n],+1)) # failed trt Th to Ic
  transitions[traind[59,n]]<-ssa.maketrans(V,rbind(varind[8,n], -1, varind[8,n],+1)) # failed trt Th to Th
  
  #Vivax indep flows (59-165)
  transitions[traind[60,n]]<-ssa.maketrans(V,rbind(varind[11,n],0,varind[11,n], +1)) # birth humans
  transitions[traind[61,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[11,n],0)) # death S=11 
  transitions[traind[62,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1,varind[11,n],0)) # death In=12 
  transitions[traind[63,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1,varind[11,n],0)) # death Ia=13
  transitions[traind[64,n]]<-ssa.maketrans(V,rbind(varind[14,n], -1,varind[11,n],0)) # death Ic=14
  transitions[traind[65,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1,varind[11,n],0)) # death Is=15  
  transitions[traind[66,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1,varind[11,n],0)) # death To=16 
  transitions[traind[67,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1,varind[11,n],0)) # death Tv=17
  transitions[traind[68,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1,varind[11,n],0)) # death Th=18
  transitions[traind[69,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1,varind[11,n],0)) # death R=19
  transitions[traind[70,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1,varind[11,n],0)) # death L=20
  transitions[traind[71,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1,varind[11,n],0)) # death Togd=21
  transitions[traind[72,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1,varind[11,n],0)) # death Tvgd=22
  transitions[traind[73,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1,varind[11,n],0)) # death Thgd=23
  transitions[traind[74,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[12,n],+1)) # incidence S to In
  transitions[traind[75,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[13,n],+1)) # incidence S to Ia
  transitions[traind[76,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[14,n],+1)) # incidence S to Ic 
  transitions[traind[77,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[16,n],+1)) # incidence S to To
  transitions[traind[78,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[17,n],+1)) # incidence S to Tv
  transitions[traind[79,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[18,n],+1)) # incidence S to Th
  transitions[traind[80,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[21,n],+1)) # incidence S to Togd
  transitions[traind[81,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[22,n],+1)) # incidence S to Tvgd
  transitions[traind[82,n]]<-ssa.maketrans(V,rbind(varind[11,n], -1,varind[23,n],+1)) # incidence S to Thgd
  transitions[traind[83,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1, varind[14,n],+1)) # recovery Is to Ic
  transitions[traind[84,n]]<-ssa.maketrans(V,rbind(varind[14,n], -1, varind[13,n],+1)) # recovery Ic to Ia
  transitions[traind[85,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[12,n],+1)) # recovery Ia to In
  transitions[traind[86,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[19,n],+1)) # recovery In to R
  transitions[traind[87,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[20,n],+1)) # recovery In to L
  transitions[traind[88,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[19,n], +1)) # recovery To to R
  transitions[traind[89,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[19,n], +1)) # recovery Tv to R
  transitions[traind[90,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[19,n], +1)) # recovery Th to R
  transitions[traind[91,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[19,n], +1)) # recovery Togd to R
  transitions[traind[92,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[19,n], +1)) # recovery Tvgd to R
  transitions[traind[93,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[19,n], +1)) # recovery Thgd to R
  transitions[traind[94,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1, varind[19,n], +1)) # recovery Is to R
  transitions[traind[95,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[20,n], +1)) # recovery To to L
  transitions[traind[96,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[20,n], +1)) # recovery Tv to L
  transitions[traind[97,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[20,n], +1)) # recovery Th to L
  transitions[traind[98,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[20,n], +1)) # recovery Togd to L
  transitions[traind[99,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[20,n], +1)) # recovery Tvgd to L
  transitions[traind[100,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[20,n], +1)) # recovery Thgd to L
  transitions[traind[101,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1, varind[20,n], +1)) # recovery Is to L
  transitions[traind[102,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[13,n], +1)) # incidence In to Ia
  transitions[traind[103,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[14,n], +1)) # incidence In to Ic
  transitions[traind[104,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[16,n], +1)) # incidence In to To
  transitions[traind[105,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[17,n], +1)) # incidence In to Tv
  transitions[traind[106,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[18,n], +1)) # incidence In to Th
  transitions[traind[107,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[21,n], +1)) # incidence In to Togd
  transitions[traind[108,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[22,n], +1)) # incidence In to Tvgd
  transitions[traind[109,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[23,n], +1)) # incidence In to Thgd
  transitions[traind[110,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[14,n], +1)) # incidence Ia to Ic
  transitions[traind[111,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[16,n], +1)) # incidence Ia to To
  transitions[traind[112,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[17,n], +1)) # incidence Ia to Tv
  transitions[traind[113,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[18,n], +1)) # incidence Ia to Th
  transitions[traind[114,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[21,n], +1)) # incidence Ia to Togd
  transitions[traind[115,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[22,n], +1)) # incidence Ia to Tvgd
  transitions[traind[116,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[23,n], +1)) # incidence Ia to Thgd
  transitions[traind[117,n]]<-ssa.maketrans(V,rbind(varind[14,n], -1, varind[15,n], +1)) # incidence Ic to Is
  transitions[traind[118,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[12,n], +1)) # incidence R to In 
  transitions[traind[119,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[13,n], +1)) # incidence R to Ia 
  transitions[traind[120,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[14,n], +1)) # incidence R to Ic 
  transitions[traind[121,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[16,n], +1)) # incidence R to To 
  transitions[traind[122,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[17,n], +1)) # incidence R to Tv
  transitions[traind[123,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[18,n], +1)) # incidence R to Th
  transitions[traind[124,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[21,n], +1)) # incidence R to Togd
  transitions[traind[125,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[22,n], +1)) # incidence R to Tvgd
  transitions[traind[126,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[23,n], +1)) # incidence R to Thgd
  transitions[traind[127,n]]<-ssa.maketrans(V,rbind(varind[19,n], -1, varind[11,n], +1)) # loss imm R to S 
  transitions[traind[128,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[12,n], +1)) # incidence L to In 
  transitions[traind[129,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[13,n], +1)) # incidence L to Ia
  transitions[traind[130,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[14,n], +1)) # incidence L to Ic
  transitions[traind[131,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[16,n], +1)) # incidence L to To
  transitions[traind[132,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[17,n], +1)) # incidence L to Tv
  transitions[traind[133,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[18,n], +1)) # incidence L to Th
  transitions[traind[134,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[21,n], +1)) # incidence L to Togd
  transitions[traind[135,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[22,n], +1)) # incidence L to Tvgd
  transitions[traind[136,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[23,n], +1)) # incidence L to Thgd
  transitions[traind[137,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[12,n], +1)) # relapse L to In 
  transitions[traind[138,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[13,n], +1)) # relapse L to Ia
  transitions[traind[139,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[14,n], +1)) # relapse L to Ic
  transitions[traind[140,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[16,n], +1)) # relapse L to To
  transitions[traind[141,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[17,n], +1)) # relapse L to Tv
  transitions[traind[142,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[18,n], +1)) # relapse L to Th
  transitions[traind[143,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[21,n], +1)) # relapse L to Togd
  transitions[traind[144,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[22,n], +1)) # relapse L to Tvgd
  transitions[traind[145,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[23,n], +1)) # relapse L to Thgd
  transitions[traind[146,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[11,n], +1)) # death of hypnozoites L to S
  transitions[traind[147,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1, varind[11,n],0)) # death Is fatal malaria
  
  transitions[traind[148,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[13,n], +1)) # failed trt To to Ia
  transitions[traind[149,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[14,n], +1)) # failed trt To to Ic
  transitions[traind[150,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[18,n], +1)) # failed trt To to Th
  
  transitions[traind[151,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[13,n], +1)) # failed trt Tv to Ia
  transitions[traind[152,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[14,n], +1)) # failed trt Tv to Ic
  transitions[traind[153,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[17,n], +1)) # failed trt Tv to Tv
  
  transitions[traind[154,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[13,n], +1)) # failed trt Th to Ia
  transitions[traind[155,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[14,n], +1)) # failed trt Th to Ic
  transitions[traind[156,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[18,n], +1)) # failed trt Th to Th
  
  transitions[traind[157,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[13,n], +1)) # failed trt Togd to Ia
  transitions[traind[158,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[14,n], +1)) # failed trt Togd to Ic
  transitions[traind[159,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[23,n], +1)) # failed trt Togd to Thgd
  
  transitions[traind[160,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[13,n], +1)) # failed trt Tvgd to Ia
  transitions[traind[161,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[14,n], +1)) # failed trt Tvgd to Ic
  transitions[traind[162,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[22,n], +1)) # failed trt Tvgd to Tvgd

  transitions[traind[163,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[13,n], +1)) # failed trt Thgd to Ia
  transitions[traind[164,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[14,n], +1)) # failed trt Thgd to Ic
  transitions[traind[165,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[23,n], +1)) # failed trt Thgd to Thgd
  
  #Entanglements (166-216)
  transitions[traind[166,n]]<-ssa.maketrans(V,rbind(varind[2,n], -1, varind[10,n],+1)) # dual trt Fal In to H
  transitions[traind[167,n]]<-ssa.maketrans(V,rbind(varind[3,n], -1, varind[10,n],+1)) # dual trt Fal Ia to H
  transitions[traind[168,n]]<-ssa.maketrans(V,rbind(varind[4,n], -1, varind[8,n],+1)) # dual trt Fal Ic to H
  transitions[traind[169,n]]<-ssa.maketrans(V,rbind(varind[5,n], -1, varind[6,n],+1)) # dual trt Fal Is to H
  transitions[traind[170,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[19,n],+1)) # dual trt Viv In to R
  transitions[traind[171,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[19,n],+1)) # dual trt Viv Ia to R
  transitions[traind[172,n]]<-ssa.maketrans(V,rbind(varind[14,n], -1, varind[19,n],+1)) # dual trt Viv Ic to R
  transitions[traind[173,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1, varind[19,n],+1)) # dual trt Viv Is to R
  transitions[traind[174,n]]<-ssa.maketrans(V,rbind(varind[12,n], -1, varind[20,n],+1)) # dual trt Viv In to L
  transitions[traind[175,n]]<-ssa.maketrans(V,rbind(varind[13,n], -1, varind[20,n],+1)) # dual trt Viv Ia to L
  transitions[traind[176,n]]<-ssa.maketrans(V,rbind(varind[14,n], -1, varind[20,n],+1)) # dual trt Viv Ic to L
  transitions[traind[177,n]]<-ssa.maketrans(V,rbind(varind[15,n], -1, varind[20,n],+1)) # dual trt Viv Is to L
  transitions[traind[178,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[19,n],+1)) # primaquine recovery To to R
  transitions[traind[179,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[19,n],+1)) # primaquine recovery Tv to R
  transitions[traind[180,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[19,n],+1)) # primaquine recovery Th to R
  
  transitions[traind[181,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[20,n],+1)) # primaquine recovery To to L
  transitions[traind[182,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[20,n],+1)) # primaquine recovery Tv to L
  transitions[traind[183,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[20,n],+1)) # primaquine recovery Th to L
  transitions[traind[184,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[19,n],+1)) # primaquine recovery Togd to R
  transitions[traind[185,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[19,n],+1)) # primaquine recovery Tvgd to R
  transitions[traind[186,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[19,n],+1)) # primaquine recovery Thgd to R
  transitions[traind[187,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[20,n],+1)) # primaquine recovery Togd for Gdef to L
  transitions[traind[188,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[20,n],+1)) # primaquine recovery for Gdef Tvgd to L
  transitions[traind[189,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[20,n],+1)) # primaquine recovery for Gdef Thgd to L
  transitions[traind[190,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[19,n],+1)) # primaquine recovery Togd to R
  transitions[traind[191,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[19,n],+1)) # primaquine recovery Tvgd to R
  transitions[traind[192,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[19,n],+1)) # primaquine recovery Thgd to R
  transitions[traind[193,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[20,n],+1)) # primaquine recovery Togd for Gdef to L
  transitions[traind[194,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[20,n],+1)) # primaquine recovery for Gdef Tvgd to L
  transitions[traind[195,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[20,n],+1)) # primaquine recovery for Gdef Thgd to L
  transitions[traind[196,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[19,n],+1)) # masked ACT recovery To to R
  transitions[traind[197,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[19,n],+1)) # masked ACT recovery Tv to R
  transitions[traind[198,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[19,n],+1)) # masked ACT recovery Th to R
  transitions[traind[199,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[20,n],+1)) # masked ACT recovery To to L
  transitions[traind[200,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[20,n],+1)) # masked ACT recovery Tv to L
  transitions[traind[201,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[20,n],+1)) # masked ACT recovery Th to L
  transitions[traind[202,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[19,n],+1)) # masked ACT recovery Togd to R
  transitions[traind[203,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[19,n],+1)) # masked ACT recovery Tvgd to R
  transitions[traind[204,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[19,n],+1)) # masked ACT recovery Thgd to R
  transitions[traind[205,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[20,n],+1)) # masked ACT recovery Togd for Gdef to L
  transitions[traind[206,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[20,n],+1)) # masked ACT recovery for Gdef Tvgd to L
  transitions[traind[207,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[20,n],+1)) # masked ACT recovery for Gdef Thgd to L
  
  transitions[traind[208,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[12,n], +1)) # relapse L to In + triggering
  transitions[traind[209,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[13,n], +1)) # relapse L to Ia+ triggering
  transitions[traind[210,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[14,n], +1)) # relapse L to Ic+ triggering
  transitions[traind[211,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[16,n], +1)) # relapse L to To+ triggering
  transitions[traind[212,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[17,n], +1)) # relapse L to Tv+ triggering
  transitions[traind[213,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[18,n], +1)) # relapse L to Th+ triggering
  transitions[traind[214,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[21,n], +1)) # relapse L to Togd+ triggering
  transitions[traind[215,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[22,n], +1)) # relapse L to Tvgd+ triggering
  transitions[traind[216,n]]<-ssa.maketrans(V,rbind(varind[20,n], -1, varind[23,n], +1)) # relapse L to Thgd+ triggering

  #failed treatments [217:261]
  transitions[traind[217,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[13,n], +1)) # primaquine failed trt To to Ia
  transitions[traind[218,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[14,n], +1)) # primaquine failed trt To to Ic
  transitions[traind[219,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[18,n], +1)) # primaquine failed trt To to Th
  
  transitions[traind[220,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[13,n], +1)) # primaquine failed trt Tv to Ia
  transitions[traind[221,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[14,n], +1)) # primaquine failed trt Tv to Ic
  transitions[traind[222,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[17,n], +1)) # primaquine failed trt Tv to Tv
  
  transitions[traind[223,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[13,n], +1)) # primaquine failed trt Th to Ia
  transitions[traind[224,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[14,n], +1)) # primaquine failed trt Th to Ic
  transitions[traind[225,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[18,n], +1)) # primaquine failed trt Th to Th

  transitions[traind[226,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[13,n], +1)) # masked act failed trt To to Ia
  transitions[traind[227,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[14,n], +1)) # masked act failed trt To to Ic
  transitions[traind[228,n]]<-ssa.maketrans(V,rbind(varind[16,n], -1, varind[18,n], +1)) # masked act failed trt To to Th
  
  transitions[traind[229,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[13,n], +1)) # masked act failed trt Tv to Ia
  transitions[traind[230,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[14,n], +1)) # masked act failed trt Tv to Ic
  transitions[traind[231,n]]<-ssa.maketrans(V,rbind(varind[17,n], -1, varind[17,n], +1)) # masked act failed trt Tv to Tv
  
  transitions[traind[232,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[13,n], +1)) # masked act failed trt Th to Ia
  transitions[traind[233,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[14,n], +1)) # masked act failed trt Th to Ic
  transitions[traind[234,n]]<-ssa.maketrans(V,rbind(varind[18,n], -1, varind[18,n], +1)) # masked act failed trt Th to Th
  
  transitions[traind[235,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[13,n], +1)) # masked failed trt Togd to Ia
  transitions[traind[236,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[14,n], +1)) # masked failed trt Togd to Ic
  transitions[traind[237,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[23,n], +1)) # masked failed trt Togd to Thgd
  
  transitions[traind[238,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[13,n], +1)) # masked failed trt Tvgd to Ia
  transitions[traind[239,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[14,n], +1)) # masked failed trt Tvgd to Ic
  transitions[traind[240,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[22,n], +1)) # masked failed trt Tvgd to Tvgd
  
  transitions[traind[241,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[13,n], +1)) # masked failed trt Thgd to Ia
  transitions[traind[242,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[14,n], +1)) # masked failed trt Thgd to Ic
  transitions[traind[243,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[23,n], +1)) # masked failed trt Thgd to Thgd

  transitions[traind[244,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[13,n], +1)) # test + ACT failed trt Togd to Ia
  transitions[traind[245,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[14,n], +1)) # test + ACT failed trt Togd to Ic
  transitions[traind[246,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[23,n], +1)) # test + ACT failed trt Togd to Thgd
  
  transitions[traind[247,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[13,n], +1)) # test + ACT failed trt Tvgd to Ia
  transitions[traind[248,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[14,n], +1)) # test + ACT failed trt Tvgd to Ic
  transitions[traind[249,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[22,n], +1)) # test + ACT failed trt Tvgd to Tvgd
  
  transitions[traind[250,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[13,n], +1)) # test + ACT failed trt Thgd to Ia
  transitions[traind[251,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[14,n], +1)) # test + ACT failed trt Thgd to Ic
  transitions[traind[252,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[23,n], +1)) # test + ACT failed trt Thgd to Thgd
  
  transitions[traind[253,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[13,n], +1)) # test - primaquine failed trt Togd to Ia
  transitions[traind[254,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[14,n], +1)) # test - primaquine failed trt Togd to Ic
  transitions[traind[255,n]]<-ssa.maketrans(V,rbind(varind[21,n], -1, varind[23,n], +1)) # test - primaquine failed trt Togd to Thgd
  
  transitions[traind[256,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[13,n], +1)) # test - primaquine failed trt Tvgd to Ia
  transitions[traind[257,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[14,n], +1)) # test - primaquine failed trt Tvgd to Ic
  transitions[traind[258,n]]<-ssa.maketrans(V,rbind(varind[22,n], -1, varind[22,n], +1)) # test - primaquine failed trt Tvgd to Tvgd
  
  transitions[traind[259,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[13,n], +1)) # test - primaquine failed trt Thgd to Ia
  transitions[traind[260,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[14,n], +1)) # test - primaquine failed trt Thgd to Ic
  transitions[traind[261,n]]<-ssa.maketrans(V,rbind(varind[23,n], -1, varind[23,n], +1)) # test - primaquine failed trt Thgd to Thgd
  
}
#Alternate formulation of transitions matrix (used in epimodel function)
transitions2<-NULL
for (i in 1: length(transitions)){
  transitions2<-rbind(transitions2,cbind(as.integer(names(transitions[[i]]))[1],as.integer(names(transitions[[i]]))[2], transitions[[i]][1], transitions[[i]][2]))
}
row.names(transitions2)<-NULL
transitionsiu1<-transitions2[,1]
transitionsiu2<-transitions2[,2]
transitionsiv1<-transitions2[,3]
transitionsiv2<-transitions2[,4]


# ************************************************************************************* #
# Set the parameters (vivax-specific parameters prefixed with v)
# ************************************************************************************* #
pars = list(
  # Common parameters to both species  
  Pmaxf=runif(N)*pvxy[,6], # the maximum populations of mosquitoes in each patch
  Pmaxv=runif(N)*pvxy[,6], # the maximum populations of mosquitoes in each patch
  delta_m=365.25/14, # death rate of mosquitoes
  amp=rep(1,N), # amplitude of seasonal forcing
  phi=pvxy[,8], # phase angle of seasonal forcing   
  bites=365.25/3, # biting rate
  prob_h=0.5, # probability that a bite will result in infection (mosquito to human)
  indexhet=c(4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,3,4.5,4.5,3,4.5,12,12,4.5,4.5,12,12,4.5), # parameter for connectivity about 1
  demog=1/50, # birth/death rate (mu)
  eff_irs=0.15, #efficacy of IRS
  eff_itn=0.15, #efficacy of bednets
  hl_net=1.5, # half-life of bednets
  eff_vmw=0.6, # treatment-seeking*effectiveness of VMW
  eff_his=0.3, # treatment-seeking*effectiveness of HIS
  eff_oth=0.1, # treatment-seeking*effectiveness of other systems e.g. Private
  elneff=21, # smoothing parameter of el nino patterns. 
  propgd=pvxy[,9], #Proportion of patch that is G6PDdef
  pgd=0.3, #probability of clinical if G6PDdef
  sensgd=0.97, #Sensitivity of test to detect G6PDdef
  mask=0.05, #probability of vivax mix cases masked as falciparum
  # Humans
  #falciparum parameters 
  gamma_m=365.25/10, # latent period
  prob_m=0.5, # probability that a bite will result in infection (human to mosquito)
  ps=0.9, # probability of clinical if non-immune
  psn=0.1, # probability of sub-patent given not clin, non-immune
  pr=0.1, # probability of clinical if immune
  prn=0.5, # probability of sub-patent given not clin, immune
  nuc=365.25/10, # recovery from clinical symptoms untreated (r_c)
  nua=365.25/130, # recovery from asym untreated (r_a)
  nus=365.25/5, # recovery from severe (r_s)
  psev=0.03, # probability of severe disease given clinical
  pmort=0.2, # proportional of all untreated severe cases that die (theta1)
  pmortt=0.15, # proportional of all treated severe cases that die (theta2)
  tausev=0.8, # probability that a severe infection is treated
  gamma_h=365.25/21, # latent period in humans
  zeta_a=12.6/27, # relative infectiousness of asym to clinical
  zeta_n=3.9/27, # relative infectiousness of submicro to clinical
  chi=365.25/28, # rate of loss of detectable HRP2 by RDT 
  # this will change if the RDT detection limit changes chi=chi_old*(mn_c-dl_RDTold)/(mn_c-dl_RDTnew)
  omega=1, # loss of immunity
  #control
  nut=365.25/3, # recovery rate under treatment (r_t)
  nuq=365.25/6, # recovery rate under quinine treatment (r_q)
  ptf=0.05, #Probability of treatment failure
  ptfc=0.75, # Probability of being clinical after treatment failure
  ptftr=0.27, #Probability of seeking trt if clinical, after treatment failure
  
  # diagnostics
  dl_RDT=log10(200), # standard RDT detection limit
  dl_micro=log10(100), # micro detection limit
  dl_qPCR=log10(2), # standard PCR detection limit
  mn_n=log10(5), # mean parasiteamia for sub micro
  mn_a=log10(1000), # mean parasiteamia for asym
  mn_c=log10(25000), # mean parasiteamia for clinical
  mn_s=log10(350000), # mean parasiteamia for severe
  sd_n=0.75, # sd parasiteamia for sub micro
  sd_a=1.5, # sd parasiteamia for asym
  sd_c=1.3, # sd parasiteamia for clinical
  sd_s=0.26, # sd parasiteamia for severe           
  
  #vivax parameters
  vgamma_m=365.25/12, # latent period
  vprob_h=0.23, # probability that a bite will result in infection (human to mosquito)
  vps=0.9, # probability of clinical if non-immune
  vpsn=0.1, # probability of sub-patent given not clin, non-immune
  vpr=0.1, # probability of clinical if immune
  vprn=0.17, # probability of sub-patent given not clin, immune
  vnuc=365.25/20, # recovery from clinical symptoms untreated (r_c)
  vnua=365.25/365.25, # recovery from asym untreated (r_a)
  vnus=365.25/3, # recovery from severe (r_s)
  vpsev=0.03, # probability of severe disease given clinical
  vpmort=0.2, # proportional of all untreated severe cases that die (theta)
  vpmortt=0.02, # proportional of all treated severe cases that die (theta)
  vtausev=0.8, # probability that a severe infection is treated
  vgamma_h=365.25/17, # latent period in humans
  vzeta_a=1, # relative infectiousness of asym to clinical
  vzeta_n=1, # relative infectiousness of submicro to clinical
  vomega=1, # loss of immunity
  vprel=0.25, # probability of relapse
  vincprel=0.30, #Probabily of relapse due to triggering
  vrel=365.25/100, # rate of relapse
  vph=0.68 , # probability of recovering with hypnozoites under ACT
  vphprim=0.13, # probability of recovering with hypnozoites under primaquine
  vkappa=365.25/400 , #  hypnozoite death rate 
  vnut=365.25/3, # recovery rate under treatment (r_t)
  vnuq=365.25/6, # recovery rate under quinine treatment (r_q)
  vnup=365/14, # recovery rate under primaquine
  vptf=0.05, #Probability of treatment failure on FLT
  vptfp=1-(0.9*0.85), #Probability of treatment failure on Primaquine( clinical failure and adherance)
  vptfc=0.5, # Probability of being clinical after treatment failure
  vptftr=0.27, #Probability of seeking trt if clinical, after treatment failure
  
  vdl_RDT=log10(200), # standard RDT detection limit
  vdl_micro=log10(100), # micro detection limit
  vdl_qPCR=log10(2), # standard PCR detection limit
  vmn_n=log10(5), # mean parasiteamia for sub micro
  vmn_a=log10(750), # mean parasiteamia for asym
  vmn_c=log10(5000), # mean parasiteamia for clinical
  vmn_s=log10(20000), # mean parasiteamia for severe
  vsd_n=0.75, # sd parasiteamia for sub micro
  vsd_a=1.5, # sd parasiteamia for asym
  vsd_c=1.3, # sd parasiteamia for clinical
  vsd_s=log10(9900), # sd parasiteamia for severe   
  
  t1=1,                  # Entanglement 1 - dual treatment switch
  t2=1                 # Entanglement 2 - triggering relapse from Pf infection switch

);

# ************************************************************************************* #
# Function to calculate inputs to transition rates
# ************************************************************************************* #
# define the error function
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

inputs<-function(parmal){
  # sensitivities
  sens_c_micro<-1-0.5*(1+erf((parmal$dl_micro-parmal$mn_c)/(parmal$sd_c*(2^0.5))))
  sens_c_RDT<-1-0.5*(1+erf((parmal$dl_RDT-parmal$mn_c)/(parmal$sd_c*(2^0.5))))
  vsens_c_micro<-1-0.5*(1+erf((parmal$vdl_micro-parmal$vmn_c)/(parmal$vsd_c*(2^0.5))))
  vsens_c_RDT<-1-0.5*(1+erf((parmal$vdl_RDT-parmal$vmn_c)/(parmal$vsd_c*(2^0.5))))
  
  sens_vmw<-sens_c_RDT  # test used by VMW
  sens_his<-sens_c_micro # test used by HIS
  vsens_vmw<-vsens_c_RDT  # test used by VMW
  vsens_his<-vsens_c_micro # test used by HIS
  sens_oth<-1 # test used by other
  
  # detection limits
  dl_0<-parmal$dl_qPCR
  vdl_0<-parmal$vdl_qPCR
  
  # durations of infection
  nun<-1/((1/parmal$nua)*(parmal$mn_n-dl_0)/(parmal$mn_a-parmal$mn_n)) #(r_n)
  vnun<-1/((1/parmal$vnua)*(parmal$vmn_n-vdl_0)/(parmal$vmn_a-parmal$vmn_n)) #(r_n)
  
  # ************************************************************************************* #
  # set up historical coverage over time using data and parameter values
  # ************************************************************************************* #
  
  #irs coverage
  cov_irs<-matrix(0,nrow=length(irsdat[,1]),ncol=N)
  for (n in 1:N){ for(i in 1:length(irsdat[,1])){
    cov_irs[i,n]<-min(1,as.numeric(irs_dis[i,n])/pvxy[n,6])
  }}
  
  irs_time2<-seq(irs_time[1]+dtout,(irs_time[length(irs_time)]+1/12),dtout)
  if (irs_time2[1]>startyear){
    lirs_0<-seq(startyear,(irs_time2[1]-dtout),by=dtout)
    irs_time2<-c(lirs_0,irs_time2)
    cov_irs<-rbind(matrix(0,nrow=length(lirs_0),ncol=N),cov_irs)
  }
  
  firs<-length(irs_time2)
  
  if (irs_time2[firs]<(startyear+tyears+dtout)){
    irs_time2<-c(irs_time2,irs_time2[firs]+1/12,(startyear+tyears)+2/12)
    cov_irs<-rbind(cov_irs,cov_irs[firs,],cov_irs[firs,])
  }
  c_irs<-cov_irs*parmal$eff_irs
  
  #itn coverage
  cov_itn<-matrix(0,nrow=length(itndat[,1]),ncol=N)
  
  for (n in 1:N){ for(i in 1:length(itndat[,1])){
    cov_itn[i,n]<-min(1,as.numeric(itn_dis[i,n])/pvxy[n,6])
  }}
  
  itn_time2<-seq(itn_time[1],(itn_time[length(itn_time)]+1/12),dtout)
  cov_itn<-rbind(matrix(0,nrow=1,ncol=N),cov_itn)
  if (itn_time2[1]>startyear){
    litn_0<-seq(startyear,(itn_time2[1]-dtout),by=dtout)
    itn_time2<-c(litn_0,itn_time2)
    cov_itn<-rbind(matrix(0,nrow=length(litn_0),ncol=N),cov_itn)
  }  
  
  fitn<-length(itn_time2)
  
  if (itn_time2[fitn]<(startyear+tyears+dtout)){
    itn_time2<-c(itn_time2,itn_time2[fitn]+1/12,(startyear+tyears)+2/12)
    cov_itn<-rbind(cov_itn,cov_itn[fitn,],cov_itn[fitn,])
  }
  
  c_itn<-matrix(0,nrow=length(itn_time2),ncol=N)
  c_itn[1,]<-cov_itn[1,]
  for (tt in 2:length(itn_time2))
  {
    dt_itn<-itn_time2[tt]-itn_time2[tt-1]
    c_itn[tt,]<-cov_itn[tt,]+0.5*c_itn[tt-1,]*exp(-dt_itn/parmal$hl_net)  #0.33 allows 0.8 effectiveness
  }
  
  c_itn<-parmal$eff_itn*c_itn
  
  if (itn_time2[1]>startyear){
    itn_time2<-c(startyear,itn_time2[1]-1/12,itn_time2)
    c_itn<-rbind(matrix(0,nrow=2,ncol=N),c_itn)
  }
 
  #vmw coverage
  cov_vmw<-matrix(0,nrow=length(vmwdat[,1]),ncol=N)
  for (n in 1:N){ for(i in 1:length(vmwdat[,1])){
    cov_vmw[i,n]<-min(1,as.numeric(vmw_dis[i,n])/pvxy[n,3])
  }}
  
  vmw_time2<-seq(vmw_time[1]+dtout,(vmw_time[length(vmw_time)]+1/12),dtout)
  if (vmw_time2[1]>startyear){
    lvmw_0<-seq(startyear,(vmw_time2[1]-dtout),by=dtout)
    vmw_time2<-c(lvmw_0,vmw_time2)
    cov_vmw<-rbind(matrix(0,nrow=length(lvmw_0),ncol=N),cov_vmw)
  }
  
  fvmw<-length(vmw_time2)
  
  if (vmw_time2[fvmw]<(startyear+tyears+dtout)){
    vmw_time2<-c(vmw_time2,vmw_time2[fvmw]+1/12,(startyear+tyears)+2/12)
    cov_vmw<-rbind(cov_vmw,cov_vmw[fvmw,],cov_vmw[fvmw,])
  }
  c_vmw<-cov_vmw*parmal$eff_vmw
  
   
  
  #set up El nino effect
  eln_inp<-runmed(elnino_index[,4], parmal$elneff)  #smoothing effect
  eln_t<-elnino_index[,1]+elnino_index[,2]/12 
 
  
  # ************************************************************************************* #
  # set up spatial connectivity
  # ************************************************************************************* #
  # matrix of distances in km
  dist<-matrix(0,nrow=N,ncol=N)
  for (i in 1:N){
    for (j in 1:N){
      dist[i,j]<-0.001*((pvxy[i,5]-pvxy[j,5])^2+(pvxy[i,4]-pvxy[j,4])^2)^0.5
    }
  }
  
  # matrix of connectivities 
  connect<-matrix(0,nrow=N,ncol=N)
  het<-10^(parmal$indexhet)
  # distance only
  for (i in 1:N){
    for (j in 1:N){
      connect[i,j]<-(1/(1+het[i]*dist[i,j]))/(sum(1/(1+het[i]*dist[i,])))
    }
  }
  #weight distance probabilities by population
  connect2<-connect 
  for (i in 1:N){
    for (j in 1:N){
      connect2[i,j]<-pvxy[j,6]*connect[i,j]/sum(pvxy[,6])
    }
  }
  #scaling pop-dis weights back to sum to 1 to indicate probability of movement
  connect3<-connect2
  for (i in 1:N){
    for (j in 1:N){
      connect3[i,j]<-connect2[i,j]/sum(connect2[i,])
    }
  }
  
  return(list(nun=nun,vnun=vnun, connect=connect, connect3=connect3,cov_irs=cov_irs, c_irs=c_irs, irs_time=irs_time,irs_time2=irs_time2, c_itn=c_itn, itn_time=itn_time,itn_time2=itn_time2, vmw_time2=vmw_time2, c_vmw=c_vmw, sens_vmw=sens_vmw, sens_his=sens_his,vsens_vmw=vsens_vmw, vsens_his=vsens_his, sens_oth=sens_oth, eln_inp=eln_inp, eln_t=eln_t))
}
########################################################################
#Set up matrices for malrates

# ************************************************************************************* #
# Function to calculate transition rates, given variables and parameters
# ************************************************************************************* #
malrates <- function(x, input, parmal, t,ti) {
  t_internal<-(ti-1)*dtout+t+startyear
  #Set up matrices
  seas<-c(rep(1,N)) # seasonality
  popf<-c(rep(0,N))    # population sizes  
  popv<-c(rep(0,N))    # population sizes   
  foif<-c(rep(0,N))     # forces of infection falciparum
  foiv<-c(rep(0,N))     # forces of infection vivax
  tranrate<-matrix(0,nrow=N,ncol=A)   # transition rate matrix
  irs<-c(rep(0,N))  #IRS coverage
  itn<-c(rep(0,N))  #ITN coverage
  c_vmw<-c(rep(0,N)) #VMW coverage
  veff<-c(rep(1,N)) #VMW effect decline
  tau<-c(rep(0,N)) #Treatment probabilty
  tauo<-c(rep(0,N)) #Treatment with other
  tauh<-c(rep(0,N)) #Treatment with HIS
  tauv<-c(rep(0,N)) #Treatment with VMW

  
  eln<-approx(input$eln_t,input$eln_inp,xout=t_internal)$y
  
  for (n in 1:N){
    seas[n]<-1+eln*parmal$amp[n]*cos(2*pi*(t_internal-parmal$phi[n])) # seasonal forcing signal
    popf[n]<-sum(x[varind[falpop,n]]) # list of the variable indices for the human population
    popv[n]<-sum(x[varind[vivpop,n]]) # list of the variable indices for the human population
    itn[n]<-approx(input$itn_time2,input$c_itn[,n],xout=t_internal)$y
    irs[n]<-approx(input$irs_time2,input$c_irs[,n],xout=t_internal)$y
    foif[n]<-sum(input$connect[n,]*seas[n]*((((1-itn[n])*parmal$bites)^2*parmal$prob_m*parmal$prob_h*((1-irs[n])*parmal$Pmaxf[n])/popf[n]*(parmal$zeta_n*x[varind[2,n]]+parmal$zeta_a*x[varind[3,n]]+x[varind[4,n]]+x[varind[5,n]])/sum(input$connect[n,]*popf))/((1-itn[n])*parmal$bites*parmal$prob_h*(1-irs[n])*parmal$Pmaxf[n]/popf[n]+parmal$delta_m)*(parmal$gamma_m/(parmal$gamma_m+parmal$delta_m))))
    foiv[n]<-sum(input$connect[n,]*seas[n]*((((1-itn[n])*parmal$bites)^2*parmal$prob_m*parmal$vprob_h*((1-irs[n])*parmal$Pmaxv[n])/popv[n]*(parmal$vzeta_n*x[varind[12,n]]+parmal$vzeta_a*x[varind[13,n]]+x[varind[14,n]]+x[varind[15,n]])/sum(input$connect[n,]*popv))/((1-itn[n])*parmal$bites*parmal$vprob_h*(1-irs[n])*parmal$Pmaxv[n]/popv[n]+parmal$delta_m)*(parmal$vgamma_m/(parmal$vgamma_m+parmal$delta_m))))
    
    foifa<-(1/(1/(foif[n])+1/parmal$gamma_h+1/parmal$gamma_m))
    foiva<-(1/(1/(foiv[n])+1/parmal$vgamma_h+1/parmal$vgamma_m))
    
    c_vmw[n]<-approx(input$vmw_time2,input$c_vmw[,n],xout=t_internal)$y
    tau[n]<-c_vmw[n]*input$sens_vmw*veff[n]+(1-c_vmw[n]*veff[n])*input$sens_his*parmal$eff_his+(1-c_vmw[n]*veff[n]-(1-c_vmw[n]*veff[n])*parmal$eff_his)*input$sens_oth*parmal$eff_oth #prob of treatment
    tauo[n]<-(1-c_vmw[n]*veff[n]-(1-c_vmw[n]*veff[n])*parmal$eff_his)*input$sens_oth*parmal$eff_oth #Prob treated by Other
    tauh[n]<-(1-c_vmw[n]*veff[n])*input$sens_his*parmal$eff_his #Prob treated by HIS
    tauv[n]<-c_vmw[n]*veff[n]*input$sens_vmw  #Prob treated by VMW
    
    t3=0
    if (t_internal>yrprim[n]){t3=1} #start radical cure from policy start date
    
    tranrate[n,]<-c(     
      parmal$demog*popf[n]+(1-parmal$tausev)*parmal$pmort*parmal$nuq*x[varind[5,n]]+parmal$tausev*parmal$nuq*parmal$pmortt*x[varind[5,n]], # rate of birth 1
      parmal$demog*x[varind[1,n]],       # rate of death of S                                      2
      parmal$demog*x[varind[2,n]],       # rate of death of In                                     3
      parmal$demog*x[varind[3,n]],       # rate of death of Ia                                     4
      parmal$demog*x[varind[4,n]],       # rate of death of Ic                                     5
      parmal$demog*x[varind[5,n]],       # rate of death of Is                                     6
      parmal$demog*x[varind[6,n]],       # rate of death of To                                     7
      parmal$demog*x[varind[7,n]],       # rate of death of Tv                                     8
      parmal$demog*x[varind[8,n]],       # rate of death of Th                                     9
      parmal$demog*x[varind[9,n]],       # rate of death of R                                      10
      parmal$demog*x[varind[10,n]],       # rate of death of H                                     11
      parmal$psn*(1-parmal$ps)*foifa*x[varind[1,n]], #incidence S to In               12
      (1-parmal$psn)*(1-parmal$ps)*foifa*x[varind[1,n]],  #      incidence S to Ia    13
      (1-tau[n])*parmal$ps*foifa*x[varind[1,n]],          #      incidence S to Ic    14
      tauo[n]*parmal$ps*foifa*x[varind[1,n]],      #             incidence S to To    15
      tauv[n]*parmal$ps*foifa*x[varind[1,n]],      #             incidence S to Tv    16
      tauh[n]*parmal$ps*foifa*x[varind[1,n]],      #             incidence S to Th    17
      (1-parmal$tausev)*(1-parmal$pmort)*parmal$nuq*x[varind[5,n]],             #recovery Is to Ic    18
      (1-parmal$psev)*parmal$nuc*x[varind[4,n]],      # recovery Ic to Ia                         19
      parmal$nua*x[varind[3,n]],      # recovery Ia to In                                         20
      input$nun*x[varind[2,n]],      # recovery In to R                                           21
      (1-parmal$ptf)*parmal$nut*x[varind[6,n]],      # recovery To to H                                          22
      (1-parmal$ptf)*parmal$nut*x[varind[7,n]],      # recovery Tv to H                                          23
      (1-parmal$ptf)*parmal$nut*x[varind[8,n]],      # recovery Th to H                                          24
      parmal$tausev*parmal$nuq*(1-parmal$pmortt)*x[varind[5,n]],      # recovery Is to H                            25
      (1-parmal$prn)*(1-parmal$pr)*foifa*x[varind[2,n]], # incidence In to Ia          26
      parmal$pr*(1-tau[n])*foifa*x[varind[3,n]],         # incidence Ia to Ic          27      
      parmal$pr*tauo[n]*foifa*x[varind[3,n]],         # incidence Ia to To             28      
      parmal$pr*tauv[n]*foifa*x[varind[3,n]],         # incidence Ia to Tv             29      
      parmal$pr*tauh[n]*foifa*x[varind[3,n]],         # incidence Ia to Th             30      
      parmal$psev*parmal$nuc*x[varind[4,n]],         # incidence Ic to Is              31
      parmal$pr*(1-tau[n])*foifa*x[varind[2,n]],   # incidence In to Ic                32
      parmal$pr*tauo[n]*foifa*x[varind[2,n]],      # incidence In to To                33
      parmal$pr*tauv[n]*foifa*x[varind[2,n]],      # incidence In to Tv                34
      parmal$pr*tauh[n]*foifa*x[varind[2,n]],      # incidence In to Th                35
      parmal$prn*(1-parmal$pr)*foifa*x[varind[9,n]],     #       incidence R to In     36
      (1-parmal$prn)*(1-parmal$pr)*foifa*x[varind[9,n]],  #      incidence R to Ia     37
      (1-tau[n])*parmal$pr*foifa*x[varind[9,n]],          #      incidence R to Ic     38
      tauo[n]*parmal$pr*foifa*x[varind[9,n]],      #             incidence R to To     39
      tauv[n]*parmal$pr*foifa*x[varind[9,n]],      #             incidence R to Tv     40
      tauh[n]*parmal$pr*foifa*x[varind[9,n]],      #             incidence R to Th     41
      parmal$omega*x[varind[9,n]],                 #   loss of immunity from R to S    42
      parmal$prn*(1-parmal$pr)*foifa*x[varind[10,n]],     #       incidence H to In    43
      (1-parmal$prn)*(1-parmal$pr)*foifa*x[varind[10,n]],  #      incidence H to Ia    44
      (1-tau[n])*parmal$pr*foifa*x[varind[10,n]],          #      incidence H to Ic    45
      tauo[n]*parmal$pr*foifa*x[varind[10,n]],      #             incidence H to To    46
      tauv[n]*parmal$pr*foifa*x[varind[10,n]],      #             incidence H to Tv    47
      tauh[n]*parmal$pr*foifa*x[varind[10,n]],      #             incidence H to Th    48
      parmal$chi*x[varind[10,n]],                    #      loss HRP2 H to  R          49
      (1-parmal$tausev)*parmal$pmort*parmal$nuq*x[varind[5,n]]+parmal$tausev*parmal$nuq*parmal$pmortt*x[varind[5,n]],     # death untreated+treated Is  50
      (1-parmal$ptfc)*parmal$ptf*parmal$nut*x[varind[6,n]],      # failure To to Ia                 51
      parmal$ptfc*(1-parmal$ptftr)*parmal$ptf*parmal$nut*x[varind[6,n]],      # failure To to Ic    52
      parmal$ptfc*parmal$ptftr*parmal$ptf*parmal$nut*x[varind[6,n]],      # failure To to Th        53
      (1-parmal$ptfc)*parmal$ptf*parmal$nut*x[varind[7,n]],      # failure Tv to Ia                54
      parmal$ptfc*(1-parmal$ptftr)*parmal$ptf*parmal$nut*x[varind[7,n]],      # failure Tv to Ic   55
      parmal$ptfc*parmal$ptftr*parmal$ptf*parmal$nut*x[varind[7,n]],      # failure Tv to Tv       56
      (1-parmal$ptfc)*parmal$ptf*parmal$nut*x[varind[8,n]],          # failure Th to Ia             57
      parmal$ptfc*(1-parmal$ptftr)*parmal$ptf*parmal$nut*x[varind[8,n]],  # failure Th to Ic        58
      parmal$ptfc*parmal$ptftr*parmal$ptf*parmal$nut*x[varind[8,n]],    # failure Th to Th          59
      
      #vivax
      parmal$demog*popv[n]+(1-parmal$vtausev)*parmal$vpmort*parmal$vnuq*x[varind[15,n]]+parmal$vtausev*parmal$vpmortt*parmal$vnuq*x[varind[15,n]], # rate of birth 60
      parmal$demog*x[varind[11,n]],       # rate of death of S                                      61
      parmal$demog*x[varind[12,n]],       # rate of death of In                                     62
      parmal$demog*x[varind[13,n]],       # rate of death of Ia                                     63
      parmal$demog*x[varind[14,n]],       # rate of death of Ic                                     64
      parmal$demog*x[varind[15,n]],       # rate of death of Is                                     65
      parmal$demog*x[varind[16,n]],       # rate of death of To                                     66
      parmal$demog*x[varind[17,n]],       # rate of death of Tv                                     67
      parmal$demog*x[varind[18,n]],       # rate of death of Th                                     68
      parmal$demog*x[varind[19,n]],       # rate of death of R                                      69
      parmal$demog*x[varind[20,n]],       # rate of death of L                                      70
      parmal$demog*x[varind[21,n]],       # rate of death of Togd                                   71
      parmal$demog*x[varind[22,n]],       # rate of death of Tvgd                                   72
      parmal$demog*x[varind[23,n]],       # rate of death of Thgd                                   73
      parmal$vpsn*((1-parmal$vps)*(1-parmal$propgd[n])+(1-parmal$pgd)*parmal$propgd[n])*foiva*x[varind[11,n]], #incidence S to In          74
      (1-parmal$vpsn)*((1-parmal$vps)*(1-parmal$propgd[n])+(1-parmal$pgd)*parmal$propgd[n])*foiva*x[varind[11,n]],  # incidence S to Ia    75
      (1-tau[n])*(parmal$vps*(1-parmal$propgd[n])+parmal$pgd*parmal$propgd[n])*foiva*x[varind[11,n]],              #  incidence S to Ic    76
      (1-parmal$propgd[n])*parmal$vps*tauo[n]*foiva*x[varind[11,n]],                                     #         incidence S to To    77
      (1-parmal$propgd[n])*tauv[n]*parmal$vps*foiva*x[varind[11,n]],                                     #         incidence S to Tv    78
      (1-parmal$propgd[n])*tauh[n]*parmal$vps*foiva*x[varind[11,n]],                                     #         incidence S to Th    79
      parmal$propgd[n]*parmal$pgd*tauo[n]*foiva*x[varind[11,n]],                                       #         incidence S to Togd    80
      parmal$propgd[n]*parmal$pgd*tauv[n]*foiva*x[varind[11,n]],                                       #         incidence S to Tvgd    81
      parmal$propgd[n]*parmal$pgd*tauh[n]*foiva*x[varind[11,n]],                                       #         incidence S to Thgd    82
      (1-parmal$tausev)*(1-parmal$vpmort)*parmal$vnuq*x[varind[15,n]],            #          recovery Is to Ic    83
      (1-parmal$vpsev)*parmal$vnuc*x[varind[14,n]],                      # recovery Ic to Ia                      84
      parmal$vnua*x[varind[13,n]],                                      # recovery Ia to In                       85
      (1-parmal$vph)*input$vnun*x[varind[12,n]],                         # recovery In to R                       86
      parmal$vph*input$vnun*x[varind[12,n]],                             # recovery In to L                       87
      (1-t3)*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[16,n]],        # recovery To to R         88
      (1-t3)*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[17,n]],        # recovery Tv to R         89
      (1-t3)*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[18,n]],        # recovery Th to R         90
      (1-t3)*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[21,n]],      # recovery Togd to R         91
      (1-t3)*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[22,n]],      # recovery Tvgd to R         92
      (1-t3)*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[23,n]],      # recovery Thgd to R         93 
      (1-parmal$vph)*parmal$vtausev*(1-parmal$vpmortt)*parmal$vnuq*x[varind[15,n]],  # recovery Is to R    94
      (1-t3)*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[16,n]],            # recovery To to L         95
      (1-t3)*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[17,n]],           # recovery Tv to L          96
      (1-t3)*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[18,n]],            # recovery Th to L         97
      (1-t3)*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[21,n]],          # recovery Togd to L         98
      (1-t3)*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[22,n]],        #  recovery Tvgd to L          99
      (1-t3)*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[23,n]],          # recovery Thgd to L         100
      parmal$vph*parmal$vtausev*(1-parmal$vpmortt)*parmal$vnuq*x[varind[15,n]], # recovery Is to L         101
      (1-parmal$vprn)*(1-parmal$vpr)*foiva*x[varind[12,n]],                   # incidence In to Ia         102
      parmal$vpr*(1-tau[n])*foiva*x[varind[12,n]],                     # incidence In to Ic                103
      (1-parmal$propgd[n])*parmal$vpr*tauo[n]*foiva*x[varind[12,n]],      # incidence In to To             104
      (1-parmal$propgd[n])*parmal$vpr*tauv[n]*foiva*x[varind[12,n]],      # incidence In to Tv             105
      (1-parmal$propgd[n])*parmal$vpr*tauh[n]*foiva*x[varind[12,n]],      # incidence In to Th             106
      parmal$propgd[n]*parmal$vpr*tauo[n]*foiva*x[varind[12,n]],        # incidence In to Togd             107
      parmal$propgd[n]*parmal$vpr*tauv[n]*foiva*x[varind[12,n]],        # incidence In to Tvgd             108
      parmal$propgd[n]*parmal$vpr*tauh[n]*foiva*x[varind[12,n]],        # incidence In to Thgd             109
      parmal$vpr*(1-tau[n])*foiva*x[varind[13,n]],                        # incidence Ia to Ic             110  
      (1-parmal$propgd[n])*parmal$vpr*tauo[n]*foiva*x[varind[13,n]],         # incidence Ia to To          111     
      (1-parmal$propgd[n])*parmal$vpr*tauv[n]*foiva*x[varind[13,n]],         # incidence Ia to Tv          112     
      (1-parmal$propgd[n])*parmal$vpr*tauh[n]*foiva*x[varind[13,n]],         # incidence Ia to Th          113
      parmal$propgd[n]*parmal$vpr*tauo[n]*foiva*x[varind[13,n]],           # incidence Ia to Togd          114     
      parmal$propgd[n]*parmal$vpr*tauv[n]*foiva*x[varind[13,n]],           # incidence Ia to Tvgd          115      
      parmal$propgd[n]*parmal$vpr*tauh[n]*foiva*x[varind[13,n]],           # incidence Ia to Thgd          116 
      parmal$vpsev*parmal$vnuc*x[varind[14,n]],                           # incidence Ic to Is             117
      parmal$vprn*(1-parmal$vpr)*foiva*x[varind[19,n]],                        #  incidence R to In        118
      (1-parmal$vprn)*(1-parmal$vpr)*foiva*x[varind[19,n]],                    # incidence R to Ia         119
      (1-tau[n])*parmal$vpr*foiva*x[varind[19,n]],                            #  incidence R to Ic         120
      (1-parmal$propgd[n])*tauo[n]*parmal$vpr*foiva*x[varind[19,n]],      #         incidence R to To      121
      (1-parmal$propgd[n])*tauv[n]*parmal$vpr*foiva*x[varind[19,n]],      #         incidence R to Tv      122
      (1-parmal$propgd[n])*tauh[n]*parmal$vpr*foiva*x[varind[19,n]],      #         incidence R to Th      123
      parmal$propgd[n]*tauo[n]*parmal$vpr*foiva*x[varind[19,n]],        #         incidence R to Togd      124
      parmal$propgd[n]*tauv[n]*parmal$vpr*foiva*x[varind[19,n]],        #         incidence R to Tvgd      125
      parmal$propgd[n]*tauh[n]*parmal$vpr*foiva*x[varind[19,n]],        #         incidence R to Thgd      126
      parmal$vomega*x[varind[19,n]],                                   # loss of immunity from R to S      127
      parmal$vprn*(1-parmal$vpr)*(foiva)*x[varind[20,n]],                 #       incidence L to In        128
      (1-parmal$vprn)*(1-parmal$vpr)*(foiva)*x[varind[20,n]],              #      incidence L to Ia        129
      (1-tau[n])*parmal$vpr*(foiva)*x[varind[20,n]],                       #      incidence L to Ic        130
      (1-parmal$propgd[n])*tauo[n]*parmal$vpr*(foiva)*x[varind[20,n]], #             incidence L to To     131
      (1-parmal$propgd[n])*tauv[n]*parmal$vpr*(foiva)*x[varind[20,n]], #             incidence L to Tv     132
      (1-parmal$propgd[n])*tauh[n]*parmal$vpr*(foiva)*x[varind[20,n]], #             incidence L to Th     133
      parmal$propgd[n]*tauo[n]*parmal$vpr*(foiva)*x[varind[20,n]],   #             incidence L to Togd     134
      parmal$propgd[n]*tauv[n]*parmal$vpr*(foiva)*x[varind[20,n]],   #             incidence L to Tvgd     135
      parmal$propgd[n]*tauh[n]*parmal$vpr*(foiva)*x[varind[20,n]],   #             incidence L to Thgd     136
      
      (1-parmal$t2)*parmal$vprn*(1-parmal$vpr)*parmal$vprel*parmal$vrel*x[varind[20,n]],                 #       relapse L to In   137
      (1-parmal$t2)*(1-parmal$vprn)*(1-parmal$vpr)*parmal$vprel*parmal$vrel*x[varind[20,n]],             #       relapse L to Ia   138
      (1-parmal$t2)*(1-tau[n])*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[20,n]],                      #       relapse L to Ic   139
      (1-parmal$t2)*(1-parmal$propgd[n])*tauo[n]*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[20,n]],    #       relapse L to To   140
      (1-parmal$t2)*(1-parmal$propgd[n])*tauv[n]*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[20,n]],    #       relapse L to Tv   141
      (1-parmal$t2)*(1-parmal$propgd[n])*tauh[n]*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[20,n]],    #       relapse L to Th   142
      (1-parmal$t2)*parmal$propgd[n]*tauo[n]*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[20,n]],        #       relapse L to Togd 143
      (1-parmal$t2)*parmal$propgd[n]*tauv[n]*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[20,n]],        #       relapse L to Tvgd 144
      (1-parmal$t2)*parmal$propgd[n]*tauh[n]*parmal$vpr*parmal$vprel*parmal$vrel*x[varind[20,n]],        #       relapse L to Thgd 145  
      parmal$vkappa*x[varind[20,n]],                                                                     #       death of hypnozoites L to  S     146                                                                         
      (1-parmal$vtausev)*parmal$vpmort*parmal$vnuq*x[varind[15,n]]+parmal$vtausev*parmal$vpmortt*parmal$vnuq*x[varind[15,n]],   # death untreated+treated Is  147 
      
      (1-t3)*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[16,n]],                      # failed trt To to Ia         148
      (1-t3)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[16,n]],        # failed trt To to Ic         149
      (1-t3)*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[16,n]],            # failed trt To to Th         150
      
      (1-t3)*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[17,n]],                      # failed trt Tv to Ia         151
      (1-t3)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[17,n]],        # failed trt Tv to Ic         152
      (1-t3)*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[17,n]],            # failed trt Tv to Tv         153
      
      (1-t3)*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[18,n]],                      # failed trt Th to Ia         154
      (1-t3)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[18,n]],        # failed trt Th to Ic         155
      (1-t3)*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[18,n]],            # failed trt Th to Th         156
      
      (1-t3)*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[21,n]],                    # failed trt Togd to Ia         157
      (1-t3)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[21,n]],      # failed trt Togd to Ic         158
      (1-t3)*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[21,n]],          # failed trt Togd to Thgd       159
      
      (1-t3)*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[22,n]],                    # failed trt Tvgd to Ia         160
      (1-t3)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[22,n]],      # failed trt Tvgd to Ic         161
      (1-t3)*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[22,n]],        # failed trt Tvgd to Tvgd         162
      
      (1-t3)*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[23,n]],                    # failed trt Thgd to Ia         163 
      (1-t3)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[23,n]],      # failed trt Thgd to Ic         164
      (1-t3)*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[23,n]],        # failed trt Thgd to Thgd         165 
      
      #Entanglements (dual treat)
      parmal$t1*(1-parmal$ptf)*parmal$nut*x[varind[2,n]]*(x[varind[16,n]]+x[varind[17,n]]+x[varind[18,n]]+x[varind[21,n]]+x[varind[22,n]]+x[varind[23,n]])/popv[n],                     # dual trt Fal In to H  166
      parmal$t1*(1-parmal$ptf)*parmal$nut*x[varind[3,n]]*(x[varind[16,n]]+x[varind[17,n]]+x[varind[18,n]]+x[varind[21,n]]+x[varind[22,n]]+x[varind[23,n]])/popv[n],                     # dual trt Fal Ia to H  167
      parmal$t1*(1-parmal$ptf)*parmal$nut*x[varind[4,n]]*(x[varind[16,n]]+x[varind[17,n]]+x[varind[18,n]]+x[varind[21,n]]+x[varind[22,n]]+x[varind[23,n]])/popv[n],                     # dual trt Fal Ic to H  168
      parmal$t1*(1-parmal$ptf)*parmal$nut*x[varind[5,n]]*(x[varind[16,n]]+x[varind[17,n]]+x[varind[18,n]]+x[varind[21,n]]+x[varind[22,n]]+x[varind[23,n]])/popv[n],                     # dual trt Fal Is to H  169
      parmal$t1*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[12,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],      # dual trt Viv In to R  170
      parmal$t1*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[13,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],      # dual trt Viv Ia to R  171
      parmal$t1*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[14,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],      # dual trt Viv Ic to R  172
      parmal$t1*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[15,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],      # dual trt Viv Is to R  173
      parmal$t1*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[12,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],          # dual trt Viv In to L  174
      parmal$t1*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[13,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],          # dual trt Viv Ia to L  175
      parmal$t1*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[14,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],          # dual trt Viv Ic to L  176
      parmal$t1*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[15,n]]*(x[varind[6,n]]+x[varind[7,n]]+x[varind[8,n]])/popf[n],           # dual trt Viv Is to L 177
  
      t3*(1-parmal$mask)*(1-parmal$vphprim)*(1-parmal$vptfp)*parmal$vnup*x[varind[16,n]],                               # primaquine recovery To to R     178
      t3*(1-parmal$mask)*(1-parmal$vphprim)*(1-parmal$vptfp)*parmal$vnup*x[varind[17,n]],                              #primaquine recovery Tv to R       179
      t3*(1-parmal$mask)*(1-parmal$vphprim)*(1-parmal$vptfp)*parmal$vnup*x[varind[18,n]],                             #primaquine recovery Th to R        180
      t3*(1-parmal$mask)*parmal$vphprim*(1-parmal$vptfp)*parmal$vnup*x[varind[16,n]],                               # primaquine recovery To to L         181
      t3*(1-parmal$mask)*parmal$vphprim*(1-parmal$vptfp)*parmal$vnup*x[varind[17,n]],                              #primaquine recovery Tv to L           182
      t3*(1-parmal$mask)*parmal$vphprim*(1-parmal$vptfp)*parmal$vnup*x[varind[18,n]],                             #primaquine recovery Th to L            183
      
      t3*(1-parmal$mask)*(parmal$sensgd*parmal$vnut*(1-parmal$vph)*(1-parmal$vptf))*x[varind[21,n]],                  # test + act recovery Togd to R     184     
      t3*(1-parmal$mask)*(parmal$sensgd*parmal$vnut*(1-parmal$vph)*(1-parmal$vptf))*x[varind[22,n]],                 #test + act recovery Tvgd to R       185
      t3*(1-parmal$mask)*(parmal$sensgd*parmal$vnut*(1-parmal$vph)*(1-parmal$vptf))*x[varind[23,n]],                #test + act recovery Thgd to R        186
      t3*(1-parmal$mask)*(parmal$sensgd*parmal$vnut*parmal$vph*(1-parmal$vptf))*x[varind[21,n]],             # test + act recovery of GDef Togd to L      187
      t3*(1-parmal$mask)*(parmal$sensgd*parmal$vnut*parmal$vph*(1-parmal$vptf))*x[varind[22,n]],            #test + act recovery of GDef Tvgd to L        188
      t3*(1-parmal$mask)*(parmal$sensgd*parmal$vnut*parmal$vph*(1-parmal$vptf))*x[varind[23,n]],           #test + act recovery of GDef Thgd to L         189
      
      t3*(1-parmal$mask)*((1-parmal$sensgd)*(1-parmal$vptfp)*(1-parmal$vphprim)*parmal$vnup)*x[varind[21,n]],        # test - primaquine recovery Togd to R     190
      t3*(1-parmal$mask)*((1-parmal$sensgd)*(1-parmal$vptfp)*(1-parmal$vphprim)*parmal$vnup)*x[varind[22,n]],      # test - primaquine recovery Tvgd to R       191
      t3*(1-parmal$mask)*((1-parmal$sensgd)*(1-parmal$vptfp)*(1-parmal$vphprim)*parmal$vnup)*x[varind[23,n]],      #test - primaquine recovery Thgd to R        192
      t3*(1-parmal$mask)*((1-parmal$sensgd)*parmal$vphprim*(1-parmal$vptfp)*parmal$vnup)*x[varind[21,n]],    # test - primaquine recovery of GDef Togd to L     193
      t3*(1-parmal$mask)*((1-parmal$sensgd)*parmal$vphprim*(1-parmal$vptfp)*parmal$vnup)*x[varind[22,n]],   #test - primaquine recovery of GDef Tvgd to L       194
      t3*(1-parmal$mask)*((1-parmal$sensgd)*parmal$vphprim*(1-parmal$vptfp)*parmal$vnup)*x[varind[23,n]],  #test - primaquine recovery of GDef Thgd to L        195
      
      t3*parmal$mask*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[16,n]],                             # masked act recovery To to R     196
      t3*parmal$mask*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[17,n]],                            #masked act recovery Tv to R       197
      t3*parmal$mask*(1-parmal$vph)*(1-parmal$vptf)*parmal$vnut*x[varind[18,n]],                           #masked act recovery Th to R        198
      t3*parmal$mask*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[16,n]],                             # masked act recovery To to L         199
      t3*parmal$mask*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[17,n]],                            #masked act recovery Tv to L           200
      t3*parmal$mask*parmal$vph*(1-parmal$vptf)*parmal$vnut*x[varind[18,n]],                           #masked act recovery Th to L            201
      t3*parmal$mask*parmal$vnut*(1-parmal$vph)*(1-parmal$vptf)*x[varind[21,n]],                           # masked act recovery Togd to R     202
      t3*parmal$mask*parmal$vnut*(1-parmal$vph)*(1-parmal$vptf)*x[varind[22,n]],                          #masked act recovery Tvgd to R       203
      t3*parmal$mask*parmal$vnut*(1-parmal$vph)*(1-parmal$vptf)*x[varind[23,n]],                         #masked act recovery Thgd to R        204
      t3*parmal$mask*parmal$vnut*parmal$vph*(1-parmal$vptf)*x[varind[21,n]],                       # masked act recovery of GDef Togd to L     205
      t3*parmal$mask*parmal$vnut*parmal$vph*(1-parmal$vptf)*x[varind[22,n]],                     #masked act recovery of GDef Tvgd to L        206
      t3*parmal$mask*parmal$vnut*parmal$vph*(1-parmal$vptf)*x[varind[23,n]],                   #masked act recovery of GDef Thgd to L          207
      
      parmal$t2*parmal$vprn*(1-parmal$vpr)*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],                       #       relapse L to In +triggering  208
      parmal$t2*(1-parmal$vprn)*(1-parmal$vpr)*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],                   #      relapse L to Ia  +triggering  209
      parmal$t2*(1-tau[n])*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],                           #      relapse L to Ic+triggering     210
      parmal$t2*(1-parmal$propgd[n])*tauo[n]*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],     #          relapse L to To +triggering    211
      parmal$t2*(1-parmal$propgd[n])*tauv[n]*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],   #           relapse L to Tv+triggering      212
      parmal$t2*(1-parmal$propgd[n])*tauh[n]*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],     #          relapse L to Th  +triggering   213
      parmal$t2*parmal$propgd[n]*tauo[n]*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],       #          relapse L to Togd +triggering    214
      parmal$t2*parmal$propgd[n]*tauv[n]*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],      #          relapse L to Tvgd+triggering      215
      parmal$t2*parmal$propgd[n]*tauh[n]*parmal$vpr*(parmal$vprel*((popf[n]-x[varind[9,n]]-x[varind[10,n]])/popf[n])+parmal$vincprel*((x[varind[9,n]]+x[varind[10,n]])/popf[n]))*parmal$vrel*x[varind[20,n]],       #           relapse L to Thgd  +triggering  216  
      
      #Failed treatments
      
      t3*(1-parmal$mask)*(1-parmal$vptfc)*parmal$vptfp*parmal$vnup*x[varind[16,n]],                            # primaquine failed trt To to Ia     217
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptfp*parmal$vnup*x[varind[16,n]],              # primaquine failed trt To to Ic     218
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*parmal$vptfp*parmal$vnup*x[varind[16,n]],                  # primaquine failed trt To to Th     219
    
      t3*(1-parmal$mask)*(1-parmal$vptfc)*parmal$vptfp*parmal$vnup*x[varind[17,n]],                           #primaquine failed trt Tv to Ia       220
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptfp*parmal$vnup*x[varind[17,n]],             #primaquine failed trt Tv to Ic       221
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*parmal$vptfp*parmal$vnup*x[varind[17,n]],                 #primaquine failed trt Tv to Tv       222
      
      t3*(1-parmal$mask)*(1-parmal$vptfc)*parmal$vptfp*parmal$vnup*x[varind[18,n]],                           #primaquine failed trt Th to Ia       223
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*parmal$vptfp*parmal$vnup*x[varind[18,n]],             #primaquine failed trt Th to Ic       224
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*parmal$vptfp*parmal$vnup*x[varind[18,n]],                 #primaquine failed trt Th to Th       225
      
      t3*parmal$mask*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[16,n]],                                 # masked act failed trt To to Ia     226
      t3*parmal$mask*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[16,n]],                   # masked act failed trt To to Ic     227
      t3*parmal$mask*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[16,n]],                       # masked act failed trt To to Th     228
      
      t3*parmal$mask*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[17,n]],                                 #masked act failed trt Tv to Ia      229
      t3*parmal$mask*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[17,n]],                   #masked act failed trt Tv to Ic      230
      t3*parmal$mask*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[17,n]],                       #masked act failed trt Tv to Tv      231
      
      t3*parmal$mask*(1-parmal$vptfc)*parmal$vptf*parmal$vnut*x[varind[18,n]],                                 #masked act failed trt Th to Ia      232
      t3*parmal$mask*(1-parmal$vptftr)*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[18,n]],                   #masked act failed trt Th to Ic      233
      t3*parmal$mask*parmal$vptftr*parmal$vptfc*parmal$vptf*parmal$vnut*x[varind[18,n]],                       #masked act failed trt Th to Th      234
      
      t3*parmal$mask*(1-parmal$vptfc)*parmal$vnut*parmal$vptf*x[varind[21,n]],                                 # masked act failed trt Togd to Ia   235
      t3*parmal$mask*(1-parmal$vptftr)*parmal$vptfc*parmal$vnut*parmal$vptf*x[varind[21,n]],                   # masked act failed trt Togd to Ic   236
      t3*parmal$mask*parmal$vptftr*parmal$vptfc*parmal$vnut*parmal$vptf*x[varind[21,n]],                     # masked act failed trt Togd to Thgd   237
      
      t3*parmal$mask*(1-parmal$vptfc)*parmal$vnut*parmal$vptf*x[varind[22,n]],                                #masked act failed trt Tvgd to Ia     238
      t3*parmal$mask*(1-parmal$vptftr)*parmal$vptfc*parmal$vnut*parmal$vptf*x[varind[22,n]],                  #masked act failed trt Tvgd to Ic     239
      t3*parmal$mask*parmal$vptftr*parmal$vptfc*parmal$vnut*parmal$vptf*x[varind[22,n]],                    #masked act failed trt Tvgd to Tvgd     240
      
      t3*parmal$mask*(1-parmal$vptfc)*parmal$vnut*parmal$vptf*x[varind[23,n]],                               #masked act failed trt Thgd to Ia      241
      t3*parmal$mask*(1-parmal$vptftr)*parmal$vptfc*parmal$vnut*parmal$vptf*x[varind[23,n]],                 #masked act failed trt Thgd to Ic      242
      t3*parmal$mask*parmal$vptftr*parmal$vptfc*parmal$vnut*parmal$vptf*x[varind[23,n]],                   #masked act failed trt Thgd to Thgd      243
      
      t3*(1-parmal$mask)*(1-parmal$vptfc)*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[21,n]],                  # test + act failed trt Togd to Ia     244     
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[21,n]],    # test + act failed trt Togd to Ic     245     
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[21,n]],        # test + act failed trt Togd to Thgd   246     
      
      t3*(1-parmal$mask)*(1-parmal$vptfc)*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[22,n]],                  #test + act failed trt Tvgd to Ia      247
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[22,n]],    #test + act failed trt Tvgd to Ic      248
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[22,n]],        #test + act failed trt Tvgd to Tvgd    249
      
      t3*(1-parmal$mask)*(1-parmal$vptfc)*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[23,n]],                  #test + act failed trt Thgd to Ia      250
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[23,n]],    #test + act failed trt Thgd to Ic      251
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*parmal$sensgd*parmal$vnut*parmal$vptf*x[varind[23,n]],        #test + act failed trt Thgd to Thgd    252
      
      t3*(1-parmal$mask)*(1-parmal$vptfc)*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[21,n]],               # test - primaquine failed trt Togd to Ia       253
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[21,n]], # test - primaquine failed trt Togd to Ic       254
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[21,n]],     # test - primaquine failed trt Togd to Thgd     255
      
      t3*(1-parmal$mask)*(1-parmal$vptfc)*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[22,n]],               # test - primaquine failed trt Tvgd to Ia       256
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[22,n]], # test - primaquine failed trt Tvgd to Ic       257
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[22,n]],     # test - primaquine failed trt Tvgd to Tvgd     258
      
      t3*(1-parmal$mask)*(1-parmal$vptfc)*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[23,n]],               #test - primaquine failed trt Thgd to Ia        259
      t3*(1-parmal$mask)*(1-parmal$vptftr)*parmal$vptfc*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[23,n]], #test - primaquine failed trt Thgd to Ic        260
      t3*(1-parmal$mask)*parmal$vptftr*parmal$vptfc*(1-parmal$sensgd)*parmal$vptfp*parmal$vnup*x[varind[23,n]]      #test - primaquine failed trt Thgd to Thgd      261
      
      )
  }
  return(c(t(tranrate)))
}

# POST PROCESSING function
postproc <- function(parpro,out,tran) {
  # sensitivities
  sens_n_micro<-1-0.5*(1+erf((parpro$dl_micro-parpro$mn_n)/(parpro$sd_n*(2^0.5))))
  sens_n_RDT<-1-0.5*(1+erf((parpro$dl_RDT-parpro$mn_n)/(parpro$sd_n*(2^0.5))))
  sens_n_qPCR<-1-0.5*(1+erf((parpro$dl_qPCR-parpro$mn_n)/(parpro$sd_n*(2^0.5))))
  sens_a_micro<-1-0.5*(1+erf((parpro$dl_micro-parpro$mn_a)/(parpro$sd_a*(2^0.5))))
  sens_a_RDT<-1-0.5*(1+erf((parpro$dl_RDT-parpro$mn_a)/(parpro$sd_a*(2^0.5))))
  sens_a_qPCR<-1-0.5*(1+erf((parpro$dl_qPCR-parpro$mn_a)/(parpro$sd_a*(2^0.5))))
  sens_c_micro<-1-0.5*(1+erf((parpro$dl_micro-parpro$mn_c)/(parpro$sd_c*(2^0.5))))
  sens_c_RDT<-1-0.5*(1+erf((parpro$dl_RDT-parpro$mn_c)/(parpro$sd_c*(2^0.5))))
  sens_c_qPCR<-1-0.5*(1+erf((parpro$dl_qPCR-parpro$mn_c)/(parpro$sd_c*(2^0.5))))
  sens_s_micro<-1-0.5*(1+erf((parpro$dl_micro-parpro$mn_s)/(parpro$sd_s*(2^0.5))))
  sens_s_RDT<-1-0.5*(1+erf((parpro$dl_RDT-parpro$mn_s)/(parpro$sd_s*(2^0.5))))
  sens_s_qPCR<-1-0.5*(1+erf((parpro$dl_qPCR-parpro$mn_s)/(parpro$sd_s*(2^0.5))))
  sens_H_micro<-0
  sens_H_RDT<-1
  sens_H_qPCR<-0
  vsens_n_micro<-1-0.5*(1+erf((parpro$vdl_micro-parpro$vmn_n)/(parpro$vsd_n*(2^0.5))))
  vsens_n_RDT<-1-0.5*(1+erf((parpro$vdl_RDT-parpro$vmn_n)/(parpro$vsd_n*(2^0.5))))
  vsens_n_qPCR<-1-0.5*(1+erf((parpro$vdl_qPCR-parpro$vmn_n)/(parpro$vsd_n*(2^0.5))))
  vsens_a_micro<-1-0.5*(1+erf((parpro$vdl_micro-parpro$vmn_a)/(parpro$vsd_a*(2^0.5))))
  vsens_a_RDT<-1-0.5*(1+erf((parpro$vdl_RDT-parpro$vmn_a)/(parpro$vsd_a*(2^0.5))))
  vsens_a_qPCR<-1-0.5*(1+erf((parpro$vdl_qPCR-parpro$vmn_a)/(parpro$vsd_a*(2^0.5))))
  vsens_c_micro<-1-0.5*(1+erf((parpro$vdl_micro-parpro$vmn_c)/(parpro$vsd_c*(2^0.5))))
  vsens_c_RDT<-1-0.5*(1+erf((parpro$vdl_RDT-parpro$vmn_c)/(parpro$vsd_c*(2^0.5))))
  vsens_c_qPCR<-1-0.5*(1+erf((parpro$vdl_qPCR-parpro$vmn_c)/(parpro$vsd_c*(2^0.5))))
  vsens_s_micro<-1-0.5*(1+erf((parpro$vdl_micro-parpro$vmn_s)/(parpro$vsd_s*(2^0.5))))
  vsens_s_RDT<-1-0.5*(1+erf((parpro$vdl_RDT-parpro$vmn_s)/(parpro$vsd_s*(2^0.5))))
  vsens_s_qPCR<-1-0.5*(1+erf((parpro$vdl_qPCR-parpro$vmn_s)/(parpro$vsd_s*(2^0.5))))
  
  sens_vmw<-sens_c_RDT  # test used by VMW
  sens_his<-sens_c_micro # test used by HIS
  vsens_vmw<-vsens_c_RDT  # test used by VMW
  vsens_his<-vsens_c_micro # test used by HIS
  sens_oth<-1 # test used by other
  
  
  # ************************************************************************************* #
  # for outputting the  time series for each patch
  # ************************************************************************************* #
  
  # VMW cases
  vmw_predv<-vmw_predf<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    vmw_predf[,n]<-rowSums(tran[,c(traind[c(16,29,34,40,47,56),n])])/12
    vmw_predv[,n]<-rowSums(tran[,c(traind[c(78,81,105,108,112,115,122,125,132,135,141,144,153,162,212,215,222,231,240,249,258),n])])/12
  }
  
  # HIS cases (Uncomplicated cases accessing Health system (incidence +relapses+failures))
  his_predv<-his_predf<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    his_predf[,n]<-rowSums(tran[,c(traind[c(17,30,35,41,48,53,59),n])])/12
    his_predv[,n]<-rowSums(tran[,c(traind[c(79,82,106,109,113,116,123,126,133,136,142,145,150,156,159,165,213,216,219,225,228,234,237,243,246,252,255,261),n])])/12
  }
  
  his_predmix<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    his_predmix[,n]<-his_predv[,n]*(0.5*out[,1+(3+(n-1)*B)]+out[,1+(4+(n-1)*B)]+out[,1+(5+(n-1)*B)]+out[,1+(6+(n-1)*B)]+out[,1+(7+(n-1)*B)]+out[,1+(8+(n-1)*B)])/sum(out[length(out[,1]),(1+falpop+(n-1)*B)]) + his_predf[,n]*(0.5*out[,1+(13+(n-1)*B)]+out[,1+(14+(n-1)*B)]+out[,1+(15+(n-1)*B)]+out[,1+(16+(n-1)*B)]+out[,1+(17+(n-1)*B)]+out[,1+(21+(n-1)*B)]+out[,1+(22+(n-1)*B)])/sum(out[length(out[,1]),(1+(vivpop+(n-1)*B))])
  }
  
  his_predv1<-his_predf1<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    his_predf1[,n]<-his_predf[,n]-his_predmix[,n]
    his_predv1[,n]<-his_predv[,n]-his_predmix[,n]
  }
  
  # Treated cases (His+other)(Uncomplicated cases accessing Health system (incidence +relapses+failures))
  trt_predf<-trt_predv<-trt_predv1<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    trt_predf[,n]<-rowSums(tran[,c(traind[c(15:17,28:30,33:35,39:41,46:48,53,56,59),n])])/12
    trt_predv[,n]<-rowSums(tran[,c(traind[c(77:82,104:109,111:116,121:126,131:136,140:145,150,153,156,159,162,165,211:216, 219,222, 225,228,231,234,237,240,243,246,249,252,255,258,261),n])])/12
  }
  
  trt_predmix<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    trt_predmix[,n]<-trt_predv[,n]*(0.5*out[,1+(3+(n-1)*B)]+out[,1+(4+(n-1)*B)]+out[,1+(5+(n-1)*B)]+out[,1+(6+(n-1)*B)]+out[,1+(7+(n-1)*B)]+out[,1+(8+(n-1)*B)])/sum(out[length(out[,1]),(1+falpop+(n-1)*B)]) + trt_predf[,n]*(0.5*out[,1+(13+(n-1)*B)])/sum(out[length(out[,1]),(1+(vivpop+(n-1)*B))])
  }
  
  trt_predv1<- trt_predv-trt_predmix
  
  # Reported severe cases (successfully treated +deaths (trt+untrt))
  severe_predf<-severe_predv<-severe_predv1<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    severe_predf[,n]<-(tran[,traind[25,n]]+tran[,traind[50,n]])/12
    severe_predv[,n]<-(tran[,traind[94,n]]+tran[,traind[101,n]]+tran[,traind[147,n]])/12
  }

  severe_predmix<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    severe_predmix[,n]<-severe_predv[,n]*(0.5*out[,1+(3+(n-1)*B)]+out[,1+(4+(n-1)*B)]+out[,1+(5+(n-1)*B)]+out[,1+(6+(n-1)*B)]+out[,1+(7+(n-1)*B)]+out[,1+(8+(n-1)*B)])/sum(out[length(out[,1]),(1+falpop+(n-1)*B)]) + severe_predf[,n]*(0.5*out[,1+(13+(n-1)*B)]+out[,1+(14+(n-1)*B)]+out[,1+(16+(n-1)*B)]+out[,1+(17+(n-1)*B)]+out[,1+(18+(n-1)*B)]+out[,1+(21+(n-1)*B)]+out[,1+(22+(n-1)*B)]+out[,1+(23+(n-1)*B)])/sum(out[length(out[,1]),(1+(vivpop+(n-1)*B))])
  }
  
  severe_predv1=severe_predv - severe_predmix

  # Fatalities
  fatal_pred<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    fatal_pred[,n]<-(tran[,traind[50,n]]+tran[,traind[147,n]])/12
  }
  
  # True clinical burden (all clinical cases+failures - uncomp+severe; trt+untreated)
  totclininc_predf<-totclininc_predv<-totclininc_predv1<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    totclininc_predf[,n]<-rowSums(tran[,c(traind[c(14:17,27:35,38:41,45:48,52,53,55,56,58,59),n])])/12
    totclininc_predv[,n]<-rowSums(tran[,c(traind[c(76:82,103:117,120:126,130:136,139:145,149,150,152,153,155,156,158,159,161,162,164,165,210:216,218,219,221,222,224,225,227,228,230,231,233,234,236,237,239,240,242,243,245,246,248,249,251,252,254,255,257,258,260,261),n])])/12
  }
  
  totclininc_predmix<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    totclininc_predmix[,n]<-totclininc_predv[,n]*(0.5*out[,1+(3+(n-1)*B)]+out[,1+(4+(n-1)*B)]+out[,1+(5+(n-1)*B)]+out[,1+(6+(n-1)*B)]+out[,1+(7+(n-1)*B)]+out[,1+(8+(n-1)*B)])/sum(out[length(out[,1]),(1+falpop+(n-1)*B)]) + totclininc_predf[,n]*(0.5*out[,1+(13+(n-1)*B)])/sum(out[length(out[,1]),(1+(vivpop+(n-1)*B))])
  }  
  
  totclininc_predv1=totclininc_predv- totclininc_predmix
  
  # True incidence (all uncomplicated +severe+failures)
  totinc_predf<-totinc_predv<-totinc_predv1<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    totinc_predf[,n]<-rowSums(tran[,c(traind[c(12:17,26:41,43:48,51:59),n])])/12
    totinc_predv[,n]<-rowSums(tran[,c(traind[c(74:82,102:126,128:145,148:165,208:261),n])])/12
  }
 
  totinc_predmix<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    totinc_predmix[,n]<-totinc_predv[,n]*(1-(out[,1+(1+(n-1)*B)]+out[,1+(9+(n-1)*B)]+out[,1+(10+(n-1)*B)])/sum(out[length(out[,1]),(1+falpop+(n-1)*B)])) 
  }  
 
  totinc_predv1=totinc_predv - totinc_predmix
    
  # ************************************************************************************* #
  # for predicting prevalence with different tests
  # ************************************************************************************* #
  # true prevalence
  totalinf_pred<-rowSums(out[,c(varind[2:5,],varind[12:15,])+1])
  
  prevalence_predf<-prevalence_predv<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    prevalence_predf[,n]<-100*rowSums(out[,c(varind[2:5,n])+1])/rowSums(out[,(varind[falpop,n]+1)])
    prevalence_predv[,n]<-100*rowSums(out[,c(varind[12:15,n])+1])/rowSums(out[,(varind[vivpop,n]+1)])
  }
  prev_micro_predf<-prev_micro_predv<-matrix(0,nrow=length(out[,1]),ncol=N)
  prev_RDT_predf<- prev_RDT_predv<-matrix(0,nrow=length(out[,1]),ncol=N)
  prev_qPCR_predf<-prev_qPCR_predv<-matrix(0,nrow=length(out[,1]),ncol=N)
  for (n in 1:N){
    # prevalence by Micro
    prev_micro_predf[,n]<-100*(sens_n_micro*out[,c(varind[2,n])+1]+sens_a_micro*out[,c(varind[3,n])+1]+sens_c_micro*out[,c(varind[4,n])+1]+sens_s_micro*out[,c(varind[5,n])+1]+sens_H_micro*out[,c(varind[10,n])+1])/rowSums(out[,(varind[falpop,n]+1)])
    prev_micro_predv[,n]<-100*(vsens_n_micro*out[,c(varind[12,n])+1]+vsens_a_micro*out[,c(varind[13,n])+1]+vsens_c_micro*out[,c(varind[14,n])+1]+vsens_s_micro*out[,c(varind[15,n])+1])/rowSums(out[,(varind[vivpop,n]+1)])
    
    # prevalence by RDT
    prev_RDT_predf[,n]<-100*(sens_n_RDT*out[,c(varind[2,n])+1]+sens_a_RDT*out[,c(varind[3,n])+1]+sens_c_RDT*out[,c(varind[4,n])+1]+sens_s_RDT*out[,c(varind[5,n])+1]+sens_H_RDT*out[,c(varind[10,n])+1])/rowSums(out[,(varind[falpop,n]+1)])
    prev_RDT_predv[,n]<-100*(vsens_n_RDT*out[,c(varind[12,n])+1]+vsens_a_RDT*out[,c(varind[13,n])+1]+vsens_c_RDT*out[,c(varind[14,n])+1]+vsens_s_RDT*out[,c(varind[15,n])+1])/rowSums(out[,(varind[vivpop,n]+1)])
    
    # prevalence by qPCR
    prev_qPCR_predf[,n]<-100*(sens_n_qPCR*out[,c(varind[2,n])+1]+sens_a_qPCR*out[,c(varind[3,n])+1]+sens_c_qPCR*out[,c(varind[4,n])+1]+sens_s_qPCR*out[,c(varind[5,n])+1]+sens_H_qPCR*out[,c(varind[10,n])+1])/rowSums(out[,(varind[falpop,n]+1)])
    prev_qPCR_predv[,n]<-100*(vsens_n_qPCR*out[,c(varind[12,n])+1]+vsens_a_qPCR*out[,c(varind[13,n])+1]+vsens_c_qPCR*out[,c(varind[14,n])+1]+vsens_s_qPCR*out[,c(varind[15,n])+1])/rowSums(out[,(varind[vivpop,n]+1)])
  }
  
  return(cbind(vmw_predf,    #1
               vmw_predv,     #2
               his_predf,    #3
               his_predv,    #4
               his_predf1,  #5
               his_predv1,   #6
               his_predmix,    #7
               fatal_pred,    #8
               prevalence_predf,    #9
               prevalence_predv,    #10
               prev_micro_predf,    #11
               prev_micro_predv,    #12
               prev_RDT_predf,    #13
               prev_RDT_predv,    #14
               prev_qPCR_predf,    #15
               prev_qPCR_predv,    #16
               totclininc_predf,    #17
               totclininc_predv1,    #18
               trt_predf, #19
               trt_predv1, #20
               severe_predf, #21
               severe_predv, #22
               totinc_predf, #23
               totinc_predv1, #24
               totalinf_pred
  ))
  
}


ti<-1

# ************************************************************************************* #
# ************************************************************************************* #
# ************************************************************************************* #
# ODE SOLVER
# ************************************************************************************* #
# ************************************************************************************* #

epiModel<-function(t,state, parode,input) 
{
  with(as.list(c(state, parode)),
       {
         
         #   # ************************************************************************************* #
         #   # define variables
         #   # ************************************************************************************* #
         
         Z=state
         
         # rates of change
         ti<-1
         transit<-malrates(Z[1:V],input,parode,Z[V+1],ti)  
         
         if (sum(is.na(transit))>0)  {
           stop("transit NA   ",Z[V+1], "                                      ", 
                as.data.frame(transit))
         }
         
         eq<-rep(0.0, V)
         
         eq<-EQ(L, N, eq, transit,transitionsiu1,transitionsiu2,transitionsiv1,transitionsiv2)
         
         eq[V+1]<-1
         
         dZ<-eq
         
         # return the rate of change
         list(c(dZ))
       }
  ) 
}

# ************************************************************************************* #
# RUN FUNCTION
# ************************************************************************************* #
# ************************************************************************************* #
# Adjust starting conditions to best available data


run<-function(parrun){
  
  # ************************************************************************************* #
  # define initial conditions
  inp<-inputs(parrun)
  initcondrun<-NULL
  for(n in 1:N){
    initcondrun<-c(initcondrun,c(S=1/10*pvxy[n,6], In=1/10*pvxy[n,6], Ia=1/10*pvxy[n,6], Ic=1/10*pvxy[n,6], Is=1/10*pvxy[n,6], To=1/10*pvxy[n,6], Tv=1/10*pvxy[n,6], Th=1/10*pvxy[n,6], R=1/10*pvxy[n,6], H=1/10*pvxy[n,6],Sv=1/13*pvxy[n,6], Inv=1/13*pvxy[n,6], Iav=1/13*pvxy[n,6], Icv=1/13*pvxy[n,6], Isv=1/13*pvxy[n,6], Tov=1/13*pvxy[n,6], Tvv=1/13*pvxy[n,6],Thv=1/13*pvxy[n,6], Rv=1/13*pvxy[n,6], Lv=1/13*pvxy[n,6],Togd=1/13*pvxy[n,6], Tvgd=1/13*pvxy[n,6], Thgd=1/13*pvxy[n,6]))
  }
  initoderun<-round(initcondrun)
  staterun <- c(initoderun,0)
  ti<-1  
  inp<-inputs(parrun)
  transitrun <- malrates(initoderun,inp,parrun,0,ti)
  
  
   # SOLVE THE ODEs and get output
  timesrun <- seq(0, tyears, by = dtout) # Model run time
  #Solve ODE
  outoderun <- ode(y = staterun, times = timesrun, func = epiModel, parms = parrun, method  = "lsoda", input=inp)
  # Compute transitions at each time step
  tranoderun<-matrix(0,nrow=length(outoderun[,1]),ncol=length(transitions))
  for (ti in 1:(tsteps+1)){
    tranoderun[ti,]<-t(malrates(outoderun[ti,2:(1+V)],inp, parrun,0,ti))
  }
  #Compute outputs
  ppoutrun<-postproc(parrun,outoderun,tranoderun)
  modeltimes<-outoderun[,1]+startyear
  
  vmw_predf_run<-ppoutrun[,1:N]
  vmw_predv_run<-ppoutrun[,(N+1):(2*N)]
  his_predf_run<-ppoutrun[,(2*N+1):(3*N)]
  his_predv_run<-ppoutrun[,(3*N+1):(4*N)]
  his_predf1_run<-ppoutrun[,(4*N+1):(5*N)]
  his_predv1_run<-ppoutrun[,(5*N+1):(6*N)]
  his_predmix_run<-ppoutrun[,(6*N+1):(7*N)]
  fatal_pred_run<-ppoutrun[,(7*N+1):(8*N)]
  totclininc_predf_run<-ppoutrun[,(16*N+1):(17*N)]
  totclininc_predv1_run<-ppoutrun[,(17*N+1):(18*N)]
  trt_predf_run<-ppoutrun[,(18*N+1):(19*N)]
  trt_predv1_run<-ppoutrun[,(19*N+1):(20*N)]
  severe_predf_run<-ppoutrun[,(20*N+1):(21*N)]
  severe_predv1_run<-ppoutrun[,(21*N+1):(22*N)]
  
  #Monthly Output object
  outputM<-list(vmw_predf_run, 
                vmw_predv_run,
                his_predf_run,
                his_predv_run,
                his_predf1_run,
                his_predv1_run,
                his_predmix_run,
                fatal_pred_run,
                totclininc_predf_run,
                totclininc_predv1_run,
                trt_predf_run,
                trt_predv1_run,
                severe_predf_run,
                severe_predv1_run)
  
  #Annual Summation
  sumtimes<-startyear:(startyear+tyears)
  failedf_run_year<-failedv_run_year<-relapse_run_year<-severef_run_year<-severev_run_year<-trtf_run_year<-trtv_run_year<-hisf_run_year<-hisv_run_year<-totclinincf_run_year<-totclinincv_run_year<-fatal_run_year<-matrix(0,nrow=length(sumtimes)-1,ncol=N)
  for (i in 1:(length(sumtimes)-1)){
    for (n in 1:N){
      hisf_run_year[i,n]<-sum(his_predf_run[which(floor(modeltimes)==((startyear-1)+i)),n])  # runting to Pf+mixed
      hisv_run_year[i,n]<-sum(his_predv1_run[which(floor(modeltimes)==((startyear-1)+i)),n]) #runting to Pv mono
      trtf_run_year[i,n]<-sum(trt_predf_run[which(floor(modeltimes)==((startyear-1)+i)),n])  # runting to Pf+mixed all treated
      trtv_run_year[i,n]<-sum(trt_predv1_run[which(floor(modeltimes)==((startyear-1)+i)),n])  # runting to Pf mono
      fatal_run_year[i,n]<-sum(fatal_pred_run[which(floor(modeltimes)==((startyear-1)+i)),n])
      totclinincf_run_year[i,n]<-sum(totclininc_predf_run[which(floor(modeltimes)==((startyear-1)+i)),n])
      totclinincv_run_year[i,n]<-sum(totclininc_predv1_run[which(floor(modeltimes)==((startyear-1)+i)),n])
    }  }
  
  #Annual Output object
  outputA<-list(hisf_run_year,
                hisv_run_year,
                trtf_run_year,
                trtv_run_year,
                fatal_run_year,
                totclinincf_run_year,
                totclinincv_run_year)
  
  return(list(outputM,outputA)) 
}

Sys.time()
run1<-run(pars)
Sys.time()

#Time parameters
t_yr<-startyear:(startyear+tyears-1)
t_mon<-seq(startyear,(startyear+tyears), 1/12)
#Plot Viewing Window
y<-6:20
m<-62:241



#Plot 1: Annual HIS Incidence(Pf and Pv)
layout(matrix(c(1:24), 6, 4, byrow = TRUE))
par(oma = c(2, 3, 2, 2)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(3, 2, 2, 2)) # make the plots be closer together
for (i in 1:N){
  plot(t_yr[y],run1[[2]][[3]][y,i],type="l",col="red",ylim=c(0,max(run1[[2]][[3]][y,i],run1[[2]][[4]][y,i])),main=pvxy[i,1])
  lines(t_yr[y],run1[[2]][[4]][y,i],type="l",col="blue")
 }
plot.new()
plot.new()
mtext("Annual HIS Incidence(Pf (red) and Pv(blue))", side=3, line=1, outer=T)

#Plot 2: Monthly HIS Incidence(Pf and Pv)
layout(matrix(c(1:24), 6, 4, byrow = TRUE))
par(oma = c(2, 3, 2, 2)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(3, 2, 2, 2)) # make the plots be closer together
for (i in 1:N){
  plot(t_mon[m],run1[[1]][[3]][m,i],type="l",col="red",ylim=c(0,max(run1[[1]][[3]][m,i],run1[[1]][[4]][m,i])),main=pvxy[i,1])
  lines(t_mon[m],run1[[1]][[4]][m,i],type="l",col="blue")
}
plot.new()
plot.new()
mtext("Monthly HIS Incidence(Pf (red) and Pv(blue))", side=3, line=1, outer=T)

#Plot 3: Monthly Reported Fatalities(Pf and Pv)
layout(matrix(c(1:24), 6, 4, byrow = TRUE))
par(oma = c(2, 3, 2, 2)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(3, 2, 2, 2)) # make the plots be closer together
for (i in 1:N){
  plot(t_mon[m],run1[[1]][[8]][m,i],type="l",col="green4",ylim=c(0,max(run1[[1]][[8]][m,i])),main=pvxy[i,1])
}
plot.new()
plot.new()
mtext("Monthly Fatalities", side=3, line=1, outer=T)

#Plot 4: Annual Reported Fatalities(Pf and Pv)
layout(matrix(c(1:24), 6, 4, byrow = TRUE))
par(oma = c(2, 3, 2, 2)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(3, 2, 2, 2)) # make the plots be closer together
for (i in 1:N){
  plot(t_yr[y],run1[[2]][[5]][y,i],type="l",col="green4",ylim=c(0,max(run1[[2]][[5]][y,i])),main=pvxy[i,1])
}
plot.new()
plot.new()
mtext("Annual Fatalities", side=3, line=1, outer=T)

#Plot 5: Clinical and Reported Treated (P. falciparum)
layout(matrix(c(1:24), 6, 4, byrow = TRUE))
par(oma = c(2, 3, 2, 2)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(3, 2, 2, 2)) # make the plots be closer together
for (i in 1:N){
  plot(t_mon[m],run1[[1]][[9]][m,i],type="l",col="red",ylim=c(0,max(run1[[1]][[9]][m,i],run1[[1]][[11]][m,i])),main=pvxy[i,1])
  polygon(c(t_mon[m],rev(t_mon[m])),c(run1[[1]][[9]][m,i],rev(run1[[1]][[11]][m,i])),col="LightPink",border="red")
}
plot.new()
plot.new()
mtext("P. falciparum: Clinical (upper) and Treated (lower) Cases", side=3, line=1, outer=T)

#Plot 6: Clinical and Reported Treated (P. vivax)
layout(matrix(c(1:24), 6, 4, byrow = TRUE))
par(oma = c(2, 3, 2, 2)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(3, 2, 2, 2)) # make the plots be closer together
for (i in 1:N){
  plot(t_mon[m],run1[[1]][[10]][m,i],type="l",col="blue",ylim=c(0,max(run1[[1]][[10]][m,i],run1[[1]][[12]][m,i])),main=pvxy[i,1])
  polygon(c(t_mon[m],rev(t_mon[m])),c(run1[[1]][[10]][m,i],rev(run1[[1]][[12]][m,i])),col="LightCyan",border="blue")
}
plot.new()
plot.new()
mtext("P. vivax: Clinical (upper) and Treated (lower) Cases", side=3, line=1, outer=T)


