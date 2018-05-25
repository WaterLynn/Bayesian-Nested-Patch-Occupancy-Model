#May 22, 2018
#Author: Lynn Waterhouse
#Contact information: lwaterhouse@ucsd.edu
##################################################

#The purpose of this code is to generate tag detection data for a system
#You must have JAGS installed in order to run this code


#set working directory
setwd("~/PROJECT_SalmonPitPatch/ToyExample")

#load in packages
#install.packages("VGAM")
#install.packages("MASS")
#install.packages("gtools")
#install.packages("rjags")
#install.packages("runjags")
#install.packages("R2jags")
#install.packages("gdata")
library(VGAM) # for fitting dirichlet
library(MASS) # for fitting beta
library(gtools)
library(gdata)
#library(rjags) #for running jags
library(runjags) #for running jags
library(R2jags) #for running jags


###########################################################
######### Set Up Simulation for Fitting Simulated Data ####
###########################################################
###########################################################
#Establish how many final populations are in the system
#This system has 3 main branches (A,B,C) past the entry point (with the option to stay at the entry point)
npops <- 4 #for A, B, C, stay (top level) #main level population

#With branch A where are the stopping points
npops_ALow <-3 #for lower A- A1, up, or stay
npops_Aup <- 2 # for upper A- A2, stay

#Within branch B where are the stopping points
npops_BLow <- 2 # for B- up, stay
npops_Bup <- 3 #for upper B- B1, B2, stay

#Within branch C where are the stopping points
npops_CLow <- 3 # for C- C1, up, stay
npops_Cup <- 3 #for upper C- C2, C3 or stay

####################################################
####################################################
#In addition to knowing the stopping points you will need to know which locations 
#have double versus single arrays
#In this example all mainstems are single arrays (1 at each of lowA, lowB, lowC)
#All other locations are doubles

#detection arrays in A (7)
#"lowA","A1_down","A1_up","Aup_down","Aup_up","A2_down","A2_up"

#detection arrays in B (7)
#"lowB","Bup_down","Bup_up","B1_down","B1_up","B2_down","B2_up"

#detection arrays in C (9)
#"lowC","C1_down","C1_up","Cup_down","Cup_up","C2_down","C2_up",
#                       "C3_down","C3_up"

####################################################
####################################################

set.seed(12345) #make this reproducible
nsim<- 10   #How many simulations are we running

#Number of fish
Nfish<-1000  	#Total fish in the system
nfish<-100 	#Number of tagged fish in system


###############################
#Jags parameters for later on
my.chains = 4		  #number of chains (4)
my.iter = 2000 	#number of iterations in each chain (25,000)
my.burnin = 500 	#burnin per chain (5000)
my.thin=4		    #thining (40)
#At the end how many MCMC results will we have
my.mcmc.count<-my.chains*(my.iter-my.burnin)/my.thin 
################################

#Create result matrix for population estimates to go into from simulation code
#Here number of fish is out of the # of fish tagged
#Number of fish winding up in each of the main branches for each MCMC run of each simulation
#4 main choices (A, B, C, stay)
res.mat<- array(rep(NA, nsim*4*my.mcmc.count), dim=c(nsim, 4,my.mcmc.count))

#create 3 result matrices for all the pop estimates
#Matrix of where the fish go within A (by location)
#A LOW--proportion going to mainstem (1), A1 (2), up(3)
#A UP--proportion to mainstem (1),  A2 (2)
res.A.mat<-array(rep(NA, nsim*5*my.mcmc.count), dim=c(nsim, 5,my.mcmc.count))

#Matrix of where the fish go within B (by location)
# B low- There are 3 bins, Mainstem (1), up(2)
# B Up-- mainstem(1), B1 (2), B3 (3)
res.B.mat<-array(rep(NA, nsim*5*my.mcmc.count), dim=c(nsim, 5,my.mcmc.count))

#Matrix of where the fish go within C (by location)
#C low--- mainstem (1), c1 (2), up (3)
#C up-- mainstem (1), C2 (2), and C3 (3)
res.C.mat<-array(rep(NA, nsim*6*my.mcmc.count), dim=c(nsim, 6,my.mcmc.count))

#note to expand from #tagged to # in system, you would need to: 
# divide these results by nfish (to get proportions)
# multiply by Nfish (to get total number is system)

##############################################################
##############################################################
# These matrices will be overwritten each simulation run
#create matrix to store original data - 
#movement of fish in system
#by detection array
orig.A.mat<- array(rep(NA, nsim*7), dim=c(nsim, 7))
orig.B.mat<- array(rep(NA, nsim*7), dim=c(nsim, 7))
orig.C.mat<- array(rep(NA, nsim*9), dim=c(nsim, 9))

#create matrix to store original data with observation error- 
#observed movement of fish in system
#by detection array
obs.A.mat<- array(rep(NA, nsim*7), dim=c(nsim, 7))
obs.B.mat<- array(rep(NA, nsim*7), dim=c(nsim, 7))
obs.C.mat<- array(rep(NA, nsim*9), dim=c(nsim, 9))

#create matrix to store individual fish movement in the system
#by detection array
Abase<-as.data.frame(matrix(rep(NA,nfish*7),nrow=nfish,ncol=7))
Bbase<-as.data.frame(matrix(rep(NA,nfish*7),nrow=nfish,ncol=7))
Cbase<-as.data.frame(matrix(rep(NA,nfish*9),nrow=nfish,ncol=9))

#################################################################################
#################################################################################
#Run the simulation code
for (j in 1:nsim) {
        #get the header titles
  Abase<-as.data.frame(matrix(rep(NA,nfish*7),nrow=nfish,ncol=7))
  Bbase<-as.data.frame(matrix(rep(NA,nfish*7),nrow=nfish,ncol=7))
  Cbase<-as.data.frame(matrix(rep(NA,nfish*9),nrow=nfish,ncol=9))
    colnames(Abase)<-c("lowA","A1_down","A1_up","Aup_down","Aup_up","A2_down","A2_up") 
    colnames(Bbase)<-c("lowB","Bup_down","Bup_up","B1_down","B1_up","B2_down","B2_up")
    colnames(Cbase)<-c("lowC","C1_down","C1_up","Cup_down","Cup_up","C2_down","C2_up",
                       "C3_down","C3_up")
    
    # Set up simulation parameters
    N.simfish <- nfish #how many tagged fish?
    
    #create matrices with the fish locations by array WITHOUT obs. error
    A.sim <- as.data.frame(matrix(data = 0, nrow = N.simfish, ncol = 7)) # matrix of final destinations for all tagged fish in A
    colnames(A.sim)=colnames(Abase) # match the column names of real and simulated data
    B.sim <- as.data.frame(matrix(data = 0, nrow = N.simfish, ncol = 7)) # matrix of final destinations for all tagged fish in B
    colnames(B.sim)=colnames(Bbase) # match the column names of real and simulated data
    C.sim <- as.data.frame(matrix(data = 0, nrow = N.simfish, ncol = 9)) # matrix of final destinations for all tagged fish in C
    colnames(C.sim)=colnames(Cbase) # match the column names of real and simulated data
    
    # dulplicate above matrices to add in obs. error based on detection efficiency estimates
    A.sim.obs <- A.sim
    B.sim.obs <- B.sim
    C.sim.obs <- C.sim
    
    #Detection probabilities for our system
    #Here all detection probabilities are uniform between a value and 1
    #There are three options for that value
    a<-.7 #high detection (.7 to 1)
    b<-.2 #low .2 to 1 is possible
    c<-0  #low 0 to 1 is possible
    #Instead of using a,b,c you could replace these with actual values, or other distributions
 
    lowA_median<-runif(1,a,1)
    A1_up_median<-runif(1,b,1)
    A1_down_median<-runif(1,b,1)
    Aup_up_median<-runif(1,c,1)
    Aup_down_median<-runif(1,a,1)
    A2_up_median<-runif(1,a,1)
    A2_down_median<-runif(1,a,1)
    lowB_median<-runif(1,b,1)
    Bup_up_median<-runif(1,b,1)
    Bup_down_median<-runif(1,a,1)
    B1_up_median<-runif(1,a,1)
    B1_down_median<-runif(1,c,1)
    B2_up_median<-runif(1,a,1)
    B2_down_median<-runif(1,b,1)
    lowC_median<-runif(1,a,1)
    C1_up_median<-runif(1,a,1)
    C1_down_median<-runif(1,a,1)
    Cup_up_median<-runif(1,c,1)
    Cup_down_median<-runif(1,c,1)
    C2_up_median<-runif(1,a,1)
    C2_down_median<-runif(1,a,1)
    C3_up_median<-runif(1,b,1)
    C3_down_median<-runif(1,b,1)
    
    #Movement probabillities to main branches (A, B, C, stay)
    main_p_now<-c(.7,.02,.13,.15) #needs to be 4 probabilities, summing to 1 (A, B, C, stay)
    
    #MAIN TRIB SPLITS ################################################################################
    # Wenache vs. entiat vs. wells  vs. nowhere
    main_p_a <- main_p_now
    
    # Now put fish into the system
    main_p_n <- rmultinom(1, N.simfish, main_p_a)
    
    #now fill in fish location matrix without & with observation error
    A.sim$lowA[1:main_p_n[1]] <- 1 #fish going to A
    A.sim.obs$lowA<-rbinom(c(1:N.simfish),1,A.sim$lowA*lowA_median) #observation error
    
    B.sim$lowB[(main_p_n[1]+1):(main_p_n[1]+main_p_n[2])] <- 1 #fish going to B
    B.sim.obs$lowB<-rbinom(c(1:N.simfish),1,B.sim$lowB*lowB_median) #observation error
    
    C.sim$lowC[(main_p_n[1]+main_p_n[2]+1):(main_p_n[1]+main_p_n[2]+main_p_n[3])] <- 1 #fish going to C
    C.sim.obs$lowC<-rbinom(c(1:N.simfish),1,C.sim$lowC*lowC_median) #observation error
    
    #A LOW #####################################################################################
    # proportion going to mainstem (1), A2 (2), upstream (3)
    lowA_a <- c(.1,.2,.7) #3 probs, sum to 1
    
    # Now put fish into the system
    lowA_n <- rmultinom(1, main_p_n[1], lowA_a)
    
    #now fill in fish location matrix without & with observation error
    # note that random binomial draws are done with rbinom(n, size, prob)
    
    A.sim$A1_down[1:lowA_n[2]] 			<- 1 #A1_down
    A.sim.obs$A1_down<-rbinom(c(1:N.simfish),1,A.sim$A1_down*A1_down_median) #observation error
    A.sim$A1_up[1:lowA_n[2]] 			<- 1 #A1_down
    A.sim.obs$A1_up<-rbinom(c(1:N.simfish),1,A.sim$A1_up*A1_up_median) #observation error
    
    A.sim$Aup_down[(lowA_n[2]+1):(lowA_n[2]+lowA_n[3])] 	<- 1 #Aup_down
    A.sim.obs$Aup_down<-rbinom(c(1:N.simfish),1,A.sim$Aup_down*Aup_down_median) #observation error
    A.sim$Aup_up[(lowA_n[2]+1):(lowA_n[2]+lowA_n[3])] 		<- 1 #Aup_up
    A.sim.obs$Aup_up<-rbinom(c(1:N.simfish),1,A.sim$Aup_up*Aup_up_median) #observation error
    
    #the ID num of the first fish going to upper A
    num.Aup<-(lowA_n[2]+1)
    
    #Upper A#####################################################################################
    # proportion of fish going to mainstem (1),  A2 (2)
    AUp_a <- c(.3,.7) #2 probabilities, sum to 1
    
    # Now put fish into the system
    AUp_n <- rmultinom(1, lowA_n[3], AUp_a)
    
    #now fill in fish location matrix without & with observation error
    if(AUp_n[2]>0){
      A.sim$A2_down[num.Aup:(num.Aup+AUp_n[2]-1)] 			<- 1 #A2_down
      A.sim.obs$A2_down<-rbinom(c(1:N.simfish),1,A.sim$A2_down*A2_down_median) #observation error
      A.sim$A2_up[num.Aup:(num.Aup+AUp_n[2]-1)] 			<- 1 #A2_up
      A.sim.obs$A2_up<-rbinom(c(1:N.simfish),1,A.sim$A2_up*A2_up_median) #observation error
    }
    
    #B Branch #################################################################
    # (1) mainstem, (2) up river
    
    lowB_a <- c(.2,.8) #2 probs, sum 1
    
    # Now put fish into the system
    lowB_n <- rmultinom(1, main_p_n[2], lowB_a)
    
    #now fill in fish location matrix without & with observation error
    B.sim$Bup_down[(main_p_n[1]+1):(main_p_n[1]+lowB_n[2])] 			<- 1 #Bup_down
    B.sim.obs$Bup_down<-rbinom(c(1:N.simfish),1,B.sim$Bup_down*Bup_down_median) #observation error
    B.sim$Bup_up[(main_p_n[1]+1):(main_p_n[1]+lowB_n[2])] 				<- 1 #Bup_up
    B.sim.obs$Bup_up<-rbinom(c(1:N.simfish),1,B.sim$Bup_up*Bup_up_median) #observation error
    
    #the ID num of the first fish going to upper B
    num.Bup<-(main_p_n[1]+1)
    
    #Upper B#####################################################################################
    # proportion of fish going to mainstem (1),  B1 (2), and B2(3)
    BUp_a <- c(.2,.5,.3) #3 probabilities, sum to 1
    
    # Now put fish into the system
    BUp_n <- rmultinom(1, lowB_n[2], BUp_a)
    
    #now fill in fish location matrix without & with observation error
    #the if statements are to avoid counting errors when no fish go places
    if(BUp_n[2]>0){
      B.sim$B1_down[num.Bup:(num.Bup+BUp_n[2]-1)] 			<- 1 #B1_down
      B.sim.obs$B1_down<-rbinom(c(1:N.simfish),1,B.sim$B1_down*B1_down_median) #observation error
      B.sim$B1_up[num.Bup:(num.Bup+BUp_n[2]-1)] 			<- 1 #B1_up
      B.sim.obs$B1_up<-rbinom(c(1:N.simfish),1,B.sim$B1_up*B1_up_median) #observation error
    }
    if(BUp_n[3]>0){
      if(BUp_n[2]>0){
        B.sim$B2_down[(num.Bup+BUp_n[2]):(num.Bup+BUp_n[2]+BUp_n[3]-1)] 			<- 1 #B2_down
        B.sim.obs$B2_down<-rbinom(c(1:N.simfish),1,B.sim$B2_down*B2_down_median) #observation error
        B.sim$B2_up[(num.Bup+BUp_n[2]):(num.Bup+BUp_n[2]+BUp_n[3]-1)] 			<- 1 #B2_up
        B.sim.obs$B2_up<-rbinom(c(1:N.simfish),1,B.sim$B2_up*B2_up_median) #observation error
      }else if(BUp_n[2]==0){
        B.sim$B2_down[(num.Bup):(num.Bup+BUp_n[3]-1)] 			<- 1 #B2_down
        B.sim.obs$B2_down<-rbinom(c(1:N.simfish),1,B.sim$B2_down*B2_down_median) #observation error
        B.sim$B2_up[(num.Bup):(num.Bup+BUp_n[3]-1)] 			<- 1 #B2_up
        B.sim.obs$B2_up<-rbinom(c(1:N.simfish),1,B.sim$B2_up*B2_up_median) #observation error
      }
      
    }
    
    #C LOW #####################################################################################
    # proportion going to mainstem (1), C1 (2), upstream (3)
    lowC_a <- c(.1,.2,.7) #3 probs, sum to 1
    
    # Now put fish into the system
    lowC_n <- rmultinom(1, main_p_n[3], lowC_a)
    
    #now fill in fish location matrix without & with observation error
    # note that random binomial draws are done with rbinom(n, size, prob)
    C.sim$C1_down[(main_p_n[1]+main_p_n[2]+1):(main_p_n[1]+main_p_n[2]+lowC_n[2])] <- 1 #C1_down
    C.sim.obs$C1_down<-rbinom(c(1:N.simfish),1,C.sim$C1_down*C1_down_median) #observation error
    C.sim$C1_up[(main_p_n[1]+main_p_n[2]+1):(main_p_n[1]+main_p_n[2]+lowC_n[2])]<- 1 #C1_down
    C.sim.obs$C1_up<-rbinom(c(1:N.simfish),1,C.sim$C1_up*C1_up_median) #observation error
    
    C.sim$Cup_down[(main_p_n[1]+main_p_n[2]+lowC_n[2]+1):
                     (main_p_n[1]+main_p_n[2]+lowC_n[2]+lowC_n[3])]	<- 1 #Cup_down
    C.sim.obs$Cup_down<-rbinom(c(1:N.simfish),1,C.sim$Cup_down*Cup_down_median) #observation error
    C.sim$Cup_up[(main_p_n[1]+main_p_n[2]+lowC_n[2]+1):
                   (main_p_n[1]+main_p_n[2]+lowC_n[2]+lowC_n[3])]		<- 1 #Cup_up
    C.sim.obs$Cup_up<-rbinom(c(1:N.simfish),1,C.sim$Cup_up*Cup_up_median) #observation error
    
    #the ID num of the first fish going to upper C
    num.Cup<-(main_p_n[1]+main_p_n[2]+lowC_n[2]+1)
    
    #Upper C#####################################################################################
    # proportion of fish going to mainstem (1),  C2 (2), and C3(3)
    CUp_a <- c(.2,.5,.3) #3 probabilities, sum to 1
    
    # Now put fish into the system
    CUp_n <- rmultinom(1, lowC_n[3], CUp_a)
    
    #now fill in fish location matrix without & with observation error
    # the if statements are to avoid counting errors when no fish move in the system
    if(CUp_n[2]>0){
      C.sim$C2_down[num.Cup:(num.Cup+CUp_n[2]-1)] 			<- 1 #C2_down
      C.sim.obs$C2_down<-rbinom(c(1:N.simfish),1,C.sim$C2_down*C2_down_median) #observation error
      C.sim$C2_up[num.Cup:(num.Cup+CUp_n[2]-1)] 			<- 1 #C2_up
      C.sim.obs$C2_up<-rbinom(c(1:N.simfish),1,C.sim$C2_up*C2_up_median) #observation error
    }
    if(CUp_n[3]>0){
      C.sim$C3_down[(num.Cup+CUp_n[2]):(num.Cup+CUp_n[2]+CUp_n[3]-1)] 			<- 1 #C3_down
      C.sim.obs$C3_down<-rbinom(c(1:N.simfish),1,C.sim$C3_down*C3_down_median) #observation error
      C.sim$C3_up[(num.Cup+CUp_n[2]):(num.Cup+CUp_n[2]+CUp_n[3]-1)] 			<- 1 #C3_up
      C.sim.obs$C3_up<-rbinom(c(1:N.simfish),1,C.sim$C3_up*C3_up_median) #observation error
    }
    
    #### End of Data Simulation #################
    #############################################
    #with no observation error
    #Bbase<-as.matrix(B.sim) 
    #Cbase<-as.matrix(C.sim)
    #Abase<-as.matrix(A.sim)
    
    #with observation error
    Bbase<-as.matrix(B.sim.obs) 
    Cbase<-as.matrix(C.sim.obs)
    Abase<-as.matrix(A.sim.obs)
    
    #####################################
    ## Save the Simulated Data results###
    #####################################
    
    #sum over all the fish - to get branch totals
    #for original data
    orig.B.mat[j,] <-colSums(B.sim)
    orig.A.mat[j,] <-colSums(A.sim)
    orig.C.mat[j,] <-colSums(C.sim)
    #data with observation error
    obs.B.mat[j,] <-colSums(B.sim.obs)
    obs.A.mat[j,] <-colSums(A.sim.obs)
    obs.C.mat[j,] <-colSums(C.sim.obs)
    
    #file names for storing data- for original data
    filenameA.fake <- paste(format(Sys.time(), "%Y.%m.%d"),"simdata.A.mat",j, sep=" ")
    filenameB.fake <- paste(format(Sys.time(), "%Y.%m.%d"),"simdata.B.mat",j, sep=" ")
    filenameC.fake <- paste(format(Sys.time(), "%Y.%m.%d"),"simdata.C.mat",j, sep=" ")
    
    #file names for storing data- for original data with observation error
    filenameA.fake.obs <- paste(format(Sys.time(), "%Y.%m.%d"),"simdata.obserror.A.mat",j, sep=" ")
    filenameB.fake.obs <- paste(format(Sys.time(), "%Y.%m.%d"),"simdata.obserror.B.mat",j, sep=" ")
    filenameC.fake.obs <- paste(format(Sys.time(), "%Y.%m.%d"),"simdata.obserror.C.mat",j, sep=" ")
    
    write.csv(orig.A.mat[j,],file=filenameA.fake)
    write.csv(orig.B.mat[j,],file=filenameB.fake)
    write.csv(orig.C.mat[j,],file=filenameC.fake)
    
    write.csv(obs.A.mat[j,],file=filenameA.fake.obs)
    write.csv(obs.B.mat[j,],file=filenameB.fake.obs)
    write.csv(obs.C.mat[j,],file=filenameC.fake.obs)
    
    ###########################################################
    ###########################################################
    ######### Set Up Simulation for Fitting Simulated Data ####
    ###########################################################
    ###########################################################
    
    
    npops <- 4 #for A, B, C, stay (top level)
    npops_ALow <-3 #for lower A - mainstem, A1, upstream
    npops_Aup <- 2 # for A up-- mainstem, A2
    npops_BLow <- 2 # for B --mainstem, Bup
    npops_Bup <- 3 #for B up-- mainstem, B1, B2
    npops_CLow <- 3 # for C-- Mainstem, C1, up
    npops_Cup <- 3 #for Cup--- mainstem, C2, C3
    
    alpha<-rep(1,10) #general purpose vector of ones (arbitrarily large, then pull subsets in model for Dirichlet priors)
    zero_vec<-rep(0,10) #general purpose vector of zeros (arbitrarily large, then pull subsets in model for "switching" matrices)
    
    ###################################################################
    # Inits for JAGS
    ###################################################################
    #first lets create inits matrices
    a_init = array(0, dim=c(nfish, npops))
    a_ALow_init = array(0, dim=c(nfish, npops_ALow+1))
    a_Aup_init = array(0, dim=c(nfish, npops_Aup+1))
    a_BLow_init = array(0,dim=c(nfish, npops_BLow+1))
    a_Bup_init = array(0,dim=c(nfish, npops_Bup+1))
    a_CLow_init = array(0, dim=c(nfish, npops_CLow+1))
    a_Cup_init= array(0, dim=c(nfish, npops_Cup+1))
    
    a_init[,1]=apply(Abase,1,max) #A detects
    a_init[,2]=apply(Bbase,1,max) #B detects
    a_init[,3]=apply(Cbase,1,max) #C detects
    a_init[,4]=abs(apply(a_init,1,max)-1) #not seen anywhere
    
    ##colnames(Abase)<-c("lowA","A1_down","A1_up","Aup_down","Aup_up","A2_down","A2_up") 
    #mainstem (1), A1 (2), Aup (3), and not there (4)
    a_ALow_init[,4]= abs(a_init[,1] -1) #not in A
    a_ALow_init[,2]= apply(Abase[,2:3],1,max) #heard in A1
    a_ALow_init[,3]= apply(Abase[,4:7],1,max) #heard upstream
    a_ALow_init[,1]= a_init[,1] - apply(a_ALow_init[,2:3],1,max)#heard in A Mainstem only
    
    #mainstem (1),  A2 (2), and not there (3)
    a_Aup_init[,2]= apply(Abase[,6:7],1,max) #heard in A2
    a_Aup_init[,1]= a_ALow_init[,3] - apply(Abase[,6:7],1,max) #heard in Mainstem only
    a_Aup_init[,3]= abs(a_ALow_init[,3] -1) # not in upper A
    
    #colnames(Bbase)<-c("lowB","Bup_down","Bup_up","B1_down","B1_up","B2_down","B2_up")
    #mainstem (1), Bup (2), and not there (3)
    a_BLow_init[,3]=abs(a_init[,2] -1) #not in B
    a_BLow_init[,2]=apply(Bbase[,2:7],1,max) #heard upstream
    a_BLow_init[,1]= a_init[,2] - apply(Bbase[,2:7],1,max)#heard in B Mainstem only
    
    #mainstem (1),  B1 (2), B2 (3), and not there (4)
    a_Bup_init[,2]=apply(Bbase[,4:5],1,max) #heard in B1
    a_Bup_init[,3]=apply(Bbase[,6:7],1,max) #heard in B2
    a_Bup_init[,1] =a_BLow_init[,2] - apply(Bbase[,4:7],1,max) #heard in Mainstem only
    a_Bup_init[,4]= abs(a_BLow_init[,2] -1) # not in upper B
    
    #colnames(Cbase)<-c("lowC","C1_down","C1_up","Cup_down","Cup_up","C2_down","C2_up",
    #                   "C3_down","C3_up")
    #mainstem (1), C1 (2), Cup (3), and not there (4)
    a_CLow_init[,4]= abs(a_init[,3] -1) #not in C
    a_CLow_init[,2]= apply(Cbase[,2:3],1,max) #heard in C1
    a_CLow_init[,3]= apply(Cbase[,4:9],1,max) #heard upstream
    a_CLow_init[,1]= a_init[,3] - apply(a_CLow_init[,2:3],1,max)#heard in C Mainstem only
    
    #mainstem (1),  C2 (2), C3 (3), and not there (4)
    a_Cup_init[,2]=apply(Cbase[,6:7],1,max) #heard in C2
    a_Cup_init[,3]=apply(Cbase[,8:9],1,max) #heard in C3
    a_Cup_init[,1] =a_CLow_init[,3] - apply(Cbase[,6:9],1,max) #heard in Mainstem only
    a_Cup_init[,4]= abs(a_CLow_init[,3] -1) # not in upper B
    
    ########################################################################
    #now fit using the posteriors
    #+1 bin is the stay bin
    jags.inits <- function(){ 
      list(
        "a"=c(a_init%*%seq(1,npops)), "a_ALow"=c(a_ALow_init%*%seq(1,npops_ALow+1)), 
        "a_Aup"=c(a_Aup_init%*%seq(1,npops_Aup+1)), 
        "a_BLow"=c(a_BLow_init%*%seq(1,npops_BLow+1)), 
        "a_Bup"=c(a_Bup_init%*%seq(1,npops_Bup+1)),
        "a_CLow"=c(a_CLow_init%*%seq(1,npops_CLow+1)), 
        "a_Cup"=c(a_Cup_init%*%seq(1,npops_Cup+1))
      )
    }
    
    
    ###################################################################
    # below builds the model file:
    ###################################################################
    cat("
        model{
        
        # Set up array detection efficiency priors
        #here we use all beta(1,1) priors, but this could be swapped
        #can be set as a constant
        #or provided with another distribution or information
        #This is where estimates of detection efficiency would be used
        lowA_p~dbeta(1,1)
        A1_up_p~dbeta(1,1)
        A1_down_p~dbeta(1,1)
        Aup_up_p~dbeta(1,1)
        Aup_down_p~dbeta(1,1)
        A2_up_p~dbeta(1,1)
        A2_down_p~dbeta(1,1)
        lowB_p~dbeta(1,1)
        Bup_up_p~dbeta(1,1)
        Bup_down_p~dbeta(1,1)
        B1_up_p~dbeta(1,1)
        B1_down_p~dbeta(1,1)
        B2_up_p~dbeta(1,1)
        B2_down_p~dbeta(1,1)
        lowC_p~dbeta(1,1)
        C1_up_p~dbeta(1,1)
        C1_down_p~dbeta(1,1)
        Cup_up_p~dbeta(1,1)
        Cup_down_p~dbeta(1,1)
        C2_up_p~dbeta(1,1)
        C2_down_p~dbeta(1,1)
        C3_up_p~dbeta(1,1)
        C3_down_p~dbeta(1,1)
        
# First we step the yes/no bin for going to A vs. B vs. C  vs. stay ############
    main_p[1:(npops)]~ddirch(alpha[1:(npops)]); #Dirichlet for probs for going to npops bins
      for (ii in 1:nfish) {
        a[ii] ~ dcat( main_p[1:(npops)] ) 
      }
    # expand the dcat variable into a matrix of zeros and ones
    for (ii in 1:nfish){
      for (jj in 1:(npops))	{
        catexp[ii,jj] <- equals(a[ii],jj) #equals(x,y) is a test for equality, returns [1,0]
              # 1 = A, 2 = B, 3 = C, 4 = not there
      }
    }

################################################################################################
    # NOW WE DEAL WITH A POOLS###############################################################
    # LOWER A- between mouth of A and upper A #----------------------------------------------------
    # There are 4 bins ---  mainstem (1), A1 (2), upstream (3), and not there (4)
    
    p_pop_ALow[1:(npops_ALow)]~ddirch(alpha[1:(npops_ALow)]); #Dirichlet for probs for going to bins
        
    # set up a matrix that deals with yes/no in the tributary or not
    pMatALow[1,1:npops_ALow]<-zero_vec[1:(npops_ALow)] # when not in trib, 0 prob of being in sub areas
    pMatALow[1,(npops_ALow+1)]<-1 #set the 'not there' bin to prob = 1
    pMatALow[2,1:npops_ALow]<-p_pop_ALow # when in trib, >0 probs of being in sub areas
    pMatALow[2,(npops_ALow+1)]<-0 #set the 'not there' bin to prob = 0
        
    for (ii in 1:nfish) {
    # the row number acts as switch between rows 1&2 using stochastic node
    a_ALow[ii] ~ dcat( pMatALow[(catexp[ii,1]+1),1:(npops_ALow+1)] ) # the row number acts as on/off switch
    for (jj in 1:(npops_ALow+1))	{ #now expand the dcat into matrix of zeros and ones
    catexp_ALow[ii,jj] <- equals(a_ALow[ii],jj) #equals(x,y) is a test for equality, returns [1,0]
    }
        
    # OBSERTVATION PART FOR LOWER A
    #first do mainstem (if it's seen anywhere in mainstem OR tribs in lower A -- thus the max statement)
    #first array...
    Abase[ii,1]~dbern(lowA_p*max(catexp_ALow[ii,1],catexp_ALow[ii,2],catexp_ALow[ii,3]))

    # next do A1
    #first array...
    Abase[ii,2]~dbern(A1_down_p*catexp_ALow[ii,2])
    #second array...   
    Abase[ii,3]~dbern(A1_up_p*catexp_ALow[ii,2])
    
   } #ends the ifish loop started at the top of this section

########### UPPER A  ###########################
    #There are 3 bins --- mainstem (1),  A2 (2), and not there (3)
    p_pop_Aup[1:(npops_Aup)]~ddirch(alpha[1:(npops_Aup)]); #Dirichlet for probs for going to bins
        
    # set up a matrix that deals with yes/no in the tributary or not
    pMatAup[1,1:npops_Aup]<-zero_vec[1:(npops_Aup)] # when not in trib, 0 prob of being in sub areas
    pMatAup[1,(npops_Aup+1)]<-1 #set the 'not there' bin to prob = 1
    pMatAup[2,1:npops_Aup]<-p_pop_Aup # when in trib, >0 probs of being in sub areas
    pMatAup[2,(npops_Aup+1)]<-0 #set the 'not there' bin to prob = 0
        
    for (ii in 1:nfish) {
    # the row number acts as switch between rows 1&2 using stochastic node
    a_Aup[ii] ~ dcat( pMatAup[(catexp_ALow[ii,3]+1),1:(npops_Aup+1)] ) # the row number acts as on/off switch (6= Lower-->upper)
    for (jj in 1:(npops_Aup+1))	{ #now expand the dcat into matrix of zeros and ones
    catexp_Aup[ii,jj] <- equals(a_Aup[ii],jj) #equals(x,y) is a test for equality, returns [1,0]
    }
        
    # OBSERTVATION PART FOR UPPER A
    #first do main stem (if it's seen anywhere in mainstem OR tribs in upper A -- thus the max statement)
    #first array...
    Abase[ii,4]~dbern(Aup_down_p*max(catexp_Aup[ii,1],catexp_Aup[ii,2]))
    Abase[ii,5]~dbern(Aup_up_p*max(catexp_Aup[ii,1],catexp_Aup[ii,2]))
    
    # next do A2
    #first array...
    Abase[ii,6]~dbern(A2_down_p*catexp_Aup[ii,2])
    #second array...
    Abase[ii,7]~dbern(A2_up_p*catexp_Aup[ii,2])	
    
    } #ends the iifish loop started at the top of this section   
   
##########
# NOW WE DEAL WITH B#######################################################
# There are 3 bins, Mainstem (1), upper B (2), not there (3)
    p_pop_BLow[1:(npops_BLow)]~ddirch(alpha[1:(npops_BLow)]); #Dirichlet for probs for going to bins
        
    # set up a matrix that deals with yes/no in the tributary or not
    pMatBLow[1,1:npops_BLow]<-zero_vec[1:(npops_BLow)] # when not in trib, 0 prob of being in sub areas
    pMatBLow[1,(npops_BLow+1)]<-1 #set the 'not there' bin to prob = 1
    pMatBLow[2,1:npops_BLow]<-p_pop_BLow # when in trib, >0 probs of being in sub areas
    pMatBLow[2,(npops_BLow+1)]<-0 #set the 'not there' bin to prob = 0
    
    for (ii in 1:nfish) {
    # the row number acts as switch between rows 1&2 using stochastic node
    a_BLow[ii] ~ dcat( pMatBLow[(catexp[ii,2]+1),1:(npops_BLow+1)] ) # the row number acts as on/off switch
    for (jj in 1:(npops_BLow+1))	{ #now expand the dcat into matrix of zeros and ones
    catexp_BLow[ii,jj] <- equals(a_BLow[ii],jj) #equals(x,y) is a test for equality, returns [1,0]
    }
        
    # OBSERTVATION PART FOR LOWER B
    #first do mainstem (if it's seen anywhere in mainstem OR tribs in lower B -- thus the max statement)
    #first array...
    Bbase[ii,1]~dbern(lowB_p*max(catexp_BLow[ii,1],catexp_BLow[ii,2]))
    
    } #ends nfish loop in this section

###### UPPER B  ###############
# 	There are 4 bins --- mainstem (1),  B1 (2), B2 (3), and not there (4)
    p_pop_Bup[1:(npops_Bup)]~ddirch(alpha[1:(npops_Bup)]); #Dirichlet for probs for going to bins
    
    # set up a matrix that deals with yes/no in the tributary or not
    pMatBup[1,1:npops_Bup]<-zero_vec[1:(npops_Bup)] # when not in trib, 0 prob of being in sub areas
    pMatBup[1,(npops_Bup+1)]<-1 #set the 'not there' bin to prob = 1
    pMatBup[2,1:npops_Bup]<-p_pop_Bup # when in trib, >0 probs of being in sub areas
    pMatBup[2,(npops_Bup+1)]<-0 #set the 'not there' bin to prob = 0
    
    for (ii in 1:nfish) {
    # the row number acts as switch between rows 1&2 using stochastic node
    a_Bup[ii] ~ dcat( pMatBup[(catexp_BLow[ii,2]+1),1:(npops_Bup+1)] ) # the row number acts as on/off switch (6= Lower-->upper)
    for (jj in 1:(npops_Bup+1))	{ #now expand the dcat into matrix of zeros and ones
    catexp_Bup[ii,jj] <- equals(a_Bup[ii],jj) #equals(x,y) is a test for equality, returns [1,0]
    }
        
    # OBSERTVATION PART FOR UPPER B
    #first do main stem (if it's seen anywhere in mainstem OR tribs in upper B -- thus the max statement)
    #first array...
    Bbase[ii,2]~dbern(Bup_down_p*max(catexp_Bup[ii,1],catexp_Bup[ii,2],catexp_Bup[ii,3]))
    Bbase[ii,3]~dbern(Bup_up_p*max(catexp_Bup[ii,1],catexp_Bup[ii,2],catexp_Bup[ii,3]))
        
    # next do B1
    #first array...
    Bbase[ii,4]~dbern(B1_down_p*catexp_Bup[ii,2])
    #second array...
    Bbase[ii,5]~dbern(B1_up_p*catexp_Bup[ii,2])	

    # next do B2
    #first array...
    Bbase[ii,6]~dbern(B2_down_p*catexp_Bup[ii,3])
    #second array...
    Bbase[ii,7]~dbern(B2_up_p*catexp_Bup[ii,3])
      
    } #ends the iifish loop started at the top of this section 

################################################################################################
    # NOW WE DEAL WITH C POOLS###############################################################
        # LOWER C- between mouth of C and upper C #----------------------------------------------------
        # There are 4 bins ---  mainstem (1), C1 (2), upstream (3), and not there (4)
        
        p_pop_CLow[1:(npops_CLow)]~ddirch(alpha[1:(npops_CLow)]); #Dirichlet for probs for going to bins
        
        # set up a matrix that deals with yes/no in the tributary or not
        pMatCLow[1,1:npops_CLow]<-zero_vec[1:(npops_CLow)] # when not in trib, 0 prob of being in sub areas
        pMatCLow[1,(npops_CLow+1)]<-1 #set the 'not there' bin to prob = 1
        pMatCLow[2,1:npops_CLow]<-p_pop_CLow # when in trib, >0 probs of being in sub areas
        pMatCLow[2,(npops_CLow+1)]<-0 #set the 'not there' bin to prob = 0
        
        for (ii in 1:nfish) {
        # the row number acts as switch between rows 1&2 using stochastic node
        a_CLow[ii] ~ dcat( pMatCLow[(catexp[ii,3]+1),1:(npops_CLow+1)] ) # the row number acts as on/off switch
        for (jj in 1:(npops_CLow+1))	{ #now expand the dcat into matrix of zeros and ones
        catexp_CLow[ii,jj] <- equals(a_CLow[ii],jj) #equals(x,y) is a test for equality, returns [1,0]
        }
        
        # OBSERTVATION PART FOR LOWER C
        #first do mainstem (if it's seen anywhere in mainstem OR tribs in lower C -- thus the max statement)
        #first array...
        Cbase[ii,1]~dbern(lowC_p*max(catexp_CLow[ii,1],catexp_CLow[ii,2],catexp_CLow[ii,3]))
        
        # next do C1
        #first array...
        Cbase[ii,2]~dbern(C1_down_p*catexp_CLow[ii,2])
        #second array...   
        Cbase[ii,3]~dbern(C1_up_p*catexp_CLow[ii,2])
        
        } #ends the ifish loop started at the top of this section
###### UPPER C  ###############
# 	There are 4 bins --- mainstem (1),  C2 (2), C3 (3), and not there (4)
        p_pop_Cup[1:(npops_Cup)]~ddirch(alpha[1:(npops_Cup)]); #Dirichlet for probs for going to bins
        
        # set up a matrix that deals with yes/no in the tributary or not
        pMatCup[1,1:npops_Cup]<-zero_vec[1:(npops_Cup)] # when not in trib, 0 prob of being in sub areas
        pMatCup[1,(npops_Cup+1)]<-1 #set the 'not there' bin to prob = 1
        pMatCup[2,1:npops_Cup]<-p_pop_Cup # when in trib, >0 probs of being in sub areas
        pMatCup[2,(npops_Cup+1)]<-0 #set the 'not there' bin to prob = 0
        
        for (ii in 1:nfish) {
        # the row number acts as switch between rows 1&2 using stochastic node
        a_Cup[ii] ~ dcat( pMatCup[(catexp_CLow[ii,3]+1),1:(npops_Cup+1)] ) # the row number acts as on/off switch (3= Lower-->upper)
        for (jj in 1:(npops_Cup+1))	{ #now expand the dcat into matrix of zeros and ones
        catexp_Cup[ii,jj] <- equals(a_Cup[ii],jj) #equals(x,y) is a test for equality, returns [1,0]
        }
        
        # OBSERTVATION PART FOR UPPER C
        #first do main stem (if it's seen anywhere in mainstem OR tribs in upper C -- thus the max statement)
        #first array...
        Cbase[ii,4]~dbern(Cup_down_p*max(catexp_Cup[ii,1],catexp_Cup[ii,2],catexp_Cup[ii,3]))
        Cbase[ii,5]~dbern(Cup_up_p*max(catexp_Cup[ii,1],catexp_Cup[ii,2],catexp_Cup[ii,3]))
        
        # next do C2
        #first array...
        Cbase[ii,6]~dbern(C2_down_p*catexp_Cup[ii,2])
        #second array...
        Cbase[ii,7]~dbern(C2_up_p*catexp_Cup[ii,2])	
        
        # next do C3
        #first array...
        Cbase[ii,8]~dbern(C3_down_p*catexp_Cup[ii,3])
        #second array...
        Cbase[ii,9]~dbern(C3_up_p*catexp_Cup[ii,3])
        
        } #ends the iifish loop started at the top of this section 
        
    
    #end model
  }",file="toy_model_basic.txt")
    
    
    ###################################################################
    # now send the data and code to JAGS:
    ###################################################################
    
    #Data that Jags needs in order to run    
    dat = list( "Abase", "Bbase", "Cbase", "npops","npops_ALow", "npops_Aup",
                "npops_BLow", "npops_Bup","npops_CLow", "npops_Cup",
                "alpha", "zero_vec","nfish")

    #parameters to be monitored and saved by JAGS
    parameters <- c("main_p","p_pop_ALow","p_pop_Aup","p_pop_BLow","p_pop_Bup","p_pop_CLow","p_pop_Cup",
                    "lowA_p","A1_up_p", "A1_down_p","Aup_up_p","Aup_down_p","A2_up_p","A2_down_p",
                    "lowB_p","Bup_up_p","Bup_down_p","B1_up_p","B1_down_p","B2_up_p","B2_down_p",
                    "lowC_p","C1_up_p","C1_down_p","Cup_up_p","Cup_down_p","C2_up_p","C2_down_p",
                    "C3_up_p","C3_down_p")
    
    model = "toy_model_basic.txt"
    
    library(runjags)
    library(R2jags)
    library(gtools)
    library(gdata)
    # run the model in JAGS, using default settings
    
    toy_basic <- jags(data=dat, inits=jags.inits, parameters.to.save=parameters, 
                              model.file=model, n.chains = my.chains, n.iter = my.iter, 
                              n.burnin = my.burnin,n.thin=my.thin, DIC = TRUE, 
                              progress.bar="text")

    attach.jags(toy_basic)
    
    save.image(,file=paste(format(Sys.time(), "%Y-%m-%d %I-%p"),"Simulation.Data.toy_basic",j,sep=" "))

    #colnames(Abase)<-c("lowA","A1_down","A1_up","Aup_down","Aup_up","A2_down","A2_up") 
    #colnames(Bbase)<-c("lowB","Bup_down","Bup_up","B1_down","B1_up","B2_down","B2_up")
    #colnames(Cbase)<-c("lowC","C1_down","C1_up","Cup_down","Cup_up","C2_down","C2_up",
    #                  "C3_down","C3_up")
    
    #save output - estimates of fish in each stream (out of nfish tagged)
    #4 main choices (A, B, C, stay)
    res.mat[j,1,]<-apply(cbind(main_p[,1]),1,prod)*nfish  #A
    res.mat[j,2,]<-apply(cbind(main_p[,2]),1,prod)*nfish   #B
    res.mat[j,3,]<-apply(cbind(main_p[,3]),1,prod)*nfish   #C
    res.mat[j,4,]<-apply(cbind(main_p[,4]),1,prod)*nfish   #stay at dam
    
    #for A stream
    #A LOW--proportion going to mainstem (1), A1 (2), up(3)
    #A UP--proportion to mainstem (1),  A2 (2)
    res.A.mat[j,1,]<-apply(cbind(p_pop_ALow[,1],main_p[,1]),1,prod)*nfish #stay A main
    res.A.mat[j,2,]<-apply(cbind(p_pop_ALow[,2],main_p[,1]),1,prod)*nfish #A1
    res.A.mat[j,3,]<-apply(cbind(p_pop_ALow[,3],main_p[,1]),1,prod)*nfish #pass through upstream
    res.A.mat[j,4,]<-apply(cbind(p_pop_ALow[,3],main_p[,1],p_pop_Aup[,1]),1,prod)*nfish #stay up
    res.A.mat[j,5,]<-apply(cbind(p_pop_ALow[,3],main_p[,1],p_pop_Aup[,2]),1,prod)*nfish #A2
    
    # B low- There are 3 bins, Mainstem (1), up(2)
    # B Up-- mainstem(1), B1 (2), B2 (3)
    res.B.mat[j,1,]<-apply(cbind(p_pop_BLow[,1],main_p[,2]),1,prod)*nfish #stay B main
    res.B.mat[j,2,]<-apply(cbind(p_pop_BLow[,2],main_p[,2]),1,prod)*nfish #head B up
    res.B.mat[j,3,]<-apply(cbind(p_pop_BLow[,1],main_p[,2],p_pop_Bup[,1]),1,prod)*nfish #stay B up
    res.B.mat[j,4,]<-apply(cbind(p_pop_BLow[,1],main_p[,2],p_pop_Bup[,2]),1,prod)*nfish #B1
    res.B.mat[j,5,]<-apply(cbind(p_pop_BLow[,1],main_p[,2],p_pop_Bup[,3]),1,prod)*nfish #B2
    
    #ABOVE C--- mainstem (1), c1 (2), up (3)
    #upper C-- mainstem (1), C2 (2), and C3 (3)
    res.C.mat[j,1,]<-apply(cbind(p_pop_CLow[,1],main_p[,3]),1,prod)*nfish #stay C main
    res.C.mat[j,2,]<-apply(cbind(p_pop_CLow[,2],main_p[,3]),1,prod)*nfish #head to C1
    res.C.mat[j,3,]<-apply(cbind(p_pop_CLow[,3],main_p[,3]),1,prod)*nfish #head to C up
    res.C.mat[j,4,]<-apply(cbind(p_pop_CLow[,3],main_p[,3],p_pop_Cup[,1]),1,prod)*nfish #stay at C up main
    res.C.mat[j,5,]<-apply(cbind(p_pop_CLow[,3],main_p[,3],p_pop_Cup[,2]),1,prod)*nfish #C2
    res.C.mat[j,6,]<-apply(cbind(p_pop_CLow[,3],main_p[,3],p_pop_Cup[,3]),1,prod)*nfish #C3
    
    #file names for storing data- model estimates 
    filenameres <- paste(format(Sys.time(), "%Y.%m.%d"),"res.mat.basic",j, sep=" ")
    filenameA <- paste(format(Sys.time(), "%Y.%m.%d"),"res.A.mat.basic",j, sep=" ")
    filenameB <- paste(format(Sys.time(), "%Y.%m.%d"),"res.B.mat.basic",j, sep=" ")
    filenameC <- paste(format(Sys.time(), "%Y.%m.%d"),"res.C.mat.basic",j, sep=" ")

    #save as CSV files
    write(res.mat[j,,], file=filenameres, sep=" ")
    write(res.A.mat[j,,], file=filenameA, sep=" ")
    write(res.B.mat[j,,], file=filenameB, sep=" ")
    write(res.C.mat[j,,], file=filenameC, sep=" ")
    
    rm(toy_basic)
    rm(C.sim)
    rm(C.sim.obs)
    rm(A.sim)
    rm(A.sim.obs)
    rm(B.sim)
    rm(B.sim.obs)
}

save.image("NestedPatchOccSimMod.RDATA") #save entire data set

#note- while this code saves every output from each simulation as a file- 
#it is easier to work with the 3 dimension model output array files
#otherwise you need to read back in the files into R and work with them

#############################################################################
#############################################################################
# How to get Comparisons

#simulated data are stored in [simulation #, stream #]
#A mat columns-"lowA","A1_down","A1_up","Aup_down","Aup_up","A2_down","A2_up" 
#B mat columns- "lowB","Bup_down","Bup_up","B1_down","B1_up","B2_down","B2_up"
#C mat columns-"lowC","C1_down","C1_up","Cup_down","Cup_up","C2_down","C2_up",
#                       "C3_down","C3_up"
#note that in data with no observation error the columns with be the same
   #in A: 2 and 3, 4 and 5, 6 and 7
   #in B: 2 and 3, 4 and 5, 6 and 7
   #in C: 2 and 3, 4 and 5, 6 and 7, 8 and 9
#orig.B.mat
#orig.A.mat
#orig.C.mat 
#simulation data with observation error are stored in
#obs.B.mat 
#obs.A.mat
#obs.C.mat 

#results are all stored in 3d array [simulation #, stream, mcmc result #]
#res.mat #(A, B, C, stay)
#res.A.mat #(A main, A1, A upstream, stay A up, A2)
#res.B.mat #(B main, B up, stay B up, B1, B2)
#res.C.mat #(C main, C1, C up, stay C up, C2, C3)

#############################################################
#Comparing mainstem-passage
orig.A.mat[,1] 	#true #'s passing through A
res.mat[,1,]    #model estimates of A passage (all simulations, all MCMC results)

orig.B.mat[,1]	#true #'s passing through B
res.mat[,2,]	#model estimates of B passage (all simulations, all MCMC results)

orig.C.mat[,1]	#true #'s passing through C
res.mat[,3,]	#model estimates of C passage (all simulations, all MCMC results)

#comparing what stays in mainstem
orig.A.mat[,1]-orig.A.mat[,2]-orig.A.mat[,4] 	#true #'s stay at A [enter A-A1-Aup]
res.A.mat[,1,]    #model estimates of A stay (all simulations, all MCMC results)

orig.B.mat[,1]-orig.B.mat[,2]	#true #'s stay at B [enter B- Bup]
res.B.mat[,2,]	#model estimates of B stay (all simulations, all MCMC results)

orig.C.mat[,1]-orig.C.mat[,2]-orig.C.mat[,4]	#true #'s stay at C [enter C- C1-Cup]
res.C.mat[,3,]	#model estimates of C stay (all simulations, all MCMC results)

#############################################################
#comparing within A branch
#orig.A.mat columns-"lowA","A1_down","A1_up","Aup_down","Aup_up","A2_down","A2_up" 
#res.A.mat #(A main, A1, A upstream, stay A up, A2)

orig.A.mat[,2]	#true #'s to A1
res.A.mat[,2,] 	#model estimate of A1

orig.A.mat[,4]	#true #'s passing Aup
res.A.mat[,3,]	#model estimate passing Aup

orig.A.mat[,4]-orig.A.mat[,6]	#true #'s staying Aup
res.A.mat[,4,]			#model estimate stay Aup

orig.A.mat[,6]	#true #'s to A2
res.A.mat[,5,]	#model estimate to A2

#############################################################
#comparing within B branch
#orig.B.mat columns- "lowB","Bup_down","Bup_up","B1_down","B1_up","B2_down","B2_up"
#res.B.mat #(B main, B up, stay B up, B1, B2)

orig.A.mat[,2]	#true #'s passing Bup
res.B.mat[,2,]	#model estimate passing Bup

orig.A.mat[,2]-orig.A.mat[,4]-orig.A.mat[,6] #true #'s staying Bup
res.B.mat[,3,]	#model estimate staying at Bup

orig.A.mat[,4]	#true # stay at B1
res.B.mat[,4,]	#model estimate of B1

orig.A.mat[,6]	#true # stay at B2
res.B.mat[,5,] 	#model estimate of B2

#############################################################
#comparing within A branch
#orig.C.mat columns-"lowC","C1_down","C1_up","Cup_down","Cup_up","C2_down","C2_up",
#                       "C3_down","C3_up"
#res.C.mat #(C main, C1, C up, stay C up, C2, C3)

orig.C.mat[,2]	#true # C1
res.C.mat[,2,]	#model estimate C1

orig.C.mat[,4]	#true # to Cup
res.C.mat[,3,]	#model estimate to Cup

orig.C.mat[,4]-orig.C.mat[,6]-orig.C.mat[,8]	#true # stay Cup
res.C.mat[,4,]	#model estimate stay Cup

orig.C.mat[,6]	#true # C2
res.C.mat[,5,]	#model estimate to C2

orig.C.mat[,8]	#true # C3
res.C.mat[,6,]	#model estimate to C3

################################################################
# Example of how to compare
#find the median and 95CI for each of nsim simulations of model results
main.est.vals<-apply(res.mat,c(1,2),function(x) quantile(x, c(.025,.5,.975)))

main.est.vals[,,1] #this is the nsim values of 2.5%, median, and 97.5% for A
main.est.vals[,,2] #this is the nsim value of 2.5%, median, and 97.5% for B
main.est.vals[,,3] #this is the nsim value of 2.5%, median, and 97.5% for C

#to calculate how many times the truth is in the 95% CI
A.inside<-rep(NA,length=nsim)
for(i in 1:nsim){
  if(orig.A.mat[i,1]>=main.est.median[1,i,1]&orig.A.mat[i,1]<=main.est.median[3,i,1]){
    A.inside[i]<-1
  }else{A.inside[i]<-0}
}
table(A.inside) #how many times was the true value for # fish entering A inside the 
                #95% Bayes CI

#to calculate % bias of median of MCMC estimate
A.percbias<-100*(main.est.median[2,,1]-orig.A.mat[,1])/orig.A.mat[,1]
A.percbias
boxplot(A.percbias) #boxplot of % bias of median

#violin plot of % vias of median
install.packages("vioplot")
library(vioplot)
vioplot(A.percbias)
