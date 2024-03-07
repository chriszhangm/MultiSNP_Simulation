###Continuous Traits Simulation###
library(pwr)
#ABF
#output:Approximate Bayes factor
abf = function(es0,sd0,esr,sdr,heter=F,log=T){
  #es0:estimated effect size of the original study
  #sd0:estimated standard error of the original study
  #esr:estimated effect size(s) of the replication(s)
  #sdr:estimated standard error(s) of the replication(s)
  #heter: T->assume the existence of the heterogeneity, F -> homogeneity
  #log: T -> return log(BF), otherwise return BF
  abf_w = sd0^2
  z0=es0/sd0
  z_rep = esr/sdr
  if (heter) {
    tauS=tSDSL(esr,sdr)
  }
  else{tauS=0}
  abf_value = sqrt(abf_w*(1/abf_w+sum(1/(sdr+sqrt(tauS))^2)))*exp(1/2*(z0^2*sum(1/(sdr+sqrt(tauS))^2)-2*(z0/sqrt(abf_w))*sum(z_rep/(sdr+sqrt(tauS)))-(sum(z_rep/(sdr+sqrt(tauS)))^2))*(1/abf_w+sum(1/(sdr+sqrt(tauS))^2))^(-1))
  if(log){return(log(1/abf_value))}
  else{return(1/abf_value)}
}

#DSL-estimate tauS #Assume we observe theta_hat and se from the results using logistic regression
tSDSL<-function(theta_hat,se_theta) 
{
  se_inv = 1/(se_theta)^2
  Q = sum(se_inv*(theta_hat-sum(theta_hat*se_inv)/sum(se_inv))^2)
  tSDSL = max(0,(Q-length(theta_hat)-1)/(sum(se_inv)-sum(se_inv*se_inv)/sum(se_inv)))
  return(tSDSL)
}
##Quantative Traits
multipleSNPsimu = function(num_snp=20,num_rep=3,n0=5000,nratio=0.5,rr=0.1,tauS=0.1){
  #num_snp: the number of significant SNPs
  #num_rep: the number of replications
  #alpha: type I error rate
  #power: statistical power
  #n0: the sample size for the original study
  #nratio: the ratio of large sample size and small sample size for replications, where large sample size is defined as a range from 0.5*n0~1.5n0 and small size from 0.2n0~0.5n0.
  #heritability for all locus
  #coefficients
  if(tauS==0){
    beta_all = runif(num_snp,1.5,2)
  }
  else{
    beta_all = runif(num_snp,1.5,2)+rnorm(num_snp,0,tauS)
  }
  #Minor allele frequency - 5 groups.
  q1n=q2n=q3n=q4n=floor(num_snp/5)
  q5n=num_snp-sum(c(q1n,q2n,q3n,q4n))
  q_all = c(runif(q1n,0.05,0.1),
            runif(q2n,0.1,0.15),
            runif(q3n,0.15,0.2),
            runif(q4n,0.2,0.3),
            runif(q5n,0.3,0.4))
  #calculate the residual variance
  ve = ((1-rr)*sum(2*beta_all^2*(q_all)*(1-q_all)))/rr
  
  #record the power for all loci
  power_loci8 = numeric(num_snp)
  power_loci6 = numeric(num_snp)
  rr2i = beta_all^2*2*q_all*(1-q_all)/(sum(beta_all^2*2*q_all*(1-q_all))+ve)
  ff2i = rr2i/(1-rr2i)
  for (i in 1:num_snp) {
    power_loci8[i] = pwr.f2.test(u=1,v=n0-2,f2 = ff2i[i],sig.level = 5e-8)$power
    power_loci6[i] = pwr.f2.test(u=1,v=n0-2,f2 = ff2i[i],sig.level = 5e-6)$power
  }
  #generate genotypes for the original study
  genotype = matrix(NA,nrow = n0,ncol=num_snp)
  for (i in 1:num_snp) {
    genotype[,i] = sample(c(0,1,2),n0,replace = T,prob=c((1-q_all[i])^2,2*(1-q_all[i])*q_all[i],q_all[i]^2))
  }
  beta_matrix = matrix(rep(beta_all,n0),nrow = n0,byrow = T)
  y0 = rowSums(beta_matrix*genotype) + rnorm(n0,0,sd = sqrt(ve))
  #fit linear regression to get estimated effect size, standard error, t0, and p values.
  #y0info: row -> SNP1,2,...
  y0info = matrix(NA,nrow = num_snp,ncol=4)
  for (i in 1:num_snp) {
    lmfit = summary(lm(y0~genotype[,i]))
    y0info[i,] = c(lmfit$coefficients[c(2,4,6,8)])
  }
  #colnames(y0info) = c('es0','sd0','z0','p0')
  
  #replications
  #sample size
  n_rep_large = floor(nratio*num_rep)
  n_rep_small = num_rep - n_rep_large
  if (n_rep_large !=0) {
    nr = floor(c(runif(n_rep_large,min = 0.5*n0,max = 1.5*n0),
           runif(n_rep_small,min = 0.1*n0,max = 0.5*n0)))
  }
  else{nr=floor(runif(num_rep,0.1*n0,0.5*n0))}
  #generate genotypes for the replications using 3d-array
  yrinfo_total = array(NA,dim = c(num_snp,4,num_rep))
  for (j in 1:num_rep) {
    genotype = matrix(NA,nrow = nr[j],ncol=num_snp)
    for (i in 1:num_snp) {
      genotype[,i] = sample(c(0,1,2),nr[j],replace = T,prob=c((1-q_all[i])^2,2*(1-q_all[i])*q_all[i],q_all[i]^2))
    }
    beta_matrix = matrix(rep(beta_all,nr[j]),nrow = nr[j],byrow = T)
    yr = rowSums(beta_matrix*genotype) + rnorm(nr[j],0,sd = sqrt(ve))
    yrinfo = matrix(NA,nrow = num_snp,ncol=4)
    for (i in 1:num_snp) {
      lmfit = summary(lm(yr~genotype[,i]))
      yrinfo[i,] = c(lmfit$coefficients[c(2,4,6,8)])
    }
    #colnames(y0info) = c('es','sd','z','p')
    yrinfo_total[,,j] = yrinfo
  }
  return(list(y0=y0info,yr=yrinfo_total,power_orig8 = power_loci8,power_orig6 = power_loci6))
  
} 
##generate a set of K=200 simulations##
##50 SNPs
##num_rep = {3,5,10}
##rr=0.2,0.25,0.3
##tauS=0,0.1,0.2
##prior = 1e-5

Simu_Qtl = function(M=200,rr=0.3,tauS=0,num_rep=5,nratio=0.5){
  abf_total = numeric(M)
  for (m in 1:M) {
    tt1 = multipleSNPsimu(num_snp = 50,rr = rr,tauS=tauS,num_rep =num_rep,nratio = nratio)
    abf1 = numeric(50)
    for (i in 1:50) {
      esr1=numeric(num_rep)
      sdr1=numeric(num_rep)
      for (j in 1:num_rep) {
        esr1[j] = tt1$yr[i,1,j]
        sdr1[j] = tt1$yr[i,2,j]
      }
      abf1[i] = abf(es0=tt1$y0[i,1],sd0 = tt1$y0[i,2],esr=esr1,sdr=sdr1,heter = T)
    }
    abf_total[m] = sum(exp(abf1)*1e-5>=3)/50
  }
  return(abf_total)
}

rr1 = c(0.2,0.25,0.3)
tauS1 = c(0,0.2,0.4)
numrep1 = c(3,5,10)

data_total_mix=numeric(0)
for (iii in 1:3) {
    for (jjj in 1:3) {
      for (kkk in 1:3) {
        hsnp=rr1[iii]
        hetero=tauS1[jjj]
        mm=numrep1[kkk]
        data_total_mix = rbind(data_total_mix,Simu_Qtl(M=200,rr=hsnp,tauS = hetero,num_rep = mm,nratio = 0.5))
      }
    }
}

data_total_ss=numeric(0)
for (iii in 1:3) {
  for (jjj in 1:3) {
    for (kkk in 1:3) {
      hsnp=rr1[iii]
      hetero=tauS1[jjj]
      mm=numrep1[kkk]
      data_total_ss = rbind(data_total_ss,Simu_Qtl(M=200,rr=hsnp,tauS = hetero,num_rep = mm,nratio = 1))
    }
  }
}

data_total = cbind(data_total_ss,data_total_mix)
write.csv(data_total,"df1.csv")

multipleSNPsimu(num_snp = 20,n0 = 2500,rr =.1,tauS=0,num_rep =3,nratio = 0)
