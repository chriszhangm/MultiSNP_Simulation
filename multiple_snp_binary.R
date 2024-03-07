#Case-Control

#OABF
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
  #abf_value = sqrt(abf_w*(1/abf_w+sum(1/(sdr+sqrt(tauS))^2)))*exp(1/2*(z0^2*sum(1/(sdr+sqrt(tauS))^2)-2*(z0/sqrt(abf_w))*sum(z_rep/(sdr+sqrt(tauS)))-(sum(z_rep/(sdr+sqrt(tauS)))^2))*(1/abf_w+sum(1/(sdr+sqrt(tauS))^2))^(-1))
  
  temp1 = sd0^2*sum(1/(sdr^2+tauS))
  abf_value = 1/sqrt(1+temp1) * exp(-1/2*(z0^2*sum(1/(sdr^2+tauS))-2*z0/sd0*sum(z_rep/sqrt(sdr^2+tauS))-sum(z_rep/sqrt(sdr^2+tauS))^2) * (sd0^-2+sum(1/(sdr^2+tauS)))^-1)
  #if(log){return(log(1/abf_value))}
  #else{return(1/abf_value)}
  if(log){return(log(abf_value))}
}

#DSL-estimate tauS #Assume we observe theta_hat and se from the results using logistic regression
tSDSL<-function(theta_hat,se_theta) 
{
  se_inv = 1/(se_theta)^2
  Q = sum(se_inv*(theta_hat-sum(theta_hat*se_inv)/sum(se_inv))^2)
  tSDSL = max(0,(Q-length(theta_hat)-1)/(sum(se_inv)-sum(se_inv*se_inv)/sum(se_inv)))
  return(tSDSL)
}

multipleSNPsimu = function(num_snp=50,num_rep=3,n0=2500,nratio=0.5,tauS=0,prev=0.2){
  #num_snp: the number of significant SNPs
  #num_rep: the number of replications
  #alpha: type I error rate
  #power: statistical power
  #n0: the sample size for the original study
  #nratio: the ratio of large sample size and small sample size for replications, where large sample size is defined as a range from 0.5*n0~1.5n0 and small size from 0.2n0~0.5n0.
  #heritability for all locus
  #prev: prevalence rate
  #coefficients
  numsnp_ll = ceiling(num_snp*0) #large OR ~ 1.2-1.5
  numsnp_ss = num_snp - numsnp_ll #small OR ~ 1.01-1.2
  if(tauS==0){
    beta_all = c(runif(numsnp_ss,1.05,1.5))
  }
  else{
    beta_all = abs(c(runif(numsnp_ss,1.05,1.5)) + rnorm(num_snp,0,tauS))
  }
  beta_all = log(beta_all)
  #Minor allele frequency - 5 groups.
  q1n=q2n=q3n=q4n=floor(num_snp/5)
  q5n=num_snp-sum(c(q1n,q2n,q3n,q4n))
  q_all = c(runif(q1n,0.05,0.1),
            runif(q2n,0.1,0.15),
            runif(q3n,0.15,0.2),
            runif(q4n,0.2,0.3),
            runif(q5n,0.3,0.5))
  #get the intercept in the disease model such that we get our expected prevalence rate.
  #get the design matrix
  X_design = matrix(NA,nrow = 5000,ncol=num_snp)
  Xbeta = matrix(NA,nrow = 5000,ncol=num_snp)
  for (i in 1:num_snp) {
    qi = q_all[i]
    X_design[,i] = sample(c(0,1,2),size = 5000,replace = T,prob = c((1-qi)^2,2*qi*(1-qi),qi^2))
    Xbeta[,i] = X_design[,i]*beta_all[i]
  }
  
  #alpha_cand = seq(-10,-15,-0.1)
  #PrevalenceY = numeric(length(alpha_cand))
  # if(prev==0.01){alpha=-15}
  # else if(prev==0.05){alpha=-12}
  # else{alpha=-9}
  alpha=-1
  PrevalenceY=0.5
  while (abs(PrevalenceY-prev)>0.01) {
    ProbLogistic = plogis(rowSums(Xbeta) + alpha)
    Ydisease = rbinom(n = 5000,size = 1,prob = ProbLogistic)
    PrevalenceY = sum(Ydisease)/5000
    alpha = alpha - 0.01
  }
  #alpha
  
  #generate the original study
  ns=1.1*n0/prev
  X_design = matrix(NA,nrow = ns,ncol=num_snp)
  Xbeta = matrix(NA,nrow = ns,ncol=num_snp)
  for (i in 1:num_snp) {
    qi = q_all[i]
    X_design[,i] = sample(c(0,1,2),size = ns,replace = T,prob = c((1-qi)^2,2*qi*(1-qi),qi^2))
    Xbeta[,i] = X_design[,i]*beta_all[i]
  }
  ProbLogistic = plogis(rowSums(Xbeta) + alpha) #alpha has been found through the procedure above to make sure the prevalence we want.
  Ydisease = numeric(ns)
  for (j in 1:ns) {
    Ydisease[j] = rbinom(n = 1,size = 1,prob = ProbLogistic[j])
  }
  #pick n0/2 patients for the control group and n0/2 for the case group.
  ind_case = sample(which(Ydisease==1),n0/2,replace = T)
  ind_control = sample(which(Ydisease==0),n0/2,replace = T)
  Ydisease_truncate = Ydisease[c(ind_case,ind_control)]
  X_design_truncate = X_design[c(ind_case,ind_control),]
  #fit logistic regression to get estimated effect size, standard error, z0, and p values.
  #y0info: row -> SNP1,2,...
  y0info = matrix(NA,nrow = num_snp,ncol=4)
  for (i in 1:num_snp) {
    glmfit = summary(glm(Ydisease_truncate~X_design_truncate[,i],family = binomial(link = "logit")))
    y0info[i,] = c(glmfit$coefficients[c(2,4,6,8)])
  }
  
  #replications
  #sample size
  n_rep_large = floor(nratio*num_rep)
  n_rep_small = num_rep - n_rep_large
  if (n_rep_large !=0) {
    nr = floor(c(runif(n_rep_large,min = 0.5*n0,max = 1.5*n0),
                 runif(n_rep_small,min = 0.2*n0,max = 0.5*n0)))
  }
  else{nr=floor(runif(num_rep,0.2*n0,0.5*n0))}
  nrs = 1.1*nr/prev
  #generate genotypes for the replications using 3d-array
  yrinfo_total = array(NA,dim = c(num_snp,4,num_rep))
  for (j in 1:num_rep) {
    X_design = matrix(NA,nrow = nrs[j],ncol=num_snp)
    Xbeta = matrix(NA,nrow = nrs[j],ncol=num_snp)
    for (i in 1:num_snp) {
      qi = q_all[i]
      X_design[,i] = sample(c(0,1,2),size = nrs[j],replace = T,prob = c((1-qi)^2,2*qi*(1-qi),qi^2))
      Xbeta[,i] = X_design[,i]*beta_all[i]
    }
    ProbLogistic = plogis(rowSums(Xbeta) + alpha) #alpha has been found through the procedure above to make sure the prevalence we want.
    Yrdisease = numeric(nrs[j])
    for (k in 1:nrs[j]) {
      Yrdisease[k] = rbinom(n = 1,size = 1,prob = ProbLogistic[k])
    }
    #pick floor(nrs[j]/2) patients for the control group and floor(nrs[j]/2) for the case group.
    ind_case = sample(which(Yrdisease==1),floor(nr[j]/2),replace = T)
    ind_control = sample(which(Yrdisease==0),floor(nr[j]/2),replace = T)
    Yrdisease_truncate = Yrdisease[c(ind_case,ind_control)]
    X_design_truncate = X_design[c(ind_case,ind_control),]
    
    yrinfo = matrix(NA,nrow = num_snp,ncol=4)
    for (i in 1:num_snp) {
      glmfit = summary(glm(Yrdisease_truncate~X_design_truncate[,i],family = binomial(link = "logit")))
      yrinfo[i,] = c(glmfit$coefficients[c(2,4,6,8)])
    }
    #colnames(y0info) = c('es','sd','z','p')
    yrinfo_total[,,j] = yrinfo
  }
  return(list(y0=y0info,yr=yrinfo_total))
}
#kkk = multipleSNPsimu()

Simu_Binary = function(M=200,tauS=0,num_rep=5,nratio=0.5,prev=0.1){
  abf_total = numeric(M)
  for (m in 1:M) {
    tt1 = multipleSNPsimu(num_snp = 50,tauS=tauS,n0 = 5000,num_rep =num_rep,nratio = nratio,prev = prev)
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
    #print(m)
  }
  return(abf_total)
}


prev1 = c(0.05,0.1,0.2)
tauS1 = c(0,0.1)
numrep1 = c(3,5,10)

data_total_mix=numeric(0)
for (iii in 1:3) {
  for (jjj in 1:2) {
    for (kkk in 1:3) {
      prev=prev1[iii]
      hetero=tauS1[jjj]
      mm=numrep1[kkk]
      data_total_mix = rbind(data_total_mix,Simu_Binary(M=300,tauS=hetero,num_rep=mm,nratio=0.5,prev=prev))
    }
  }
}

data_total_ss=numeric(0)
for (iii in 1:3) {
  for (jjj in 1:2) {
    for (kkk in 1:3) {
      prev=prev1[iii]
      hetero=tauS1[jjj]
      mm=numrep1[kkk]
      data_total_ss = rbind(data_total_ss,Simu_Binary(M=300,tauS=hetero,num_rep=mm,nratio=0,prev=prev))
    }
  }
}

data_total = cbind(data_total_ss,data_total_mix)
write.csv(data_total,"df_binary.csv")
