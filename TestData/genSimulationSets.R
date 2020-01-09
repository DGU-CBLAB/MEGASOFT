#### Generated on 2019-12-11
#### Created by Cue Hyunkyu Lee
#### Last updated on 2019-12-12

## The package mvtnorm can computes multivariate normal and t probabilities, quantiles, random deviates and densities.
require(mvtnorm)

tau2.zero.prob=c(0.8415,0.7709,0.7413,0.7135,0.6989,0.6742,0.663,0.6542,0.6506,0.646,0.6411,0.6281,0.6229,0.6308,0.6299,0.6152,0.6103,0.6121,0.618,0.6049,0.6068,0.597,0.5924,0.5946,0.5954,0.5881,0.5899,0.5809,0.5922,0.5854,0.583,0.5794,0.5768,0.5817,0.5773,0.5848,0.5764,0.5821,0.5777,0.5804,0.5735,0.5736,0.5691,0.5674,0.5681,0.5739,0.5705,0.568,0.5587) ## From nstudy 2 to nstudy 50 (49 entries)

RE2 <- function(beta, stders, ...) {
  ##------------------------------
  ## Simplified RE2
  ##------------------------------
  n=length(beta)
  vars <- stders**2
  sigma=diag(vars)
  ws <- 1/vars
  sumBetaProdW <- as.vector(beta %*% ws)
  sumW <- sum(ws)
  sumW2 <- sum(ws**2)
  meanBeta <- sumBetaProdW / sumW
  Q <- (beta - meanBeta)**2 %*% ws
  meanW <- mean(ws)
  Sw2 <- 1/(n-1) * (sumW2 - n*meanW*meanW)
  U = (n-1)*(meanW - Sw2/(n*meanW))
  tau2 <- max( (Q-(n-1))/U, 0 )
  ##----------------------------------------------
  ## Eigen-decomposition optimization (EMMA- style)
  ##----------------------------------------------
  K <- sigma
  eig.L <- eigen(K,symmetric=TRUE)
  L.values <- eig.L$values
  L.vectors <- eig.L$vectors
  S <- diag(n)-matrix(1/n,n,n)
  eig.R<- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  R.values <- eig.R$values[1:(n-1)]
  R.vectors <- eig.R$vectors[,1:(n-1)]
  etas <- crossprod(R.vectors,beta)
  etasqs <- etas^2
  xis <- L.values
  lambdas <- R.values
  
  LL.fun <- function(x) {
    0.5*(n*log(2*pi)+sum(log(xis+x))+sum(etasqs/(lambdas+x)))
  }
  optim.rst <- optim(par=c(tau2),fn=LL.fun, method = "L-BFGS-B",
                     lower = 0, upper=Inf)
  mle.tau2 <- optim.rst$par[1]
  Hinv <- solve(sigma+mle.tau2*diag(n))
  mle.beta <- colSums(Hinv) %*% beta / sum(Hinv)
  mle.ll <- -optim.rst$value
  
  ll.fun <- function(x)
    -sum(dmvnorm(beta, mean=rep(x[1],length(beta)),
                 sigma=sigma+diag(x[2],length(beta)),log=T))
  null.ll=-ll.fun(c(0,0))
  lrt.mle <- -2*(null.ll-mle.ll)
  ##p.RE2 <- 0.5*pchisq(lrt.mle,1,lower.tail=F)+0.5*pchisq(lrt.mle,2,lower.tail=F)
  stat1=-2*(null.ll+ll.fun(c(meanBeta,0)))
  stat2=-2*(-ll.fun(c(meanBeta,0))+ll.fun(c(mle.beta,tau2)))
  
  p.RE2 <- tau2.zero.prob[n-1]*pchisq(lrt.mle,1,lower.tail=F)+(1-tau2.zero.prob[n-1])*pchisq(lrt.mle,2,lower.tail=F)
  return(p.RE2)
}

generate_studies = function(N, MAF, gamma, p = 1e-10){
  # We assume that the equal number of cases and controls are collected and refers to the total number of individuals as sample size N
  N_case = floor(N/2)
  N_cont = N - N_case

  # ### 
  # MAF : minor allele frequency (In simulation, the value is fixed), 
  # N: the number of total samples, 
  # M : the number of total simulations 
   
  MAF_cont = MAF # If we assume a small disease prevalence, the expected frequency in controls is the population MAF
  # and the exprected frequency in cases is 
  MAF_case = (gamma * MAF) / ((gamma - 1) * MAF + 1)
  
  # 2 by 2 contingency table 
  #              cases | controls
  # exposed    | a     | b
  # reference  | c     | d 
  
  while ( p < 1e-8 ){
    case_genotypes = rbinom(n = N_case, size = 2, prob = MAF_case)
    cont_genotypes = rbinom(n = N_cont, size = 2, prob = MAF_cont)
    a = sum(case_genotypes)
    c = N_case * 2 - a
    b = sum(cont_genotypes)
    d = N_case * 2 - b
    
    beta = log( (a * d) / (b * c) )
    var = (1/a + 1/b + 1/c + 1/d)
    se = sqrt(var)
    
    p = pchisq( (beta / se)^2, df=1, lower.tail = F)
  }

  return(c(beta, se))
}

generate_a_simulation = function(RE2_p = 1)
{
  
  while ( RE2_p > 1e-8 )
  {
    betas = rep(0, 8)
    stders = rep(0, 8)
    
    ### In this simulation example, we assume four different types of studies.
    ## the first type 
    N = 2000
    gamma = 1.3 # relative risk 
    MAF = 0.3 # population minor allele frequency (fixed)
    
    for ( j in 1:2 ){
      res = generate_studies(N, MAF, gamma, p = 1e-10)
      betas[j] = res[1]
      stders[j] = res[2]
    }
    
    ## the second type 
    N = 100
    gamma = 1.3
    MAF = 0.3 # population minor allele frequency (fixed)
    
    for ( j in 3:4 ){
      res = generate_studies(N, MAF, gamma, p = 1e-10)
      betas[j] = res[1]
      stders[j] = res[2]
    }
    
    ## the thirt type 
    N = 2000
    gamma = 1.0
    MAF = 0.3 # population minor allele frequency (fixed)
    
    for ( j in 5:6 ){
      res = generate_studies(N, MAF, gamma, p = 1e-10)
      betas[j] = res[1]
      stders[j] = res[2]
    }
    
    ## the fourth type 
    N = 100
    gamma = 1.0
    MAF = 0.3 # population minor allele frequency (fixed)
    
    for ( j in 7:8 ){
      res = generate_studies(N, MAF, gamma, p = 1e-10)
      betas[j] = res[1]
      stders[j] = res[2]
    }
    
    RE2_p = RE2(betas, stders)
  }

  a_simulation_set = rbind(betas, stders)
  return(a_simulation_set)
}

main = function(){
  M = 1000
  betas = matrix(rep(0, M * 8), ncol = 8)
  stders = matrix(rep(0, M * 8), ncol = 8)
  for (iter in 1:M)
  {
    res = generate_a_simulation()
    betas[iter,] = res['betas',]
    stders[iter,] = res['stders',]
  }

  METASOFT_input = matrix(rep(0,(8*2+1)*M), nrow = M)
  colnames(METASOFT_input) = c('SNP',paste(rep(paste('Study_',1:8, sep=''),each=2), rep(c('_beta','_se'),8),sep='')  )
  SNP = paste(1:M,sep='')
  for (i in 1:M){
    METASOFT_input[i,] = c(SNP[i],as.vector(rbind(betas[i,],stders[i,])))
  }
  
  write.table(METASOFT_input,'./inputMS.txt',quote=F,row.names=F,col.names=F)
  
  
  }

cat('The process will begin generating 1000 sets of METASOFT inputs \n\n')
begin = Sys.time()
main()
cat('Done \n')



# write.table(METASOFT_input,'[out.txt]',quote=F,row.names=F,col.names=F)

# MetaSOFT command line arguments
# java -jar Metasoft.jar -input [out.txt] -output [out_meta.txt] -pvalue_table [HanEskinPvalueTable.txt] -mvalue

# Replicate the result in Figure of the the paper, 'Interpreting Meta-Analyses of Genome-Wide Association Studies'.
# We should remove the column index of out_meta.txt (the length of header of METASOFT is not the same as the length of output columns)
# meta_out = read.table('out_meta.txt',header = F)
# 
# mvlues = meta_out[,25:32]
# type1 = as.vector(as.matrix(mvlues[,1:2]))
# type2 = as.vector(as.matrix(mvlues[,3:4]))
# type3 = as.vector(as.matrix(mvlues[,5:6]))
# type4 = as.vector(as.matrix(mvlues[,7:8]))
# 
# 
# hist(type1[!is.na(type1) ], breaks=100)
# hist(type2[!is.na(type2) ], breaks=100)
# hist(type3[!is.na(type3) ], breaks=100)
# hist(type4[!is.na(type4) ], breaks=100)

passed = Sys.time()-begin
cat('Analysis time: ',passed,'secs')
