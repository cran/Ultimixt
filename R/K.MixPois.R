K.MixPois <-
function(xobs, k, alpha0, alpha, Nsim){
 alpha0=rep(alpha0, k)
 alpha1=rep(alpha, k)
 n=length(xobs)
 p0=rep(1/k, k)
 gamma0=rep(1/k, k)
    
 sd_lambda=.1
 eps1=500
 eps2=500
 betta=rep(1, 2)
 T=30000
 m=length(betta)
 naccept1=rep(0,m)
 naccept2=rep(0,m)
 naccept3=rep(0,m)
 naccept5=rep(0,m)
 na1=0
 na2=0
 na3=0
 na5=0
 adapt_lambda_ra=rep(0,T)
 adapt_sd_lambda=rep(0,T)
 adapt_p_ra=rep(0,T)
 adapt_eps_p=rep(0,T)
 adapt_gamma_ra=rep(0,T)
 adapt_eps_gamma=rep(0,T)
    
 n=length(xobs)
 meanobs=mean(xobs)
    
 mui=rep(mean(xobs), m)
 output_lambda_s=matrix(mean(xobs), nrow=T, ncol=m)
    
 pp=matrix(p0, nrow=m, ncol=k)
 output_p_s=array(p0, dim=c(T,m,k))
 
 gama=matrix(gamma0, nrow=m, ncol=k)
 output_gamma_s=array(gamma0, dim=c(T,m,k))
 
 component_mean=array(mean(xobs), dim=c(T,m,k))
    
 a_adapt=0 #batch number of 100 iterations
 a=1
    
den_dirichlet<-function (x, alpha)
 {
  dirichlet1 <- function(x, alpha) {
   logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
   s <- sum((alpha - 1) * log(x))
   exp(sum(s) - logD)
 }
  if (!is.matrix(x))
  if (is.data.frame(x))
  x <- as.matrix(x)
  else x <- t(x)
  if (!is.matrix(alpha))
  alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x),byrow = TRUE)
  if (any(dim(x) != dim(alpha)))
  stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
  pd <- vector(length = nrow(x))
  for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i,])
  pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
  pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
  return(pd)
}
rdirichlet=function(n, alpha){
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}

rept=1
adapt_condition_lambda=TRUE
adapt_condition_p=TRUE
adapt_condition_gamma=TRUE
adapt_condition=TRUE
T1=1
T_stop=1000
# T_stop=15000
criter_condition=FALSE
a_lambda=0
a_p=0
a_gamma=0
A=1
repeat{
#################################
for(t in 1:(T-1)){
T1[rept==1]=T1[rept==1]+1

################# Hasting independent proposal for lambda:
comp_mean=sweep(output_gamma_s[t,,],MARGIN=1,output_lambda_s[t,],`*`)/output_p_s[t,,]

v=1
mix_prob=matrix(0,ncol=m, nrow=n)
 repeat{
mixprob=apply(cbind(comp_mean[v,], output_p_s[t,,][v,]), 1, function(u)u[2]*dpois(xobs, lambda=u[1], log=FALSE))
mixprob=rowSums(mixprob)
mix_prob[,v]=mixprob
if(v==m)
break
v=v+1
}
cur=mix_prob
lcur=log(cur)
slcur= colSums(lcur, na.rm=FALSE, dims=1)
curlik=slcur-log(output_lambda_s[t,])
if(adapt_condition==TRUE){
prop1=exp(log(output_lambda_s[t,])+rnorm(m,0, sd_lambda))

log_mix_lambdat=-log(output_lambda_s[t,])+dnorm(log(output_lambda_s[t,]),log(prop1),sd_lambda, log=TRUE)
log_mix_lambda_prop=-log(prop1)+dnorm(log(prop1),log(output_lambda_s[t,]),sd_lambda, log=TRUE)
}else{
    meanobs=mean(xobs)
prop1=exp(log(meanobs)+rnorm(m,0, sd_lambda))
log_mix_lambdat=-log(output_lambda_s[t,])+dnorm(log(output_lambda_s[t,]),log(meanobs),sd_lambda, log=TRUE)
log_mix_lambda_prop=-log(prop1)+dnorm(log(prop1),log(meanobs),sd_lambda, log=TRUE)
}

comp_mean=sweep(output_gamma_s[t,,],MARGIN=1,prop1,`*`)/output_p_s[t,,]
v=1
mix_prob=matrix(0,ncol=m, nrow=n)
repeat{
mixprob=apply(cbind(comp_mean[v,], output_p_s[t,,][v,]), 1, function(u)u[2]*dpois(xobs, lambda=u[1], log=FALSE))
mixprob=rowSums(mixprob)
mix_prob[,v]=mixprob
if(v==m)
break
v=v+1
}
propl=mix_prob
lpropl=log(propl)
slpropl=colSums(lpropl, na.rm=FALSE, dims=1)
proplik=slpropl-log(prop1)
                
log_ro_lambda=betta*proplik-betta*curlik+log_mix_lambdat-log_mix_lambda_prop
Accept_lambda=(log(runif(length(betta)))<(log_ro_lambda))

mui[Accept_lambda]=prop1[Accept_lambda]
output_lambda_s[t+1,]=mui
                
slcur[Accept_lambda]=slpropl[Accept_lambda]
naccept1[Accept_lambda]=naccept1[Accept_lambda]+1
na1[Accept_lambda[1]]=na1[Accept_lambda[1]]+1
a_adapt[a_adapt<50]=a_adapt[a_adapt<50]+1

if(a>300){men_adapt_lambda=median(adapt_lambda_ra[(a-10):a])
cond=(men_adapt_lambda>.35&men_adapt_lambda<.55)
adapt_condition_lambda[cond==TRUE]=FALSE
}
if(adapt_condition_lambda==TRUE){
 ac_lambda_ra=na1/50
 adapt_lambda_ra[a][a_adapt==50]=ac_lambda_ra
 ac_cond_lambda1=(ac_lambda_ra<.44)
 ac_cond_lambda2=(ac_lambda_ra>.44)
 ls_lambda=log(sd_lambda[1])
 ls_lambda[ac_cond_lambda1&a_adapt==50]=ls_lambda[ac_cond_lambda1&a_adapt==50]-min(.01,1/sqrt(t))
 ls_lambda[ac_cond_lambda2&a_adapt==50]=ls_lambda[ac_cond_lambda2&a_adapt==50]+min(.01,1/sqrt(t))
 sd_lambda[1]=exp(ls_lambda)
 adapt_sd_lambda[a][a_adapt==50]=sd_lambda[1]
 a_lambda[a_adapt==50]=a_lambda[a_adapt==50]+1}else{adapt_lambda_ra[a][a_adapt==50]=adapt_lambda_ra[a-1][a_adapt==50]
 adapt_sd_lambda[a][a_adapt==50]=adapt_sd_lambda[a-1][a_adapt==50]}
################# Hasting independent proposal for gammas:
curlik=slcur+log(apply(output_gamma_s[t,,], 1, function(u)den_dirichlet(u, alpha1)))

alp1=output_gamma_s[t,,]*eps2+1
prop4=matrix(0, nrow=m, ncol=k)
v=1
repeat{
    prop4[v,]=rdirichlet(1, alp1[v,])
    if(v==m)
    break
    v=v+1
}
alppt=alp1
alpp=prop4*eps2+1
log_dirgamma1=rep(0,m)
log_dirgamma2=rep(0,m)
for(l in 1:m){
    log_dirgamma1[l]=log(apply(matrix(output_gamma_s[t,,][l,], nrow=1), 1, function(u)den_dirichlet(u, alpp[l,])))
    log_dirgamma2[l]=log(apply(matrix(prop4[l,], nrow=1), 1, function(u)den_dirichlet(u, alppt[l,])))
}
###############
comp_mean=sweep(prop4,MARGIN=1,output_lambda_s[t+1,],`*`)/output_p_s[t,,]

v=1
mix_prob=matrix(0,ncol=m, nrow=n)
repeat{
    mixprob=apply(cbind(comp_mean[v,], output_p_s[t,,][v,]), 1, function(u)u[2]*dpois(xobs, lambda=u[1], log=FALSE))
    mixprob=rowSums(mixprob)
    mix_prob[,v]=mixprob
    if(v==m)
    break
    v=v+1
}

propl=mix_prob
lpropl=log(propl)
slpropl=colSums(lpropl, na.rm=FALSE, dims=1)
proplik=slpropl+log(apply(prop4, 1, function(u)den_dirichlet(u, alpha1)))

log_ro_gamma=betta*proplik-betta*curlik+log_dirgamma1-log_dirgamma2
Accept_gamma=(log(runif(length(betta)))<(log_ro_gamma))

gama[Accept_gamma,]=prop4[Accept_gamma,]
output_gamma_s[t+1,,]=gama

slcur[Accept_gamma]=slpropl[Accept_gamma]
naccept5[Accept_gamma]=naccept5[Accept_gamma]+1
na5[Accept_gamma[1]]=na5[Accept_gamma[1]]+1

if(a>300){men_adapt_gamma=median(adapt_gamma_ra[(a-10):a])
        cond=(men_adapt_gamma>.2&men_adapt_gamma<.28)
    adapt_condition_gamma[cond==TRUE]=FALSE
}

if(adapt_condition_gamma==TRUE){
    ac_gamma_ra=na5/50
    adapt_gamma_ra[a][a_adapt==50]=ac_gamma_ra

            ac_cond_gamma1=(ac_gamma_ra<.234)
            ac_cond_gamma2=(ac_gamma_ra>.234)
            
            coco1=(eps2<=5)
            coco2=(eps2>5)
            eps2[ac_cond_gamma1&coco2&a_adapt==50]=eps2[ac_cond_gamma1&coco2&a_adapt==50]+5
            eps2[ac_cond_gamma2&coco2&a_adapt==50]=eps2[ac_cond_gamma2&coco2&a_adapt==50]-5
            eps2[ac_cond_gamma1&coco1&a_adapt==50]=eps2[ac_cond_gamma1&coco1&a_adapt==50]+5
            eps2[ac_cond_gamma2&coco1&a_adapt==50]=eps2[ac_cond_gamma2&coco1&a_adapt==50]-eps2[ac_cond_gamma2&coco1&a_adapt==50]/100

        adapt_eps_gamma[a][a_adapt==50]=eps2
    a_gamma[a_adapt==50]=a_gamma[a_adapt==50]+1}else{adapt_gamma_ra[a][a_adapt==50]=adapt_gamma_ra[a-1][a_adapt==50]
        adapt_eps_gamma[a][a_adapt==50]=adapt_eps_gamma[a-1][a_adapt==50]}

################### Hasting independent proposal for p:
curlik=slcur+log(apply(output_p_s[t,,], 1, function(u)den_dirichlet(u, alpha0)))
                            
alp1=output_p_s[t,,]*eps1+1
prop4=matrix(0, nrow=m, ncol=k)
v=1
repeat{
prop4[v,]=rdirichlet(1, alp1[v,])
if(v==m)
break
v=v+1
}
alppt=alp1
alpp=prop4*eps1+1
log_dirp1=rep(0,m)
log_dirp2=rep(0,m)
for(l in 1:m){
log_dirp1[l]=log(apply(matrix(output_p_s[t,,][l,], nrow=1), 1, function(u)den_dirichlet(u, alpp[l,])))
log_dirp2[l]=log(apply(matrix(prop4[l,], nrow=1), 1, function(u)den_dirichlet(u, alppt[l,])))
}
###############
comp_mean=sweep(output_gamma_s[t+1,,],MARGIN=1,output_lambda_s[t+1,],`*`)/prop4
                            
v=1
mix_prob=matrix(0,ncol=m, nrow=n)
repeat{
mixprob=apply(cbind(comp_mean[v,], prop4[v,]), 1, function(u)u[2]*dpois(xobs, lambda=u[1], log=FALSE))
mixprob=rowSums(mixprob)
mix_prob[,v]=mixprob
if(v==m)
break
v=v+1
}
                            
propl=mix_prob
lpropl=log(propl)
slpropl=colSums(lpropl, na.rm=FALSE, dims=1)
proplik=slpropl+log(apply(prop4, 1, function(u)den_dirichlet(u, alpha0)))
                            
log_ro_p=betta*proplik-betta*curlik+log_dirp1-log_dirp2
Accept_p=(log(runif(length(betta)))<(log_ro_p))
                            
pp[Accept_p,]=prop4[Accept_p,]
output_p_s[t+1,,]=pp
                            
slcur[Accept_p]=slpropl[Accept_p]
naccept3[Accept_p]=naccept3[Accept_p]+1
na3[Accept_p[1]]=na3[Accept_p[1]]+1
                            
if(a>300){men_adapt_p=median(adapt_p_ra[(a-10):a])
cond=(men_adapt_p>.2&men_adapt_p<.28)
 adapt_condition_p[cond==TRUE]=FALSE
}
                            
if(adapt_condition_p==TRUE){
ac_p_ra=na3/50
adapt_p_ra[a][a_adapt==50]=ac_p_ra
ac_cond_p1=(ac_p_ra<.234)
ac_cond_p2=(ac_p_ra>.234)
                                        
coco1=(eps1<=5)
coco2=(eps1>5)
eps1[ac_cond_p1&coco2&a_adapt==50]=eps1[ac_cond_p1&coco2&a_adapt==50]+5
eps1[ac_cond_p2&coco2&a_adapt==50]=eps1[ac_cond_p2&coco2&a_adapt==50]-5
eps1[ac_cond_p1&coco1&a_adapt==50]=eps1[ac_cond_p1&coco1&a_adapt==50]+5
eps1[ac_cond_p2&coco1&a_adapt==50]=eps1[ac_cond_p2&coco1&a_adapt==50]-eps1[ac_cond_p2&coco1&a_adapt==50]/100
adapt_eps_p[a][a_adapt==50]=eps1
a_p[a_adapt==50]=a_p[a_adapt==50]+1}else{adapt_p_ra[a][a_adapt==50]=adapt_p_ra[a-1][a_adapt==50]
adapt_eps_p[a][a_adapt==50]=adapt_eps_p[a-1][a_adapt==50]}

############## component means:
comp_mean=sweep(output_gamma_s[t+1,,],MARGIN=1,output_lambda_s[t+1,],`*`)/output_p_s[t+1,,]
component_mean[t+1,,]=comp_mean
 # Permutation of the labeling:
if(rept==2){
	perm <- sample(k)
output_p_s[t,1,perm]=output_p_s[t,1,]
output_p_s[t,2,perm]=output_p_s[t,2,]
output_gamma_s[t,1,perm]=output_gamma_s[t,1,]
output_gamma_s[t,2,perm]=output_gamma_s[t,2,]
component_mean[t,1,perm]=component_mean[t,1,]
component_mean[t,2,perm]=component_mean[t,2,]
}                               
na1[a_adapt==50]=0
na3[a_adapt==50]=0
na5[a_adapt==50]=0
a[a_adapt==50]=a[a_adapt==50]+1
a_adapt[a_adapt==50]=0
stop_condition1=(adapt_condition_lambda==FALSE & adapt_condition_gamma==FALSE & adapt_condition_p==FALSE)
if(rept==1 & stop_condition1==TRUE)
break
                                    
}
        
if(rept==2)
break
        
adapt_condition=FALSE
adapt_condition_lambda=adapt_condition_gamma=adapt_condition_p=FALSE
betta=rep(1, 2)
T=Nsim
m=length(betta)
p0=output_p_s[t+1,,][1:m,]
gamma0=output_gamma_s[t+1,,][1:m,]

naccept1=rep(0,m)
naccept2=rep(0,m)
naccept3=rep(0,m)
naccept5=rep(0,m)
        
n=length(xobs)
meanobs=mean(xobs)
        
mui=rep(mean(xobs), m)
output_lambda_s=matrix(mean(xobs), nrow=T, ncol=m)
        
pp=matrix(p0, nrow=m, ncol=k)
output_p_s=array(p0, dim=c(T,m,k))
for(i in 1:m){
   for(j in 1:k){
   output_p_s[,i,j]=p0[i,j]
}
}

gama=matrix(gamma0, nrow=m, ncol=k)
output_gamma_s=array(gamma0, dim=c(T,m,k))
for(i in 1:m){
    for(j in 1:k){
        output_gamma_s[,i,j]=gamma0[i,j]
    }
}

 component_mean=array(0, dim=c(T,m,k))
 rept=rept+1
}
adapt_lambda_ra=adapt_lambda_ra[1:a_lambda]
adapt_p_ra=adapt_p_ra[1:a_p]
adapt_gamma_ra=adapt_gamma_ra[1:a_gamma]
adapt_sd_lambda=adapt_sd_lambda[1:a_lambda]
adapt_eps_p=adapt_eps_p[1:a_p]
adapt_eps_gamma=adapt_eps_gamma[1:a_gamma]
adapt_rat=list(adapt_lambda_ra, adapt_p_ra, adapt_gamma_ra)
adapt_scale=list(adapt_sd_lambda, adapt_eps_p, adapt_eps_gamma)
output_lambda_s=output_lambda_s[1:t,]
output_p_s=output_p_s[1:t,,]
output_gamma_s=output_gamma_s[1:t,,]
component_mean=component_mean[1:t,,]
n1=naccept1/t
n3=naccept3/t
n5=naccept5/t
l=which.max(n1)
mean_global=output_lambda_s[,l] #[!output_mu_s[,l] %in% boxplot.stats(output_mu_s[,l])$out]
weights=output_p_s[,l,]
gammas=output_gamma_s[,l,]
epsilon_p=eps1
epsilon_gamma=eps2

component_means=component_mean[,l,]
 optimal_para=c(sd_lambda, epsilon_p, epsilon_gamma)
 accept_rat=c(n1[l], n3[l], n5[l])
 estimate=list(mean_global, weights, gammas, accept_rat, optimal_para, adapt_rat, adapt_scale, component_means)
 names(estimate)=c("mean_global", "weights", "gammas", "accept_rat", "optimal_para", "adapt_rat", "adapt_scale", "component_means")

    return(estimate)
}
