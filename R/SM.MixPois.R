SM.MixPois <-
function(estimate, xobs){
    k=dim(estimate[[3]])[2]
    mean_global=estimate[[1]][!estimate[[1]] %in% boxplot.stats(estimate[[1]])$out]
    # Global lambda:
    Mixture_mean=list(c(summary(mean_global)[3], summary(mean_global)[4]))
    Mixture_mean=data.frame(Mixture_mean)
    nameM=c("lambda: Mean of mixture distribution")
    names(Mixture_mean)=nameM
    
    # component means, gammas and weights:
    T=dim(estimate[[8]])[1]
        mus=estimate[[8]]
        gammas=estimate[[3]]
        ps=estimate[[2]]
        
        perm_mus=matrix(0, nrow=T, ncol=k)
        perm_gammas=matrix(0, nrow=T, ncol=k)
        perm_p=matrix(0, nrow=T, ncol=k)
        
        v=1
        t=1
        repeat{
            x=rbind(mus[t,], gammas[t,], ps[t,])
            y=matrix(0, ncol=k, nrow=3)
            l=1
            repeat{
                MIN=which.min(x[3,])
                y[,l]=x[,MIN]
                x=x[,-MIN]
                cond=(is.vector(x)==TRUE)
                if(cond==TRUE)
                break
                l=l+1
            }
            y[,k]=x
            perm_mus[t,]=y[1,]
            perm_gammas[t,]=y[2,]
            perm_p[t,]=y[3,]
            
            if(t==dim(mus)[1]-1)
            break
            t=t+1
        }
        
        mean_p=colSums(perm_p)/T
        f=1
        repeat{
            cond=(mean_p[f]>.06)
            if(cond==TRUE |f==k)
            break
            f=f+1
        }
        f=f-1
        
        k_prim=k-f
        new_perm_mus=perm_mus[1:t,(f+1):k]
        new_perm_gammas=perm_gammas[1:t,(f+1):k]
        new_perm_ps=perm_p[1:t,(f+1):k]
        
        dada=cbind(c(new_perm_mus), c(new_perm_ps))
        results=kmeans(dada, k_prim)
        
        A=cbind(c(new_perm_mus), c(new_perm_gammas), results$cluster, results$cluster)
        S=length(c(new_perm_mus))
        m1=s1=matrix(0, ncol=k_prim, nrow=S)
        means_cluster=gammas_cluster=p_cluster=list(0)
        l2=rep(1,k_prim)
        l3=rep(1,k_prim)
        
        for(i in 1:S){
            j=1
            repeat{
                if(A[i,3]==j){m1[l2[j],j]=A[i,1];l2[j]=l2[j]+1}
                if(A[i,4]==j){s1[l3[j],j]=A[i,2];l3[j]=l3[j]+1}
                cond=(j==k_prim)
                if(cond==TRUE)
                break
                j=j+1
            }
        }
        l2=l2-1
        l3=l3-1
        for(j in 1:k_prim){
            means_cluster[[j]]=m1[1:l2[j],j]
            gammas_cluster[[j]]=s1[1:l3[j],j]
        }
        
        summary_means_cluster=summary_p_cluster=summary_gammas_cluster=list(0)
        if(f>0){
            for(j in 1:f){
                summary_means_cluster[[j]]=summary(perm_mus[1:t,j])[3:4]
                summary_gamma_cluster[[j]]=summary(perm_gammas[1:t,j])[3:4]
                summary_p_cluster[[j]]=summary(perm_p[1:t,j])[3:4]
            }
        }
        for(j in 1:k_prim){
            summary_means_cluster[[(f+j)]]=summary(means_cluster[[j]])[3:4]
            summary_gammas_cluster[[(f+j)]]=summary(gammas_cluster[[j]])[3:4]
            summary_p_cluster[[(f+j)]]=summary(perm_p[1:t,(f+j)])[3:4]
        }
        
        # Max likelihood:
        if(k>10){n_permutation=permutations(n=10,r=10,v=1:10)}else{n_permutation=permutations(n=k,r=k,v=1:k)}
        PERM=dim(n_permutation)[1]+1
        like_mean=rep(0, PERM)
        mean_weights=matrix(0, ncol=k, nrow=PERM)
        
        median_weights=matrix(0, ncol=k, nrow=PERM)
 
        mean_weights[1,]=pps_mean=matrix(unlist(summary_p_cluster), ncol=2, byrow=TRUE)[,2]
        means_mean=matrix(unlist(summary_means_cluster), ncol=2, byrow=TRUE)[,2]
        median_weights[1,]=pps_median=matrix(unlist(summary_p_cluster), ncol=2, byrow=TRUE)[,1]
        means_median=matrix(unlist(summary_means_cluster), ncol=2, byrow=TRUE)[,1]
        gamas_mean=matrix(unlist(summary_gammas_cluster), ncol=2, byrow=TRUE)[,2]
        gamas_median=matrix(unlist(summary_gammas_cluster), ncol=2, byrow=TRUE)[,1]
        
        f_mean = f_median=rep(0,length(xobs),1)
            
            for (i in 1:length(pps_mean)){f_mean = f_mean + pps_mean[i]*dpois(xobs,lambda=means_mean[i],log=FALSE)}
            like_mean[1]=sum(f_mean)
        
        v=2
        repeat{
            f_mean = f_median=rep(0,length(xobs),1)
            y=n_permutation[(v-1),]
            pps_mean=pps_mean[y]
            pps_median=pps_median[y]
            
            for (i in 1:length(pps_mean)){f_mean = f_mean + pps_mean[i]*dpois(xobs,lambda=means_mean[i],log=FALSE)}
            like_mean[v]=sum(f_mean)
            mean_weights[v,]=pps_mean
            
            median_weights[v,]=pps_median
            if(v==PERM)
            break
            v=v+1
        }
        
        MT= which.max(like_mean)
        median_weights[MT,1]=1-sum(median_weights[MT,2:k])
        mean_weights[MT,1]=1-sum(mean_weights[MT,2:k])
for(w in 1:k){
            summary_p_cluster[[w]]=c(median_weights[MT,w], mean_weights[MT,w])
            summary_means_cluster[[w]]=c(means_median[w], means_mean[w])
            summary_gammas_cluster[[w]]=c(gamas_median[w], gamas_mean[w])
            
            names(summary_means_cluster[[w]])=c("Median", "mean")
            names(summary_p_cluster[[w]])=c("Median", "mean")
            names(summary_gammas_cluster[[w]])=c("Median", "mean")
        }
        
        
        names(summary_means_cluster)=rep("lambda",k)
        names(summary_gammas_cluster)=rep("gamma",k)
        names(summary_p_cluster)=rep("weight",k)
        
        # Acceptance rate of proposals and optimal scales:
            Acc_rat=estimate[[4]]
            nameA=c("lambda", "p", "gamma")
            names(Acc_rat)=nameA
            Acc_rat=list(Acc_rat)
            nameA=c("Acceptance rate of proposals")
            names(Acc_rat)=nameA
            
            Opt_scale=estimate[[5]]
            nameO=c("s_lambda", "epsilon_p", "epsilon_gamma")
            names(Opt_scale)=nameO
            Opt_scale=list(Opt_scale)
            nameO=c("Optimal proposal scales")
            names(Opt_scale)=nameO
        
        print(Mixture_mean)
        cat("", iter <- c("##############################"), "\n")
        cat("Component means, standard deviations and weights:", iter <- c(""), "\n")
        cat("", iter <- c(""), "\n")
        print(data.frame(summary_p_cluster))
        cat("", iter <- c(""), "\n")
        print(data.frame(summary_means_cluster))
        cat("", iter <- c(""), "\n")
        print(data.frame(summary_gammas_cluster))
        cat("", iter <- c("##############################"), "\n")
        print(Acc_rat)
        cat("", iter <- c("##############################"), "\n")
        print(Opt_scale) 
        
        return(list(data.frame(summary_p_cluster), data.frame(summary_means_cluster)))
}
