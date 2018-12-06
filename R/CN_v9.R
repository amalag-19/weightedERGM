#########################################################################################################
########################################## All Weighted code ############################################
#########################################################################################################
## Loading the required packages
library(foreach)
library(doParallel)
library(locfit)
library(Rcpp)
library(RcppArmadillo)

#########################################################################################################
######################################### All estimation code ########################################### 
#########################################################################################################
## Defining a wrapper function for 'Non-Parametric' weighted static undirected density case with block densities estimated by local likelihood with kernel approach (Loader). This is the generalized code that differentiates between adjacency matrix and weight matrix 
wrapper_weighted_stat_undir<-function(net_adjacency,net_weight,nclust,thres=10^(-6),theta_init,sim_indicator,theta_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(MASS)
  library(quadprog)
  library(combinat)
  library(locfit)
  library(Rcpp)
  library(RcppArmadillo)
  sourceCpp(file = "/Users/Amal/Box Sync/PSU/Fall 2017/Main_Research/Network Models/Project 4 (Weighted)/code/CONTundirv6/src/rcpparma_hello_world.cpp")
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,block_dens_mat,net_adjacency,N,K){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    
    quad_lin_coeff<-gamma_update_weighted_stat_undir(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, block_dens_mat=block_dens_mat, net_adjacency=net_adjacency, N=N, K=K)
    
    solve_QP_wrapper<-function(quad_lin_coeff,constraint_matrix,constraint_vector,node_ID,gamma.curr){
      try_QP<-try(gamma_next_vec<-solve.QP((-2*diag(as.vector(quad_lin_coeff[node_ID,,1]))),as.vector(quad_lin_coeff[node_ID,,2]),t(constraint_matrix),constraint_vector)$solution,silent = TRUE)
      if('numeric' %in% class(try_QP)){
        gamma_next_vec<-try_QP
      }
      else{
        gamma_next_vec<-gamma.curr[node_ID,]
        print("red flag 1")
      }
      return(gamma_next_vec)
    }
    
    #print(quad_lin_coeff)
    for (i in 1:N){
      if(sum(is.nan(quad_lin_coeff[i,,1]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,1])),1]<--Inf
      }
      if(sum(is.nan(quad_lin_coeff[i,,2]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,2])),2]<-Inf
      }
      if(sum(quad_lin_coeff[i,,1]<=-10^8)==0){
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K){
          gamma.next[i,]<-gamma.curr[i,]
        }
      }else if((sum(quad_lin_coeff[i,,1]<=-10^8)>0)&&(sum(quad_lin_coeff[i,,1]<=-10^8)<K)){
        quad_lin_coeff[i,which(quad_lin_coeff[i,,1]<=-10^8),1]<--(10^10)
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K){
          gamma.next[i,]<-gamma.curr[i,]
        }
      }else if(sum(quad_lin_coeff[i,,1]<=-10^8)==K){
        gamma.next[i,]<-gamma.curr[i,]
      }
      #print(i)
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,gamma,net_adjacency,N,K){
    theta.next<-rep(NA_real_,K)
    gradient<-grad_theta_weighted_stat_undir(theta=theta.curr, gamma=gamma, net_adjacency=net_adjacency, N=N, K=K)
    hess<-hess_theta_weighted_stat_undir(theta=theta.curr, gamma=gamma, N=N, K=K)
    try_hess_inv<-try(hess_inv<-solve(hess),silent = TRUE)
    if('matrix' %in% class(try_hess_inv)){
      theta.next<-theta.curr-(try_hess_inv%*%gradient)
    }
    else{
      theta.next<-theta.curr-(ginv(hess)%*%gradient)
      print("red flag 2")
    }
    
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,net_adjacency,net_weight,K,n_iter,thres){
    N<-dim(net_adjacency)[1]
    
    ## storing the non-zero weights in a single vector
    weights<-c()
    weights_ids<-data.frame("row_id"=NA_integer_,"col_id"=NA_integer_)
    k<-1
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        if(net_adjacency[i,j]!=0){
          weights<-c(weights,net_weight[i,j])
          weights_ids[k,]<-c(i,j)
          k<-k+1
        }
      }
    }
    
    ## Initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-matrix(NA_real_,K,n_iter)
    block_dens<-array(NA_real_,dim=c(N,N,K,K,n_iter))
    
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,1]<-start[[3]]
    block_dens[,,,,1]<-start[[4]]
    block_dens_mat<-matrix(NA_real_,N*K,N*K)
    for (k in 1:K){
      for (l in k:K){
        for (i in 1:(N-1)){
          for (j in (i+1):N){
            block_dens_mat[(k-1)*N+i,(l-1)*N+j]<-block_dens[i,j,k,l,1]
          }
        }
      }
    }
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<150)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      
      ## Variational E step:
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr = theta[,iter_index-1],block_dens_mat = block_dens_mat,net_adjacency=net_adjacency,N=N,K=K)
      #print("gamma")
      
      ## M step:
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      
      ## Updating the theta vector
      theta[,iter_index]<-theta.update(theta.curr = theta[,iter_index-1],gamma = gamma[,,iter_index],net_adjacency=net_adjacency,N = N,K = K)
      print(theta[,iter_index])
      
      ## Getting the cluster memberships of all nodes from variational E step
      clust_est<-as.vector(apply(X = gamma[,,iter_index],MARGIN = 1,FUN = which.max))
      
      ## parttioning the ties into pairs of clusters and string in a list (output of a field in R)
      weights_list<-tie_clust_partition(clust_est = clust_est,net_adjacency = net_adjacency,net_weight = net_weight,N = N,K = K)
      
      ## Updating the block densities
      for (l in 1:K){
        for (k in 1:l){
          weights_kl<-weights_list[[(l-1)*K+k]]
          #dens_est_non_zero<-density.lf(x =  weights_non_zero_kl, ev=weights_non_zero)$y
          if(sum(!is.na(weights_kl))>1){
            dens_est<-density.lf(x =  weights_kl, ev=weights)$y
            for (m in 1:length(dens_est)){
              block_dens[weights_ids[m,1],weights_ids[m,2],k,l,iter_index]<-dens_est[m]
            }
          }else{
            for (m in 1:nrow(weights_ids)){
              block_dens[weights_ids[m,1],weights_ids[m,2],k,l,iter_index]<-1/nrow(weights_ids)
            }
          }
          ##weights=gamma[weights_non_zero_ids[,1],k,iter_index]*gamma[weights_non_zero_ids[,2],l,iter_index]
        }
      }
      
      ## Converting the 4 dimensional array block dens into a large super matrix with outer matrix dimensions K*K, whose each       element is a submatrix of N*N.
      block_dens_mat<-matrix(NA_real_,N*K,N*K)
      for (k in 1:K){
        for (l in k:K){
          for (i in 1:(N-1)){
            for (j in (i+1):N){
              block_dens_mat[(k-1)*N+i,(l-1)*N+j]<-block_dens[i,j,k,l,iter_index]
            }
          }
        }
      }
      
      #print("block_dens")
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_weighted_stat_undir(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta=theta[,iter_index],block_dens_mat=block_dens_mat, net_adjacency=net_adjacency, N=N, K=K)
      #print(ELBO_grid.curr)
      if((ELBO_grid.curr==-Inf)|(ELBO_grid.prev==-Inf)){
        error<-Inf
      }else{
        error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      }
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta,block_dens[,,,,iter_index-1]))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,net_adjacency,N){
    gradient<-grad_theta_weighted_stat_undir_K1(theta=theta.curr, net_adjacency=net_adjacency, N=N)
    hess<-hess_theta_weighted_stat_undir_K1(theta=theta.curr, N=N)
    theta.next<-theta.curr-(gradient/hess)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,net_adjacency,net_weight,n_iter,thres){
    N<-dim(net_adjacency)[1]
    
    ## storing the non-zero weights in a single vector
    weights<-c()
    weights_ids<-data.frame("row_id"=NA_integer_,"col_id"=NA_integer_)
    k<-1
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        if(net_adjacency[i,j]!=0){
          weights<-c(weights,net_weight[i,j])
          weights_ids[k,]<-c(i,j)
          k<-k+1
        }
      }
    }
    
    ## initializing the arrays for parameters
    theta<-rep(NA_real_,n_iter)
    block_dens<-array(NA_real_,dim=c(N,N,n_iter))
    
    theta[1]<-start[[1]]
    block_dens[,,1]<-start[[2]]
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<150)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[iter_index]<-theta.update_K1(theta.curr=theta[iter_index-1], net_adjacency=net_adjacency, N=N)
      print(theta[iter_index])
      
      ## Updating the block density
      dens_est<-density.lf(x =  weights, ev=weights)$y
      ##weights=gamma[weights_non_zero_ids[,1],k,iter_index]*gamma[weights_non_zero_ids[,2],l,iter_index]
      for (m in 1:nrow(weights_ids)){
        block_dens[weights_ids[m,1],weights_ids[m,2],iter_index]<-dens_est[m]
      }
      
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_weighted_stat_undir_K1(theta = theta[iter_index], block_dens_mat=block_dens[,,iter_index], net_adjacency=net_adjacency, N=N)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(theta,block_dens))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  #################################################
  ## Defining a function to calculate the estimate of complete log likelihood
  comp_loglik<-function(gamma=NA,pi=NA,theta,net_adjacency,N,K,block_dens,cluster_ids_est=NA){
    if(K!=1){
      ## 1st term
      t1<-0
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          cluster_id_i<-cluster_ids_est[i]
          cluster_id_j<-cluster_ids_est[j]
          exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
          if(net_adjacency[i,j]!=0){
            t1<-t1+((theta[cluster_id_i]+theta[cluster_id_j])-(log(1+exp_val))+log(block_dens[i,j,min(cluster_id_i,cluster_id_j),max(cluster_id_i,cluster_id_j)]))
          }else if(net_adjacency[i,j]==0){
            t1<-t1-(log(1+exp_val))
          }
        }
      }
      ## 2nd term
      t2<-0
      for (i in 1:N){
        t2<-t2+log(pi[cluster_ids_est[i]])
      }
      comp_val<-t1+t2
    }else if(K==1){
      comp_val<-0
      exp_val<-exp(2*theta)
      log_exp_val<-log(1+exp_val)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          if(net_adjacency[i,j]!=0){
            comp_val<-comp_val+((2*theta)-log_exp_val+log(block_dens[i,j]))
          }else if(net_adjacency[i,j]==0){
            comp_val<-comp_val-log_exp_val
          }
        }
      }
    }
    return(comp_val)
  }
  
  #################################################
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,pi=NA,theta,net_adjacency,N,K,block_dens,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, pi = pi, theta = theta,net_adjacency = net_adjacency,N = N,K = K, block_dens=block_dens, cluster_ids_est = cluster_ids_est)
      t2<-(K-1)*log(N)
      t3<-K*log((N*(N-1))/2)
      ICL_val<-t1-t2-t3
    }else if(K==1){
      t1<-comp_loglik(theta = theta,net_adjacency = net_adjacency,N = N,K = K,block_dens=block_dens)
      t2<-log((N*(N-1))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(net_adjacency)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  K<-nclust ## Defining the number of clusters
  ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
  set.seed((2))
  MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = net_adjacency,alpha = 1/2,num.iterations = 100)
  gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
  ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1.
  for(i in 1:N){
    gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
    gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
  }
  #gamma.start<-matrix(rep(1/K,N*K),N,K)
  #################################################
  ## initializing the block densities
  if(K==1){
    block_dens_init<-matrix(NA_real_,N,N)
    ## Calculating non-zero weights
    weights<-c()
    weights_ids<-data.frame("row_id"=NA_integer_,"col_id"=NA_integer_)
    k<-1
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        if(net_adjacency[i,j]!=0){
          weights<-c(weights,net_weight[i,j])
          weights_ids[k,]<-c(i,j)
          k<-k+1
        }
      }
    }
    for(m in 1:nrow(weights_ids)){
      block_dens_init[weights_ids[m,1],weights_ids[m,2]]<-1/length(weights)
      block_dens_init[weights_ids[m,2],weights_ids[m,1]]<-block_dens_init[weights_ids[m,1],weights_ids[m,2]]
    }
  }else{
    block_dens_init<-array(NA_real_,dim=c(N,N,K,K))
    ## Calculating non-zero weights
    weights<-c()
    weights_ids<-data.frame("row_id"=NA_integer_,"col_id"=NA_integer_)
    k<-1
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        if(net_adjacency[i,j]!=0){
          weights<-c(weights,net_weight[i,j])
          weights_ids[k,]<-c(i,j)
          k<-k+1
        }
      }
    }
    
    for (k in 1:K){
      for(l in 1:K){
        for(m in 1:nrow(weights_ids)){
          block_dens_init[weights_ids[m,1],weights_ids[m,2],k,l]<-1/length(weights)
          block_dens_init[weights_ids[m,2],weights_ids[m,1],k,l]<-block_dens_init[weights_ids[m,1],weights_ids[m,2],k,l]
        }
      }
    }
  }
  
  
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start[[1]]<-theta_init
    start[[2]]<-block_dens_init
    param<-iterator_K1(start=start, net_adjacency = net_adjacency,net_weight = net_weight, n_iter=200, thres=thres)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## pi (mixture distribution)
    start[[3]]<-theta_init ## theta
    start[[4]]<-block_dens_init ## block_dens
    #debug(iterator)
    param<-iterator(start=start, net_adjacency=net_adjacency, net_weight = net_weight, K=K, n_iter=200, thres=thres)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[[1]][n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge[[1]]<-param[[1]][n_last]
    param_converge[[2]]<-param[[2]][,,n_last]
    ICL_val<-ICL(theta = param_converge[[1]],net_adjacency = net_adjacency,N = N,K = K,block_dens = param_converge[[2]])
    if(sim_indicator==1){
      if(K==K_true){
        debug(RASE_theta)
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,ICL_val,RASE_theta)
      }else{
        output_list<-list(param_converge,ICL_val)
      }
    }else{output_list<-list(param_converge,ICL_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    param_converge[[4]]<-param[[4]]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    ICL_val<-ICL(gamma = param_converge[[1]], pi=param_converge[[2]] ,theta = param_converge[[3]],net_adjacency = net_adjacency, N = N,K = K, block_dens = param_converge[[4]], cluster_ids_est = cluster_ids_est)
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
      print(RI_val)
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[3]][K_permute_mat[k,]]
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id_theta<-which.min(RASE_theta_vec)
        RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
        output_list<-list(param_converge,cluster_ids_est,ICL_val,RI_val,RASE_theta_val)
      }else{
        output_list<-list(param_converge,cluster_ids_est,ICL_val,RI_val)
      }
    }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
  }
  return(output_list)
}

#########################################################################################################
## Defining a wrapper function for 'Parametric Normal' weighted static undirected density case. This is the generalized code that differentiates between adjacency matrix and weight matrix 
wrapper_weighted_stat_undir_parametric_normal<-function(net_adjacency, net_weight, nclust,thres=10^(-6), theta_init, mu_init, sig2_init, sim_indicator, theta_true=NA, K_true=NA, cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(Rcpp)
  library(RcppArmadillo)
  sourceCpp(file = "/Users/Amal/Box Sync/PSU/Fall 2017/Main_Research/Network Models/Project 4 (Weighted)/code/CONTundirv6/src/rcpparma_hello_world.cpp")
  sourceCpp(file = "/Users/Amal/Box Sync/PSU/Fall 2017/Main_Research/Network Models/Project 4 (Weighted)/code/CONTundirv7/src/rcpparma_hello_world.cpp")
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,mu_mat.curr,sig2_mat.curr,net_adjacency,net_weight,N,K){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    
    quad_lin_coeff<-gamma_update_weighted_stat_undir_parametric_normal(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, mu_mat=mu_mat.curr, sig2_mat=sig2_mat.curr, net_adjacency=net_adjacency, net_weight=net_weight, N=N, K=K)
    
    solve_QP_wrapper<-function(quad_lin_coeff,constraint_matrix,constraint_vector,node_ID,gamma.curr){
      try_QP<-try(gamma_next_vec<-solve.QP((-2*diag(as.vector(quad_lin_coeff[node_ID,,1]))),as.vector(quad_lin_coeff[node_ID,,2]),t(constraint_matrix),constraint_vector)$solution,silent = TRUE)
      if('numeric' %in% class(try_QP)){
        gamma_next_vec<-try_QP
      }
      else{
        gamma_next_vec<-gamma.curr[node_ID,]
        print("red flag 1")
      }
      return(gamma_next_vec)
    }
    
    #print(quad_lin_coeff)
    for (i in 1:N){
      if(sum(is.nan(quad_lin_coeff[i,,1]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,1])),1]<--Inf
      }
      if(sum(is.nan(quad_lin_coeff[i,,2]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,2])),2]<-Inf
      }
      if(sum(quad_lin_coeff[i,,1]<=-10^8)==0){
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K){
          gamma.next[i,]<-gamma.curr[i,]
        }
      }else if((sum(quad_lin_coeff[i,,1]<=-10^8)>0)&&(sum(quad_lin_coeff[i,,1]<=-10^8)<K)){
        quad_lin_coeff[i,which(quad_lin_coeff[i,,1]<=-10^8),1]<--(10^10)
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K){
          gamma.next[i,]<-gamma.curr[i,]
        }
      }else if(sum(quad_lin_coeff[i,,1]<=-10^8)==K){
        gamma.next[i,]<-gamma.curr[i,]
      }
      #print(i)
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,gamma.curr,net_adjacency,net_weight,N,K){
    theta.next<-matrix(NA_real_,K,K)
    gradient<-grad_theta_weighted_stat_undir(theta=theta.curr, gamma=gamma.curr, net_adjacency=net_adjacency, N=N, K=K)
    hess<-hess_theta_weighted_stat_undir(theta=theta.curr, gamma=gamma.curr, N=N, K=K)
    theta.next<-theta.curr-(solve(hess)%*%gradient)
    return(theta.next)
  }
  
  mu_mat.update<-function(mu_mat.curr,sig2_mat.curr,gamma.curr,net_adjacency,net_weight,N,K){
    mu_mat.next<-matrix(NA_real_,K,K)
    gradient<-grad_mu_mat_weighted_stat_undir_parametric_normal(mu_mat=mu_mat.curr, sig2_mat=sig2_mat.curr, gamma=gamma.curr, net_adjacency=net_adjacency, net_weight=net_weight, N=N, K=K)
    hess<-hess_mu_mat_weighted_stat_undir_parametric_normal(sig2_mat=sig2_mat.curr, gamma=gamma.curr, net_adjacency=net_adjacency, N=N, K=K)
    for(k in 1:K){
      for(l in k:K){
        if(l==k){
          mu_mat.next[k,l]<-mu_mat.curr[k,l]-gradient[k,l]/hess[k,l]
        }else if(l!=k){
          mu_mat.next[k,l]<-mu_mat.curr[k,l]-((gradient[k,l]+gradient[l,k])/(hess[k,l]+hess[l,k]))
        }
      }
    }
    return(mu_mat.next)
  }
  
  sig2_mat.update<-function(mu_mat.curr,sig2_mat.curr,gamma.curr,net_adjacency,net_weight,N,K){
    sig2_mat.next<-matrix(NA_real_,K,K)
    gradient<-grad_sig2_mat_weighted_stat_undir_parametric_normal(mu_mat=mu_mat.curr, sig2_mat=sig2_mat.curr, gamma=gamma.curr, net_adjacency=net_adjacency, net_weight=net_weight, N=N, K=K)
    hess<-hess_sig2_mat_weighted_stat_undir_parametric_normal(mu_mat=mu_mat.curr, sig2_mat=sig2_mat.curr, gamma=gamma.curr, net_adjacency=net_adjacency, net_weight=net_weight, N=N, K=K)
    for(k in 1:K){
      for(l in k:K){
        if(l==k){
          sig2_mat.next[k,l]<-sig2_mat.curr[k,l]-gradient[k,l]/hess[k,l]
        }else if(l!=k){
          sig2_mat.next[k,l]<-sig2_mat.curr[k,l]-((gradient[k,l]+gradient[l,k])/(hess[k,l]+hess[l,k]))
        }
      }
    }
    return(sig2_mat.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,net_adjacency,net_weight,K,n_iter,thres){
    N<-dim(net_adjacency)[1]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-matrix(NA_real_,K,n_iter)
    mu<-array(NA_real_,dim=c(K,K,n_iter))
    sig2<-array(NA_real_,dim=c(K,K,n_iter))
    
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,1]<-start[[3]]
    mu[,,1]<-start[[4]]
    sig2[,,1]<-start[[5]]
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<150)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr = theta[,iter_index-1], mu_mat.curr = mu[,,iter_index-1], sig2_mat.curr = sig2[,,iter_index-1], net_adjacency=net_adjacency, net_weight = net_weight, N=N, K=K)
      
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      
      ## Updating the network parameters
      theta[,iter_index]<-theta.update(theta.curr = theta[,iter_index-1], gamma=gamma[,,iter_index], net_adjacency = net_adjacency, net_weight = net_weight, N = N, K = K)
      #print(theta[,iter_index])
      
      ## Updating the alpha matrix
      mu[,,iter_index]<-mu_mat.update(mu_mat.curr = mu[,,iter_index-1], sig2_mat.curr = sig2[,,iter_index-1], gamma = gamma[,,iter_index], net_adjacency = net_adjacency, net_weight = net_weight, N = N, K = K)
      #print(mu[,,iter_index])
      
      ## Updating the beta matrix
      sig2[,,iter_index]<-sig2_mat.update(mu_mat.curr = mu[,,iter_index], sig2_mat.curr = sig2[,,iter_index-1], gamma = gamma[,,iter_index], net_adjacency = net_adjacency, net_weight = net_weight, N = N, K = K)
      #print(sig2[,,iter_index])
      
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_weighted_stat_undir_parametric_normal(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta= theta[,iter_index], mu_mat=mu[,,iter_index], sig2_mat=sig2[,,iter_index], net_adjacency=net_adjacency, net_weight=net_weight, N=N, K=K)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
      #Sys.sleep(time = 2)
    }
    #print(theta[,iter_index-1])
    #print(mu[,,iter_index-1])
    #print(sig2[,,iter_index-1])
    return(list(gamma,pi,theta,mu,sig2))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,network,N){
    gradient<-grad_HMM_stat_undir_K1(theta=theta.curr, network=network, N=N)
    hess<-hess_HMM_stat_undir_K1(theta=theta.curr, N=N)
    theta.next<-theta.curr-(gradient/hess)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    ## initializing the arrays for parameters
    theta<-rep(NA_real_,n_iter)
    theta[1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<200)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[iter_index]<-theta.update_K1(theta.curr=theta[iter_index-1], network=network, N=N)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_stat_undir_K1(theta = theta[iter_index], network=network, N=N)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  #################################################
  ## Defining a function to calculate the estimate of complete log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,theta,network,N,K,cluster_ids_est=NA){
    if(K!=1){
      ## 1st term
      t1<-0
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          cluster_id_i<-cluster_ids_est[i]
          cluster_id_j<-cluster_ids_est[j]
          exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
          t1<-t1+((network[i,j]*(theta[cluster_id_i]+theta[cluster_id_j]))-(log(1+exp_val)))
        }
      }
      ## 2nd term
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i]])
      }
      comp_val<-t1+t2
    }else if(K==1){
      comp_val<-0
      exp_val<-exp(2*theta)
      log_exp_val<-log(1+exp_val)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          comp_val<-comp_val+((network[i,j]*(2*theta))-log_exp_val)
        }
      }
    }
    return(comp_val)
  }
  
  #################################################
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,theta,network,N,K,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, theta = theta,network = network,N = N,K = K, cluster_ids_est = cluster_ids_est)
      t2<-K*log((N*(N-1))/2)
      ICL_val<-t1-t2
    }else if(K==1){
      t1<-comp_loglik(theta = theta,network = network,N = N,K = K)
      t2<-log((N*(N-1))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(net_adjacency)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  K<-nclust ## Defining the number of clusters
  ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
  set.seed((2))
  MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = net_adjacency,alpha = 1/2,num.iterations = 100)
  gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
  ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1.
  for(i in 1:N){
    gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
    gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-theta_init
    param<-iterator_K1(start=start, net_adjacency = net_adjacency, net_weight = net_weight, K=K, n_iter=1000, thres=thres)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## alpha (initial distribution)
    start[[3]]<-theta_init ## weight
    start[[4]]<-mu_init ## alpha
    start[[5]]<-sig2_init ## beta
    #debug(iterator)
    param<-iterator(start=start, net_adjacency = net_adjacency, net_weight = net_weight, K=K, n_iter=1000, thres=thres)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[n_last]
    ICL_val<-ICL(theta = param_converge,network = sim.net,N = N,K = K)
    if(sim_indicator==1){
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,ICL_val,RASE_theta)
      }else{
        output_list<-list(param_converge,ICL_val)
      }
    }else{output_list<-list(param_converge,ICL_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    param_converge[[4]]<-param[[4]][,,n_last]
    param_converge[[5]]<-param[[5]][,,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    #ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]] ,theta = param_converge[[3]],network = sim.net,N = N,K = K, cluster_ids_est = cluster_ids_est)
    ICL_val<-0
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
      print(RI_val)
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        RASE_pi_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[3]][K_permute_mat[k,]]
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id_theta<-which.min(RASE_theta_vec)
        RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
        output_list<-list(param_converge,cluster_ids_est,ICL_val,RI_val,RASE_theta_val)
      }else{
        output_list<-list(param_converge,cluster_ids_est,ICL_val,RI_val)
      }
    }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
  }
  return(output_list)
}

#########################################################################################################
## Defining a wrapper function for 'Parametric Gamma' weighted static undirected density case. This is the generalized code that differentiates between adjacency matrix and weight matrix 
wrapper_weighted_stat_undir_parametric_gamma<-function(net_adjacency, net_weight, nclust, thres=10^(-6),theta_init,alpha_init,beta_init,sim_indicator,theta_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  library(Rcpp)
  library(RcppArmadillo)
  sourceCpp(file = "/Users/Amal/Box Sync/PSU/Fall 2017/Main_Research/Network Models/Project 4 (Weighted)/code/CONTundirv6/src/rcpparma_hello_world.cpp")
  sourceCpp(file = "/Users/Amal/Box Sync/PSU/Fall 2017/Main_Research/Network Models/Project 4 (Weighted)/code/CONTundirv8/src/rcpparma_hello_world.cpp")
  
  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input: 
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,alpha_mat.curr,beta_mat.curr,net_adjacency,net_weight,N,K){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))
    
    quad_lin_coeff<-gamma_update_weighted_stat_undir_parametric_gamma(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, alpha_mat=alpha_mat.curr, beta_mat=beta_mat.curr, log_Gamma_alpha_mat=log(gamma(alpha_mat.curr)), net_adjacency=net_adjacency, net_weight=net_weight, N=N, K=K)
    solve_QP_wrapper<-function(quad_lin_coeff,constraint_matrix,constraint_vector,node_ID,gamma.curr){
      try_QP<-try(gamma_next_vec<-solve.QP((-2*diag(as.vector(quad_lin_coeff[node_ID,,1]))),as.vector(quad_lin_coeff[node_ID,,2]),t(constraint_matrix),constraint_vector)$solution,silent = TRUE)
      if('numeric' %in% class(try_QP)){
        gamma_next_vec<-try_QP
      }
      else{
        gamma_next_vec<-gamma.curr[node_ID,]
        print("red flag 1")
      }
      return(gamma_next_vec)
    }
    
    #print(quad_lin_coeff)
    for (i in 1:N){
      if(sum(is.nan(quad_lin_coeff[i,,1]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,1])),1]<--Inf
      }
      if(sum(is.nan(quad_lin_coeff[i,,2]))>0){
        quad_lin_coeff[i,which(is.nan(quad_lin_coeff[i,,2])),2]<-Inf
      }
      if(sum(quad_lin_coeff[i,,1]<=-10^8)==0){
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K){
          gamma.next[i,]<-gamma.curr[i,]
        }
      }else if((sum(quad_lin_coeff[i,,1]<=-10^8)>0)&&(sum(quad_lin_coeff[i,,1]<=-10^8)<K)){
        quad_lin_coeff[i,which(quad_lin_coeff[i,,1]<=-10^8),1]<--(10^10)
        if(sum(quad_lin_coeff[i,,2]>=10^8)==0){
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if((sum(quad_lin_coeff[i,,2]>=10^8)>0)&&(sum(quad_lin_coeff[i,,2]>=10^8)<K)){
          quad_lin_coeff[i,which(quad_lin_coeff[i,,2]>=10^8),2]<-10^10
          gamma.next[i,]<-solve_QP_wrapper(quad_lin_coeff = quad_lin_coeff,constraint_matrix = constraint_matrix,constraint_vector = constraint_vector,node_ID = i,gamma.curr = gamma.curr)
          #gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
        }else if(sum(quad_lin_coeff[i,,2]>=10^8)==K){
          gamma.next[i,]<-gamma.curr[i,]
        }
      }else if(sum(quad_lin_coeff[i,,1]<=-10^8)==K){
        gamma.next[i,]<-gamma.curr[i,]
      }
      #print(i)
    }
    
    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }
  
  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,gamma.curr,net_adjacency,net_weight,N,K){
    theta.next<-matrix(NA_real_,K,K)
    gradient<-grad_theta_weighted_stat_undir(theta=theta.curr, gamma=gamma.curr, net_adjacency=net_adjacency, N=N, K=K)
    hess<-hess_theta_weighted_stat_undir(theta=theta.curr, gamma=gamma.curr, N=N, K=K)
    theta.next<-theta.curr-(solve(hess)%*%gradient)
    return(theta.next)
  }
  
  alpha_mat.update<-function(alpha_mat.curr,gamma.curr,beta_mat.curr,net_adjacency,net_weight,N,K){
    alpha_mat.next<-matrix(NA_real_,K,K)
    gradient<-grad_alpha_mat_weighted_stat_undir_parametric_gamma(gamma=gamma.curr,beta_mat=beta_mat.curr, Digamma_alpha_mat=digamma(alpha_mat.curr), net_adjacency=net_adjacency, net_weight=net_weight, N=N, K=K)
    hess<-hess_alpha_mat_weighted_stat_undir_parametric_gamma(gamma=gamma.curr, Trigamma_alpha_mat=trigamma(alpha_mat.curr), net_adjacency=net_adjacency, N=N, K=K)
    for(k in 1:K){
      for(l in k:K){
        if(l==k){
          alpha_mat.next[k,l]<-alpha_mat.curr[k,l]-gradient[k,l]/hess[k,l]
        }else if(l!=k){
          alpha_mat.next[k,l]<-alpha_mat.curr[k,l]-((gradient[k,l]+gradient[l,k])/(hess[k,l]+hess[l,k]))
        }
      }
    }
    return(alpha_mat.next)
  }
  
  beta_mat.update<-function(beta_mat.curr,gamma.curr,alpha_mat.curr,net_adjacency,net_weight,N,K){
    beta_mat.next<-matrix(NA_real_,K,K)
    gradient<-grad_beta_mat_weighted_stat_undir_parametric_gamma(gamma=gamma.curr, alpha_mat=alpha_mat.curr,  beta_mat=beta_mat.curr, net_adjacency=net_adjacency, net_weight=net_weight, N=N, K=K)
    hess<-hess_beta_mat_weighted_stat_undir_parametric_gamma(gamma=gamma.curr, alpha_mat=alpha_mat.curr,  beta_mat=beta_mat.curr, net_adjacency=net_adjacency, N=N, K=K)
    for(k in 1:K){
      for(l in k:K){
        if(l==k){
          beta_mat.next[k,l]<-beta_mat.curr[k,l]-gradient[k,l]/hess[k,l]
        }else if(l!=k){
          beta_mat.next[k,l]<-beta_mat.curr[k,l]-((gradient[k,l]+gradient[l,k])/(hess[k,l]+hess[l,k]))
        }
      }
    }
    return(beta_mat.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,net_adjacency,net_weight,K,n_iter,thres){
    N<-dim(net_adjacency)[1]
    
    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-matrix(NA_real_,K,n_iter)
    alpha<-array(NA_real_,dim=c(K,K,n_iter))
    beta<-array(NA_real_,dim=c(K,K,n_iter))
    
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,1]<-start[[3]]
    alpha[,,1]<-start[[4]]
    beta[,,1]<-start[[5]]
    
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<150)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr = theta[,iter_index-1], alpha_mat.curr = alpha[,,iter_index-1], beta_mat.curr = beta[,,iter_index-1], net_adjacency=net_adjacency, net_weight = net_weight, N=N, K=K)
      
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      
      ## Updating the network parameters
      theta[,iter_index]<-theta.update(theta.curr = theta[,iter_index-1], gamma=gamma[,,iter_index], net_adjacency = net_adjacency, net_weight = net_weight, N = N, K = K)
      #print(theta[,iter_index])
      
      ## Updating the alpha matrix
      alpha[,,iter_index]<-alpha_mat.update(alpha_mat.curr = alpha[,,iter_index-1], gamma.curr = gamma[,,iter_index], beta_mat.curr = beta[,,iter_index-1], net_adjacency = net_adjacency, net_weight = net_weight, N = N, K = K)
      #print(alpha[,,iter_index])
      
      ## Updating the beta matrix
      beta[,,iter_index]<-beta_mat.update(beta_mat.curr = beta[,,iter_index-1], gamma.curr = gamma[,,iter_index], alpha_mat.curr =  alpha[,,iter_index], net_adjacency = net_adjacency, net_weight = net_weight, N = N, K = K)
      
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_weighted_stat_undir_parametric_gamma(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta= theta[,iter_index], alpha_mat=alpha[,,iter_index], beta_mat=beta[,,iter_index], log_Gamma_alpha_mat=log(gamma(alpha[,,iter_index])), net_adjacency=net_adjacency, net_weight=net_weight, N=N, K=K)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    print(theta[,iter_index-1])
    print(alpha[,,iter_index-1])
    print(beta[,,iter_index-1])
    return(list(gamma,pi,theta,alpha,beta))
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  ## Defining the functions for K=1
  
  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,network,N){
    gradient<-grad_HMM_stat_undir_K1(theta=theta.curr, network=network, N=N)
    hess<-hess_HMM_stat_undir_K1(theta=theta.curr, N=N)
    theta.next<-theta.curr-(gradient/hess)
    return(theta.next)
  }
  
  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    ## initializing the arrays for parameters
    theta<-rep(NA_real_,n_iter)
    theta[1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<200)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[iter_index]<-theta.update_K1(theta.curr=theta[iter_index-1], network=network, N=N)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_stat_undir_K1(theta = theta[iter_index], network=network, N=N)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }
  
  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters
  
  #################################################
  ## Defining a function to calculate the estimate of complete log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,theta,network,N,K,cluster_ids_est=NA){
    if(K!=1){
      ## 1st term
      t1<-0
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          cluster_id_i<-cluster_ids_est[i]
          cluster_id_j<-cluster_ids_est[j]
          exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
          t1<-t1+((network[i,j]*(theta[cluster_id_i]+theta[cluster_id_j]))-(log(1+exp_val)))
        }
      }
      ## 2nd term
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i]])
      }
      comp_val<-t1+t2
    }else if(K==1){
      comp_val<-0
      exp_val<-exp(2*theta)
      log_exp_val<-log(1+exp_val)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          comp_val<-comp_val+((network[i,j]*(2*theta))-log_exp_val)
        }
      }
    }
    return(comp_val)
  }
  
  #################################################
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,theta,network,N,K,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, theta = theta,network = network,N = N,K = K, cluster_ids_est = cluster_ids_est)
      t2<-K*log((N*(N-1))/2)
      ICL_val<-t1-t2
    }else if(K==1){
      t1<-comp_loglik(theta = theta,network = network,N = N,K = K)
      t2<-log((N*(N-1))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }
  
  
  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }
  
  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }
  
  #########################################################################################################
  ## Defining the parameters
  N<-dim(net_adjacency)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  K<-nclust ## Defining the number of clusters
  ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
  set.seed((2))
  MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = net_adjacency,alpha = 1/2,num.iterations = 100)
  gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
  ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1.
  for(i in 1:N){
    gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
    gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
  }
  
  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()
  
  if(K==1){
    start<-theta_init
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## alpha (initial distribution)
    start[[3]]<-theta_init ## weight
    start[[4]]<-alpha_init ## alpha
    start[[5]]<-beta_init ## beta
    #debug(iterator)
    param<-iterator(start=start, net_adjacency = net_adjacency, net_weight = net_weight, K=K, n_iter=1000, thres=thres)
  }
  
  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[n_last]
    ICL_val<-ICL(theta = param_converge,network = sim.net,N = N,K = K)
    if(sim_indicator==1){
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,ICL_val,RASE_theta)
      }else{
        output_list<-list(param_converge,ICL_val)
      }
    }else{output_list<-list(param_converge,ICL_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    param_converge[[4]]<-param[[4]][,,n_last]
    param_converge[[5]]<-param[[5]][,,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    #ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]] ,theta = param_converge[[3]],network = sim.net,N = N,K = K, cluster_ids_est = cluster_ids_est)
    ICL_val<-0
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
      print(RI_val)
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        RASE_pi_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[3]][K_permute_mat[k,]]
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id_theta<-which.min(RASE_theta_vec)
        RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
        output_list<-list(param_converge,cluster_ids_est,ICL_val,RI_val,RASE_theta_val)
      }else{
        output_list<-list(param_converge,cluster_ids_est,ICL_val,RI_val)
      }
    }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
  }
  return(output_list)
}

#########################################################################################################
######################################### All simulation code ########################################### 
#########################################################################################################
## Defining a function to simulate network adjacency and weight matrices along with cluster memberships of nodes for 'Parametric Normal' weighted static undirected density case
simulate_network_weighted_stat_undir_normal<-function(N,pi,theta,mu_mat,sig2_mat){
  net_adjacency<-matrix(NA_integer_,N,N)
  net_weight<-matrix(NA_integer_,N,N)
  clust_id_sampling<-rmultinom(n = N, size = 1, prob = pi)
  clust_ids<-apply(X = clust_id_sampling,MARGIN = 2,FUN = function(x){
    y<-which(x==1)
    return(y)
  })
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      edge_ij<-exp(theta[clust_ids[i]]+theta[clust_ids[j]])/(1+exp(theta[clust_ids[i]]+theta[clust_ids[j]]))
      edge_indicator_0_sample<-rbinom(n = 1, size = 1, prob = edge_ij)
      if(edge_indicator_0_sample==0){
        net_adjacency[i,j]<-0
        net_adjacency[j,i]<-net_adjacency[i,j]
      }else if(edge_indicator_0_sample==1){
        net_adjacency[i,j]<-1
        net_adjacency[j,i]<-net_adjacency[i,j]
        net_weight[i,j]<-rnorm(n = 1,mean =  mu_mat[min(clust_ids[i],clust_ids[j]),max(clust_ids[i],clust_ids[j])],sd = sqrt(sig2_mat[min(clust_ids[i],clust_ids[j]),max(clust_ids[i],clust_ids[j])]))
        net_weight[j,i]<-net_weight[i,j]
      }
    }
  }
  diag(net_adjacency)<-0
  return(list(net_adjacency,net_weight,clust_ids))
}

#########################################################################################################
## Defining a function to simulate network adjacency and weight matrices along with cluster memberships of nodes for 'Parametric Gamma' weighted static undirected density case
simulate_network_weighted_stat_undir_gamma<-function(N,pi,theta,alpha_mat,beta_mat){
  net_adjacency<-matrix(NA_integer_,N,N)
  net_weight<-matrix(NA_integer_,N,N)
  clust_id_sampling<-rmultinom(n = N, size = 1, prob = pi)
  clust_ids<-apply(X = clust_id_sampling,MARGIN = 2,FUN = function(x){
    y<-which(x==1)
    return(y)
  })
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      edge_ij<-exp(theta[clust_ids[i]]+theta[clust_ids[j]])/(1+exp(theta[clust_ids[i]]+theta[clust_ids[j]]))
      edge_indicator_0_sample<-rbinom(n = 1, size = 1, prob = edge_ij)
      if(edge_indicator_0_sample==0){
        net_adjacency[i,j]<-0
        net_adjacency[j,i]<-net_adjacency[i,j]
      }else if(edge_indicator_0_sample==1){
        net_adjacency[i,j]<-1
        net_adjacency[j,i]<-net_adjacency[i,j]
        net_weight[i,j]<-rgamma(n = 1,shape = alpha_mat[min(clust_ids[i],clust_ids[j]),max(clust_ids[i],clust_ids[j])],rate = beta_mat[min(clust_ids[i],clust_ids[j]),max(clust_ids[i],clust_ids[j])])
        net_weight[j,i]<-net_weight[i,j]
      }
    }
  }
  diag(net_adjacency)<-0
  return(list(net_adjacency,net_weight,clust_ids))
}

#########################################################################################################
## Defining a function to simulate network adjacency and weight matrices along with cluster memberships of nodes for 'Parametric t (Non central)' weighted static undirected density case
simulate_network_weighted_stat_undir_t<-function(N,pi,theta,df_mat,ncp_mat){
  net_adjacency<-matrix(NA_integer_,N,N)
  net_weight<-matrix(NA_integer_,N,N)
  clust_id_sampling<-rmultinom(n = N, size = 1, prob = pi)
  clust_ids<-apply(X = clust_id_sampling,MARGIN = 2,FUN = function(x){
    y<-which(x==1)
    return(y)
  })
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      edge_ij<-exp(theta[clust_ids[i]]+theta[clust_ids[j]])/(1+exp(theta[clust_ids[i]]+theta[clust_ids[j]]))
      edge_indicator_0_sample<-rbinom(n = 1, size = 1, prob = edge_ij)
      if(edge_indicator_0_sample==0){
        net_adjacency[i,j]<-0
        net_adjacency[j,i]<-net_adjacency[i,j]
      }else if(edge_indicator_0_sample==1){
        net_adjacency[i,j]<-1
        net_adjacency[j,i]<-net_adjacency[i,j]
        net_weight[i,j]<-rt(n = 1,df = df_mat[min(clust_ids[i],clust_ids[j]),max(clust_ids[i],clust_ids[j])],ncp = ncp_mat[min(clust_ids[i],clust_ids[j]),max(clust_ids[i],clust_ids[j])])
        net_weight[j,i]<-net_weight[i,j]
      }
    }
  }
  diag(net_adjacency)<-0
  return(list(net_adjacency,net_weight,clust_ids))
}

#########################################################################################################
#########################################################################################################
#########################################################################################################
