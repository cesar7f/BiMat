#This functions get the modularity matrix
modularity.modularity_matrix <- function(matrix) {
  
  obj <- list()

  obj$adjacency <- matrix
  obj$n_nodes <- nrow(obj$adjacency)
  obj$kk <- apply(obj$adjacency,2,sum)
  obj$n_edges <-  sum(obj$kk)/2
  
  if(obj$n_edges == 0){obj$bb_matrix = obj$adjacency}
  else{
    obj$bb_matrix <- obj$adjacency - (obj$kk %*% t(obj$kk))/(2 * obj$n_edges)
  }
  
  return(obj)
}

#This function calculates the leading eigenvector algorithm of a adjacency matrix or igraph object
modularity.leading.eigenvector.commmunities <- function(adjacency_or_igraph_object)
{
  
  if(class(adjacency_or_igraph_object) == 'igraph'){
    adjacency_matrix <- as.matrix(get.adjacency(adjacency_or_igraph_object,attr='weight'))
  }
  else{
    adjacency_matrix <- as.matrix(adjacency_or_igraph_object)
  }
  
  obj <- modularity.modularity_matrix(adjacency_matrix)
  obj <- modularity.main_algorithm(obj)
  
  obj$membership <- as.vector(obj$membership)
  obj$modularity <- as.double(obj$Q)
  obj$vcount <- as.double(obj$n_nodes)

  
  return(obj)
  
}


#Main algorithm method. It start with a single module. Then, for each existing module, it finds the biggest
#local eigenvalue, which eigenvector is used to subdivide such module only if the total global modularity increase
modularity.main_algorithm <- function(obj){
  
  obj$N <- 1
  obj$Q <- 0
  
  obj$membership <- matrix(1,obj$n_nodes, 1)
  
  qinc <- 1
  
  divisible = TRUE
  
  while(qinc == 1){
    
    qinc <- 0
    nn <- obj$N
    
    #For each module, try to divide in two
    for (i in 1:nn){
      
      #Avoid check something that is not divisible
      if(divisible[i] == FALSE){ next;}
      
      nindices <- c(which(obj$membership == i))
      nsize <- length(nindices)
      blocal <- obj$bb_matrix[nindices,nindices]
      #browser()
      diag(blocal) <- diag(blocal) - apply(blocal,1,sum)
      
      #browser()
      #Take the vector with the biggest eigenvalue
      max_eigenvalue_vector <- eigen(blocal)$vectors[,1]
      
      #Use the vector to divide the module in two
      membership_local <- matrix(1,nsize,1)
      membership_local[which(max_eigenvalue_vector <= 0)] <- -1
      
      #Increase of the total global modularity by the division of the module
      deltaQ <- (t(membership_local) %*% blocal %*% membership_local)/(4*obj$n_edges)
      
      #Divide only if increase in modularity is detected
      if(deltaQ > 0 && length(which(membership_local==-1))>0 ){
        qinc <- 1
        newind <- which(membership_local == 1) 
        obj$membership[nindices[newind]] <- obj$N + 1
        obj$N <- obj$N + 1; 
        obj$Q <- obj$Q + deltaQ
        divisible = c(divisible,TRUE)
        print(sprintf('%3d %3d %.16f ', i, nn, obj$Q))
        
        #if(obj$N == 18) browser();
        #break;
      }
      #Otherwise the module is not divisible and will not be tested again
      else
      {
        divisible[i] <- FALSE
      }
      
    }
  }
  return (obj)
}

#Function to make heat map of the adjacency matrix in modularity sorting
modularity.plot_modularity <- function(obj)
{
  sorted_index <- order(obj$membership)
  matrix_to_plot <- obj$adjacency[sorted_index,sorted_index]
  image(matrix_to_plot[,obj$n_nodes:1],axes=FALSE,col=topo.colors(120))
}

#This function is not used by the Newman algorithm. I just coded to test bipartite networks too.
#If you want to test bipartite matrices, this function converts them to unipartite, such thay you can apply
#the newman algorithm directly on them
bipartite_to_unipartite <- function(bip_matrix){
  
  bip_matrix = as.matrix(bip_matrix)
  n_rows <- nrow(bip_matrix)
  n_cols <- ncol(bip_matrix)
  
  uni_matrix <- rbind(  cbind( matrix(0,n_rows,n_rows),bip_matrix),  cbind( t(bip_matrix), matrix(0,n_cols,n_cols)  ))
  
  return(uni_matrix)
  
}


#matrix = read.csv('moebus_matrix.txt',sep=' ',header = FALSE)
#matrix = as.matrix(matrix)