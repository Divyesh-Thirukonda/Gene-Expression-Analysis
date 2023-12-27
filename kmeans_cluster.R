kmeans_cluster <- function(data, k){
  # We want the gene expression values to be the columns, conditions to be rows
  data <- t(data)

  # Randomize the initial gene cluster assignments
  cluster_labels <- sample(1:k, size=ncol(data), replace=TRUE)
  # initialize the center at 0
  center <- matrix(rep(0,nrow(data)*k), ncol=k)
  # initialize the distance vector to 0
  dis <- rep(0, times=k)
  cost <- 0
  # For each gene expr (col) add distance formula result to get the total cost
  for(i in 1:ncol(data)){
    cost <- cost + sqrt(sum((data[,i]-center[,cluster_labels[i]])^2,na.rm=T))
  }
  print(paste('Initial Cost ',cost,sep=''))
  count <- 0

  # loop until break cmd
  while(TRUE){

    # reset cost
    cost_old <- cost
    
    for(i in 1:k){

      # if there are no values of i
      if(sum(cluster_labels == i) == 0){
        # keep same
        center[,i] <- center[,i]
      }
      # otherwise if all are i
      else if(sum(cluster_labels == i) == 1){
        # set the center columns to the data columns where it's i
        center[,i] <- data[,cluster_labels == i]
      }
      else{
        # otherwise, set the cetner columns to the condition's average where it's i
        center[,i] <- rowMeans(data[,cluster_labels == i])
      }
    }
    # for every gene
    for(i in 1:ncol(data)){
      # for every centroid
      for(j in 1:ncol(center)){
        # set matrix value's distance to the difference squared
        dis[j] <- sqrt(sum((data[,i]-center[,j])^2,na.rm=T))
      }
      # assign into closest cluster
      cluster_labels[i] <- which.min(dis)
    }
    
    cost <- 0
    for(i in 1:ncol(data)){
      # for each gene, add the costs using the distance formula and recalculate the total
      cost <- cost + sqrt(sum((data[,i]-center[,cluster_labels[i]])^2,na.rm=T))
    }
    print(paste('Iteration: ',count,' -- cost:',cost,sep=''))
    # if the difference is smaller than .0001, break out (stopping criterion)
    if(abs(cost - cost_old) < 0.0001){
      break
    }
    # otherwise, increase the count
    count <- count + 1
  }
  return(cluster_labels)
}