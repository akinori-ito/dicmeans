#
# Succesive Binary Partition k-means clustering algorithm
#
# 2021/11/14 Akinori Ito
#

renumber <- function(x) {
  u <- unique(x)
  r <- list()
  for (i in 1:length(u)) {
    r[[u[i]]] <- i
  }
  y <- rep(0,length(x))
  for (i in 1:length(x)) {
    y[i] <- r[[x[i]]]
  }
  y
}

fairest <- function(clss) {
  fairness <- rep(0,length(clss))
  #cat("Fairness:\n")
  for (i in 1:length(clss)) {
    fq <- as.numeric(table(clss[[i]]$cluster))
    #cat("  ",fq,"->")
    fq <- fq/sum(fq)
    fairness[i] <- -sum(fq*log(fq))
    #cat("fairness ",fairness[i],"\n")
  }
  which.max(fairness)
}

kmeans2 <- function(x,centers=NULL,iter.max=10,ineq=0) {
  if (ineq < 0 | 1 < ineq) {
    stop("ineq should be between 0 and 1")
  }
  if (nrow(x) < 2) {
    stop("x should be more than one point")
  }
  if (nrow(x) == 2) {
    return(list(
      cluster=c(1,2),
      center=x
    ))
  }
  if (is.null(centers)) {
    g <- colMeans(x)
    delta <- rnorm(ncol(x))
    center <- rbind(g+delta,g-delta)
  }
  cluster <- rep(0,nrow(x))
  d <- matrix(0,nrow=nrow(x),ncol=2)
  half <- nrow(x)/2
  for (k in 1:iter.max) {
    for (i in 1:2) {
      d[,i] <- apply(x,1,function(z){sum((z-center[i,])^2)})
    }
    nd <- apply(d,1,function(x){x[1]/sum(x)})
    ind <- order(nd,decreasing=TRUE)
    n1 <- sum(nd>=0.5)
    h <- floor(half*(1-ineq)+n1*ineq)
    ncenter <- matrix(0,nrow=2,ncol=ncol(x))
    ncenter[1,] <- colMeans(x[ind[1:h],])
    ncenter[2,] <- colMeans(x[ind[(h+1):length(ind)],])
    if (mean((center-ncenter)^2) == 0) break
    center <- ncenter
  }
  cluster[ind[1:h]] <- 1
  cluster[ind[(h+1):length(ind)]] <- 2
  list(cluster=cluster,center=center)
}

#' The successive binary partition k-means algorithm.
#' @param x A matrix or data frame. Each row corresponds to the each data.
#' @param n Number of clusters.
#' @param iter.max Number of maximum iteration in each split.
#' @param iter.final Number of iteration at the final centroid refinement.
#' @param algorithm Gigen to kmeans().
#' @return List of the clustering result.
#'   cluster: a vector indicating the cluster numbers of the samples.
#'   centers: a matrix of the centroids
#'   size: a vector of the cluster sizes
#' @export
sbp_kmeans <- function(x,n,iter.final=0,
                     iter.max=20,algorithm="Hartigan-Wong") {
  x.max <- x
  org_ind.max <- 1:nrow(x)
  current_cluster_num <- 2
  current_class <- rep(1,nrow(x))
  while (current_cluster_num <= n) {
    # clss <- list()
    # for (i in 1:3) {
    #   clss[[i]] <- stats::kmeans(x.max,2,
    #                       iter.max=iter.max,
    #                       algorithm=algorithm)
    # }
    # cls <- clss[[fairest(clss)]]
    cls <- kmeans2(x.max,iter.max=iter.max)
    current_class[org_ind.max] <- cls$cluster+current_cluster_num-1
    current_class <- renumber(current_class)
    tbl <- base::table(current_class)
    max_class <- as.integer(names(tbl)[which.max(tbl)])
    org_ind.max <- which(current_class==max_class)
    x.max <- x[org_ind.max,]
    current_cluster_num <- current_cluster_num+1
  }
  centers <- matrix(0,nrow=n,ncol=ncol(x))
  for (i in 1:n) {
    ind <- which(current_class==i)
    if (length(ind) == 1) {
      centers[i,] <- x[ind,]
    } else {
      centers[i,] <- colMeans(x[ind,])
    }
  }
  if (iter.final > 0) {
    cls <- kmeans(x,centers,iter.max=iter.final,algorithm=algorithm)
    return(cls)
  }
  list(
    cluster=current_class,
    centers=centers,
    size=as.numeric(tbl)
  )
}
