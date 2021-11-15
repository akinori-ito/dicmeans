#
# Dichotomic k-means clustering algorithm
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

#' The dichotomic k-means algorithm.
#' @param x A matrix or data frame. Each row corresponds to the each data.
#' @param n Number of clusters.
#' @param iter.max Given to kmeans().
#' @param algorithm Gigen to kmeans().
#' @return The clustering result.
#' @export
dicmeans <- function(x,n,iter.max=20,algorithm="Hartigan-Wong") {
  x.max <- x
  org_ind.max <- 1:nrow(x)
  current_cluster_num <- 2
  current_class <- rep(1,nrow(x))
  while (current_cluster_num <= n) {
    clss <- list()
    for (i in 1:3) {
      clss[[i]] <- stats::kmeans(x.max,2,
                           iter.max=iter.max,
                           algorithm=algorithm)
      
    }
    cls <- clss[[fairest(clss)]]
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
  list(
    cluster=current_class,
    centers=centers,
    size=as.numeric(tbl)
  )
}
