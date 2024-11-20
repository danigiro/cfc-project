# Code extract from https://github.com/nielsaka/modelconf

#' @param loss A matrix of size (n x m). The columns contain the estimated losses
#'   for each of the \code{m} models.
#' @param test A character string. It specifies the test statistic to be used.
#'   Available tests are "t.max".
#' @param B A scalar, the number of bootstrap samples.
#' @param l A scalar, the block length used in the moving-block bootstrap.
estMCS <- function(loss, test="t.max", B=1000, l=2){
    if (!any(test %in% c("t.range", "t.max", "t.min"))) stop("MCS: Unknown test statistic.")
    if (NCOL(loss) == 1) stop("MCS: Need more data. Only one model entered.")
    n <- nrow(loss)
    m <- ncol(loss)
    model.names <- dimnames(loss)[[2]]
    if(is.null(model.names)) model.names <- 1:m
    blocks <- makeBlocks(n, l)
    boot.index <- makeIndex(B, blocks)
    mcs <- matrix(NA, nrow=m, ncol=3)
    colnames(mcs) <- c("model", "p-val", "MCS p-val")
    models <- 1:m
    stats <- makeStats(loss, boot.index)
    for(i in 1:(m-1)){
      rejection <- do.call(test, list(stats))  
      mcs[i, ] <- c(models[rejection$candidate], 
                    rejection$p.value,         
                    max(mcs[, "p-val"], rejection$p.value, na.rm=T))
      models <- models[-rejection$candidate]
      stats <- list(data.mean = stats$data.mean[-rejection$candidate],
                    boot.data.mean = stats$boot.data.mean[-rejection$candidate,])
    }
    mcs[m, ] <- c(models, 1, 1)
    mcs <- mcs[order(mcs[,"model"]), ]
    dimnames(mcs)[[1]] <- model.names
    return(mcs)
    
}

makeIndex <- function(B, blocks){
    n <- ncol(blocks)
    l <- nrow(blocks)
    z <- ceiling(n/l)
    start.points <- sample.int(n, z * B, replace = TRUE) 
    index <- blocks[, start.points]
    keep <- c(rep(TRUE, n), rep(FALSE, z*l - n))
    boot.index <- matrix(as.vector(index)[keep], nrow = n, ncol = B)
    return(boot.index)
  }

makeStats <- function(data, boot.index) {
    data.mean <- colMeans(data)
    n <- nrow(data)
    B <- ncol(boot.index)
    weights <- vapply(1L:B,
                      function(j) tabulate(boot.index[, j], nbins = n),
                      integer(n)) * (1 / n)
    if (all(apply(weights, 1, function(x) length(unique(x)) == 1))){
      warning("makeStats: Insufficient observations or block length too large.")
    } 
    boot.data.mean <- t(data) %*% weights
    return(list(data.mean = data.mean,
                boot.data.mean = boot.data.mean))
  }

makeBlocks <- function(n, l){
    blocks <- matrix(NA, nrow=l, ncol=n)
    blocks[1,] <- 1:n
    if (l > 1){
      for (i in 2:l) {
        blocks[i,] <- c(i:n, 1:(i - 1))
      }
    }
    return(blocks)
  }

t.max <- function(stats){ 
    data.mean      <- stats$data.mean
    boot.data.mean <- stats$boot.data.mean
    data.mean.all <- mean(data.mean)
    boot.data.mean <- boot.data.mean - data.mean
    boot.data.mean.all <- colMeans(boot.data.mean)
    boot.data.mean.shift <- t(t(boot.data.mean) - boot.data.mean.all)
    #data.sd <- boot.data.mean.shift*boot.data.mean.shift
    #data.sd <- sqrt(rowMeans(data.sd))
    # https://github.com/nielsaka/modelconf/issues/2
    data.sd <- apply(boot.data.mean.shift,1,sd)
    t.stat <- (data.mean - data.mean.all) / data.sd
    maxT <- max(t.stat)
    boot.t.stat <- boot.data.mean.shift / data.sd
    boot.maxT <-  apply(boot.t.stat, 2, max)
    p.value <- mean((boot.maxT - maxT) > 0)
    if(is.na(p.value)) stop(paste("p.value=", p.value))
    return(list(candidate=which(t.stat==maxT), p.value=p.value))
  }
