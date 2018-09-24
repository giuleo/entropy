
# =============================================================================
# Function to call on sparseMatrix in order to compute categorical entropy.
# It also return the data frame with the coordinates of the blocks
# =============================================================================

catEntropy <- function(X = sparseMatrix()) {
  if (class(X) != "ngCMatrix") {
    stop("Input must be a Sparse Matrix")
  }
  blocks <- blocks(X)
  
  bkarea <- blocks[blocks[,7] > 1, 7] 
  if (length(bkarea) > 0) {
    tabarea <- as.data.frame(table(bkarea))
    parea <- tabarea$Freq / sum(tabarea$Freq)
    catentropy = -sum(parea * log(parea))
  }
  else {
    catentropy <- NA
  }
  return(list(Entropy = catentropy, Blocks = blocks))
}

# =============================================================================
# Second attempt: Function to call on sparseMatrix in order to compute 
# categorical entropy. 
# =============================================================================

catEnt <- function(X = sparseMatrix()){
  
  if (class(X) != "ngCMatrix") {stop("Input must be a Sparse Matrix")}
  points <- findBlocks(X)
  
  nbk <- length(unique(points$x))
  areas <- sapply(1:nbk, function(i) length(points[points$x == i, 3]))
  bkarea <- areas[areas > 1]
  if (length(bkarea) > 0) {
    tabarea <- as.data.frame(table(bkarea))
    parea <- tabarea$Freq / sum(tabarea$Freq)
    catentropy = -sum(parea * log(parea))
  }
  else {catentropy <- NA}
  return(catentropy)
}

# ============================================================================
# It takes in the block named sparse matrix and returns a data frame 
# with the coordinates of the blocks/lines/points present in it.
# In the data frame are also three additional column for width, height 
# and area of the blocks.
# ============================================================================

blocks <- function(X) {
  if (class(X) != "ngCMatrix") {stop("Input must be a Sparse Matrix")}
  ind <- findBlocks(X)
  nbk <- length(unique(ind$x))
#  blk <- as.data.frame(matrix(nrow = nbk, ncol = 7))
#  colnames(blk) <- c("x1", "y1", "x2", "y2", "W", "H", "A")
#  for (loop in 1:nbk) {
    #a <- ind[ind$x == loop,]
    #x1 <- min(a[,1])
    #y1 <- min(a[,2])
    #x2 <- max(a[,1])
    #y2 <- max(a[,2])
    #W <- (x2 + 1) - x1
    #H <- (y2 + 1) - y1
    #A <- W*H
    #blk[loop,] <- c(x1, x2, y1, y2, W, H, A)
#  }
  x1 <- sapply(1:nbk, function(i) min(ind[ind$x == i, 1]))
  y1 <- sapply(1:nbk, function(i) min(ind[ind$x == i, 2]))
  x2 <- sapply(1:nbk, function(i) max(ind[ind$x == i, 1]))
  y2 <- sapply(1:nbk, function(i) max(ind[ind$x == i, 2]))
  blk <- data.frame("x1" = x1, "y1" = y1, "x2" = x2, "y2" = y2) 
  blk$W <- (blk[,3]+1) - blk[,1]
  blk$H <- (blk[,4]+1) - blk[,2]
  blk$A <- blk[,5]*blk[,6]
  return(blk)
}

# ===============================================================
# Matrix transformation by @alexis_laz (Stackexchange)
# Assign to non-zero elements of a sparse Matrix
# (i.e. categorical recurrence plot) a numerical code identifying 
# its block membership
# ===============================================================

findBlocks = function(x)
{
  sm = as.matrix(summary(x))
  
  gr = integer(nrow(sm))
  ngr = 0L 
  gr[1] = ngr 
  
  lastSeenRow = integer(nrow(x))
  lastSeenCol = integer(ncol(x))
  
  for(k in 1:nrow(sm)) {
    kr = sm[k, 1] 
    kc = sm[k, 2]
    i = lastSeenRow[kr]
    j = lastSeenCol[kc]
    if (i && (abs(kc - sm[i, 2]) == 1)) gr[k] = gr[i]
    else if (j && (abs(kr - sm[j, 1]) == 1)) gr[k] = gr[j]  
    else { 
      ngr = ngr + 1L; gr[k] = ngr 
    }
    lastSeenRow[kr] = k
    lastSeenCol[kc] = k        
  }
  return(data.frame(i = sm[, "i"], j = sm[, "j"], x = gr))
}

# ===============================================================
# Matrix transformation with distance matrix and clustering 
# ===============================================================

indx <- function(X){
  
  require(Rclusterpp)
  
  sm <- summary(X)
  sm$i <- as.numeric(sm$i)
  sm$j <- as.numeric(sm$j)
  # d <- as.dist(sapply(1:nrow(sm), function(i) rowSums(abs(sm[i, col(sm)] - sm))))
  # d <- dist(sm, method = "manhattan")
  b <- data.frame(i = sm[, "i"], j = sm[, "j"], 
                  x = cutree(Rclusterpp.hclust(sm, method = "single", distance = "manhattan"), 
                             h = 1))
  # b <- sparseMatrix(i = sm[, "i"], j = sm[, "j"], x = cutree(hclust(d, "single"), h = 1))
  return(b)
}