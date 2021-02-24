## load data
dta <- readRDS("data set/dta.rds")
dim(dta)

## Define functions
## Normalization
new_Normalization <- function(X) {
  X_sum <- sum(X)
  X_rowsum <- rowSums(X)
  X_colSum <- colSums(X)
  u_hat <- X_rowsum %*% t(X_colSum) / X_sum
  Z <- (X - u_hat) / sqrt(u_hat + u_hat^2 / 100)
  return(Z)
}
## kendall correlation
kendall_fun <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  n <- length(x)
  index_x0 <- which(x == 0)
  index_y0 <- which(y == 0)
  index_x1 <- setdiff(seq_len(n), index_x0)
  index_y1 <- setdiff(seq_len(n), index_y0)
  index11 <- intersect(index_x1, index_y1)
  p11 <- length(intersect(index_x1, index_y1)) / n
  p00 <- length(intersect(index_x0, index_y0)) / n
  p10 <- length(intersect(index_x1, index_y0)) / n
  p01 <- length(intersect(index_x0, index_y1)) / n
  if (identical(x[index11], y[index11])){
    tau <- 1
  } else{
    tau <- cor(x[index11], y[index11], method = "kendall")
  }
  return(p11 ^ 2 * tau + 2 * (p00 * p11 - p01 * p10))
}

spearman_fun <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  n <- length(x)
  index_x0 <- which(x == 0)
  index_y0 <- which(y == 0)
  index_x1 <- setdiff(seq_len(n), index_x0)
  index_y1 <- setdiff(seq_len(n), index_y0)
  index11 <- intersect(index_x1, index_y1)
  p11 <- length(intersect(index_x1, index_y1)) / n
  p00 <- length(intersect(index_x0, index_y0)) / n
  p10 <- length(intersect(index_x1, index_y0)) / n
  p01 <- length(intersect(index_x0, index_y1)) / n
  px1 <- length(index_x1) / n
  py1 <- length(index_y1) / n
  if (identical(x[index11], y[index11])){
    rho <- 1
  } else{
    rho <- cor(x[index11], y[index11], method = "spearman")
  }
  return(p11 * px1 * py1 * rho + 3 * (p00 * p11 - p01 * p10))
}

## Three correlation method
cor_pearson <- cor(dta, method = "pearson")
cor_kendall <- cor(dta, method = "kendall")
cor_spearman <- cor(dta, method = "spearman")
print("Finish three basic methods.")

# adjustment
cor_kendall_adj <- matrix(1, nrow = n, ncol = n)
for (i in seq_len(n)){
  for (j in i:n){
    cor_kendall_adj[i, j] <- kendall_fun(dta[i, ], dta[j, ])
    cor_kendall_adj[j, i] <- cor_kendall_adj[i, j]
  }
}
print("Finish adjustment kendall.")

cor_spearman_adj <- matrix(1, nrow = n, ncol = n)
for (i in seq_len(n)){
  for (j in i:n){
    cor_spearman_adj[i, j] <- spearman_fun(dta[i, ], dta[j, ])
    cor_spearman_adj[j, i] <- cor_spearman_adj[i, j]
  }
}
print("Finish adjustment spearman.")

## normalization data
dta_norm <- new_Normalization(dta)
cor_pearson_norm <- cor(dta_norm, method = "pearson")
cor_kendall_norm <- cor(dta_norm, method = "kendall")
cor_spearman_norm <- cor(dta_norm, method = "spearman")
print("Finish three basic methods for normalized data.")

# adjustment
cor_kendall_norm_adj <- matrix(1, nrow = n, ncol = n)
for (i in seq_len(n)){
  for (j in i:n){
    cor_kendall_norm_adj[i, j] <- kendall_fun(dta_norm[i, ], dta_norm[j, ])
    cor_kendall_norm_adj[j, i] <- cor_kendall_norm_adj[i, j]
  }
}
print("Finish adjustment kendall for normalized data.")

cor_spearman_norm_adj <- matrix(1, nrow = n, ncol = n)
for (i in seq_len(n)){
  for (j in i:n){
    cor_spearman_norm_adj[i, j] <- spearman_fun(dta_norm[i, ], dta_norm[j, ])
    cor_spearman_norm_adj[j, i] <- cor_spearman_norm_adj[i, j]
  }
}
print("Finish adjustment spearman for normalized data.")

## return results
res <- list()
res$cor_pearson <- cor_pearson
res$cor_kendall <- cor_kendall
res$cor_spearman <- cor_spearman
res$cor_kendall_adj <- cor_kendall_adj
res$cor_spearman_adj <- cor_spearman_adj

res$cor_pearson_norm <- cor_pearson_norm
res$cor_kendall_norm <- cor_kendall_norm
res$cor_spearman_norm <- cor_spearman_norm
res$cor_kendall_norm_adj <- cor_kendall_norm_adj
res$cor_spearman_norm_adj <- cor_spearman_norm_adj

## Calculate the difference
diff_mat <- matrix(0, nrow = 10, ncol = 10)
for (i in seq_len(9)) {
  for (j in (i + 1):10) {
    mat1 <- res[[i]]
    mat1[which(mat1 > 0)] <- 1
    mat1[which(mat1 < 0)] <- -1
    mat2 <- res[[j]]
    mat2[which(mat2 > 0)] <- 1
    mat2[which(mat2 < 0)] <- -1
    diff_mat <- sum((mat1 - mat2)^2)
  }
}
res$diff_mat <- diff_mat
print("Finish the first part of calculating difference.")

diff_mat_round <- matrix(0, nrow = 10, ncol = 10)
for (i in seq_len(9)) {
  for (j in (i + 1):10) {
    mat1 <- round(res[[i]], 2)
    mat1[which(mat1 > 0)] <- 1
    mat1[which(mat1 < 0)] <- -1
    mat2 <- round(res[[j]], 2)
    mat2[which(mat2 > 0)] <- 1
    mat2[which(mat2 < 0)] <- -1
    diff_mat_round <- sum((mat1 - mat2)^2)
  }
}
res$diff_mat_round <- diff_mat_round
saveRDS(res, "cor_res.rds")