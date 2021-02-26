## load data
dta <- readRDS("data set/dta.rds")
p <- ncol(dta)
n <- nrow(dta)

index0 <- apply(dta, 1, function(x){which(as.numeric(x) == 0)})
index1 <- lapply(index0, function(x){setdiff(seq_len(p), x)})
print("Begin to calculate the correlation.")

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
kendall_fun <- function(dta) {
  n <- nrow(dta)
  p <- ncol(dta)
  kendall_mat <- matrix(1, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in i:n) {
      index11 <- intersect(index1[[i]], index1[[j]])
      p11 <- length(intersect(index1[[i]],  index1[[j]])) / p
      p00 <- length(intersect(index0[[i]], index0[[j]])) / p
      p10 <- length(intersect(index1[[i]], index0[[j]])) / p
      p01 <- length(intersect(index0[[i]], index1[[j]])) / p
      if (length(index11) == 0){
        tau <- 0
      } else{
        x <- as.numeric(dta[i, index11[1]])
        y <- as.numeric(dta[j, index11[1]])
        if (identical(x, y)){
          tau <- 1
        } else{
          tau <- suppressWarnings(cor(x, y, method = "kendall"))
          if (is.na(tau)) {
            tau <- 0
          }
        }
      }
      kendall_mat[i, j] <- p11 ^ 2 * tau + 2 * (p00 * p11 - p01 * p10)
      kendall_mat[j, i] <- kendall_mat[i, j] 
    }
  }
  return(kendall_mat)
}

## spearman correlation
spearman_fun <- function(dta) {
  n <- nrow(dta)
  p <- ncol(dta)
  spearman_mat <- matrix(1, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in i:n){
      index11 <- intersect(index1[[i]], index1[[j]])
      p11 <- length(intersect(index1[[i]],  index1[[j]])) / p
      p00 <- length(intersect(index0[[i]], index0[[j]])) / p
      p10 <- length(intersect(index1[[i]], index0[[j]])) / p
      p01 <- length(intersect(index0[[i]], index1[[j]])) / p
      px1 <- length(index1[i]) / n
      py1 <- length(index1[j]) / n
      if (length(index11) == 0) {
        rho <- 0
      } else{
        x <- as.numeric(dta[i, index11[1]])
        y <- as.numeric(dta[j, index11[1]])
        if (identical(x, y)) {
          rho <- 1
        } else{
          rho <- suppressWarnings(cor(x, y, method = "spearman"))
          if (is.na(rho)) {
            rho <- 0
          }
        }
      }
      spearman_mat[i, j] <- p11 * px1 * py1 * rho + 3 * (p00 * p11 - p01 * p10)
      spearman_mat[j, i] <- spearman_mat[i, j] 
    }
  }
  return(spearman_mat)
}

## Three correlation method
cor_pearson <- cor(t(dta), method = "pearson")
print("Finish pearson method.")
# cor_kendall <- cor(t(dta), method = "kendall")
# print("Finish kendall method.")
cor_spearman <- cor(t(dta), method = "spearman")
print("Finish three basic methods.")
res <- list()
res$cor_pearson <- cor_pearson
# res$cor_kendall <- cor_kendall
res$cor_spearman <- cor_spearman

# adjustment
cor_spearman_adj <- spearman_fun(dta = dta)
print("Finish adjustment spearman.")

cor_kendall_adj <- kendall_fun(dta = dta)
print("Finish adjustment kendall.")

res$cor_kendall_adj <- cor_kendall_adj
res$cor_spearman_adj <- cor_spearman_adj
saveRDS(res, "results/cor_res.rds")

## normalization data
dta_norm <- new_Normalization(dta)
cor_pearson_norm <- cor(t(dta_norm), method = "pearson")
# cor_kendall_norm <- cor(t(dta_norm), method = "kendall")
cor_spearman_norm <- cor(t(dta_norm), method = "spearman")
print("Finish three basic methods for normalized data.")

## return results
res$cor_pearson_norm <- cor_pearson_norm
# res$cor_kendall_norm <- cor_kendall_norm
res$cor_spearman_norm <- cor_spearman_norm
saveRDS(res, "results/cor_res.rds")

## Calculate the difference
diff_mat <- matrix(0, nrow = 6, ncol = 6)
for (i in seq_len(5)) {
  for (j in (i + 1):6) {
    mat1 <- res[[i]]
    mat1[which(mat1 > 0)] <- 1
    mat1[which(mat1 < 0)] <- -1
    mat2 <- res[[j]]
    mat2[which(mat2 > 0)] <- 1
    mat2[which(mat2 < 0)] <- -1
    diff_mat[i, j] <- sum((mat1 - mat2)^2)
    diff_mat[j, i] <- diff_mat[i, j]
  }
}
diff_mat <- diff_mat / (n^2)
rownames(diff_mat) <- colnames(diff_mat) <- c("pearson", "spearman", "kendall_adj", "spearman_adj",
                                              "pearson_norm", "spearman_norm")
res$diff_mat <- diff_mat
print("Finish the first part of calculating difference.")

diff_mat_round <- matrix(0, nrow = 6, ncol = 6)
for (i in seq_len(5)) {
  for (j in (i + 1):6) {
    mat1 <- round(res[[i]], 2)
    mat1[which(mat1 > 0)] <- 1
    mat1[which(mat1 < 0)] <- -1
    mat2 <- round(res[[j]], 2)
    mat2[which(mat2 > 0)] <- 1
    mat2[which(mat2 < 0)] <- -1
    diff_mat_round[i, j] <- sum((mat1 - mat2)^2)
    diff_mat_round[j, i] <- diff_mat_round[i, j]
  }
}
diff_mat_round <- diff_mat_round / (n^2)
rownames(diff_mat_round) <- colnames(diff_mat_round) <- c("pearson", "spearman", "kendall_adj", "spearman_adj",
                                                          "pearson_norm", "spearman_norm")
res$diff_mat_round <- diff_mat_round
saveRDS(res, "results/cor_res.rds")
