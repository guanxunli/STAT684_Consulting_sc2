library(Matrix)

## load data
dta <- as.matrix(readRDS("data set/dta.rds"))
p <- ncol(dta)
n <- nrow(dta)
gene_name <- rownames(dta)
res <- list()

index0 <- apply(dta, 1, function(x){which(as.numeric(x) == 0)})
index1 <- lapply(index0, function(x){setdiff(seq_len(p), x)})

dta_01 <- matrix(0, nrow = n, ncol = p)
dta_01[which(abs(dta) > 0)] <- 1
rownames(dta_01) <- rownames(dta)

## spearman correlation
spearman_mat <- matrix(NA, nrow = n, ncol = n)
cor_01 <- matrix(NA, nrow = n, ncol = n)
cor_common <- matrix(NA, nrow = n, ncol = n)

## begin calculation
for (i in seq_len(n)) {
  for (j in i:n) {
    index11 <- intersect(index1[[i]], index1[[j]])
    p11 <- length(intersect(index1[[i]],  index1[[j]])) / p
    p00 <- length(intersect(index0[[i]], index0[[j]])) / p
    p10 <- length(intersect(index1[[i]], index0[[j]])) / p
    p01 <- length(intersect(index0[[i]], index1[[j]])) / p
    px1 <- length(index1[[i]]) / p
    py1 <- length(index1[[j]]) / p
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
    tmp <- p00 * p11 - p01 * p10
    cor_common[i, j] <- cor_common[j, i] <- rho
    cor_01[i, j] <- cor_01[j, i] <- tmp
    spearman_mat[i, j] <- spearman_mat[j, i] <- p11 * px1 * py1 * rho + 3 * tmp
  }
}
cor_01[which(cor_01 > 0)] <- 1
cor_01[which(cor_01 < 0)] <- -1

spearman_01 <- matrix(NA, nrow = n, ncol = n)
spearman_01[which(spearman_mat > 0)] <- 1
spearman_01[which(spearman_mat < 0)] <- -1

rownames(cor_01) <- colnames(cor_01) <- gene_name
rownames(cor_common) <- colnames(cor_common) <- gene_name
rownames(spearman_mat) <- colnames(spearman_mat) <- gene_name
rownames(spearman_01) <- colnames(spearman_mat) <- gene_name

res$cor_01 <- cor_01
res$cor_common <- cor_common
res$spearman_mat <- spearman_mat

## calculate the difference pair with adjusted spearman
index <- which(abs(cor_01 - spearman_01) > 0)
index0 <- which(cor_01 == 0)
index <- setdiff(index, index0)

pair_mat <- matrix(NA, nrow = 2, ncol = length(index))
pair_mat[1, ] <- index %/% n + 1
pair_mat[2, ] <- index %% n
pair_mat[2, ][which(pair_mat[2, ] == 0)] <- n

gene_pair <- matrix(NA, nrow = length(index) / 2, ncol = 2)
iter <- 1
for (i in seq_len(length(index))) {
  if (pair_mat[1, i] < pair_mat[2, i]) {
    gene_pair[iter, ] <- gene_name[c(pair_mat[1, i],  pair_mat[2, i])]
    iter <- iter + 1
  }
}
res$gene_pair <- gene_pair

## For each pair do fisher exact test
pair_p <- rep(NA, nrow(gene_pair))
for (i in seq_len(length(pair_p))) {
  tmp <- table(dta_01[gene_pair[i, 1], ], dta_01[gene_pair[i, 2], ])
  if (identical(as.numeric(dim(tmp)), c(2,2))) {
    pair_p[i] <- fisher.test(dta_01[res$gene_pair[i, 1], ], dta_01[res$gene_pair[i, 2], ])$p.value
  } else{
    pair_p[i] <- 1
  }
}
res$pair_p <- pair_p

## calculate the difference pair with common part
cor_common_01 <- matrix(0, nrow = n, ncol = n)
cor_common_01[which(cor_common > 0)] <- 1
cor_common_01[which(cor_common < 0)] <- -1

index <- which(abs(cor_01 - cor_common_01) > 0)
index <- setdiff(index, union(which(cor_01 == 0), which(cor_common == 0)))

pair_mat <- matrix(NA, nrow = 2, ncol = length(index))
pair_mat[1, ] <- index %/% n + 1
pair_mat[2, ] <- index %% n
pair_mat[2, ][which(pair_mat[2, ] == 0)] <- n

gene_pair <- matrix(NA, nrow = length(index) / 2, ncol = 2)
iter <- 1
for (i in seq_len(length(index))) {
  if (pair_mat[1, i] < pair_mat[2, i]) {
    gene_pair[iter, ] <- gene_name[c(pair_mat[1, i],  pair_mat[2, i])]
    iter <- iter + 1
  }
}
res$gene_pair_common <- gene_pair

## For each pair do fisher exact test
pair_p <- rep(NA, nrow(gene_pair))
dta_01 <- matrix(0, nrow = n, ncol = p)
dta_01[which(dta > 0)] <- 1
rownames(dta_01) <- rownames(dta)
for (i in seq_len(length(pair_p))) {
  tmp <- table(dta_01[gene_pair[i, 1], ], dta_01[gene_pair[i, 2], ])
  if (identical(as.numeric(dim(tmp)), c(2,2))) {
    pair_p[i] <- fisher.test(dta_01[res$gene_pair[i, 1], ], dta_01[res$gene_pair[i, 2], ])$p.value
  } else{
    pair_p[i] <- 1
  }
}

res$pair_p_common <- pair_p
saveRDS(res, "results/res_correlation.rds")
