
# 加载必要包
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("Matrix"))

## 命令行参数定义
spec <- matrix(c(
  "grm",  "g", 1, "character", "[Required] Input MTG2 .grm file (3 columns: ID1 ID2 Value)",
  "out",  "o", 1, "character", "[Required] Output prefix",
  "tol",  "t", 1, "double",    "[Optional] Tolerance for positive definite (default 1e-6)",
  "help", "h", 0, "logical",   "Show help"
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

if(!is.null(opt$help) || is.null(opt$grm) || is.null(opt$out)){
  cat(getopt(spec, usage = TRUE))
  quit(status = -1)
}

if (is.null(opt$tol)) opt$tol = 1e-6

# 1. 读取数据
cat(paste0("Loading GRM: ", opt$grm, "\n"))
grm_data <- fread(opt$grm, header = FALSE)
setnames(grm_data, c("V1", "V2", "V3"))

# 2. 重建矩阵
# 获取所有不重复的 ID 并排序，确保矩阵行列对应
ids <- unique(c(grm_data$V1, grm_data$V2))
n <- length(ids)
id_map <- setNames(seq_along(ids), ids)

cat(paste0("Matrix size: ", n, " x ", n, "\n"))

# 使用稀疏矩阵结构填充，然后转为标准矩阵格式
G <- matrix(0, nrow = n, ncol = n, dimnames = list(ids, ids))

# 填充对称矩阵（MTG2 输出通常是下三角）
# 利用 match 快速定位索引
idx1 <- id_map[as.character(grm_data$V1)]
idx2 <- id_map[as.character(grm_data$V2)]
vals <- as.numeric(grm_data$V3)

for(i in 1:length(vals)) {
  G[idx1[i], idx2[i]] <- vals[i]
  G[idx2[i], idx1[i]] <- vals[i]
}

# 3. 保证矩阵正定 (参考用户提供逻辑)
make_positive_definite <- function(M, tol=1e-6) {
  cat("Checking positive definiteness...\n")
  # 计算特征值
  eig <- eigen(M, symmetric = TRUE)
  # 定义最小阈值：最大特征值 * tol
  min_eig <- max(eig$values) * tol
  
  if(min(eig$values) < min_eig) {
    cat("Warning: Matrix is not positive definite. Adjusting eigenvalues...\n")
    vals <- eig$values
    vals[vals < min_eig] <- min_eig
    # 重构矩阵: Q * Lambda * Q^T
    M_new <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
    dimnames(M_new) <- dimnames(M)
    return(M_new)
  } else {
    cat("Matrix is positive definite.\n")
    return(M)
  }
}

G_pd <- make_positive_definite(G, tol = opt$tol)

# 4. 求逆
cat("Inverting matrix...\n")
# 使用 solve() 函数求逆。如果仍报错，可以考虑加入 tryCatch 并在对角线加极小值
G_inv <- tryCatch({
  solve(G_pd)
}, error = function(e) {
  cat("Solve failed, adding extra epsilon to diagonal...\n")
  diag(G_pd) <- diag(G_pd) + 1e-4
  solve(G_pd)
})

# 5. 转换为三列下三角格式输出
cat("Converting to 3-column lower triangular format...\n")
# 将上三角设为 NA
G_inv[upper.tri(G_inv)] <- NA

# 转换为长表格格式
# 使用 data.table 的效率优势
G_inv_dt <- as.data.table(as.table(G_inv))
G_inv_dt <- G_inv_dt[!is.na(N)] # 剔除上三角的 NA

# 6. 写出文件
out_file <- paste0(opt$out, ".grm.inv")
fwrite(G_inv_dt, out_file, sep = " ", col.names = FALSE)

cat(paste0("Success: Inverse GRM saved to ", out_file, "\n"))
