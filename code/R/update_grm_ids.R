
# 加载必要包
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))

## 命令行参数定义
spec <- matrix(c(
  "fam",  "f", 1, "character", "[Required] Input .fam file to get Individual IDs",
  "grm",  "g", 1, "character", "[Required] Input .grm file (3 columns: Index1 Index2 Value)",
  "out",  "o", 1, "character", "[Required] Output file path/name",
  "help", "h", 0, "logical",   "Show help"
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

# 帮助信息校验
if(!is.null(opt$help) || is.null(opt$fam) || is.null(opt$grm) || is.null(opt$out)){
  cat(getopt(spec, usage = TRUE))
  quit(status = -1)
}

# --- 1. 读取 FAM 文件 ---
# 只读取第二列 (V2)，通常是 IID
if (!file.exists(opt$fam)) stop("FAM file not found.")
cat("Reading FAM file and extracting IDs...\n")
fam <- fread(opt$fam, header = FALSE, select = 2)
id_map <- fam[[1]] # 转为向量

# --- 2. 读取 GRM 文件 ---
if (!file.exists(opt$grm)) stop("GRM file not found.")
cat("Reading GRM file...\n")
# 使用 fread 自动识别分隔符，通常 GRM 是空格或制表符
grm <- fread(opt$grm, header = FALSE)

# 校验索引范围是否合法
max_idx <- max(max(grm[[1]]), max(grm[[2]]))
if (max_idx > length(id_map)) {
  stop(paste0("Error: GRM contains index (", max_idx, ") exceeding FAM row count (", length(id_map), ")."))
}

# --- 3. 替换 ID ---
cat("Mapping indices to IDs...\n")
# 将第1、2列从数字索引替换为 IID
grm[[1]] <- id_map[grm[[1]]]
grm[[2]] <- id_map[grm[[2]]]

# --- 4. 保存结果 ---
cat("Saving updated GRM to:", opt$out, "\n")
fwrite(grm, file = opt$out, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

cat("Successfully finished.\n")
