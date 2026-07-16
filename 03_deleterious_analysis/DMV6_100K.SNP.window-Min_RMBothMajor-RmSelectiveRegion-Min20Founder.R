All_Dele_Genome2=read.table("DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_RmSelectiveRegion.txt",header=T)
All_Dele_Genome_CNDLessThanA626=All_Dele_Genome2[(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$A626),]
write.table(All_Dele_Genome_CNDLessThanA626,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_CNDlessThanA626_RmSelectiveRegion.bed",quote=F,sep="\t",row.name=F)

index_CND=match(S.can_ID,names(All_Dele_Genome2))
 All_Dele_Genome_CNDLessThanA626_CND= All_Dele_Genome_CNDLessThanA626[,c(1:12,index_CND)]
 write.table(All_Dele_Genome_CNDLessThanA626_CND,"DMRef_100KSNP.window_AllDeleBurden_AtLeastOne_CNDlessThanA626_RmSelectiveRegion_CND.deleInfo.bed",quote=F,sep="\t",row.name=F)
 CND_min=apply(All_Dele_Genome_CNDLessThanA626[,index_CND],1,min)
Reduced_A626=(All_Dele_Genome_CNDLessThanA626$A626 -All_Dele_Genome_CNDLessThanA626$Hotspot_AllDele_CNDMin )

Reduced_A626_100KSNP=(All_Dele_Genome_CNDLessThanA626$A626[All_Dele_Genome_CNDLessThanA626$SNP_num==100000] -All_Dele_Genome_CNDLessThanA626$Hotspot_AllDele_CNDMin[All_Dele_Genome_CNDLessThanA626$SNP_num==100000]  )


 CND_min_col=apply(All_Dele_Genome_CNDLessThanA626[,index_CND],1,which.min)
CND_min_col2=apply(All_Dele_Genome_CNDLessThanA626[All_Dele_Genome_CNDLessThanA626_CND$SNP_num!=100000,index_CND],1,which.min)


###################################
mat=All_Dele_Genome_CNDLessThanA626_CND[,c(5,13:ncol(All_Dele_Genome_CNDLessThanA626_CND))]
prefix="A626_CND_"

#################################
All_Dele_Genome_CNDLessThanE463=All_Dele_Genome2[(All_Dele_Genome2$Hotspot_AllDele_CNDMin < All_Dele_Genome2$E463),]
index_CND=match(S.can_ID,names(All_Dele_Genome2))
All_Dele_Genome_CNDLessThanE463_CND= All_Dele_Genome_CNDLessThanE463[,c(1:12,index_CND)]
Reduced_E463=(All_Dele_Genome_CNDLessThanE463$E463 -All_Dele_Genome_CNDLessThanE463$Hotspot_AllDele_CNDMin )
sum(Reduced_E463) ### 52875.09
mat=All_Dele_Genome_CNDLessThanE463_CND[,c(6,13:ncol(All_Dele_Genome_CNDLessThanE463_CND))]
prefix="E463_CND_"
################################################################################
# 假设您的矩阵名为 mat，800行 x 200列
# 若尚未生成，可模拟随机数据测试：
# set.seed(123)
# mat <- matrix(runif(800*200), nrow=800, ncol=200)


#k <- 11
# 总共要选10列（含第一列）
greedy_select <- function(mat, k) {
n <- ncol(mat)
if (k < 1 || k > n) stop("k 必须在 1 到 ncol(mat) 之间")
# ---------- 初始化 ----------
first_col <- 1               # 以第一列为基准（您也可改为其他列）
sel <- first_col
cur <- mat[, first_col]      # 当前每行的最小值
total <- sum(cur)            # 当前总和
reductions <- numeric(k)     # 存储每一步的减少量（第1步为0）

rem <- setdiff(1:n, sel)     # 待选列

# ---------- 贪心选择其余9列 ----------
for (step in 2:k) {
  # 计算加入剩余每一列后，总和的减少量（只计下降）
  gain <- sapply(rem, function(j) {
    total - sum(pmin(cur, mat[, j]))
  })

  best_idx <- rem[which.max(gain)]  # 使减少量最大的列
  best_gain <- max(gain)

  # 更新状态
  sel <- c(sel, best_idx)
  cur <- pmin(cur, mat[, best_idx])
  total <- sum(cur)
  reductions[step] <- best_gain
  rem <- setdiff(rem, best_idx)
}
list(selected = sel, total = total, reductions = reductions)
}

detail_df <- data.frame()   # 存储每一步的细节
summary_df <- data.frame()  # 存储每个 k 的汇总

for (k in 2:51) {
  res <- greedy_select(mat, k)

  # 汇总信息
  summary_df <- rbind(summary_df, data.frame(
    k = k,
    final_sum = res$total,
    total_reduction = sum(res$reductions)
  ))

  # 详细步骤（从第2步开始）
  for (step in 2:k) {
    detail_df <- rbind(detail_df, data.frame(
      k = k,
      step = step,
      added_col = res$selected[step],
      reduction = res$reductions[step]
    ))
  }
}

write.csv(detail_df, file = paste0(prefix,"reduction_details.csv"), row.names = FALSE)
write.csv(summary_df, file =  paste0(prefix,"reduction_summary.csv"), row.names = FALSE)

# 也可输出一个易读的文本文件（带格式）
sink( paste0(prefix,"reduction_report.txt"))
cat("========== 各 k 的详细减少量 ==========\n\n")
for (k in 2:51) {
  sub <- detail_df[detail_df$k == k, ]
  cat(sprintf("k = %d\n", k))
  cat("步骤 | 新增列 | 减少量\n")
  for (i in 1:nrow(sub)) {
    cat(sprintf("%4d | %7d | %10.4f\n", sub$step[i], sub$added_col[i], sub$reduction[i]))
  }
  cat("\n")
}
sink()

# ============================================
# 5. 打印完成信息
# ============================================
cat("已生成文件：\n")
cat(paste0(prefix,"reduction_details.csv  (每一步的详细减少量)\n"))
cat(paste0(prefix,"reduction_summary.csv   (每个 k 的最终总和与总减少量)\n"))
cat(paste0(prefix,"reduction_report.txt    (格式化的文本报告)\n"))
