# 安装和加载必要的库
install.packages("stringr")
# 安装和加载必要的库
install.packages(c("data.table", "stringi", "foreach", "doParallel"))
library(data.table)
library(stringi)
library(foreach)
library(doParallel)
library(stringr)
# 设置并行计算
cores <- parallel::detectCores()
cl <- makeCluster(cores - 1)  # 使用除主核外的所有核心
registerDoParallel(cl)

# 在并行环境中加载必要的库
clusterEvalQ(cl, {
  library(data.table)
  library(stringi)
})

# 获取文件夹路径
folder_path <- "C:\\Users\\ZHH\\Desktop\\data1"

# 读取所有txt文件内容
file_paths <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)
all_text <- sapply(file_paths, function(file) {
  paste(readLines(file, warn = FALSE), collapse = " ")
}, simplify = TRUE)

# 定义关键词和正则表达式
keywords_group1 <- c("供应商", "采购")
keywords_group2 <- c("绿色", "环境", "环保", "能源", "排放", "节能", "ESG", "绿链", "绿名单", "白名单", "社会责任", "负责任", "有害", "减排", "可持续")
keywords_group3 <- c("调查", "准入", "规定", "制度", "资质", "要求", "管理", "考核", "评估", "评审", "评价", "考核", "审核", "审查", "稽核", "查核", "考察", "选择", "依据", "确认", "监督")

# 遍历文本，匹配符合条件的句子
matched_sentences <- data.frame(Text_File = character(), Extracted_Sentence = character())

for (i in seq_along(all_text)) {
  for (word1 in keywords_group1) {
    for (word2 in keywords_group2) {
      for (word3 in keywords_group3) {
        # 构建正则表达式模式
        pattern <- paste0("(.{0,20}", word1, ").{0,20}(", word2, ").{0,20}(", word3, ".{0,20})")
        matches <- str_extract_all(all_text[i], pattern)
        
        # 提取匹配的句子并添加到数据框
        for (match in unlist(matches)) {
          matched_sentences <- rbind(matched_sentences, data.frame(Text_File = file_paths[i], Extracted_Sentence = match))
        }
      }
    }
  }
}

# 合并来自同一个文本的句子到一行，之间用|隔开
merged_sentences <- aggregate(Extracted_Sentence ~ Text_File, data = matched_sentences, paste, collapse = "|")

# 保存结果到csv文件
write.csv(merged_sentences, file = "C:\\Users\\ZHH\\Desktop\\data1\\extracted_sentences.csv", row.names = FALSE)
