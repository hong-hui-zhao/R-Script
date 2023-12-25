# 安装和加载必要的库
install.packages("stringr")
library(stringr)

# 获取文件夹路径
folder_path <- "C:\\Users\\ZHH\\Desktop\\data1"

# 读取所有txt文件内容
file_paths <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)
all_text <- sapply(file_paths, function(file) {
  paste(readLines(file, warn = FALSE), collapse = " ")
}, simplify = TRUE)

# 定义关键词和正则表达式
keywords_group1 <- c("环境", "环保", "资源", "能源", "循环", "利用率")
keywords_group2 <- c("薪酬", "考核", "晋升", "奖")

# 遍历文本，匹配符合条件的句子
matched_sentences <- data.frame(Text_File = character(), Extracted_Sentence = character())

for (i in seq_along(all_text)) {
  for (word1 in keywords_group1) {
    for (word2 in keywords_group2) {
      pattern <- paste0("(.{0,20}", word1, ").{0,20}(", word2, ".{0,20})")
      matches <- str_extract_all(all_text[i], pattern)
      
      for (match in unlist(matches)) {
        matched_sentences <- rbind(matched_sentences, data.frame(Text_File = file_paths[i], Extracted_Sentence = match))
      }
    }
  }
}

# 合并来自同一个文本的句子到一行，之间用|隔开
merged_sentences <- aggregate(Extracted_Sentence ~ Text_File, data = matched_sentences, paste, collapse = "|")

# 保存结果到csv文件
write.csv(merged_sentences, file = "C:\\Users\\ZHH\\Desktop\\data1\\extracted_sentences.csv", row.names = FALSE,fileEncoding = 'GB18030')
