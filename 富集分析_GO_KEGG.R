#【GO富集分析】
#【KEGG富集分析】
# 1.下载/导入包并设置当前工作路径
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
BiocManager::install(c("DOSE","topGO","clusterProfiler","pathview"))
BiocManager::install("GOplot")
BiocManager::install("goProfiles")
BiocManager::install('KEGG.db')
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(GOplot)
library(goProfiles)
library(ggplot2)
library(KEGG.db)
setwd('D:/')

# 2.人类基因数据库编号的类型及导入TP53等一类含有的856个基因的基因集数据
keytypes(org.Hs.eg.db)
MyGeneSet_table<-read.table(file.choose())
MyGeneSet<-as.character(MyGeneSet_table$V1)
typeof(MyGeneSet)

# 3.编号转换SYMBOL->ENTREZID
MyGeneIDSet=bitr(MyGeneSet,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
MyGeneIDSet

# 4.GO富集分析
data(geneList,package = "DOSE")  
names(geneList)
display_number = c(15, 10, 15)    #将BP、MF、CC绘制于一个图中时选择的术语个数
#（1）参数设置为ALL，对MF，CC，BP进行富集分析
ego_all<-enrichGO(gene=MyGeneIDSet$ENTREZID,
                  universe = names(geneList),
                  OrgDb = org.Hs.eg.db,
                  ont="ALL",     #该参数也可以为ALL/CC/BP/MF
                  pAdjustMethod = "BH",    #该参数还可以为holm/hochberg/hommel/bonferroni/BH/BY/fdr/none
                  pvalueCutoff = 0.5,
                  qvalueCutoff = 1,
                  readable = TRUE)
#（2）参数设置为MF，对分子功能进行富集分析
ego_MF<-enrichGO(gene=MyGeneIDSet$ENTREZID,
                 universe = names(geneList),
                 OrgDb = org.Hs.eg.db,
                 ont="MF",     #该参数也可以为ALL/CC/BP/MF
                 pAdjustMethod = "BH",    #该参数还可以为holm/hochberg/hommel/bonferroni/BH/BY/fdr/none
                 pvalueCutoff = 0.5,
                 qvalueCutoff = 1,
                 readable = TRUE)
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[1], ]
#（3）参数设置为CC，对细胞组分进行富集分析
ego_CC<-enrichGO(gene=MyGeneIDSet$ENTREZID,
                 universe = names(geneList),
                 OrgDb = org.Hs.eg.db,
                 ont="CC",     #该参数也可以为ALL/CC/BP/MF
                 pAdjustMethod = "BH",    #该参数还可以为holm/hochberg/hommel/bonferroni/BH/BY/fdr/none
                 pvalueCutoff = 0.5,
                 qvalueCutoff = 1,
                 readable = TRUE)
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
#（4）参数设置为BP，对生物学过程进行富集分析
ego_BP<-enrichGO(gene=MyGeneIDSet$ENTREZID,
                 universe = names(geneList),
                 OrgDb = org.Hs.eg.db,
                 ont="BP",     #该参数也可以为ALL/CC/BP/MF
                 pAdjustMethod = "BH",    #该参数还可以为holm/hochberg/hommel/bonferroni/BH/BY/fdr/none
                 pvalueCutoff = 0.5,
                 qvalueCutoff = 1,
                 readable = TRUE)
ego_result_BP <- na.omit(as.data.frame(ego_BP)[1:display_number[3], ])
# 富集结果的文件的输出
write.csv(as.data.frame(ego_all),"ALL-enrich.csv",row.names = FALSE)
write.csv(as.data.frame(ego_MF),"MF-enrich.csv",row.names = FALSE)
write.csv(as.data.frame(ego_BP),"BP-enrich.csv",row.names = FALSE)
write.csv(as.data.frame(ego_CC),"CC-enrich.csv",row.names = FALSE)

# 5.GO富集分析结果可视化
#（1）绘制柱状图
barplot(ego_all,showCategory=20,title="EnrichmentGO_all")
barplot(ego_MF,showCategory=20,title="EnrichmentGO_MF")
barplot(ego_BP,showCategory=20,title="EnrichmentGO_BP")
barplot(ego_CC,showCategory=20,title="EnrichmentGO_CC")
#（2）用于只能绘制BP、MF、CC的网络图
goplot(ego_BP)
goplot(ego_MF)
goplot(ego_CC)
#（3）点图
dotplot(ego_all,showCategory=20,title="GO_ALL")
dotplot(ego_MF,showCategory=20,title="GO_MF")
dotplot(ego_BP,showCategory=20,title="GO_BP")
dotplot(ego_CC,showCategory=20,title="GO_CC")
#（4）使用plotProfile包绘制柱状图
plot_count<-function(input,title)
{  
  GO_ALL<-as.data.frame(input@result)
  data<-GO_ALL[1:40,c("ID","Description","Count")]
  draw_data_all<-data.frame(GOID=data$ID,Description=data$Description,Frequency=data$Count)
  plotProfiles(list(draw_data_all),aTitle=title,labelWidth = 27)
}
plot_count(ego_all,"GO_ALL")
plot_count(ego_MF,"GO_MF")
plot_count(ego_BP,"GO_BP")
plot_count(ego_CC,"GO_CC")
#（5）绘制BP，MF，CC的分块图
go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                           Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                           GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                           type=factor(c(rep("biological process", display_number[1]), rep("cellular component", display_number[2]),
                                         rep("molecular function", display_number[3])), levels=c("molecular function", "cellular component", "biological process")))
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}
labels=(sapply(
  levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
  shorten_names))
names(labels) = rev(1:nrow(go_enrich_df))
# 绘图
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_bw() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")
p

# 5.代谢通路的富集分析
ego_kegg<-enrichKEGG(gene=MyGeneIDSet$ENTREZID,use_internal_data = TRUE)
kegg<-as.data.frame(ego_kegg@result)
write.csv(kegg,"KEGG-enrich.csv",row.names = FALSE)
#（1）点图
dotplot(ego_kegg,showCategory=20)
#（2）柱状图
barplot(ego_kegg,showCategory=20,title="ego_KEGG")
#（3）差异基因关联图
emapplot(ego_kegg,showCategory = 30)
#（4）通路图
browseKEGG(ego_kegg,"hsa04934")

