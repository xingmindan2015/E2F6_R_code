rm(list = ls())
#install.packages('DAAG')
library(DAAG)
load('figure06_two_group/tumor_e2f6_top2a_deg_relationship_above5.Rdata')

e2f6_genes = b2$row[b2$column == 'E2F6']
top2a_genes = b2$row[b2$column == 'TOP2A']
common_genes = intersect(e2f6_genes,top2a_genes)
e2f6_genes = e2f6_genes[!e2f6_genes%in%common_genes]
top2a_genes = top2a_genes[!top2a_genes%in%common_genes]
hub_genes = c('E2F6','TOP2A')

type_df = data.frame(nodes = c(hub_genes,e2f6_genes,top2a_genes,common_genes),
                     type = c(rep('hub_genes',length(hub_genes)),
                              rep('e2f6_genes',length(e2f6_genes)),
                              rep('top2a_genes',length(top2a_genes)),
                              rep('common_genes',length(common_genes))))
type_df = type_df[!duplicated(type_df$nodes),]
save(type_df,file = 'figure06_two_group/e2f6_top2a_two_group_deg_correlation.Rdata')
write.csv(b2,file = 'figure06_two_group/correlation_genes_e2f6_top2a_network.csv',row.names = F,quote = F)
write.csv(type_df,file = 'figure06_two_group/correlation_genes_e2f6_top2a_network_node_type.csv',row.names = F,quote = F)




# GO富集分析 ------------------------------------------------------------------

rm(list=ls())
load('figure06_two_group/e2f6_top2a_two_group_deg_correlation.Rdata')
library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)#画图包
type_df = type_df[-(1:2),]
genes = type_df$nodes
gene.df <- bitr(genes,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)#TCGA数据框如果没有进行基因注释，那么fromType应该是Ensembl，各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求                     
gene <- gene.df$ENTREZID
ego_ALL <- enrichGO(gene = gene,#我们上面定义了
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 0.05,#P值可以取0.05
                    qvalueCutoff = 0.25,
                    readable = TRUE)
ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.25,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.25,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.25,
                   readable = TRUE)
#4、将结果保存到当前路径
ego_ALL <- as.data.frame(ego_ALL)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)
ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)
write.csv(ego_ALL,file = "figure06_two_group/two_deg_ego_ALL1.csv",row.names = T)
save(ego,ego_ALL,ego_result_BP,ego_result_CC,ego_BP,ego_CC,ego_MF,
     file = "figure06_two_group/two_deg_ego_ALL1.Rdata")
load('figure06_two_group/two_deg_ego_ALL1.Rdata')
ego[1:4,]
table(ego$ID)
ego_result_CC= arrange(ego_result_CC,Count)
ego_result_MF= arrange(ego_result_MF,Count)

