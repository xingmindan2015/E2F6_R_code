
# TCGA_ssGSEA_GSVA_biaoxing heatmap-----------------------------------------------

rm(list=ls())
load('figure08_zhengligeneset/tcga_phenotype_ssgsea.Rdata')
re_phenotype_ssgsea=re
load('figure08_zhengligeneset/tcga_pathway_ssgsea.Rdata')
re_pathway_ssgsea = re
# load('figure08_zhengligeneset/tcga_phenotype_gsva.Rdata')
# re_phenotype_gsva = re
# load('figure08_zhengligeneset/tcga_pathway_gsva.Rdata')
# re_pathway_gsva = re
load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
up_phenotype=up
down_phenotype = down
cg_phenotype = c(up_phenotype,down_phenotype)
load('figure09_tongluzaifenxi/pathway_deg.Rdata')
up_pathway = up
down_pathway = down
cg_pathway = c(up_pathway,down_pathway)
re_phenotype_ssgsea = re_phenotype_ssgsea[cg_phenotype,]
re_pathway_ssgsea = re_pathway_ssgsea[cg_pathway,]
# re_phenotype_gsva = re_phenotype_gsva[cg_phenotype,]
# re_pathway_gsva = re_pathway_gsva[cg_pathway,]
#rownames(re_phenotype_ssgsea) = paste(rownames(re_phenotype_ssgsea),'ssGSEA')
#rownames(re_phenotype_gsva)=paste(rownames(re_phenotype_gsva),'GSVA')
#rownames(re_pathway_ssgsea)=paste(rownames(re_pathway_ssgsea),'ssGSEA')
#rownames(re_pathway_gsva)=paste(rownames(re_pathway_gsva),'GSVA')
identical(colnames(re_phenotype_ssgsea),colnames(re_pathway_ssgsea))
# identical(colnames(re_pathway_ssgsea),colnames(re_phenotype_gsva))
# identical(colnames(re_pathway_ssgsea),colnames(re_pathway_gsva))
re = rbind(re_phenotype_ssgsea,re_pathway_ssgsea)
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
identical(rownames(meta),colnames(re))
Group1_name=rownames(meta)[meta$Group=='Others']
Group2_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']
re1=re[,Group1_name]
re2=re[,Group2_name]
re=cbind(re1,re2)
meta=meta[match(colnames(re),rownames(meta)),]
identical(rownames(meta),colnames(re))
an_col = data.frame(Group = meta$Group,
                    Cluster = meta$Cluster,
                    row.names = colnames(re))
an_col = arrange(an_col,Group,Cluster)
re = re[,match(rownames(an_col),colnames(re))]
identical(rownames(an_col),colnames(re))
an = data.frame(Group = an_col$Group,
                row.names = colnames(re))
#method_legend = ifelse(str_detect(rownames(re),'ssGSEA'),'ssGSEA','GSVA')
#names(method_legend)=rownames(re)
way_type_legend = ifelse(str_detect(rownames(re),'pathway'),'Signaling Pathway','Phenotype')
names(way_type_legend)=rownames(re)
table(way_type_legend)
an_row = data.frame(Type = way_type_legend,
                    #Method = method_legend,
                    row.names = rownames(re))
an2 = data.frame(Group = an_row$Type,
                 row.names = rownames(re))
dev.off()
class(re)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
table(meta$Group)
legend_col = list(Group = c('Others'='#6d9cc5','E2F6_H&TOP2A_H'='#e7857b'),
                  Cluster=c('Cluster1'='#a5d3ed','Cluster2'='#b5aad4'),
                  #Method = c('ssGSEA'='#A4CBCC','GSVA'='#F7D1E6'),
                  Type=c('Phenotype'='#dddb97','Signaling Pathway'='#b4c5e6')                  )
pdf('figure09_tupianzhengli/tcga_ssgsea_cluster_group_phenotype_pathway_heatmap.pdf',height = 6,width = 10)
pheatmap(re,scale = "row",
         show_colnames = F,
         annotation_col = an_col,
         annotation_row = an_row,
         cluster_cols = F,
         cluster_rows = F,
         column_split=an,
         row_split = an2,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         annotation_colors = legend_col)

dev.off()

# ICGC biaoxing tonglu ssGSEA GSVA heatmap --------------------------------


rm(list=ls())
load('figure08_zhengligeneset/icgc_phenotype_ssgsea.Rdata')
re_phenotype_ssgsea=re
load('figure08_zhengligeneset/icgc_pathway_ssgsea.Rdata')
re_pathway_ssgsea = re
# load('figure08_zhengligeneset/icgc_phenotype_gsva.Rdata')
# re_phenotype_gsva = re
# load('figure08_zhengligeneset/icgc_pathway_gsva.Rdata')
# re_pathway_gsva = re
load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
up_phenotype=up
down_phenotype = down
cg_phenotype = c(up_phenotype,down_phenotype)
load('figure09_tongluzaifenxi/pathway_deg.Rdata')
up_pathway = up
down_pathway = down
cg_pathway = c(up_pathway,down_pathway)
re_phenotype_ssgsea = re_phenotype_ssgsea[cg_phenotype,]
re_pathway_ssgsea = re_pathway_ssgsea[cg_pathway,]
# re_phenotype_gsva = re_phenotype_gsva[cg_phenotype,]
# re_pathway_gsva = re_pathway_gsva[cg_pathway,]
# rownames(re_phenotype_ssgsea) = paste(rownames(re_phenotype_ssgsea),'ssGSEA')
# rownames(re_phenotype_gsva)=paste(rownames(re_phenotype_gsva),'GSVA')
# rownames(re_pathway_ssgsea)=paste(rownames(re_pathway_ssgsea),'ssGSEA')
# rownames(re_pathway_gsva)=paste(rownames(re_pathway_gsva),'GSVA')
identical(colnames(re_phenotype_ssgsea),colnames(re_pathway_ssgsea))
# identical(colnames(re_pathway_ssgsea),colnames(re_phenotype_gsva))
# identical(colnames(re_pathway_ssgsea),colnames(re_pathway_gsva))
re = rbind(re_phenotype_ssgsea,re_pathway_ssgsea)
load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
identical(rownames(meta),colnames(re))
Group1_name=rownames(meta)[meta$Group=='Others']
Group2_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']
re1=re[,Group1_name]
re2=re[,Group2_name]
re=cbind(re1,re2)
meta=meta[match(colnames(re),rownames(meta)),]
identical(rownames(meta),colnames(re))
an_col = data.frame(Group = meta$Group,
                    Cluster = meta$Cluster,
                    row.names = colnames(re))
an_col = arrange(an_col,Group,Cluster)
re = re[,match(rownames(an_col),colnames(re))]
identical(rownames(an_col),colnames(re))
an = data.frame(Group = an_col$Group,
                row.names = colnames(re))
#method_legend = ifelse(str_detect(rownames(re),'ssGSEA'),'ssGSEA','GSVA')
#names(method_legend)=rownames(re)
way_type_legend = ifelse(str_detect(rownames(re),'pathway'),'Signaling Pathway','Phenotype')
names(way_type_legend)=rownames(re)
table(way_type_legend)
an_row = data.frame(Type = way_type_legend,
                    #Method = method_legend,
                    row.names = rownames(re))
an2 = data.frame(Group = an_row$Type,
                 row.names = rownames(re))
dev.off()
class(re)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
table(meta$Group)
legend_col = list(Group = c('Others'='#6d9cc5','E2F6_H&TOP2A_H'='#e7857b'),
                  Cluster=c('Cluster1'='#a5d3ed','Cluster2'='#b5aad4'),
                  #Method = c('ssGSEA'='#A4CBCC','GSVA'='#F7D1E6'),
                  Type=c('Phenotype'='#dddb97','Signaling Pathway'='#b4c5e6')                  )
pdf('figure09_tupianzhengli/icgc_ssgsea_cluster_group_phenotype_pathway_heatmap.pdf',height = 6,width = 10)
pheatmap(re,scale = "row",
         show_colnames = F,
         annotation_col = an_col,
         annotation_row = an_row,
         cluster_cols = F,
         cluster_rows = F,
         column_split=an,
         row_split = an2,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         annotation_colors = legend_col)

dev.off()

# TCGA correlation heatmap phenotype ssgsea------------------------------
# #ssgsea####
# rm(list = ls())
# load('figure02_E2F6/lihc_sur_data.Rdata')
# load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
# load('figure08_zhengligeneset/tcga_phenotype_ssgsea.Rdata')
# re_phenotype_ssgsea=re
# load('figure08_zhengligeneset/tcga_pathway_ssgsea.Rdata')
# re_pathway_ssgsea = re
# identical(colnames(re_phenotype_ssgsea),colnames(re_pathway_ssgsea))
# re = rbind(re_phenotype_ssgsea,re_pathway_ssgsea)
# load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
# up_phenotype=up
# down_phenotype = down
# cg_phenotype = c(up_phenotype,down_phenotype)
# load('figure09_tongluzaifenxi/pathway_deg.Rdata')
# up_pathway = up
# down_pathway = down
# cg_pathway = c(up_pathway,down_pathway)
# re=re[c(cg_phenotype,cg_pathway),]
# 
# exprSet[1:4,1:4]
# exp=exprSet
# identical(rownames(meta),colnames(re))
# exp=exp[,match(colnames(re),colnames(exp))]
# identical(colnames(re),colnames(exp))
# exp=exp[c('E2F6','TOP2A'),]
# identical(colnames(exp),rownames(meta))
# data=as.data.frame(t(exp))
# data$PredictScore=meta$Predict_score
# library(vegan)
# library(dplyr)
# library(ggcor)
# library(ggplot2)
# 
# quickcor(t(re), type = "lower",method = "spearman") +
#   geom_square()+
#   scale_fill_gradient2( high = 'orange', mid = 'white',low = 'navyblue')
# df_mantel <- mantel_test(data, t(re), mantel.fun = 'mantel',
#                          spec.dist.method = 'bray',
#                          env.dist.method = 'euclidean',
#                          spec.select = list(E2F6 = 1,
#                                             TOP2A = 2,
#                                             PredictScore = 3))
# df_mantel <- df_mantel %>%
#   mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.2, 0.4, Inf),
#                     labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.4", ">= 0.4")),#定义Mantel的R值范围标签
#          df_p = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
#                     labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# save(df_mantel,file='figure09_tongluzaifenxi/tcga_ssgsesa_pathway_phenotype_corrlation.Rdata')
# load('figure09_tongluzaifenxi/tcga_ssgsesa_pathway_phenotype_corrlation.Rdata')
# quickcor(t(re),method = "pearson", type = "upper", cor.test = T, cluster.type = "all") +#环境因子之间的相关性热图
#   geom_square() +#相关性显示形式
#   geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+#显著性标签
#   scale_fill_gradient2( high = 'orange', mid = 'white',low =  'navyblue') + #颜色设置
#   anno_link(df_mantel, aes(color = df_p,
#                            size = df_r))+
#   scale_size_manual(values = c(0.5, 1, 1.5, 2))+#连线粗细设置
#   scale_color_manual(values = c("orange","#6387BB","#937abf"))+#线条颜色设置
#   guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
#          size = guide_legend(title = "Mantel's r",order = 2),
#          color = guide_legend(title = "Mantel's p", order = 3),
#          linetype = "none")
# ggsave('figure09_tupianzhengli/tcga_correlation_phenotype_pathway_network_ssgsea.pdf',width = 10,height = 10)
# #GSVA####
# rm(list = ls())
# load('figure02_E2F6/lihc_sur_data.Rdata')
# load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
# load('figure08_zhengligeneset/tcga_phenotype_gsva.Rdata')
# re_phenotype_ssgsea=re
# load('figure08_zhengligeneset/tcga_pathway_gsva.Rdata')
# re_pathway_ssgsea = re
# identical(colnames(re_phenotype_ssgsea),colnames(re_pathway_ssgsea))
# re = rbind(re_phenotype_ssgsea,re_pathway_ssgsea)
# load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
# up_phenotype=up
# down_phenotype = down
# cg_phenotype = c(up_phenotype,down_phenotype)
# load('figure09_tongluzaifenxi/pathway_deg.Rdata')
# up_pathway = up
# down_pathway = down
# cg_pathway = c(up_pathway,down_pathway)
# re=re[c(cg_phenotype,cg_pathway),]
# 
# exprSet[1:4,1:4]
# exp=exprSet
# identical(rownames(meta),colnames(re))
# exp=exp[,match(colnames(re),colnames(exp))]
# identical(colnames(re),colnames(exp))
# exp=exp[c('E2F6','TOP2A'),]
# identical(colnames(exp),rownames(meta))
# data=as.data.frame(t(exp))
# data$PredictScore=meta$Predict_score
# library(vegan)
# library(dplyr)
# library(ggcor)
# library(ggplot2)
# 
# quickcor(t(re), type = "lower",method = "spearman") +
#   geom_square()+
#   scale_fill_gradient2( high = 'orange', mid = 'white',low = 'navyblue')
# df_mantel <- mantel_test(data, t(re), mantel.fun = 'mantel',
#                          spec.dist.method = 'bray', 
#                          env.dist.method = 'euclidean',
#                          spec.select = list(E2F6 = 1,
#                                             TOP2A = 2,
#                                             PredictScore = 3))
# df_mantel <- df_mantel %>%
#   mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.2, 0.4, Inf),
#                     labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.4", ">= 0.4")),#定义Mantel的R值范围标签
#          df_p = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
#                     labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# save(df_mantel,file='figure09_tongluzaifenxi/tcga_gsva_pathway_phenotype_corrlation.Rdata')
# quickcor(t(re),method = "pearson", type = "upper", cor.test = T, cluster.type = "all") +#环境因子之间的相关性热图
#   geom_square() +#相关性显示形式
#   geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+#显著性标签
#   scale_fill_gradient2( high = 'orange', mid = 'white',low =  'navyblue') + #颜色设置
#   anno_link(df_mantel, aes(color = df_p,
#                            size = df_r))+
#   scale_size_manual(values = c(0.5, 1, 1.5, 2))+#连线粗细设置
#   scale_color_manual(values = c("orange","#6387BB","#937abf"))+#线条颜色设置
#   guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
#          size = guide_legend(title = "Mantel's r",order = 2),
#          color = guide_legend(title = "Mantel's p", order = 3),
#          linetype = "none")
# ggsave('figure09_tupianzhengli/tcga_correlation_phenotype_pathway_network_gsva.pdf',width = 10,height = 10)


# ICGC phenotype pathway correlation  -------------------------------------

# #ssgsea####
# rm(list = ls())
# load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
# load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
# load('figure08_zhengligeneset/icgc_phenotype_ssgsea.Rdata')
# re_phenotype_ssgsea=re
# load('figure08_zhengligeneset/icgc_pathway_ssgsea.Rdata')
# re_pathway_ssgsea = re
# identical(colnames(re_phenotype_ssgsea),colnames(re_pathway_ssgsea))
# re = rbind(re_phenotype_ssgsea,re_pathway_ssgsea)
# load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
# up_phenotype=up
# down_phenotype = down
# cg_phenotype = c(up_phenotype,down_phenotype)
# load('figure09_tongluzaifenxi/pathway_deg.Rdata')
# up_pathway = up
# down_pathway = down
# cg_pathway = c(up_pathway,down_pathway)
# re=re[c(cg_phenotype,cg_pathway),]
# 
# exprSet[1:4,1:4]
# exp=exprSet
# identical(rownames(meta),colnames(re))
# exp=exp[,match(colnames(re),colnames(exp))]
# identical(colnames(re),colnames(exp))
# exp=exp[c('E2F6','TOP2A'),]
# identical(colnames(exp),rownames(meta))
# data=as.data.frame(t(exp))
# data$PredictScore=meta$Predict_score
# library(vegan)
# library(dplyr)
# library(ggcor)
# library(ggplot2)
# 
# quickcor(t(re), type = "lower",method = "spearman") +
#   geom_square()+
#   scale_fill_gradient2( high = 'orange', mid = 'white',low = 'navyblue')
# df_mantel <- mantel_test(data, t(re), mantel.fun = 'mantel',
#                          spec.dist.method = 'bray',
#                          env.dist.method = 'euclidean',
#                          spec.select = list(E2F6 = 1,
#                                             TOP2A = 2,
#                                             PredictScore = 3))
# df_mantel <- df_mantel %>%
#   mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.2, 0.4, Inf),
#                     labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.4", ">= 0.4")),#定义Mantel的R值范围标签
#          df_p = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
#                     labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# save(df_mantel,file='figure09_tongluzaifenxi/icgc_ssgsesa_pathway_phenotype_corrlation.Rdata')
# load('figure09_tongluzaifenxi/icgc_ssgsesa_pathway_phenotype_corrlation.Rdata')
# quickcor(t(re),method = "pearson", type = "upper", cor.test = T, cluster.type = "all") +#环境因子之间的相关性热图
#   geom_square() +#相关性显示形式
#   geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+#显著性标签
#   scale_fill_gradient2( high = 'orange', mid = 'white',low =  'navyblue') + #颜色设置
#   anno_link(df_mantel, aes(color = df_p,
#                            size = df_r))+
#   scale_size_manual(values = c(0.5, 1, 1.5, 2))+#连线粗细设置
#   scale_color_manual(values = c("orange","#6387BB","#937abf"))+#线条颜色设置
#   guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
#          size = guide_legend(title = "Mantel's r",order = 2),
#          color = guide_legend(title = "Mantel's p", order = 3),
#          linetype = "none")
# ggsave('figure09_tupianzhengli/icgc_correlation_phenotype_pathway_network_ssgsea.pdf',width = 10,height = 10)
#GSVA####
# rm(list = ls())
# load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
# load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
# load('figure08_zhengligeneset/icgc_phenotype_gsva.Rdata')
# re_phenotype_ssgsea=re
# load('figure08_zhengligeneset/icgc_pathway_gsva.Rdata')
# re_pathway_ssgsea = re
# identical(colnames(re_phenotype_ssgsea),colnames(re_pathway_ssgsea))
# re = rbind(re_phenotype_ssgsea,re_pathway_ssgsea)
# load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
# up_phenotype=up
# down_phenotype = down
# cg_phenotype = c(up_phenotype,down_phenotype)
# load('figure09_tongluzaifenxi/pathway_deg.Rdata')
# up_pathway = up
# down_pathway = down
# cg_pathway = c(up_pathway,down_pathway)
# re=re[c(cg_phenotype,cg_pathway),]
# 
# exprSet[1:4,1:4]
# exp=exprSet
# identical(rownames(meta),colnames(re))
# exp=exp[,match(colnames(re),colnames(exp))]
# identical(colnames(re),colnames(exp))
# exp=exp[c('E2F6','TOP2A'),]
# identical(colnames(exp),rownames(meta))
# data=as.data.frame(t(exp))
# data$PredictScore=meta$Predict_score
# library(vegan)
# library(dplyr)
# library(ggcor)
# library(ggplot2)
# 
# quickcor(t(re), type = "lower",method = "spearman") +
#   geom_square()+
#   scale_fill_gradient2( high = 'orange', mid = 'white',low = 'navyblue')
# df_mantel <- mantel_test(data, t(re), mantel.fun = 'mantel',
#                          spec.dist.method = 'bray', 
#                          env.dist.method = 'euclidean',
#                          spec.select = list(E2F6 = 1,
#                                             TOP2A = 2,
#                                             PredictScore = 3))
# df_mantel <- df_mantel %>%
#   mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.2, 0.4, Inf),
#                     labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.4", ">= 0.4")),#定义Mantel的R值范围标签
#          df_p = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
#                     labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
# save(df_mantel,file='figure09_tongluzaifenxi/icgc_gsva_pathway_phenotype_corrlation.Rdata')
# quickcor(t(re),method = "pearson", type = "upper", cor.test = T, cluster.type = "all") +#环境因子之间的相关性热图
#   geom_square() +#相关性显示形式
#   geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+#显著性标签
#   scale_fill_gradient2( high = 'orange', mid = 'white',low =  'navyblue') + #颜色设置
#   anno_link(df_mantel, aes(color = df_p,
#                            size = df_r))+
#   scale_size_manual(values = c(0.5, 1, 1.5, 2))+#连线粗细设置
#   scale_color_manual(values = c("orange","#6387BB","#937abf"))+#线条颜色设置
#   guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
#          size = guide_legend(title = "Mantel's r",order = 2),
#          color = guide_legend(title = "Mantel's p", order = 3),
#          linetype = "none")
# ggsave('figure09_tupianzhengli/icgc_correlation_phenotype_pathway_network_gsva.pdf',width = 10,height = 10)

# TCGA 相关性热图ssGSEA --------------------------------------------------------

#ssGSEA####
rm(list = ls())
load('figure02_E2F6/lihc_sur_data.Rdata')
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
load('figure08_zhengligeneset/tcga_phenotype_ssgsea.Rdata')
re_phenotype_ssgsea=re
load('figure08_zhengligeneset/tcga_pathway_ssgsea.Rdata')
re_pathway_ssgsea = re
identical(colnames(re_phenotype_ssgsea),colnames(re_pathway_ssgsea))
re = rbind(re_phenotype_ssgsea,re_pathway_ssgsea)
load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
up_phenotype=up
down_phenotype = down
cg_phenotype = c(up_phenotype,down_phenotype)
load('figure09_tongluzaifenxi/pathway_deg.Rdata')
up_pathway = up
down_pathway = down
cg_pathway = c(up_pathway,down_pathway)
re=re[c(cg_phenotype,cg_pathway),]

identical(rownames(meta),colnames(re))
exp=exprSet
exp=exp[,match(colnames(re),colnames(exp))]
identical(colnames(re),colnames(exp))
exp=exp[c('E2F6','TOP2A'),]
identical(colnames(re),colnames(exp))
re=rbind(re,exp)
identical(colnames(exp),rownames(meta))
identical(colnames(re),rownames(meta))
data=as.data.frame(t(re))
data$PredictScore=meta$Predict_score
M=data

library(tidyverse)
library(corrplot)
if(!require(ggcor))devtools::install_github("xukaili/ggcor",upgrade = F)
datp <- cor.mtest(M)$p[1:22,23:25] %>%
  as.data.frame()%>% 
  rownames_to_column("A") %>% 
  gather("B","p",-A)
data <- cor(M)[1:22,23:25] %>%
  as.data.frame()%>% 
  rownames_to_column("A") %>% 
  gather("B","r",-A) %>% 
  mutate(p = datp$p)
data$ps = case_when(data$p<0.01~"**",
                    data$p<0.05~"*",
                    T~"")
library(vctrs)
library(grid)

geom_rectriangle <- function(mapping = NULL, data = NULL,
                             stat = "identity", position = "identity",
                             ...,
                             linejoin = "mitre",
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRectriangle,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      linejoin = linejoin,
      na.rm = na.rm,
      ...
    )
  )
}

GeomRectriangle <- ggproto(
  "GeomRectriangle", Geom,
  default_aes = aes(r = 1, colour = "grey35", fill = NA, size = 0.25, linetype = 1,
                    alpha = NA,type = "upper"),
  required_aes = c("x", "y"),
  draw_panel = function(self, data, panel_params, coord, linejoin = "mitre",type = "upper") {
    aesthetics <- setdiff(names(data), c("x", "y"))
    
    polys <- lapply(split(data, seq_len(nrow(data))), function(row) {
      rectriangle <- point_to_rectriangle(row$x, row$y, row$r, row$type)
      aes <- new_data_frame(row[aesthetics])[rep(1, 4), ]
      GeomPolygon$draw_panel(cbind(rectriangle, aes), panel_params, coord)
    })
    
    ggplot2:::ggname("geom_rectriangle", do.call("grobTree", polys))
  },
  draw_key = draw_key_polygon
)


point_to_rectriangle <- function(x, y, r, type = type) {
  r <- 0.5 * sign(r) * sqrt(abs(r))
  #r0 = 0.5
  xmin <- - r + x
  xmax <- r + x
  ymin <- - r + y
  ymax <- r + y 
  if(type == "upper"){
    df = new_data_frame(list(
      y = c(ymax, ymax, ymin, ymax),
      x = c(xmin, xmax, xmin, xmin)
    ))
  }else if(type == "lower"){
    df = new_data_frame(list(
      y = c(ymax, ymin, ymin, ymax),
      x = c(xmax, xmax, xmin, xmax) 
    ))
  }
  df
}

data$B=factor(data$B,levels = c('E2F6','TOP2A','PredictScore'))
table(data$A)
data$A = factor(data$A,levels = c('Apoptosis','Necroptosis','Programmed Cell Death','Immunogenic cell death','Cell stemness',
                                  'Epithelial-mesenchymal transition','Invasion','Hypoxia','Fatty acid metabolism',
                                  'm5C','m6A','m7G',
                                  'Adipocytokine signaling pathway', 'AMPK signaling pathway','ErbB signaling pathway',
                                  'Glucagon signaling pathway','Hippo signaling pathway','p53 signaling pathway',
                                  'PPAR signaling pathway','Signaling pathways regulating pluripotency of stem cells','Sphingolipid signaling pathway',
                                  'VEGF signaling pathway'))
table(data$A)
table(data$B)
ggplot()+
  geom_rectriangle(data = data, aes(A, B, fill = -log10(p)),type = "upper", r = 1)+
  scale_fill_gradient(low = "white", high = "#a5d3ed",na.value = "white")+
  ggnewscale::new_scale_fill()+
  geom_rectriangle(data = data, aes(A, B, fill = r),type = "lower", r = 1)+
  scale_fill_gradient2(high = "#e7857b", mid = "white",low = "#6d9cc5")+
  labs(x = "", y = "")+
  geom_text(data = data,aes(A,B,label = ps),
            nudge_y = 0.05)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1))
ggsave('figure09_tupianzhengli/tcga_ssgsea_correlation_hubgene_predict_phenotype_pathway_deg.pdf',height = 6,width =8,dpi = 320,units = 'in')



#GSVA####
# rm(list = ls())
# load('figure02_E2F6/lihc_sur_data.Rdata')
# load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
# load('figure08_zhengligeneset/tcga_phenotype_gsva.Rdata')
# re_phenotype_ssgsea=re
# load('figure08_zhengligeneset/tcga_pathway_gsva.Rdata')
# re_pathway_ssgsea = re
# identical(colnames(re_phenotype_ssgsea),colnames(re_pathway_ssgsea))
# re = rbind(re_phenotype_ssgsea,re_pathway_ssgsea)
# load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
# up_phenotype=up
# down_phenotype = down
# cg_phenotype = c(up_phenotype,down_phenotype)
# load('figure09_tongluzaifenxi/pathway_deg.Rdata')
# up_pathway = up
# down_pathway = down
# cg_pathway = c(up_pathway,down_pathway)
# re=re[c(cg_phenotype,cg_pathway),]
# 
# identical(rownames(meta),colnames(re))
# exp=exprSet
# exp=exp[,match(colnames(re),colnames(exp))]
# identical(colnames(re),colnames(exp))
# exp=exp[c('E2F6','TOP2A'),]
# identical(colnames(re),colnames(exp))
# re=rbind(re,exp)
# identical(colnames(exp),rownames(meta))
# identical(colnames(re),rownames(meta))
# data=as.data.frame(t(re))
# data$PredictScore=meta$Predict_score
# M=data
# 
# library(tidyverse)
# library(corrplot)
# if(!require(ggcor))devtools::install_github("xukaili/ggcor",upgrade = F)
# datp <- cor.mtest(M)$p[1:15,16:18] %>%
#   as.data.frame()%>% 
#   rownames_to_column("A") %>% 
#   gather("B","p",-A)
# data <- cor(M)[1:15,16:18] %>%
#   as.data.frame()%>% 
#   rownames_to_column("A") %>% 
#   gather("B","r",-A) %>% 
#   mutate(p = datp$p)
# data$ps = case_when(data$p<0.01~"**",
#                     data$p<0.05~"*",
#                     T~"")
# library(vctrs)
# library(grid)
# 
# geom_rectriangle <- function(mapping = NULL, data = NULL,
#                              stat = "identity", position = "identity",
#                              ...,
#                              linejoin = "mitre",
#                              na.rm = FALSE,
#                              show.legend = NA,
#                              inherit.aes = TRUE) {
#   layer(
#     data = data,
#     mapping = mapping,
#     stat = stat,
#     geom = GeomRectriangle,
#     position = position,
#     show.legend = show.legend,
#     inherit.aes = inherit.aes,
#     params = list(
#       linejoin = linejoin,
#       na.rm = na.rm,
#       ...
#     )
#   )
# }
# 
# GeomRectriangle <- ggproto(
#   "GeomRectriangle", Geom,
#   default_aes = aes(r = 1, colour = "grey35", fill = NA, size = 0.25, linetype = 1,
#                     alpha = NA,type = "upper"),
#   required_aes = c("x", "y"),
#   draw_panel = function(self, data, panel_params, coord, linejoin = "mitre",type = "upper") {
#     aesthetics <- setdiff(names(data), c("x", "y"))
#     
#     polys <- lapply(split(data, seq_len(nrow(data))), function(row) {
#       rectriangle <- point_to_rectriangle(row$x, row$y, row$r, row$type)
#       aes <- new_data_frame(row[aesthetics])[rep(1, 4), ]
#       GeomPolygon$draw_panel(cbind(rectriangle, aes), panel_params, coord)
#     })
#     
#     ggplot2:::ggname("geom_rectriangle", do.call("grobTree", polys))
#   },
#   draw_key = draw_key_polygon
# )
# 
# 
# point_to_rectriangle <- function(x, y, r, type = type) {
#   r <- 0.5 * sign(r) * sqrt(abs(r))
#   #r0 = 0.5
#   xmin <- - r + x
#   xmax <- r + x
#   ymin <- - r + y
#   ymax <- r + y 
#   if(type == "upper"){
#     df = new_data_frame(list(
#       y = c(ymax, ymax, ymin, ymax),
#       x = c(xmin, xmax, xmin, xmin)
#     ))
#   }else if(type == "lower"){
#     df = new_data_frame(list(
#       y = c(ymax, ymin, ymin, ymax),
#       x = c(xmax, xmax, xmin, xmax) 
#     ))
#   }
#   df
# }
# 
# data$B=factor(data$B,levels = c('E2F6','TOP2A','PredictScore'))
# data$A = factor(data$A,levels = c('Disulfidptosis','Hypoxia',
#                                   'Cell stemness','Epithelial-mesenchymal transition','Invasion',
#                                   'm5C','m6A','m7G',
#                                   'Adipocytokine signaling pathway', 'ErbB signaling pathway',
#                                   'Hippo signaling pathway','p53 signaling pathway',
#                                   'PPAR signaling pathway','Sphingolipid signaling pathway',
#                                   'VEGF signaling pathway'))
# table(data$A)
# table(data$B)
# ggplot()+
#   geom_rectriangle(data = data, aes(A, B, fill = -log10(p)),type = "upper", r = 1)+
#   scale_fill_gradient(low = "white", high = "#66a3b1",na.value = "white")+
#   ggnewscale::new_scale_fill()+
#   geom_rectriangle(data = data, aes(A, B, fill = r),type = "lower", r = 1)+
#   scale_fill_gradient2(high = "red", mid = "white",low = "blue")+
#   labs(x = "", y = "")+
#   geom_text(data = data,aes(A,B,label = ps),
#             nudge_y = 0.05)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1))
# ggsave('figure09_tupianzhengli/tcga_gsva_correlation_hubgene_predict_phenotype_pathway_deg.pdf',height = 6,width =8,dpi = 320,units = 'in')



# ICGC相关性热图 ---------------------------------------------------------------

#ssGSEA####
rm(list = ls())
load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
load('figure08_zhengligeneset/icgc_phenotype_ssgsea.Rdata')
re_phenotype_ssgsea=re
load('figure08_zhengligeneset/icgc_pathway_ssgsea.Rdata')
re_pathway_ssgsea = re
identical(colnames(re_phenotype_ssgsea),colnames(re_pathway_ssgsea))
re = rbind(re_phenotype_ssgsea,re_pathway_ssgsea)
load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
up_phenotype=up
down_phenotype = down
cg_phenotype = c(up_phenotype,down_phenotype)
load('figure09_tongluzaifenxi/pathway_deg.Rdata')
up_pathway = up
down_pathway = down
cg_pathway = c(up_pathway,down_pathway)
re=re[c(cg_phenotype,cg_pathway),]

identical(rownames(meta),colnames(re))
exp=exprSet
exp=exp[,match(colnames(re),colnames(exp))]
identical(colnames(re),colnames(exp))
exp=exp[c('E2F6','TOP2A'),]
identical(colnames(re),colnames(exp))
re=rbind(re,exp)
identical(colnames(exp),rownames(meta))
identical(colnames(re),rownames(meta))
data=as.data.frame(t(re))
data$PredictScore=meta$Predict_score
M=data

library(tidyverse)
library(corrplot)
if(!require(ggcor))devtools::install_github("xukaili/ggcor",upgrade = F)
datp <- cor.mtest(M)$p[1:22,23:25] %>%
  as.data.frame()%>% 
  rownames_to_column("A") %>% 
  gather("B","p",-A)
data <- cor(M)[1:22,23:25] %>%
  as.data.frame()%>% 
  rownames_to_column("A") %>% 
  gather("B","r",-A) %>% 
  mutate(p = datp$p)
data$ps = case_when(data$p<0.01~"**",
                    data$p<0.05~"*",
                    T~"")
library(vctrs)
library(grid)

geom_rectriangle <- function(mapping = NULL, data = NULL,
                             stat = "identity", position = "identity",
                             ...,
                             linejoin = "mitre",
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRectriangle,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      linejoin = linejoin,
      na.rm = na.rm,
      ...
    )
  )
}

GeomRectriangle <- ggproto(
  "GeomRectriangle", Geom,
  default_aes = aes(r = 1, colour = "grey35", fill = NA, size = 0.25, linetype = 1,
                    alpha = NA,type = "upper"),
  required_aes = c("x", "y"),
  draw_panel = function(self, data, panel_params, coord, linejoin = "mitre",type = "upper") {
    aesthetics <- setdiff(names(data), c("x", "y"))
    
    polys <- lapply(split(data, seq_len(nrow(data))), function(row) {
      rectriangle <- point_to_rectriangle(row$x, row$y, row$r, row$type)
      aes <- new_data_frame(row[aesthetics])[rep(1, 4), ]
      GeomPolygon$draw_panel(cbind(rectriangle, aes), panel_params, coord)
    })
    
    ggplot2:::ggname("geom_rectriangle", do.call("grobTree", polys))
  },
  draw_key = draw_key_polygon
)


point_to_rectriangle <- function(x, y, r, type = type) {
  r <- 0.5 * sign(r) * sqrt(abs(r))
  #r0 = 0.5
  xmin <- - r + x
  xmax <- r + x
  ymin <- - r + y
  ymax <- r + y 
  if(type == "upper"){
    df = new_data_frame(list(
      y = c(ymax, ymax, ymin, ymax),
      x = c(xmin, xmax, xmin, xmin)
    ))
  }else if(type == "lower"){
    df = new_data_frame(list(
      y = c(ymax, ymin, ymin, ymax),
      x = c(xmax, xmax, xmin, xmax) 
    ))
  }
  df
}

data$B=factor(data$B,levels = c('E2F6','TOP2A','PredictScore'))
data$A = factor(data$A,levels = c('Apoptosis','Necroptosis','Programmed Cell Death','Immunogenic cell death','Cell stemness',
                                  'Epithelial-mesenchymal transition','Invasion','Hypoxia','Fatty acid metabolism',
                                  'm5C','m6A','m7G',
                                  'Adipocytokine signaling pathway', 'AMPK signaling pathway','ErbB signaling pathway',
                                  'Glucagon signaling pathway','Hippo signaling pathway','p53 signaling pathway',
                                  'PPAR signaling pathway','Signaling pathways regulating pluripotency of stem cells','Sphingolipid signaling pathway',
                                  'VEGF signaling pathway'))
table(data$A)
table(data$B)
ggplot()+
  geom_rectriangle(data = data, aes(A, B, fill = -log10(p)),type = "upper", r = 1)+
  scale_fill_gradient(low = "white", high = "#a5d3ed",na.value = "white")+
  ggnewscale::new_scale_fill()+
  geom_rectriangle(data = data, aes(A, B, fill = r),type = "lower", r = 1)+
  scale_fill_gradient2(high = "#e7857b", mid = "white",low = "#6d9cc5")+
  labs(x = "", y = "")+
  geom_text(data = data,aes(A,B,label = ps),
            nudge_y = 0.05)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1))
ggsave('figure09_tupianzhengli/icgc_ssgsea_correlation_hubgene_predict_phenotype_pathway_deg.pdf',height = 6,width =8,dpi = 320,units = 'in')


# #GSVA####
# rm(list = ls())
# load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
# load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
# load('figure08_zhengligeneset/icgc_phenotype_gsva.Rdata')
# re_phenotype_ssgsea=re
# load('figure08_zhengligeneset/icgc_pathway_gsva.Rdata')
# re_pathway_ssgsea = re
# identical(colnames(re_phenotype_ssgsea),colnames(re_pathway_ssgsea))
# re = rbind(re_phenotype_ssgsea,re_pathway_ssgsea)
# load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
# up_phenotype=up
# down_phenotype = down
# cg_phenotype = c(up_phenotype,down_phenotype)
# load('figure09_tongluzaifenxi/pathway_deg.Rdata')
# up_pathway = up
# down_pathway = down
# cg_pathway = c(up_pathway,down_pathway)
# re=re[c(cg_phenotype,cg_pathway),]
# identical(rownames(meta),colnames(re))
# exp=exprSet
# exp=exp[,match(colnames(re),colnames(exp))]
# identical(colnames(re),colnames(exp))
# exp=exp[c('E2F6','TOP2A'),]
# identical(colnames(re),colnames(exp))
# re=rbind(re,exp)
# identical(colnames(exp),rownames(meta))
# identical(colnames(re),rownames(meta))
# data=as.data.frame(t(re))
# data$PredictScore=meta$Predict_score
# M=data
# 
# library(tidyverse)
# library(corrplot)
# if(!require(ggcor))devtools::install_github("xukaili/ggcor",upgrade = F)
# datp <- cor.mtest(M)$p[1:15,16:18] %>%
#   as.data.frame()%>% 
#   rownames_to_column("A") %>% 
#   gather("B","p",-A)
# data <- cor(M)[1:15,16:18] %>%
#   as.data.frame()%>% 
#   rownames_to_column("A") %>% 
#   gather("B","r",-A) %>% 
#   mutate(p = datp$p)
# data$ps = case_when(data$p<0.01~"**",
#                     data$p<0.05~"*",
#                     T~"")
# library(vctrs)
# library(grid)
# 
# geom_rectriangle <- function(mapping = NULL, data = NULL,
#                              stat = "identity", position = "identity",
#                              ...,
#                              linejoin = "mitre",
#                              na.rm = FALSE,
#                              show.legend = NA,
#                              inherit.aes = TRUE) {
#   layer(
#     data = data,
#     mapping = mapping,
#     stat = stat,
#     geom = GeomRectriangle,
#     position = position,
#     show.legend = show.legend,
#     inherit.aes = inherit.aes,
#     params = list(
#       linejoin = linejoin,
#       na.rm = na.rm,
#       ...
#     )
#   )
# }
# 
# GeomRectriangle <- ggproto(
#   "GeomRectriangle", Geom,
#   default_aes = aes(r = 1, colour = "grey35", fill = NA, size = 0.25, linetype = 1,
#                     alpha = NA,type = "upper"),
#   required_aes = c("x", "y"),
#   draw_panel = function(self, data, panel_params, coord, linejoin = "mitre",type = "upper") {
#     aesthetics <- setdiff(names(data), c("x", "y"))
#     
#     polys <- lapply(split(data, seq_len(nrow(data))), function(row) {
#       rectriangle <- point_to_rectriangle(row$x, row$y, row$r, row$type)
#       aes <- new_data_frame(row[aesthetics])[rep(1, 4), ]
#       GeomPolygon$draw_panel(cbind(rectriangle, aes), panel_params, coord)
#     })
#     
#     ggplot2:::ggname("geom_rectriangle", do.call("grobTree", polys))
#   },
#   draw_key = draw_key_polygon
# )
# 
# 
# point_to_rectriangle <- function(x, y, r, type = type) {
#   r <- 0.5 * sign(r) * sqrt(abs(r))
#   #r0 = 0.5
#   xmin <- - r + x
#   xmax <- r + x
#   ymin <- - r + y
#   ymax <- r + y 
#   if(type == "upper"){
#     df = new_data_frame(list(
#       y = c(ymax, ymax, ymin, ymax),
#       x = c(xmin, xmax, xmin, xmin)
#     ))
#   }else if(type == "lower"){
#     df = new_data_frame(list(
#       y = c(ymax, ymin, ymin, ymax),
#       x = c(xmax, xmax, xmin, xmax) 
#     ))
#   }
#   df
# }
# 
# data$B=factor(data$B,levels = c('E2F6','TOP2A','PredictScore'))
# data$A = factor(data$A,levels = c('Disulfidptosis','Hypoxia',
#                                   'Cell stemness','Epithelial-mesenchymal transition','Invasion',
#                                   'm5C','m6A','m7G',
#                                   'Adipocytokine signaling pathway', 'ErbB signaling pathway',
#                                   'Hippo signaling pathway','p53 signaling pathway',
#                                   'PPAR signaling pathway','Sphingolipid signaling pathway',
#                                   'VEGF signaling pathway'))
# table(data$A)
# table(data$B)
# ggplot()+
#   geom_rectriangle(data = data, aes(A, B, fill = -log10(p)),type = "upper", r = 1)+
#   scale_fill_gradient(low = "white", high = "#66a3b1",na.value = "white")+
#   ggnewscale::new_scale_fill()+
#   geom_rectriangle(data = data, aes(A, B, fill = r),type = "lower", r = 1)+
#   scale_fill_gradient2(high = "red", mid = "white",low = "blue")+
#   labs(x = "", y = "")+
#   geom_text(data = data,aes(A,B,label = ps),
#             nudge_y = 0.05)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1))
# ggsave('figure09_tupianzhengli/icgc_gsva_correlation_hubgene_predict_phenotype_pathway_deg.pdf',height = 6,width =8,dpi = 320,units = 'in')

# TCGA EMT ssgsea  --------------------------------------------------------


