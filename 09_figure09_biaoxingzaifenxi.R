##1.1 tcga表型的差异分析ssgsea####
rm(list = ls())
load('figure08_zhengligeneset/tcga_phenotype_ssgsea.Rdata')
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
identical(colnames(re),rownames(meta))
library("limma")
library('dplyr')

Group=meta$Group
table(Group)
class(Group)
exp=re
class(exp)
identical(colnames(exp),rownames(meta))

design=model.matrix(~Group)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

#添加change列标记基因上调下调
logFC_t = 0
P.Value_t = 0.05
library(dplyr)
k1 = (deg$adj.P.Val < P.Value_t)&(deg$logFC < -logFC_t);table(k1)
k2 = (deg$adj.P.Val < P.Value_t)&(deg$logFC > logFC_t);table(k2)

deg$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(deg$change)
cg = rownames(deg)[deg$change !="NOT"]
library('ggbio')
pca.plot = draw_pca(exp[cg,],Group,addEllipses = F);pca.plot
library(ggplot2)
library(hexbin)
ggsave('figure09_biaoxingtongluzaifenxi/tcga_phenotype_deg_pca_ssgsea.pdf',height = 4,width = 5,units='in')

table(Group)
Other_name=rownames(meta)[meta$Group=='Others']
ET_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']
exp_Other=exp[,Other_name]
exp_ET=exp[,ET_name]
identical(rownames(exp_Other),rownames(exp_ET))
exp=cbind(exp_Other,exp_ET)
meta=meta[match(colnames(exp),rownames(meta)),]
identical(colnames(exp),rownames(meta))
Group=meta$Group
class(Group)
h = draw_heatmap(exp[cg,],Group,n_cutoff = 2,cluster_cols = F,show_rownames = T,annotation_legend = T,legend = T )
h
ggsave('figure09_biaoxingtongluzaifenxi/tcga_phenotype_deg_heatmap_ssgsea.pdf',height = 4,width = 6,units = 'in')
v = draw_volcano(deg,logFC_cutoff = logFC_t,pkg = 4,symmetry = T,adjust = T)
v
ggsave('figure09_biaoxingtongluzaifenxi/tcga_phenotype_deg_volcanno_ssgsea.pdf',height = 4,width = 5,units = 'in')
up=rownames(deg)[deg$change=='UP']
down = rownames(deg)[deg$change == 'DOWN']
save(deg,cg,up,down,file = 'figure09_biaoxingtongluzaifenxi/tcga_phenotype_deg_ssgsea.Rdata')

##1.2 tcga表型的差异分析gsva####
rm(list = ls())
load('figure08_zhengligeneset/tcga_phenotype_gsva.Rdata')
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
identical(colnames(re),rownames(meta))
library("limma")
library('dplyr')

Group=meta$Group
class(Group)
exp=re
class(exp)
identical(colnames(exp),rownames(meta))

design=model.matrix(~Group)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

#添加change列标记基因上调下调
logFC_t = 0
P.Value_t = 0.05
library(dplyr)
k1 = (deg$adj.P.Val < P.Value_t)&(deg$logFC < -logFC_t);table(k1)
k2 = (deg$adj.P.Val < P.Value_t)&(deg$logFC > logFC_t);table(k2)

deg$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(deg$change)
cg = rownames(deg)[deg$change !="NOT"]
library('ggbio')
pca.plot = draw_pca(exp[cg,],Group,addEllipses = F);pca.plot
library(ggplot2)
library(hexbin)
ggsave('figure09_biaoxingtongluzaifenxi/tcga_phenotype_deg_pca_gsva.pdf',height = 4,width = 5,units='in')

table(Group)
Other_name=rownames(meta)[meta$Group=='Others']
ET_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']
exp_Other=exp[,Other_name]
exp_ET=exp[,ET_name]
identical(rownames(exp_Other),rownames(exp_ET))
exp=cbind(exp_Other,exp_ET)
meta=meta[match(colnames(exp),rownames(meta)),]
identical(colnames(exp),rownames(meta))
Group=meta$Group
class(Group)
h = draw_heatmap(exp[cg,],Group,n_cutoff = 2,cluster_cols = F,show_rownames = T,annotation_legend = T,legend = T )
h
ggsave('figure09_biaoxingtongluzaifenxi/tcga_phenotype_deg_heatmap_gsva.pdf',height = 4,width = 6,units = 'in')
v = draw_volcano(deg,logFC_cutoff = logFC_t,pkg = 4,symmetry = T,adjust = T)
v
ggsave('figure09_biaoxingtongluzaifenxi/tcga_phenotype_deg_volcanno_gsva.pdf',height = 4,width = 5,units = 'in')
up=rownames(deg)[deg$change=='UP']
down = rownames(deg)[deg$change == 'DOWN']
save(deg,cg,up,down,file = 'figure09_biaoxingtongluzaifenxi/tcga_phenotype_deg_gsva.Rdata')

##1.3 icgc表型的差异分析ssgsea####
rm(list = ls())
load('figure08_zhengligeneset/icgc_phenotype_ssgsea.Rdata')
load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
identical(colnames(re),rownames(meta))
library("limma")
library('dplyr')

Group=meta$Group
class(Group)
exp=re
class(exp)
identical(colnames(exp),rownames(meta))

design=model.matrix(~Group)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

#添加change列标记基因上调下调
logFC_t = 0
P.Value_t = 0.05
library(dplyr)
k1 = (deg$adj.P.Val < P.Value_t)&(deg$logFC < -logFC_t);table(k1)
k2 = (deg$adj.P.Val < P.Value_t)&(deg$logFC > logFC_t);table(k2)

deg$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(deg$change)
cg = rownames(deg)[deg$change !="NOT"]
library('ggbio')
pca.plot = draw_pca(exp[cg,],Group,addEllipses = F);pca.plot
library(ggplot2)
library(hexbin)
ggsave('figure09_biaoxingtongluzaifenxi/icgc_phenotype_deg_pca_ssgsea.pdf',height = 4,width = 5,units='in')

table(Group)
Other_name=rownames(meta)[meta$Group=='Others']
ET_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']
exp_Other=exp[,Other_name]
exp_ET=exp[,ET_name]
identical(rownames(exp_Other),rownames(exp_ET))
exp=cbind(exp_Other,exp_ET)
meta=meta[match(colnames(exp),rownames(meta)),]
identical(colnames(exp),rownames(meta))
Group=meta$Group
class(Group)
h = draw_heatmap(exp[cg,],Group,n_cutoff = 2,cluster_cols = F,show_rownames = T,annotation_legend = T,legend = T )
h
ggsave('figure09_biaoxingtongluzaifenxi/icgc_phenotype_deg_heatmap_ssgsea.pdf',height = 4,width = 6,units = 'in')
v = draw_volcano(deg,logFC_cutoff = logFC_t,pkg = 4,symmetry = T,adjust = T)
v
ggsave('figure09_biaoxingtongluzaifenxi/icgc_phenotype_deg_volcanno_ssgsea.pdf',height = 4,width = 5,units = 'in')
up=rownames(deg)[deg$change=='UP']
down = rownames(deg)[deg$change == 'DOWN']
save(deg,cg,up,down,file = 'figure09_biaoxingtongluzaifenxi/icgc_phenotype_deg_ssgsea.Rdata')

##1.4 icgc表型的差异分析gsva####
rm(list = ls())
load('figure08_zhengligeneset/icgc_phenotype_gsva.Rdata')
load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
identical(colnames(re),rownames(meta))
library("limma")
library('dplyr')

Group=meta$Group
class(Group)
exp=re
class(exp)
identical(colnames(exp),rownames(meta))

design=model.matrix(~Group)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)

#添加change列标记基因上调下调
logFC_t = 0
P.Value_t = 0.05
library(dplyr)
k1 = (deg$adj.P.Val < P.Value_t)&(deg$logFC < -logFC_t);table(k1)
k2 = (deg$adj.P.Val < P.Value_t)&(deg$logFC > logFC_t);table(k2)

deg$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(deg$change)
cg = rownames(deg)[deg$change !="NOT"]
library('ggbio')
pca.plot = draw_pca(exp[cg,],Group,addEllipses = F);pca.plot
library(ggplot2)
library(hexbin)
ggsave('figure09_biaoxingtongluzaifenxi/icgc_phenotype_deg_pca_gsva.pdf',height = 4,width = 5,units='in')

table(Group)
Other_name=rownames(meta)[meta$Group=='Others']
ET_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']
exp_Other=exp[,Other_name]
exp_ET=exp[,ET_name]
identical(rownames(exp_Other),rownames(exp_ET))
exp=cbind(exp_Other,exp_ET)
meta=meta[match(colnames(exp),rownames(meta)),]
identical(colnames(exp),rownames(meta))
Group=meta$Group
class(Group)
h = draw_heatmap(exp[cg,],Group,n_cutoff = 2,cluster_cols = F,show_rownames = T,annotation_legend = T,legend = T )
h
ggsave('figure09_biaoxingtongluzaifenxi/icgc_phenotype_deg_heatmap_gsva.pdf',height = 4,width = 6,units = 'in')
v = draw_volcano(deg,logFC_cutoff = logFC_t,pkg = 4,symmetry = T,adjust = T)
v
ggsave('figure09_biaoxingtongluzaifenxi/icgc_phenotype_deg_volcanno_gsva.pdf',height = 4,width = 5,units = 'in')
up=rownames(deg)[deg$change=='UP']
down = rownames(deg)[deg$change == 'DOWN']
save(deg,cg,up,down,file = 'figure09_biaoxingtongluzaifenxi/icgc_phenotype_deg_gsva.Rdata')

#取差异交集绘制韦恩图，热图和箱线图
#绘制韦恩图
rm(list = ls())
load('figure09_biaoxingtongluzaifenxi/tcga_phenotype_deg_ssgsea.Rdata')
write.table(up,file = 'figure09_biaoxingtongluzaifenxi/tcga_phenotype_up_ssgsea.txt',row.names = F,quote = F,col.names = F)
write.table(down,file = 'figure09_biaoxingtongluzaifenxi/tcga_phenotype_down_ssgsea.txt',row.names = F,quote = F,col.names = F)
rm(list = ls())
load('figure09_biaoxingtongluzaifenxi/tcga_phenotype_deg_gsva.Rdata')
write.table(up,file = 'figure09_biaoxingtongluzaifenxi/tcga_phenotype_up_gsva.txt',row.names = F,quote = F,col.names = F)
write.table(down,file = 'figure09_biaoxingtongluzaifenxi/tcga_phenotype_down_gsva.txt',row.names = F,quote = F,col.names = F)
rm(list = ls())
load('figure09_biaoxingtongluzaifenxi/icgc_phenotype_deg_ssgsea.Rdata')
write.table(up,file = 'figure09_biaoxingtongluzaifenxi/icgc_phenotype_up_ssgsea.txt',row.names = F,quote = F,col.names = F)
write.table(down,file = 'figure09_biaoxingtongluzaifenxi/icgc_phenotype_down_ssgsea.txt',row.names = F,quote = F,col.names = F)
rm(list = ls())
load('figure09_biaoxingtongluzaifenxi/icgc_phenotype_deg_gsva.Rdata')
write.table(up,file = 'figure09_biaoxingtongluzaifenxi/icgc_phenotype_up_gsva.txt',row.names = F,quote = F,col.names = F)
write.table(down,file = 'figure09_biaoxingtongluzaifenxi/icgc_phenotype_down_gsva.txt',row.names = F,quote = F,col.names = F)

##绘制热图####
rm(list = ls())
load('figure09_biaoxingtongluzaifenxi/tcga_phenotype_deg_ssgsea.Rdata')
up_tcga_ss = up
down_tcga_ss = down
cg_tcga_ss = cg
# load('figure09_biaoxingtongluzaifenxi/tcga_phenotype_deg_gsva.Rdata')
# up_tcga_gsva = up
# down_tcga_gsva = down
# cg_tcga_gsva = cg
load('figure09_biaoxingtongluzaifenxi/icgc_phenotype_deg_ssgsea.Rdata')
up_icgc_ss = up
down_icgc_ss = down
cg_icgc_ss = cg
# load('figure09_biaoxingtongluzaifenxi/icgc_phenotype_deg_gsva.Rdata')
# up_icgc_gsva = up
# down_icgc_gsva = down
# cg_icgc_gsva = cg

up = intersect(up_tcga_ss,up_icgc_ss)
down = intersect(down_tcga_ss,down_icgc_ss)
cg =intersect(cg_tcga_ss,cg_icgc_ss)
save(up,down,file = 'figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')


#1.1tcga表型ssgsea####
rm(list = ls())
load('figure08_zhengligeneset/tcga_phenotype_ssgsea.Rdata')
load('figure02_E2F6/lihc_sur_data.Rdata')
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
re = re[c(up,down),]
library(pheatmap)
Cluster1_name=rownames(meta)[meta$Cluster=='Cluster1']
Cluster2_name=rownames(meta)[meta$Cluster=='Cluster2']
table(meta$Group)
Group1_name=rownames(meta)[meta$Group=='Others']
Group2_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']

re1=re[,Cluster1_name]
re2=re[,Cluster2_name]
re=cbind(re1,re2)
meta=meta[match(colnames(re),rownames(meta)),]
identical(rownames(meta),colnames(re))
re2=as.data.frame(re)
an = data.frame(group = meta$Cluster,
                row.names = colnames(re2))
dev.off()
re2=as.matrix(re2)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf('figure09_biaoxingtongluzaifenxi/tcga_ssgsea_phenotype_heatmap_cluster_deg.pdf',height = 4,width = 6)
pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = F,column_split=an,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()

identical(colnames(re),rownames(meta))
Group=meta$Cluster
draw_boxplot(re,Group,xlab = ' ',ylab = 'Score',sort = F)
ggsave('figure09_biaoxingtongluzaifenxi/tcga_ssgsea_phenotype_boxplot_cluster_deg.pdf',height = 6,width =5,units = 'in')

Group1_name=rownames(meta)[meta$Group=='Others']
Group2_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']
re1=re[,Group1_name]
re2=re[,Group2_name]
re=cbind(re1,re2)
meta=meta[match(colnames(re),rownames(meta)),]
identical(rownames(meta),colnames(re))
re2=as.data.frame(re)
an = data.frame(group = meta$Group,
                row.names = colnames(re2))
dev.off()
re2=as.matrix(re2)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf('figure09_biaoxingtongluzaifenxi/tcga_ssgsea_phenotype_heatmap_deg.pdf',height = 4,width = 8)
pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = F,column_split=an,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()

identical(colnames(re),rownames(meta))
Group=meta$Group
draw_boxplot(re,Group,xlab = ' ',ylab = 'Score',sort = F)
ggsave('figure09_biaoxingtongluzaifenxi/tcga_ssgsea_phenotype_boxplot_deg.pdf',height = 6,width = 5,units = 'in')

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
datp <- cor.mtest(M)$p[1:12,13:15] %>%
  as.data.frame()%>% 
  rownames_to_column("A") %>% 
  gather("B","p",-A)
data <- cor(M)[1:12,13:15] %>%
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
ggsave('figure09_biaoxingtongluzaifenxi/tcga_ssgsea_correlation_hubgene_predict_phenotype_deg.pdf',height = 5,width =6,dpi = 320,units = 'in')

#1.2tcga表型gsva####
rm(list = ls())
load('figure08_zhengligeneset/tcga_phenotype_gsva.Rdata')
load('figure02_E2F6/lihc_sur_data.Rdata')
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
re = re[c(up,down),]
library(pheatmap)
Cluster1_name=rownames(meta)[meta$Cluster=='Cluster1']
Cluster2_name=rownames(meta)[meta$Cluster=='Cluster2']
table(meta$Group)
Group1_name=rownames(meta)[meta$Group=='Others']
Group2_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']

re1=re[,Cluster1_name]
re2=re[,Cluster2_name]
re=cbind(re1,re2)
meta=meta[match(colnames(re),rownames(meta)),]
identical(rownames(meta),colnames(re))
re2=as.data.frame(re)
an = data.frame(group = meta$Cluster,
                row.names = colnames(re2))
dev.off()
re2=as.matrix(re2)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf('figure09_biaoxingtongluzaifenxi/tcga_gsva_phenotype_heatmap_cluster_deg.pdf',height = 4,width = 6)
pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = F,column_split=an,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()

identical(colnames(re),rownames(meta))
Group=meta$Cluster
draw_boxplot(re,Group,xlab = ' ',ylab = 'Score',sort = F)
ggsave('figure09_biaoxingtongluzaifenxi/tcga_gsva_phenotype_boxplot_cluster_deg.pdf',height = 6,width =5,units = 'in')

Group1_name=rownames(meta)[meta$Group=='Others']
Group2_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']
re1=re[,Group1_name]
re2=re[,Group2_name]
re=cbind(re1,re2)
meta=meta[match(colnames(re),rownames(meta)),]
identical(rownames(meta),colnames(re))
re2=as.data.frame(re)
an = data.frame(group = meta$Group,
                row.names = colnames(re2))
dev.off()
re2=as.matrix(re2)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf('figure09_biaoxingtongluzaifenxi/tcga_gsva_phenotype_heatmap_deg.pdf',height = 4,width = 8)
pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = F,column_split=an,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()

identical(colnames(re),rownames(meta))
Group=meta$Group
draw_boxplot(re,Group,xlab = ' ',ylab = 'Score',sort = F)
ggsave('figure09_biaoxingtongluzaifenxi/tcga_gsva_phenotype_boxplot_deg.pdf',height = 6,width = 5,units = 'in')

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
datp <- cor.mtest(M)$p[1:12,13:15] %>%
  as.data.frame()%>% 
  rownames_to_column("A") %>% 
  gather("B","p",-A)
data <- cor(M)[1:12,13:15] %>%
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
ggsave('figure09_biaoxingtongluzaifenxi/tcga_gsva_correlation_hubgene_predict_phenotype_deg.pdf',height = 5,width =6,dpi = 320,units = 'in')

#1.3icgc表型ssgsea####
rm(list = ls())
load('figure08_zhengligeneset/icgc_phenotype_ssgsea.Rdata')
load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
re = re[c(up,down),]
library(pheatmap)
Cluster1_name=rownames(meta)[meta$Cluster=='Cluster1']
Cluster2_name=rownames(meta)[meta$Cluster=='Cluster2']
table(meta$Group)
Group1_name=rownames(meta)[meta$Group=='Others']
Group2_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']

re1=re[,Cluster1_name]
re2=re[,Cluster2_name]
re=cbind(re1,re2)
meta=meta[match(colnames(re),rownames(meta)),]
identical(rownames(meta),colnames(re))
re2=as.data.frame(re)
an = data.frame(group = meta$Cluster,
                row.names = colnames(re2))
dev.off()
re2=as.matrix(re2)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf('figure09_biaoxingtongluzaifenxi/icgc_ssgsea_phenotype_heatmap_cluster_deg.pdf',height = 4,width = 6)
pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = F,column_split=an,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()

identical(colnames(re),rownames(meta))
Group=meta$Cluster
draw_boxplot(re,Group,xlab = ' ',ylab = 'Score',sort = F)
ggsave('figure09_biaoxingtongluzaifenxi/icgc_ssgsea_phenotype_boxplot_cluster_deg.pdf',height = 6,width =5,units = 'in')

Group1_name=rownames(meta)[meta$Group=='Others']
Group2_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']
re1=re[,Group1_name]
re2=re[,Group2_name]
re=cbind(re1,re2)
meta=meta[match(colnames(re),rownames(meta)),]
identical(rownames(meta),colnames(re))
re2=as.data.frame(re)
an = data.frame(group = meta$Group,
                row.names = colnames(re2))
dev.off()
re2=as.matrix(re2)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf('figure09_biaoxingtongluzaifenxi/icgc_ssgsea_phenotype_heatmap_deg.pdf',height = 4,width = 8)
pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = F,column_split=an,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()

identical(colnames(re),rownames(meta))
Group=meta$Group
draw_boxplot(re,Group,xlab = ' ',ylab = 'Score',sort = F)
ggsave('figure09_biaoxingtongluzaifenxi/icgc_ssgsea_phenotype_boxplot_deg.pdf',height = 6,width = 5,units = 'in')

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
datp <- cor.mtest(M)$p[1:12,13:15] %>%
  as.data.frame()%>% 
  rownames_to_column("A") %>% 
  gather("B","p",-A)
data <- cor(M)[1:12,13:15] %>%
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
ggsave('figure09_biaoxingtongluzaifenxi/icgc_ssgsea_correlation_hubgene_predict_phenotype_deg.pdf',height = 5,width =6,dpi = 320,units = 'in')

#1.4icgc表型gsva####
rm(list = ls())
load('figure08_zhengligeneset/icgc_phenotype_gsva.Rdata')
load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
load('figure09_biaoxingtongluzaifenxi/phenotype_deg.Rdata')
re = re[c(up,down),]
library(pheatmap)
Cluster1_name=rownames(meta)[meta$Cluster=='Cluster1']
Cluster2_name=rownames(meta)[meta$Cluster=='Cluster2']
table(meta$Group)
Group1_name=rownames(meta)[meta$Group=='Others']
Group2_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']

re1=re[,Cluster1_name]
re2=re[,Cluster2_name]
re=cbind(re1,re2)
meta=meta[match(colnames(re),rownames(meta)),]
identical(rownames(meta),colnames(re))
re2=as.data.frame(re)
an = data.frame(group = meta$Cluster,
                row.names = colnames(re2))
dev.off()
re2=as.matrix(re2)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf('figure09_biaoxingtongluzaifenxi/icgc_gsva_phenotype_heatmap_cluster_deg.pdf',height = 4,width = 6)
pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = F,column_split=an,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()

identical(colnames(re),rownames(meta))
Group=meta$Cluster
draw_boxplot(re,Group,xlab = ' ',ylab = 'Score',sort = F)
ggsave('figure09_biaoxingtongluzaifenxi/icgc_gsva_phenotype_boxplot_cluster_deg.pdf',height = 6,width =5,units = 'in')

Group1_name=rownames(meta)[meta$Group=='Others']
Group2_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']
re1=re[,Group1_name]
re2=re[,Group2_name]
re=cbind(re1,re2)
meta=meta[match(colnames(re),rownames(meta)),]
identical(rownames(meta),colnames(re))
re2=as.data.frame(re)
an = data.frame(group = meta$Group,
                row.names = colnames(re2))
dev.off()
re2=as.matrix(re2)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf('figure09_biaoxingtongluzaifenxi/icgc_gsva_phenotype_heatmap_deg.pdf',height = 4,width = 8)
pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = F,column_split=an,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()

identical(colnames(re),rownames(meta))
Group=meta$Group
draw_boxplot(re,Group,xlab = ' ',ylab = 'Score',sort = F)
ggsave('figure09_biaoxingtongluzaifenxi/icgc_gsva_phenotype_boxplot_deg.pdf',height = 6,width = 5,units = 'in')

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
datp <- cor.mtest(M)$p[1:12,13:15] %>%
  as.data.frame()%>% 
  rownames_to_column("A") %>% 
  gather("B","p",-A)
data <- cor(M)[1:12,13:15] %>%
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
ggsave('figure09_biaoxingtongluzaifenxi/icgc_gsva_correlation_hubgene_predict_phenotype_deg.pdf',height = 5,width =6,dpi = 320,units = 'in')

