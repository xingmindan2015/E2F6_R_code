
# TCGA-heatmap_immunecell_pathway ------------------------------------------------------------

rm(list = ls())
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
load('figure07_immune/tcga_immu_cell16_ssgsaea.Rdata')
re_cell = re
load('figure07_immune/tcga_immu_pathway13_ssgsaea.Rdata')
re_pathway = re
identical(colnames(re_cell),colnames(re_pathway))
re = rbind(re_cell,re_pathway)
table(meta$Group)
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
col_legend = ifelse(rownames(re)%in%rownames(re_cell),'Immune_cell','Immune_pathway')
names(col_legend)=rownames(re)
table(col_legend)
an_row = data.frame(col_legend,row.names = rownames(re))
dev.off()
class(re)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
table(meta$Group)
legend_col = list(Group = c('Others'='#6d9cc5','E2F6_H&TOP2A_H'='#e7857b'),
              Cluster=c('Cluster1'='#a5d3ed','Cluster2'='#b5add4'),
              col_legend=c('Immune_cell'='#dddb97','Immune_pathway'='#b4c5e6'))
pdf('figure07_immune/tcga_ssgsea_immu_cell_pathway_heatmap_cluster_group.pdf',height = 5,width = 8)
pheatmap(re,scale = "row",
         show_colnames = F,
         annotation_col = an_col,
         #annotation_row = an_row,
         cluster_cols = F,
         cluster_rows = F,
         column_split=an,
         row_split = an_row,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         annotation_colors = legend_col)

dev.off()

# TCGA_e2f6_top2a_immune_cell_pathway_correlation -----------------------------

rm(list = ls())
load('figure02_E2F6/lihc_sur_data.Rdata')
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
load('figure07_immune/tcga_immu_cell16_ssgsaea.Rdata')
re_cell = re
load('figure07_immune/tcga_immu_pathway13_ssgsaea.Rdata')
re_pathway = re
identical(colnames(re_cell),colnames(re_pathway))
re = rbind(re_cell,re_pathway)

exprSet[1:4,1:4]
exp=exprSet
identical(rownames(meta),colnames(re))
exp=exp[,match(colnames(re),colnames(exp))]
identical(colnames(re),colnames(exp))
exp=exp[c('E2F6','TOP2A'),]
identical(colnames(re),colnames(exp))
re=rbind(re,exp)
identical(colnames(exp),rownames(meta))
identical(colnames(re),rownames(meta))
data=as.data.frame(t(re))
data$PredictScore=meta$Predict_score
# library(corrplot)
# library(corrplot)
# M = cor(as.matrix(data))
# testRes = cor.mtest(as.matrix(data), conf.level = 0.95)$p
# library(tidyverse)
# g = pivot_longer(rownames_to_column(as.data.frame(M),var = "from"),
#                  cols = 2:(ncol(M)+1),
#                  names_to = "to",
#                  values_to = "cor")
# gp = pivot_longer(rownames_to_column(as.data.frame(testRes)),
#                   cols = 2:(ncol(M)+1),
#                   names_to = "gene",
#                   values_to = "p")
# g$p = gp$p
# g = g[g$from!=g$to,]
# g$group = case_when(g$cor>0.3 & g$p<0.05 ~ "positive",
#                     g$cor< -0.3 & g$p<0.05 ~ "negative",
#                     T~"not" )
# head(g)
# write.csv(g,file = 'figure07_immune/tcga_ssgsea_e2f6_top2a_predictsocre_relationship.csv')
M=data

library(tidyverse)
library(corrplot)
if(!require(ggcor))devtools::install_github("xukaili/ggcor",upgrade = F)
datp <- cor.mtest(M)$p[1:29,30:32] %>%
  as.data.frame()%>% 
  rownames_to_column("A") %>% 
  gather("B","p",-A)
data <- cor(M)[1:29,30:32] %>%
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
  theme(axis.text.x = element_text(vjust = 1,hjust = 1,angle = 80))
ggsave('figure07_immune/tcga_ssgsea_correlation_hubgene_predict_immunecell_pathway.pdf',height = 5,width =10,dpi = 320,units = 'in')

# TCGA_e2f6_top2a_immune_cell_pathway_network -----------------------------
rm(list = ls())
load('figure02_E2F6/lihc_sur_data.Rdata')
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
load('figure07_immune/tcga_immu_cell16_ssgsaea.Rdata')
re_cell = re
load('figure07_immune/tcga_immu_pathway13_ssgsaea.Rdata')
re_pathway = re
identical(colnames(re_cell),colnames(re_pathway))
re = rbind(re_cell,re_pathway)

exprSet[1:4,1:4]
exp=exprSet
identical(rownames(meta),colnames(re))
exp=exp[,match(colnames(re),colnames(exp))]
identical(colnames(re),colnames(exp))
exp=exp[c('E2F6','TOP2A'),]
identical(colnames(exp),rownames(meta))
data=as.data.frame(t(exp))
data$PredictScore=meta$Predict_score
# install.packages("ggplot2")
# install.packages("vegan")
# install.packages("dplyr")
# install.packages("devtools")
# devtools::install_github("houyunhuang/ggcor")
#加载包
library(vegan)
library(dplyr)
library(ggcor)
library(ggplot2)

quickcor(t(re), type = "lower",method = "spearman") +
  geom_square()+
  scale_fill_gradient2( high = 'orange', mid = 'white',low = 'navyblue')
df_mantel <- mantel_test(data, t(re), mantel.fun = 'mantel',
                         spec.dist.method = 'bray', 
                         env.dist.method = 'euclidean',
                         spec.select = list(E2F6 = 1,
                                            TOP2A = 2,
                                            PredictScore = 3))
df_mantel <- df_mantel %>%
  mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.2, 0.4, Inf),
                    labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.4", ">= 0.4")),#定义Mantel的R值范围标签
         df_p = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
#devtools::install_github("Hy4m/linkET", force = TRUE)
#library(linkET)
#remove.packages('linkET')

# qcorrplot(correlate(t(re)),
#           type = "lower", # 热图展示下半部分
#           diag = FALSE) + # 不展示对角线
#   geom_square() +
#   geom_couple(aes(colour = df_p, # 网络线的颜色映射mantel_test结果的相关性
#                   size=r), # 网络线的粗细映射mantel_test结果的显著性
#               data = df_mantel,
#               curvature = nice_curvature(by = 'to'),) # 基于起点或终点绘制最合适的曲线
quickcor(t(re),method = "pearson", type = "upper", cor.test = T, cluster.type = "all") +#环境因子之间的相关性热图
  geom_square() +#相关性显示形式
  geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+#显著性标签
  scale_fill_gradient2( high = 'orange', mid = 'white',low = 'navyblue') + #颜色设置
  anno_link(df_mantel, aes(color = df_p,
                           size = df_r))+
  scale_size_manual(values = c(0.5, 1, 1.5, 2))+#连线粗细设置
  scale_color_manual(values = c("orange","#6387BB","#937abf"))+#线条颜色设置
  guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "Mantel's p", order = 3),
         linetype = "none")
ggsave('figure07_immune/tcga_correlation_immunecell_pathway_network.pdf',width = 12,height = 12)

# ICGC_heatmap_immunecell_pathway -----------------------------------------
rm(list = ls())
load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
load('figure07_immune/icgc_immu_cell16_ssgsaea.Rdata')
re_cell = re
load('figure07_immune/icgc_immu_pathway13_ssgsaea.Rdata')
re_pathway = re
identical(colnames(re_cell),colnames(re_pathway))
re = rbind(re_cell,re_pathway)
table(meta$Group)
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
col_legend = ifelse(rownames(re)%in%rownames(re_cell),'Immune_cell','Immune_pathway')
names(col_legend)=rownames(re)
table(col_legend)
an_row = data.frame(col_legend,row.names = rownames(re))
dev.off()
class(re)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
table(meta$Group)
legend_col = list(Group = c('Others'='#6d9cc5','E2F6_H&TOP2A_H'='#e7857b'),
                  Cluster=c('Cluster1'='#a5d3ed','Cluster2'='#b5aad4'),
                  col_legend=c('Immune_cell'='#dddb97','Immune_pathway'='#b4c5e6'))
pdf('figure07_immune/icgc_ssgsea_immu_cell_pathway_heatmap_cluster_group.pdf',height = 5,width = 8)
pheatmap(re,scale = "row",
         show_colnames = F,
         annotation_col = an_col,
         #annotation_row = an_row,
         cluster_cols = F,
         cluster_rows = F,
         column_split=an,
         row_split = an_row,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         annotation_colors = legend_col)

dev.off()

# ICGC_e2f6_top2a_immunecell_pathway_correlation --------------------------
rm(list = ls())
load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
load('figure07_immune/icgc_immu_cell16_ssgsaea.Rdata')
re_cell = re
load('figure07_immune/icgc_immu_pathway13_ssgsaea.Rdata')
re_pathway = re
identical(colnames(re_cell),colnames(re_pathway))
re = rbind(re_cell,re_pathway)

exprSet[1:4,1:4]
exp=exprSet
identical(rownames(meta),colnames(re))
exp=exp[,match(colnames(re),colnames(exp))]
identical(colnames(re),colnames(exp))
exp=exp[c('E2F6','TOP2A'),]
identical(colnames(re),colnames(exp))
re=rbind(re,exp)
identical(colnames(exp),rownames(meta))
identical(colnames(re),rownames(meta))
data=as.data.frame(t(re))
data$PredictScore=meta$Predict_score
# library(corrplot)
# library(corrplot)
# M = cor(as.matrix(data))
# testRes = cor.mtest(as.matrix(data), conf.level = 0.95)$p
# library(tidyverse)
# g = pivot_longer(rownames_to_column(as.data.frame(M),var = "from"),
#                  cols = 2:(ncol(M)+1),
#                  names_to = "to",
#                  values_to = "cor")
# gp = pivot_longer(rownames_to_column(as.data.frame(testRes)),
#                   cols = 2:(ncol(M)+1),
#                   names_to = "gene",
#                   values_to = "p")
# g$p = gp$p
# g = g[g$from!=g$to,]
# g$group = case_when(g$cor>0.3 & g$p<0.05 ~ "positive",
#                     g$cor< -0.3 & g$p<0.05 ~ "negative",
#                     T~"not" )
# head(g)
# write.csv(g,file = 'figure07_immune/tcga_ssgsea_e2f6_top2a_predictsocre_relationship.csv')
M=data

library(tidyverse)
library(corrplot)
if(!require(ggcor))devtools::install_github("xukaili/ggcor",upgrade = F)
datp <- cor.mtest(M)$p[1:29,30:32] %>%
  as.data.frame()%>% 
  rownames_to_column("A") %>% 
  gather("B","p",-A)
data <- cor(M)[1:29,30:32] %>%
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
  theme(axis.text.x = element_text(vjust = 1,hjust = 1,angle = 80))
ggsave('figure07_immune/icgc_ssgsea_correlation_hubgene_predict_immunecell_pathway.pdf',height = 5,width =10,dpi = 320,units = 'in')


# TCGA_checkpoint_heatmap -------------------------------------------------

rm(list = ls())
gs=read.csv('figure07_immune/immunecheckpoint.csv',header = T)
load('figure02_E2F6/lihc_sur_data.Rdata')
exprSet[1:4,1:4]
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
identical(colnames(exprSet),rownames(meta))
exp=exprSet
exp[1:4,1:4]
gs=gs$Symbol
gs=gs[-24]
table(gs%in%rownames(exp))
f = exp[gs,]
re=f
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
an = data.frame(Group = meta$Group,Cluster = meta$Cluster,
                row.names = colnames(re2))
#an$Cluster = factor(an$Cluster,levels = c('Cluster2','Cluster1'))
an = arrange(an,Group,Cluster)

re2 = re2 [,match(rownames(an),colnames(re2))]
an1 = data.frame(Group = an$Group,
                row.names = colnames(re2))
identical(rownames(an),rownames(an1))
dev.off()
re2=as.matrix(re2)
legend_col = list(Group = c('Others'='#6d9cc5','E2F6_H&TOP2A_H'='#e7857b'),
                  Cluster=c('Cluster1'='#a5d3ed','Cluster2'='#b5aad4'))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf('figure07_immune/tcga_ssgsea_immu_checkpoint_heatmap_cluster_group.pdf',height =4,width = 6)

pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = F,
         column_split=an1,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         annotation_colors = legend_col)

dev.off()


# ICGC_checkpoint_heatmap -------------------------------------------------

rm(list = ls())
gs=read.csv('figure07_immune/immunecheckpoint.csv',header = T)
load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
exprSet[1:4,1:4]
load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
identical(colnames(exprSet),rownames(meta))
exp=exprSet
exp[1:4,1:4]
gs=gs$Symbol
gs=gs[-24]
table(gs%in%rownames(exp))
f = exp[gs,]
re=f
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
an = data.frame(Group = meta$Group,Cluster=meta$Cluster,
                row.names = colnames(re2))
an = arrange(an,Group,Cluster)
re2 = re2[,match(rownames(an),colnames(re2))]
identical(rownames(an),colnames(re2))
an1 = data.frame(Group = an$Group,row.names = colnames(re2))
dev.off()
re2=as.matrix(re2)
legend_col = list(Group = c('Others'='#6d9cc5','E2F6_H&TOP2A_H'='#e7857b'),
                  Cluster=c('Cluster1'='#a5d3ed','Cluster2'='#b5aad4'))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))

pdf('figure07_immune/icgc_immu_checkpoint_heatmap_cluster_group.pdf',height = 4,width = 6)
pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = F,
         column_split=an1,
         #legend_breaks=seq(-2,2,1),
         breaks=bk,
         annotation_colors = legend_col,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

dev.off()
# ICGC_e2f6_top2a_immune_cell_pathway_network -----------------------------
rm(list = ls())
load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
load('figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')
load('figure07_immune/icgc_immu_cell16_ssgsaea.Rdata')
re_cell = re
load('figure07_immune/icgc_immu_pathway13_ssgsaea.Rdata')
re_pathway = re
identical(colnames(re_cell),colnames(re_pathway))
re = rbind(re_cell,re_pathway)

exprSet[1:4,1:4]
exp=exprSet
identical(rownames(meta),colnames(re))
exp=exp[,match(colnames(re),colnames(exp))]
identical(colnames(re),colnames(exp))
exp=exp[c('E2F6','TOP2A'),]

identical(colnames(exp),rownames(meta))
data=as.data.frame(t(exp))
data$PredictScore=meta$Predict_score
# install.packages("ggplot2")
# install.packages("vegan")
# install.packages("dplyr")
# install.packages("devtools")
# devtools::install_github("houyunhuang/ggcor")
#加载包
library(vegan)
library(dplyr)
library(ggcor)
library(ggplot2)

quickcor(t(re), type = "lower",method = "spearman") +
  geom_square()+
  scale_fill_gradient2( high = 'orange', mid = 'white',low = 'navyblue')
df_mantel <- mantel_test(data, t(re), mantel.fun = 'mantel',
                         spec.dist.method = 'bray', 
                         env.dist.method = 'euclidean',
                         spec.select = list(E2F6 = 1,
                                            TOP2A = 2,
                                            PredictScore = 3))
df_mantel <- df_mantel %>%
  mutate(df_r = cut(r, breaks = c(-Inf, 0.1, 0.2, 0.4, Inf),
                    labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.4", ">= 0.4")),#定义Mantel的R值范围标签
         df_p = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
#devtools::install_github("Hy4m/linkET", force = TRUE)
#library(linkET)
#remove.packages('linkET')

# qcorrplot(correlate(t(re)),
#           type = "lower", # 热图展示下半部分
#           diag = FALSE) + # 不展示对角线
#   geom_square() +
#   geom_couple(aes(colour = df_p, # 网络线的颜色映射mantel_test结果的相关性
#                   size=r), # 网络线的粗细映射mantel_test结果的显著性
#               data = df_mantel,
#               curvature = nice_curvature(by = 'to'),) # 基于起点或终点绘制最合适的曲线
quickcor(t(re),method = "pearson", type = "upper", cor.test = T, cluster.type = "all") +#环境因子之间的相关性热图
  geom_square() +#相关性显示形式
  geom_mark(r = NA,sig.thres = 0.05, size = 3.5, colour = "black")+#显著性标签
  scale_fill_gradient2( high = 'orange', mid = 'white',low = 'navyblue') + #颜色设置
  anno_link(df_mantel, aes(color = df_p,
                           size = df_r))+
  scale_size_manual(values = c(0.5, 1, 1.5, 2))+#连线粗细设置
  scale_color_manual(values = c("orange","#6387BB","#937abf"))+#线条颜色设置
  guides(fill = guide_colorbar(title = "correlation", order = 1),#图例相关设置
         size = guide_legend(title = "Mantel's r",order = 2),
         color = guide_legend(title = "Mantel's p", order = 3),
         linetype = "none")
ggsave('figure07_immune/icgc_correlation_immunecell_pathway_network.pdf',width = 12,height = 12)

