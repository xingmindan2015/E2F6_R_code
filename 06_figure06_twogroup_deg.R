
# E2F6_TOP2A被分成两组的预后分析 -----------------------------------------------------

rm(list = ls())
load("figure02_E2F6/lihc_sur_data.Rdata")
identical(colnames(exprSet),rownames(meta))
a=c('E2F6','TOP2A')
exp=exprSet[a,]
data=as.data.frame(t(exp))
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6_TOP2A=ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='High','E2F6_H&TOP2A_H','Others')
data$Group_E2F6_TOP2A=as.factor(data$Group_E2F6_TOP2A)
table(data$Group_E2F6_TOP2A)

identical(rownames(data),rownames(meta))
meta$Group=data$Group_E2F6_TOP2A

#E2F6的生存分析
identical(rownames(meta),rownames(data))
sfit1=survfit(Surv(Time, Event)~Group, data=meta)
g1=ggsurvplot(fit = sfit1,
              palette = c('#e7857b','#6d9cc5'),
              risk.table =TRUE,
              pval =TRUE,
              conf.int =T,
              ncensor.plot = TRUE,data = meta)
library(patchwork)
p1=g1$plot
p2=g1$table
p_sur1=p1 /p2  + plot_layout(design = 'A
                         A
                         A
                         A
                         B')
ggsave('figure06_two_group/tcga_survival_E2F6_top2a_two_group.pdf',p_sur1,width = 10,height = 6,units = 'in')
save(meta,file = 'figure06_two_group/tcga_meta_two_groups.Rdata')

rm(list = ls())
load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
identical(colnames(exprSet),rownames(meta))

a=c('E2F6','TOP2A')
exp=exprSet[a,]

data=as.data.frame(t(exp))
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6_TOP2A=ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='High','E2F6_H&TOP2A_H','Others')
data$Group_E2F6_TOP2A=as.factor(data$Group_E2F6_TOP2A)
table(data$Group_E2F6_TOP2A)
identical(rownames(data),rownames(meta))
meta$Group=data$Group_E2F6_TOP2A

#E2F6的生存分析
identical(rownames(meta),rownames(data))
colnames(meta)=c("Gender","Age","Event","Stage","Time","Group")
sfit1=survfit(Surv(Time, Event)~Group, data=meta)
g1=ggsurvplot(fit = sfit1,
              palette = c('#e7857b','#6d9cc5'),
              risk.table =TRUE,
              pval =TRUE,
              conf.int =T,
              ncensor.plot = TRUE,data = meta)
library(patchwork)
p1=g1$plot
p2=g1$table
p_sur1=p1 /p2  + plot_layout(design = 'A
                         A
                         A
                         A
                         B')
ggsave('figure06_two_group/icgc_survival_E2F6_top2a_two_group.pdf',p_sur1,width = 10,height = 6,units = 'in')
save(meta,file = 'figure06_two_group/icgc_meta_two_groups.Rdata')


# E2F6_TOP2A均升高组与其他患者间的差异分析TCGA -------------------------------------------


#高水平组与其他组间的差异分析
rm(list = ls())
getwd()
####1.TCGA数据库中E2F6在正常与肝癌的差异
load('/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/data_symbols_lihc.Rdata')
load('figure06_two_group/tcga_meta_two_groups.Rdata')
exp=count1
exp[1:4,1:4]
exp=exp[,str_sub(colnames(exp),14,15)!='02']
table(str_sub(colnames(exp),14,15))
exp=exp[,str_sub(colnames(exp),14,15)=='01']
rownames(meta)[1:10]
colnames(exp)=str_sub(colnames(exp),1,12)
table(colnames(exp)%in%rownames(meta))
exp=exp[,colnames(exp)%in%rownames(meta)]
exp=exp[,match(rownames(meta),colnames(exp))]
identical(rownames(meta),colnames(exp))
meta_o=rownames(meta)[meta$Group=='Others']
meta_H=rownames(meta)[meta$Group!='Others']
exp_o=exp[,colnames(exp)%in%meta_o]
exp_H=exp[,colnames(exp)%in%meta_H]
identical(rownames(exp_o),rownames(exp_H))
exp=cbind(exp_o,exp_H)
meta=meta[match(colnames(exp),rownames(meta)),]
identical(rownames(meta),colnames(exp))
Group=meta$Group
Group=factor(Group,levels = c('Others','E2F6_H&TOP2A_H'))
table(Group)
class(Group)

a=exp
library(DESeq2)
colData <- data.frame(row.names =colnames(exp), 
                      condition=Group)
if(!file.exists("figure06_two_group/tcga_deseqdata_two_group.Rdata")){
  dds <- DESeqDataSetFromMatrix(
    countData = exp,
    colData = colData,
    design = ~ condition)
  dds <- DESeq(dds)
  save(dds,file = "figure06_two_group/tcga_deseqdata_two_group.Rdata")
}
load("figure06_two_group/tcga_deseqdata_two_group.Rdata")
class(dds)

res <- results(dds, contrast = c("condition",rev(levels(Group))))
c("condition",rev(levels(Group)))
class(res)
DEG1 <- as.data.frame(res)
DEG1 <- DEG1[order(DEG1$pvalue),] 
DEG1 = na.omit(DEG1)
head(DEG1)

#添加Change列标记基因上调下调
logFC_t = 2
p.adj = 0.05

k1 = (DEG1$padj < p.adj)&(DEG1$log2FoldChange < -logFC_t);table(k1)
k2 = (DEG1$padj < p.adj)&(DEG1$log2FoldChange > logFC_t);table(k2)

DEG1$Change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG1$Change)
head(DEG1)
write.csv(DEG1,file = 'figure06_two_group/tcga_deseqdata_two_group.csv')
save(DEG1,file='figure06_two_group/tcga_deseqdata_DEG_two_group.Rdata')
load('figure06_two_group/tcga_deseqdata_DEG_two_group.Rdata')

library(ggplot2)
library(tinyarray)
exp[1:4,1:4]
# cpm 去除文库大小的影响

dat = log2(cpm(a)+1)
dim(dat)
cg1 = rownames(DEG1)[DEG1$Change !="NOT"]
DEG1=DEG1[order(DEG1$log2FoldChange,decreasing = T),]
cg_up=rownames(DEG1)[DEG1$Change=="UP"][1:10]
cg_down=rownames(DEG1)[DEG1$Change=="DOWN"][1:10]
cg=c(cg_up,cg_down)
table(cg%in%rownames(dat))
cg=cg[cg%in%rownames(dat)]
#load('figure06_two_group/tumor_e2f6_top2a_deg_relationship_above5.Rdata')#后面的相关性分析得到的，直接跑会报错
#correlation_genes=b2$row[!duplicated(b2$row)]
h1 = draw_heatmap(dat[cg,],Group,n_cutoff = 2,cluster_cols = F,show_rownames = T,annotation_legend = T,legend = T)
#h2 = draw_heatmap(dat[correlation_genes,],Group,n_cutoff = 2,cluster_cols = F,show_rownames = T,annotation_legend = T,legend = T)
v1 = draw_volcano(DEG1[cg1,],pkg = 1,logFC_cutoff = logFC_t,adjust = T,color = c('blue','red'),symmetry = T)
library(patchwork)
h1|v1 +plot_layout(guides = 'collect') &theme(legend.position = 'right')
ggsave('figure06_two_group/heatmap_vocan_two_group.pdf',width = 12,height = 5,dpi = 320,units = 'in')
ggsave('figure06_two_group/heatmap_two_group.pdf',h1,width = 6,height = 5,dpi = 320,units = 'in')
ggsave('figure06_two_group/vocan_two_group.pdf',v1,width = 6,height = 5,dpi = 320,units = 'in')
#ggsave('figure06_two_group/correlation_heatmap_two_group.pdf',h2,width = 6,height = 5,dpi = 320,units = 'in')


#高水平组与其他组差异基因与E2F6的相关性分析，以及其在正常组与肿瘤组中的表达情况。
rm(list = ls())
#install.packages('DAAG')
library(DAAG)
load('/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/data_symbols_lihc.Rdata')
load('figure06_two_group/tcga_meta_two_groups.Rdata')
exp=count1
exp[1:4,1:4]
exp=exp[,str_sub(colnames(exp),14,15)!='02']
table(str_sub(colnames(exp),14,15))
exp=exp[,str_sub(colnames(exp),14,15)=='01']
rownames(meta)[1:10]
colnames(exp)=str_sub(colnames(exp),1,12)
table(colnames(exp)%in%rownames(meta))
exp=exp[,colnames(exp)%in%rownames(meta)]
exp=exp[,match(rownames(meta),colnames(exp))]
identical(rownames(meta),colnames(exp))
meta_o=rownames(meta)[meta$Group=='Others']
meta_H=rownames(meta)[meta$Group!='Others']
exp_o=exp[,colnames(exp)%in%meta_o]
exp_H=exp[,colnames(exp)%in%meta_H]
identical(rownames(exp_o),rownames(exp_H))
exp=cbind(exp_o,exp_H)
meta=meta[match(colnames(exp),rownames(meta)),]
identical(rownames(meta),colnames(exp))
Group=meta$Group
Group=factor(Group,levels = c('Others','E2F6_H&TOP2A_H'))
table(Group)
class(Group)
a=exp
exp=log2(cpm(exp)+1)
load('figure06_two_group/tcga_deseqdata_DEG_two_group.Rdata')
cg1 = rownames(DEG1)[DEG1$Change !="NOT"]
f=c('E2F6','TOP2A')
exp=exp[c(cg1,f),]
dim(exp)
library(Hmisc)
class(exp)
mydata=exp
res2 <- rcorr(as.matrix(t(mydata)))
res2
res2$r
res2$P
res3_r=res2$r[1:597,598:599]
res3_P=res2$P[1:597,598:599]
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat) 
  data.frame( row = rownames(cormat)[row(cormat)[ut]], 
              column = rownames(cormat)[col(cormat)[ut]], 
              cor =(cormat)[ut], 
              p = pmat[ut] )
}
b=flattenCorrMatrix(res2$r, res2$P)
write.csv(b,file = 'figure06_two_group/tumor_e2f6_top2a_deg_relationship.csv')
table(b$column%in%f)
b1=b[b$column%in%f,]
b2=b1[b1$cor>0.5,]
table(duplicated(b2$row))
max(b2$cor)
min(b2$cor)
write.csv(b2,file = 'figure06_two_group/tumor_e2f6_top2a_deg_relationship_above5.csv')
save(b2,file = 'figure06_two_group/tumor_e2f6_top2a_deg_relationship_above5.Rdata')
a = read.csv('figure06_two_group/tumor_e2f6_top2a_deg_relationship_above5.csv')
table(a$row)
# Group=ifelse(colnames(exp)%in%meta_o,0,1)
# names(Group)=colnames(exp)
# ##先进行单因素回归分析
# data=as.data.frame(t(exp))
# data$Group=Group
# data$Group <- factor(data$Group, levels = c('0','1'), labels = c('Low', 'High'))
# ncol(data)
# colnames(data)
# #install.packages('tableone')
# library(tableone)
# table2 <- CreateTableOne(vars = colnames(data)[1:597], # 指定对哪些变量进行分析
#                          data = data,
#                          strata = 'Group', # 指定分组
#                          # 是否对总样本进行分析
#                          addOverall = F)
# 
# result2 <- print(table2, 
#                  # 是否对分类变量全部展示
#                  showAllLevels = F,
#                  # 指定非参数检验变量，exact选项可以指定确切概率检验的变量
#                  nonnormal = colnames(data)[1:597])
# 
# class(result2)
# results=as.data.frame(results)
# write.csv(result2, file = 'figure06_two_group/tcga_uni_factor_analysis.csv')
# 
# #共线性分析
# library(corrplot)
# a=data[,1:597]
# corr <- cor(as.matrix(a))
# corrplot.mixed(corr)# 简单进行可视化,查看是存在多重共线性，存在才进行LASSO
# write.csv(corr, file = 'figure06_two_group/tcga_correlation_analysis.csv')
# 
# model.old <- glm(Group~., data = data, family = binomial())
# summary(model.old)$coefficients
# 
# # 逐步回归stepAIC函数在MASS包中
# library(MASS)
# model.both <- stepAIC(model.old, direction = 'both') # 'forward', 'backward'
# # 查看模型结果
# summary(model.both)$coefficients
# # 计算OR值及可信区间
# exp(cbind('OR' = coef(model.both), confint(model.both)))
# save(model.both,file = 'figure06_two_group/stepwise_deg_two_group.Rdata')
f=c('E2F6','TOP2A')
#table(colnames(data)%in%f)
library(tidyverse)
library(corrplot)
if(!require(ggcor))devtools::install_github("Github-Yilei/ggcor",upgrade = F)
exp=mydata
dim(exp)
exp[1:4,1:4]
genes_relationship=b2$row[!duplicated(b2$row)]
exp=exp[c(genes_relationship,'TOP2A'),]
e=t(exp)
dim(e)
datp <- cor.mtest(e)$p[15:16,1:14] %>%
  as.data.frame()%>% 
  rownames_to_column("A") %>% 
  gather("B","p",-A)

data <- cor(e)[15:16,1:14] %>%
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

# data <- cor(e) %>%
#   as.data.frame()%>% 
#   rownames_to_column("A") %>% 
#   gather("B","value1",-A) %>% 
#   mutate(value2 = -log2(sample(value1, length(value1))+1))
#install.packages("ggnewscale")
library(ggnewscale)
e2f6_top2a_relationship_genes=ggplot()+
  geom_rectriangle(data = data, aes(A, B, fill = -log10(p)),type = "upper", r = 1)+
  scale_fill_gradient(low = "white", high = "#a5d3ed",na.value = "white")+
  ggnewscale::new_scale_fill()+
  geom_rectriangle(data = data, aes(A, B, fill = r),type = "lower", r = 1)+
  scale_fill_gradient2(high = "#e7857b", mid = "white",low = "#6d9cc5")+
  labs(x = "", y = "")+
  geom_text(data = data,aes(A,B,label = ps),
            nudge_y = 0.05)+
  theme_bw()
ggsave('figure06_two_group/e2f6_top2a_relationship_genes.pdf',e2f6_top2a_relationship_genes,width = 5,height = 6,units = 'in')

##相关基因在肝癌与癌旁组织中的表达情况--箱式图####
##TCGA
rm(list = ls())
load('/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/data_symbols_lihc.Rdata')
load('/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/clinical_lihc_postclear.Rdata')
fpkm_lihc1[1:4,1:4]
boxplot(fpkm_lihc1[,1:20])
exp=log2(fpkm_lihc1+1)
load('figure06_two_group/tumor_e2f6_top2a_deg_relationship_above5.Rdata')
correlation_genes=b2$row[!duplicated(b2$row)]
exp=exp[correlation_genes,]
table(str_sub(colnames(exp),14,15))
exp=exp[,str_sub(colnames(exp),14,15)!='02']
Group=ifelse(str_sub(colnames(exp),14,15)=='01','Tumor','Non-tumor')
Group=factor(Group,levels = c('Non-tumor','Tumor'))
exp=exp[-15,]
rownames(exp)
table(Group)
draw_boxplot(exp,Group,xlab = ' ',ylab = 'Expression',sort = F)
ggsave('figure06_two_group/tcga_correlation_genes_deg_tumor_noraml.pdf',width = 6,height = 5,units = 'in')

##相关基因在正常组与肿瘤组中的表达情况--热图####
rm(list = ls())
#install.packages('DAAG')
library(DAAG)
load('/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/data_symbols_lihc.Rdata')
load('figure06_two_group/tcga_meta_two_groups.Rdata')
exp=fpkm_lihc1
exp[1:4,1:4]
boxplot(exp[,1:10])
exp = log2(exp+1)
exp=exp[,str_sub(colnames(exp),14,15)!='02']
normal_n = colnames(exp)[str_sub(colnames(exp),14,15)=='11']
tumor_n = colnames(exp)[str_sub(colnames(exp),14,15)=='01']
exp_normal = exp[,normal_n]
exp_tumor = exp[,tumor_n]
exp = cbind(exp_normal,exp_tumor)
Group = ifelse(str_sub(colnames(exp),14,15)=='01','Tumor','Non-tumor')
Group=factor(Group,levels = c('Non-tumor','Tumor'))
load('figure06_two_group/tumor_e2f6_top2a_deg_relationship_above5.Rdata')
correlation_genes=b2$row[!duplicated(b2$row)]
exp = exp[correlation_genes,]
rownames(exp)
exp = exp[-15,]
draw_heatmap(exp,Group,n_cutoff = 2,cluster_cols = F,show_rownames = T,annotation_legend = T,legend = T)
ggsave('figure06_two_group/tcga_correlation_genes_deg_tumor_noraml_heatmap.pdf',width = 6,height = 5,units = 'in')

##相关基因在肝癌不同分期中的表达情况####
rm(list=ls())
load('figure02_E2F6/lihc_sur_data.Rdata')
load('figure04_E2F6_TOP2A/tcga_cluster_clinialinfo.Rdata')
a = c('E2F6','TOP2A')
exp=exprSet[a,]
data=as.data.frame(t(exp))
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6_TOP2A=ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='High','E2F6_H&TOP2A_H','Others')
data$Group_E2F6_TOP2A=factor(data$Group_E2F6_TOP2A,levels = c('Others','E2F6_H&TOP2A_H'))
table(data$Group_E2F6_TOP2A)
identical(rownames(data),rownames(meta))
meta$Group = data$Group_E2F6_TOP2A

load('figure06_two_group/tumor_e2f6_top2a_deg_relationship_above5.Rdata')
correlation_genes=b2$row[!duplicated(b2$row)]
exp = exprSet[correlation_genes,]

an = data.frame(Group = meta$Group,
                Cluster = meta$Cluster,
                Stage = meta$Stage,
                row.names = rownames(meta))
an = arrange(an,Group,Cluster,Stage)
exp = exp[,match(rownames(an),colnames(exp))]
identical(rownames(an),colnames(exp))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
legend_col = list(Group = c('Others'='#6d9cc5','E2F6_H&TOP2A_H'='#e7857b'),
                  Cluster=c('Cluster1'='#a5d3ed','Cluster2'='#b5aad4'),
                  Stage = c('I'='#A4CBCC','II'='#F7D1E6','III'='#b4c5e6','IV' = '#bea095'))
exp = exp[-15,]
pdf('figure06_two_group/tcga_correlation_gene_group_cluster_stage_heatmap.pdf',width = 6,height = 4)
pheatmap(exp,scale = "row",
         show_colnames = F,
         annotation_col = an,
         cluster_cols = F,
         cluster_rows = F,
         legend_breaks=seq(-2,2,1),
         breaks=bk,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         annotation_colors = legend_col)
dev.off()
meta = meta[match(colnames(exp),rownames(meta)),]
identical(colnames(exp),rownames(meta))
Group=meta$Group
draw_boxplot(exp,Group,xlab = ' ',ylab = 'Expression',sort = F)
ggsave('figure06_two_group/tcga_correlation_genes_tumor_two_group.pdf',width = 6,height = 5,units = 'in')

Stage = meta$Stage
draw_boxplot(exp,Stage,xlab = ' ',ylab = 'Expression',sort = F)
ggsave('figure06_two_group/tcga_correlation_genes_tumor_stage.pdf',width = 8,height = 5,units = 'in')

Cluster = meta$Cluster
draw_boxplot(exp,Cluster,xlab = ' ',ylab = 'Expression',sort = F)
ggsave('figure06_two_group/tcga_correlation_genes_tumor_cluster.pdf',width = 6,height = 5,units = 'in')


##相关基因与肝癌预后的关系
logrankfile = "figure06_two_group/tcga_correlation_genes_log_rank_p.Rdata"
identical(rownames(meta),colnames(exp))
if(!file.exists(logrankfile)){
  log_rank_p <- apply(exp , 1 , function(gene){
    # gene=exprSet[1,]
    meta$group=ifelse(gene>median(gene),'high','low')  
    data.survdiff=survdiff(Surv(Time, Event)~group,data=meta)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    return(p.val)
  })
  log_rank_p=sort(log_rank_p)
  save(log_rank_p,file = logrankfile)
}
load(logrankfile)
table(log_rank_p<0.01) 
table(log_rank_p<0.05) 

##相关基因的单因素cox回归分析
coxfile = "figure06_two_group/tcga_correlation_genes_unicox.Rdata"
if(!file.exists(coxfile)){
  cox_results <-apply(exp , 1 , function(gene){
    #gene= exprSet[1,]
    meta$gene = gene
    #可直接使用连续型变量
    m = coxph(Surv(Time, Event) ~ gene, data =  meta)
    #也可使用二分类变量
    #meta$group=ifelse(gene>median(gene),'high','low') 
    #m=coxph(Surv(time, event) ~ group, data =  meta)
    beta <- coef(m)
    se <- sqrt(diag(vcov(m)))
    HR <- exp(beta)
    HRse <- HR * se
    
    #summary(m)
    tmp <- round(cbind(coef = beta, 
                       se = se, z = beta/se, 
                       p = 1 - pchisq((beta/se)^2, 1),
                       HR = HR, HRse = HRse,
                       HRz = (HR - 1) / HRse, 
                       HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                       HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                       HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
    
    return(tmp['gene',]) 
    #return(tmp['grouplow',])#二分类变量
  })
  cox_results=as.data.frame(t(cox_results))
  save(cox_results,file = coxfile)
}
load(coxfile)
table(cox_results$p<0.01)
table(cox_results$p<0.05)

lr = names(log_rank_p)[log_rank_p<0.05];length(lr)
cox = rownames(cox_results)[cox_results$p<0.05];length(cox)
length(intersect(lr,cox))
save(lr,cox,file = 'figure06_two_group/tcga_correlation_genes_log_unicox_genes.Rdata')
write.csv(cox_results,file = 'figure06_two_group/tcga_correlation_genes_log_unicox_genes.csv',quote = F)
#load('figure06_two_group/tcga_correlation_genes_log_unicox_genes.Rdata')

rm(list = ls())
a=read.csv('figure06_two_group/tcga_correlation_genes_log_unicox_genes2.csv')
tmp=a
colnames(tmp)
tmp=rbind(c("Trait","Coef" ,NA, NA, NA, "HR", "P"),tmp)
str(tmp)
tmp$HR=as.numeric(tmp$HR)
tmp$lower=as.numeric(tmp$lower)
tmp$upper=as.numeric(tmp$upper)
dev.off()
library(forestplot)
dev.off()
pdf('figure06_two_group/tcga_e2f6_top2a_uniforest_correlation_unicox.pdf',width = 6,height = 6)
forestplot(
  tmp[, c(1, 2, 6,7)],
  mean = tmp[, 3],
  lower = tmp[, 4],
  upper = tmp[, 5],
  zero = 1,
  boxsize = 0.2,
  col = fpColors(box = '#6d9cc5', lines = 'black', zero = 'grey'),
  lty.ci = "solid",
  graph.pos = 3,
  #xticks = F,
  is.summary = c(T, rep(F, 14)),
  align = "l",
  hrzl_lines = list(
    "1" = gpar(lty=1.5),
    "2" = gpar(lty=1.5),
    "16"= gpar(lty=1.5)),
  colgap = unit(5, 'mm')
)
dev.off()


#构建logistics回归模型
rm(list = ls())
load('figure02_E2F6/lihc_sur_data.Rdata')
exprSet[1:4,1:4]
load('figure04_E2F6_TOP2A/tcga_cluster_clinialinfo.Rdata')
identical(rownames(meta),colnames(exprSet))
exp=exprSet[c('E2F6','TOP2A'),]
train=as.data.frame(t(exp))
identical(rownames(train),rownames(meta))
train$Group_E2F6=ifelse(train$E2F6>median(train$E2F6),'High','Low')
train$Group_TOP2A=ifelse(train$TOP2A>median(train$TOP2A),'High','Low')
train$Group_E2F6_TOP2A=ifelse(train$Group_E2F6=='High'&train$Group_TOP2A=='High','E2F6_H&TOP2A_H','Others')
train$Group_E2F6_TOP2A=factor(train$Group_E2F6_TOP2A,levels = c('Others','E2F6_H&TOP2A_H'))
identical(rownames(meta),rownames(train))
colnames(meta)
meta$Group=train$Group_E2F6_TOP2A
model<-glm(Group_E2F6_TOP2A~E2F6+TOP2A,data=train,family = "binomial")
summary(model)
model2<-step(object = model,trace = 0)
summary(model2)
anova(object = model2,test = "Chisq")
gs=c('E2F6','TOP2A')
save(model2,gs,file = 'figure06_two_group/train_e2f6_top2a_logistic_model_genes.Rdata')
fp <- predict(model2,newdata = train)
identical(rownames(meta),rownames(train))
meta$Predict_score=fp
save(meta,file='figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
##根据构建的logistics模型画nomograph图
library(rms)
library(Hmisc);library(lattice);library(survival);library(Formula);library(ggplot2);library(SparseM)

library(regplot)
CHDfit<-glm(Group_E2F6_TOP2A ~ E2F6+TOP2A, family=binomial(link = logit), data=train)
regplot(CHDfit,
        observation=train[1,], 
        obscol = '#6d9cc5',
        plots = c("no plot","no plot"),
        points = T,
        prfail = T)


###构建验证集
load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
load('figure04_E2F6_TOP2A/icgc_cluster_clinialinfo.Rdata')
load('figure06_two_group/train_e2f6_top2a_logistic_model_genes.Rdata')
exprSet[1:4,1:4]
exp=exprSet
exp=exp[gs,]
test=as.data.frame(t(exp))
identical(rownames(test),rownames(meta))
test$Group_E2F6=ifelse(test$E2F6>median(test$E2F6),'High','Low')
test$Group_TOP2A=ifelse(test$TOP2A>median(test$TOP2A),'High','Low')
test$Group_E2F6_TOP2A=ifelse(test$Group_E2F6=='High'&test$Group_TOP2A=='High','E2F6_H&TOP2A_H','Others')
test$Group_E2F6_TOP2A=factor(test$Group_E2F6_TOP2A,levels = c('Others','E2F6_H&TOP2A_H'))
###########绘制ROC曲线
##训练集
pred_f_training<-predict(model2,train)

train$predict=pred_f_training
##验证集
pred_f_testing<-predict(model2,test)
test$predict=pred_f_testing
write.csv(train,file='figure06_two_group/train_roc_data.csv')
write.csv(test,file='figure06_two_group/test_roc_data.csv')
identical(rownames(test),rownames(meta))
colnames(test)
meta$Group=test$Group_E2F6_TOP2A
meta$Predict_score=test$predict
###########绘制校准曲线
##训练集和验证集
fit1 <- lrm(Group_E2F6_TOP2A ~ E2F6+TOP2A, data = train, x=TRUE,y=TRUE)
cal1 <- calibrate(fit1, method = "boot", B=1000)

fit2 <- lrm(Group_E2F6_TOP2A ~ E2F6+TOP2A, data = test, x=TRUE,y=TRUE)
cal2 <- calibrate(fit2, method = "boot", B=1000)
dev.off()
pdf('figure06_two_group/train_cal_curv.pdf',width = 5,height = 5)
plot(cal1,
     xlim =c(0,1),
     ylim =c(0,1),
     xlab = "Predicted Probability",
     ylab = "Observed Probability",
     legend = FALSE,
     subtitles = FALSE)
abline(0,1,col = "#224444",lty = 3,lwd = 2)
lines(cal1[,c("predy","calibrated.orig")], type = "l",lwd = 2,col="#2166AC",pch =16)
lines(cal1[,c("predy","calibrated.corrected")], type = "l",lwd = 2,col="#B2182B",pch =16)
legend("topleft",#图例的位置
       c("Apparent","Ideal","Bias-corrected"),#图例文字
       lty = c(2,2),
       lwd = c(2,2),#图例中线的粗细
       col = c("#224444","#2166AC","#B2182B"), #图例线的颜色，与文字对应
       bty = "n",#不显示图例边框,# "o"为加边框
       cex = 1.2) #图例字体大小
dev.off()
##验证集中的校准图
pdf('figure06_two_group/test_cal_curv.pdf',width = 5,height = 5)
plot(cal2,
     xlim =c(0,1),
     ylim =c(0,1),
     xlab = "Predicted Probability",
     ylab = "Observed Probability",
     legend = FALSE,
     subtitles = FALSE)
abline(0,1,col = "#de9e20",lty = 3,lwd = 2)
lines(cal2[,c("predy","calibrated.orig")], type = "l",lwd = 2,col="#0758a6",pch =16)
lines(cal2[,c("predy","calibrated.corrected")], type = "l",lwd = 2,col="#d47222",pch =16)
legend("topleft",#图例的位置
       c("Apparent","Ideal","Bias-corrected"),#图例文字
       lty = c(2,2),
       lwd = c(2,2),#图例中线的粗细
       col = c("#de9e20","#0758a6","#d47222"), #图例线的颜色，与文字对应
       bty = "n",#不显示图例边框,# "o"为加边框
       cex = 1.2) #图例字体大小
dev.off()
#devtools::install_github('yikeshu0611/ggDCA')
library(ggDCA)
library(rms)
train$Group_E2F6_TOP2A=ifelse(train$Group_E2F6_TOP2A=='Others',0,1)
PredicScore <- rms::lrm(Group_E2F6_TOP2A ~ E2F6+TOP2A, data = train)
E2F6=rms::lrm(Group_E2F6_TOP2A ~ E2F6 , data = train)
TOP2A=rms::lrm(Group_E2F6_TOP2A ~  TOP2A , data = train)
#remove.packages('dcurves')
library(ggDCA)
library(rms)
data  <- dca(PredicScore,E2F6,TOP2A)
ggplot(data)
#ggsave('figure06_two_group/train_dca.pdf',height = 5,width = 6,dpi=320,units = 'in')
test$Group_E2F6_TOP2A=ifelse(test$Group_E2F6_TOP2A=='Others',0,1)
PredicScore <- lrm(Group_E2F6_TOP2A ~ E2F6+TOP2A, data = test)
E2F6=lrm(Group_E2F6_TOP2A ~ E2F6 , data = test)
TOP2A=lrm(Group_E2F6_TOP2A ~  TOP2A , data = test)

data  <- dca(PredicScore,E2F6,TOP2A)
ggplot(data)
#ggsave('figure06_two_group/test_dca.pdf',height = 5,width =6 ,dpi=320,units = 'in')
save(meta,file='figure06_two_group/icgc_predict_factor_clinicalinformation.Rdata')

rm(list = ls())
load('figure06_two_group/train_e2f6_top2a_logistic_model_genes.Rdata')
fit.result<-summary(model2)
df1<-fit.result$coefficients
df2<-confint(model2)
df3<-cbind(df1,df2)

df4<-data.frame(df3[-1,c(1,4,5,6)])
df4$Var<-rownames(df4)
colnames(df4)<-c("OR","Pvalue","OR_1","OR_2","Var")
df5<-df4[,c(5,1,2,3,4)]
df5$OR_mean<-df5$OR
df5$OR<-paste0(round(df5$OR,2),
               "(",
               round(df5$OR_1,2),
               "~",
               round(df5$OR_2,2),
               ")")
df5$Pvalue<-round(df5$Pvalue,3)
write.csv(df5,file = "figure06_two_group/tcga_logistic_forestplot_2genes.csv",
          quote = F,row.names = F)
library(forestplot)
fp<-read.csv("figure06_two_group/tcga_logistic_forestplot_2genes2.csv",header=T)
pdf('figure06_two_group/2multilogistics_gene.pdf',width = 5,height =3)
forestplot(labeltext=as.matrix(fp[,1:3]),
           mean=fp$OR_mean,
           lower=fp$OR_1,
           upper=fp$OR_2,
           zero=1,
           boxsize=0.4,
           lineheight = unit(7,'mm'),
           colgap=unit(5,'mm'),
           col=fpColors(box='#1075BB',
                        summary='#8B008B',
                        lines = 'black',
                        zero = 'grey'),
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),
                            xlab  = gpar(cex = 0.8),
                            cex = 0.9),
           lty.ci = "solid",
           title = "", 
           line.margin = 0.08,
           graph.pos=2,
           align = "l",
           hrzl_lines = list(
             "1" = gpar(lty=1),
             "2" = gpar(lty=1),
             "4"= gpar(lty=1)))
dev.off()


# GO_KEGG富集分析 -------------------------------------------------------------

#GO富集分析####
rm(list = ls())
load('figure06_two_group/tcga_deseqdata_DEG_two_group.Rdata')
library(AnnotationDbi)
library(org.Hs.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)#画图包
diff=DEG1[DEG1$Change!='NOT',]
#2、基因id转换，kegg和go富集用的ID类型是ENTREZID）
gene.df <- bitr(rownames(diff),fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)#TCGA数据框如果没有进行基因注释，那么fromType应该是Ensembl，各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求                     
gene <- gene.df$ENTREZID
#3、GO富集
##CC表示细胞组分，MF表示分子功能，BP表示生物学过程，ALL表示同时富集三种过程，选自己需要的,我一般是做BP,MF,CC这3组再合并成一个数据框，方便后续摘取部分通路绘图。
ego_ALL <- enrichGO(gene = gene,#我们上面定义了
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，一般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 1,#P值可以取0.05
                    qvalueCutoff = 1,
                    readable = TRUE)

ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   readable = TRUE)

#4、将结果保存到当前路径
ego_ALL <- as.data.frame(ego_ALL)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)
ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)
write.csv(ego_ALL,file = "figure06_two_group/two_deg_ego_ALL.csv",row.names = T)
save(ego,ego_ALL,ego_result_BP,ego_result_CC,ego_BP,ego_CC,ego_MF,
     file = "figure06_two_group/two_deg_ego_ALL.Rdata")
load('figure06_two_group/two_deg_ego_ALL.Rdata')
ego[1:4,]
table(ego$ID)
bp = c('GO:0098742','GO:0045216','GO:0050673','GO:0048762','GO:0050678',
       'GO:0045785','GO:0050679','GO:0007162','GO:0022617','GO:0001837')
cc = c('GO:0062023','GO:0016323','GO:0045111','GO:0009925','GO:0045178',
       'GO:0098862','GO:0032838','GO:0044298','GO:0005881','GO:0043256')
mf = c('GO:0045296','GO:0005201','GO:0005516','GO:0030674','GO:0097110',
       'GO:0005200','GO:0015631','GO:0000146','GO:0003774','GO:0005178')

ego_result_BP= arrange(ego_result_BP,Count)[bp,]
ego_result_CC= arrange(ego_result_CC,Count)[cc,]
ego_result_MF= arrange(ego_result_MF,Count)[mf,]


##将以上我们摘取的部分通路重新组合成数据框
go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  Type=factor(c(rep("Biological process", 10), 
                rep("Cellular component", 10),
                rep("Molecular function", 10)), 
              levels=c("Biological process", "Cellular component","Molecular function" )))

##通路的名字太长了，我选取了通路的前五个单词作为通路的名字
for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") #这里的5就是指5个单词的意思，可以自己更改
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
GoTop10Barplot2=ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=Type)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = COLS) +
  theme_bw() +
  xlab("") +
  coord_flip()+
  ylab("Num of Genes") +
  labs(title = "The Most Enriched GO Terms")+
  theme(axis.text.x=element_text( color="black",angle = 70,vjust = 1, hjust = 1 ))+
  theme(axis.text.y=element_text( color="black" ))#angle是坐标轴字体倾斜的角度，可以自己设置

ggsave('figure06_two_group/two_deg_go.pdf',width = 8,height = 6,units = 'in')

#KEGG富集分析####
rm(list = ls())
load('figure06_two_group/tcga_deseqdata_DEG_two_group.Rdata')
#2、基因id转换，kegg和go富集用的ID类型是ENTREZID）
gene= DEG1[DEG1$Change != 'NOT',] 
gene <- bitr(rownames(gene),fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)#TCGA数据框如果没有进行基因注释，那么fromType应该是Ensembl，各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求                     


kk <- enrichKEGG(gene = gene$ENTREZID,keyType = "kegg",organism= "human", qvalueCutoff = 1, pvalueCutoff=1)
hh <- as.data.frame(kk)
rownames(hh) <- 1:nrow(hh)
#hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
write.csv(hh,file = 'figure06_two_group/two_deg_kegg.csv')
hh = arrange(hh,desc(GeneRatio))

table(hh$category)
hh1 = hh[c(1:20),]

hh1$order=factor(rev(as.integer(rownames(hh1))),labels = rev(hh1$Description))

hh1$Description

KeggTop30Pointplot2=ggplot(hh1,aes(y=order,x=Count,fill = category))+
  geom_point(aes(size=Count,color=-1*p.adjust))+# 修改点的大小
  scale_color_gradient(low="#66C3A5",high =  "#8DA1CB")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()
#（2）对上调/下调/所有差异基因进行富集分析

ggsave('figure06_two_group/two_deg_kegg.pdf',width = 6,height = 5,units = 'in')
