rm(list = ls())
getwd()
###1.TCGA数据库中E2F6在正常与肝癌的差异
load('/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/data_symbols_lihc.Rdata')
load('/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/clinical_lihc_postclear.Rdata')
fpkm_lihc1[1:4,1:4]
exp=log2(fpkm_lihc1+1)
a=c('E2F6','TOP2A')
exp=exp[a,]
data=as.data.frame(t(exp))
Group=ifelse(str_sub(rownames(data),14,15)=='01','Tumor','Non-tumor')
Group=factor(Group,levels = c('Non-tumor','Tumor'))
data$Group=Group
colnames(data)
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6_TOP2A=ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='High','E2F6_H&TOP2A_H',
                             ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='Low','E2F6_H&TOP2A_L',
                                    ifelse(data$Group_E2F6=='Low'&data$Group_TOP2A=='High','E2F6_L&TOP2A_H','E2F6_L&TOP2A_L')))
data$Group_E2F6_TOP2A=as.factor(data$Group_E2F6_TOP2A)
table(data$Group_E2F6_TOP2A)

library(dplyr)

p=ggstatsplot::ggbarstats(
  data = data,
  x = Group_E2F6_TOP2A,
  y = Group,
  ggtheme = ggplot2::theme_classic(),
  palette = "category10_d3", 
  package = "ggsci", 
  k = 3, 
  perc.k = 1 ,
  results.subtitle = T,xlab ='',
  return='plot',messages=T,
  bf.message = F
) 
ggsave('figure04_E2F6_TOP2A/tcga_e2f6_top2a_tumor_noraml.pdf',p,width = 5,height = 4,units = 'in')

rm(list = ls())
load('/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/icgc_tumor_normal_rawcount.Rdata')
colnames(exp)
exp[1:4,1:4]

Group=ifelse(str_detect(colnames(exp),'cancer'),'Tumor','Non-tumor')
Group=factor(Group,levels = c('Non-tumor','Tumor'))
table(Group)
exp_n=exp[,Group=='Non-tumor'];dim(exp_n)
exp_t = exp[,Group == "Tumor"];dim(exp_t)
identical(rownames(exp_n),rownames(exp_t))
exp=cbind(exp_n,exp_t)
Group=ifelse(str_detect(colnames(exp),'cancer'),'Tumor','Non-tumor')
Group=factor(Group,levels = c('Non-tumor','Tumor'))
names(Group)=colnames(exp)
table(Group)
exp[1:4,1:4]
exp=log2(cpm(exp)+1)
a=c('E2F6','TOP2A')
exp=exp[a,]
data=as.data.frame(t(exp))
identical(rownames(data),names(Group))
data$Group=Group
colnames(data)
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6_TOP2A=ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='High','E2F6_H&TOP2A_H',
                             ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='Low','E2F6_H&TOP2A_L',
                                    ifelse(data$Group_E2F6=='Low'&data$Group_TOP2A=='High','E2F6_L&TOP2A_H','E2F6_L&TOP2A_L')))
data$Group_E2F6_TOP2A=as.factor(data$Group_E2F6_TOP2A)
table(data$Group_E2F6_TOP2A)

library(dplyr)

p=ggstatsplot::ggbarstats(
  data = data,
  x = Group_E2F6_TOP2A,
  y = Group,
  ggtheme = ggplot2::theme_classic(),
  palette = "category10_d3", 
  package = "ggsci", 
  k = 3, 
  perc.k = 1 ,
  results.subtitle = T,xlab ='',
  return='plot',messages=T,
  bf.message = F
) 
ggsave('figure04_E2F6_TOP2A/icgc_e2f6_top2a_tumor_noraml.pdf',p,width = 5,height = 4,units = 'in')

rm(list = ls())
load("figure02_E2F6/lihc_sur_data.Rdata")
identical(colnames(exprSet),rownames(meta))
a=c('E2F6','TOP2A')
exp=exprSet[a,]
data=as.data.frame(t(exp))
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6_TOP2A=ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='High','E2F6_H&TOP2A_H',
                             ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='Low','E2F6_H&TOP2A_L',
                                    ifelse(data$Group_E2F6=='Low'&data$Group_TOP2A=='High','E2F6_L&TOP2A_H','E2F6_L&TOP2A_L')))
data$Group_E2F6_TOP2A=as.factor(data$Group_E2F6_TOP2A)
table(data$Group_E2F6_TOP2A)

identical(rownames(data),rownames(meta))
meta$Group=data$Group_E2F6_TOP2A

#E2F6的生存分析
identical(rownames(meta),rownames(data))
sfit1=survfit(Surv(Time, Event)~Group, data=meta)
g1=ggsurvplot(fit = sfit1,
              palette = c("#0f3391", "#d68d28",'#ddb419','#a11109'),
              risk.table =TRUE,
              pval =TRUE,
              conf.int =F,
              ncensor.plot = TRUE,data = meta)
library(patchwork)
p1=g1$plot
p2=g1$table
p_sur1=p1 /p2  + plot_layout(design = 'A
                         A
                         A
                         A
                         B')
ggsave('figure04_E2F6_TOP2A/tcga_survival_E2F6_top2a.pdf',p_sur1,width = 10,height = 6,units = 'in')




rm(list = ls())
load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
identical(colnames(exprSet),rownames(meta))

a=c('E2F6','TOP2A')
exp=exprSet[a,]

data=as.data.frame(t(exp))
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6_TOP2A=ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='High','E2F6_H&TOP2A_H',
                             ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='Low','E2F6_H&TOP2A_L',
                                    ifelse(data$Group_E2F6=='Low'&data$Group_TOP2A=='High','E2F6_L&TOP2A_H','E2F6_L&TOP2A_L')))
data$Group_E2F6_TOP2A=as.factor(data$Group_E2F6_TOP2A)
table(data$Group_E2F6_TOP2A)
identical(rownames(data),rownames(meta))
meta$Group=data$Group_E2F6_TOP2A

#E2F6的生存分析
identical(rownames(meta),rownames(data))
colnames(meta)=c("Gender","Age","Event","Stage","Time","Group")
sfit1=survfit(Surv(Time, Event)~Group, data=meta)
g1=ggsurvplot(fit = sfit1,
              palette = c("#0f3391", "#d68d28",'#ddb419','#a11109'),
              risk.table =TRUE,
              pval =TRUE,
              conf.int =F,
              ncensor.plot = TRUE,data = meta)
library(patchwork)
p1=g1$plot
p2=g1$table
p_sur1=p1 /p2  + plot_layout(design = 'A
                         A
                         A
                         A
                         B')
ggsave('figure04_E2F6_TOP2A/icgc_survival_E2F6_top2a.pdf',p_sur1,width = 10,height = 6,units = 'in')

rm(list = ls())
load("figure02_E2F6/lihc_sur_data.Rdata")
identical(colnames(exprSet),rownames(meta))
a=c('E2F6','TOP2A')
exp=exprSet[a,]
data=as.data.frame(t(exp))
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6_TOP2A=ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='High','E2F6_H&TOP2A_H',
                             ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='Low','E2F6_H&TOP2A_L',
                                    ifelse(data$Group_E2F6=='Low'&data$Group_TOP2A=='High','E2F6_L&TOP2A_H','E2F6_L&TOP2A_L')))
data$Group_E2F6_TOP2A=as.factor(data$Group_E2F6_TOP2A)
table(data$Group_E2F6_TOP2A)

identical(rownames(data),colnames(exp))
Group=data$Group_E2F6_TOP2A

pca.plot = draw_pca(exprSet,Group,addEllipses = F);pca.plot
ggsave(pca.plot,file='figure04_E2F6_TOP2A/tcga_e2f6_top2a_pca.pdf',width = 5,height = 4,units = 'in')

rm(list = ls())
load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
identical(colnames(exprSet),rownames(meta))

a=c('E2F6','TOP2A')
exp=exprSet[a,]

data=as.data.frame(t(exp))
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6_TOP2A=ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='High','E2F6_H&TOP2A_H',
                             ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='Low','E2F6_H&TOP2A_L',
                                    ifelse(data$Group_E2F6=='Low'&data$Group_TOP2A=='High','E2F6_L&TOP2A_H','E2F6_L&TOP2A_L')))
data$Group_E2F6_TOP2A=as.factor(data$Group_E2F6_TOP2A)
table(data$Group_E2F6_TOP2A)
identical(rownames(data),colnames(exp))
Group=data$Group_E2F6_TOP2A

pca.plot = draw_pca(exprSet,Group,addEllipses = F);pca.plot
ggsave(pca.plot,file='figure04_E2F6_TOP2A/icgc_e2f6_top2a_pca.pdf',width = 5,height = 4,units = 'in')

rm(list = ls())
load("figure02_E2F6/lihc_sur_data.Rdata")
identical(colnames(exprSet),rownames(meta))
a=c('E2F6','TOP2A')
exp=exprSet[a,]
data=as.data.frame(t(exp))
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6=factor(data$Group_E2F6,levels = c('Low','High'))
data$Group_TOP2A=factor(data$Group_TOP2A,levels = c('Low','High'))
identical(rownames(data),colnames(exp))
Group=data$Group_E2F6

pca.plot = draw_pca(exprSet,Group,addEllipses = F);pca.plot
ggsave(pca.plot,file='figure02_E2F6/tcga_e2f6_pca_tumor.pdf',width = 5,height = 4,units = 'in')

identical(rownames(data),colnames(exp))
Group=data$Group_TOP2A

pca.plot = draw_pca(exprSet,Group,addEllipses = F);pca.plot
ggsave(pca.plot,file='figure03_TOP2A/tcga_top2a_pca_tumor.pdf',width = 5,height = 4,units = 'in')



rm(list = ls())
load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
identical(colnames(exprSet),rownames(meta))

a=c('E2F6','TOP2A')
exp=exprSet[a,]

data=as.data.frame(t(exp))
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6=factor(data$Group_E2F6,levels = c('Low','High'))
data$Group_TOP2A=factor(data$Group_TOP2A,levels = c('Low','High'))
identical(rownames(data),colnames(exp))
Group=data$Group_E2F6

pca.plot = draw_pca(exprSet,Group,addEllipses = F);pca.plot
ggsave(pca.plot,file='figure02_E2F6/icgc_e2f6_pca_tumor.pdf',width = 5,height = 4,units = 'in')

identical(rownames(data),colnames(exp))
Group=data$Group_TOP2A

pca.plot = draw_pca(exprSet,Group,addEllipses = F);pca.plot
ggsave(pca.plot,file='figure03_TOP2A/icgc_top2a_pca_tumor.pdf',width = 5,height = 4,units = 'in')


# 一致性聚类 -------------------------------------------------------------------


#一致性聚类分析
rm(list = ls())
load("figure02_E2F6/lihc_sur_data.Rdata")
identical(colnames(exprSet),rownames(meta))
a=c('E2F6','TOP2A')
exp=exprSet[a,]
exp=as.matrix(exp)
df <-  sweep(exp,1, apply(exp,1,median,na.rm=T))
class(df)
library(ConsensusClusterPlus)
maxK <-  6 #最多分成几组
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 500,
                                 pItem = 0.8,
                                 pFeature = 1,
                                 clusterAlg = "km",
                                 distance = 'euclidean',
                                 seed = 1086,
                                 title="test",
                                 innerLinkage="complete",
                                 plot="pdf")

icl = calcICL(results,
              title="test",
              plot="pdf")
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK

for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK
#聚类结果
library(tidyverse)
table(results[[optK]]$consensusClass)
Cluster = results[[optK]]$consensusClass
identical(names(Cluster),rownames(meta))
meta$Cluster = Cluster
table(meta$Cluster)
meta$Cluster=ifelse(meta$Cluster==2,'Cluster1','Cluster2')
meta$Cluster=as.factor(meta$Cluster)
identical(colnames(exp),rownames(meta))
Cluster=meta$Cluster
draw_pca(exp,Cluster,addEllipses = T,color = c('#6d9cc5','#e7857b'))
ggsave('figure04_E2F6_TOP2A/tcga_e2f6_top2a_subtype_pca.pdf',height = 4,width = 5,dpi = 320,units = 'in')

draw_boxplot(exp,Cluster,xlab = ' ',ylab = 'Expression',sort = F,color = c('#6d9cc5','#e7857b'))
ggsave('figure04_E2F6_TOP2A/tcga_e2f6_top2a_subtype_boxplot.pdf',height = 4,width = 4,dpi=320,units = 'in')
colnames(meta)
library(survival)
library(survminer)
sfit1 <- survfit(Surv(Time, Event) ~ Cluster,
                data = meta)
g1=ggsurvplot(fit = sfit1,
              palette = c('#6d9cc5','#e7857b'),
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
ggsave('figure04_E2F6_TOP2A/tcga_subtype_e2f6_top2a_survive.pdf',p_sur1,width = 5,height = 4,units = 'in')


data=as.data.frame(t(exp))
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6_TOP2A=ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='High','E2F6_H&TOP2A_H',
                             ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='Low','E2F6_H&TOP2A_L',
                                    ifelse(data$Group_E2F6=='Low'&data$Group_TOP2A=='High','E2F6_L&TOP2A_H','E2F6_L&TOP2A_L')))
data$Group_E2F6_TOP2A=as.factor(data$Group_E2F6_TOP2A)
table(data$Group_E2F6_TOP2A)
identical(rownames(data),rownames(meta))
data$Cluster=meta$Cluster
library(dplyr)

p=ggstatsplot::ggbarstats(
  data = data,
  x = Group_E2F6_TOP2A,
  y = Cluster,
  ggtheme = ggplot2::theme_classic(),
  palette = "category10_d3", 
  package = "ggsci", 
  k = 3, 
  perc.k = 1 ,
  results.subtitle = T,xlab ='',
  return='plot',messages=T,
  bf.message = F
) 
ggsave('figure04_E2F6_TOP2A/tcga_e2f6_top2a_subtype_cluster.pdf',p,width = 5,height = 4,units = 'in')
save(meta,file = 'figure04_E2F6_TOP2A/tcga_cluster_clinialinfo.Rdata')



rm(list = ls())
load('figure02_E2F6/icgc_lihc_sur_data.Rdata')
identical(colnames(exprSet),rownames(meta))
a=c('E2F6','TOP2A')
exp=exprSet[a,]
exp=as.matrix(exp)
df <-  sweep(exp,1, apply(exp,1,median,na.rm=T))
library(ConsensusClusterPlus)
maxK <-  6 #最多分成几组
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 500,
                                 pItem = 0.8,
                                 pFeature = 1,
                                 clusterAlg = "hc",
                                 seed = 10086,
                                 title="test1",
                                 innerLinkage="complete",
                                 plot="pdf")

icl = calcICL(results,
              title="test1",
              plot="pdf")
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK

for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK
#聚类结果
library(tidyverse)
table(results[[optK]]$consensusClass)
Cluster = results[[optK]]$consensusClass
identical(names(Cluster),rownames(meta))
meta$Cluster = Cluster
table(meta$Cluster)
meta$Cluster=ifelse(meta$Cluster==1,'Cluster1','Cluster2')
meta$Cluster=as.factor(meta$Cluster)
identical(colnames(exp),rownames(meta))
Cluster=meta$Cluster
draw_pca(exp,Cluster,addEllipses = T,color = c('#6d9cc5','#e7857b'))
ggsave('figure04_E2F6_TOP2A/icgc_e2f6_top2a_subtype_pca.pdf',height = 4,width = 5,dpi = 320,units = 'in')

draw_boxplot(exp,Cluster,xlab = ' ',ylab = 'Expression',sort = F,color = c('#6d9cc5','#e7857b'))
ggsave('figure04_E2F6_TOP2A/icgc_e2f6_top2a_subtype_boxplot.pdf',height = 4,width = 4,dpi=320,units = 'in')
colnames(meta)=c("Gender","Age","Event","Stage","Time","Cluster")
library(survival)
library(survminer)

sfit1 <- survfit(Surv(Time, Event) ~ Cluster,
                 data = meta)
g1=ggsurvplot(fit = sfit1,
              palette = c('#6d9cc5','#e7857b'),
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
ggsave('figure04_E2F6_TOP2A/icgc_subtype_e2f6_top2a_survive.pdf',p_sur1,width = 5,height = 4,units = 'in')
data=as.data.frame(t(exp))
data$Group_E2F6=ifelse(data$E2F6>median(data$E2F6),'High','Low')
data$Group_TOP2A=ifelse(data$TOP2A>median(data$TOP2A),'High','Low')
data$Group_E2F6_TOP2A=ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='High','E2F6_H&TOP2A_H',
                             ifelse(data$Group_E2F6=='High'&data$Group_TOP2A=='Low','E2F6_H&TOP2A_L',
                                    ifelse(data$Group_E2F6=='Low'&data$Group_TOP2A=='High','E2F6_L&TOP2A_H','E2F6_L&TOP2A_L')))
data$Group_E2F6_TOP2A=as.factor(data$Group_E2F6_TOP2A)
table(data$Group_E2F6_TOP2A)
identical(rownames(data),rownames(meta))
data$Cluster=meta$Cluster
library(dplyr)

p=ggstatsplot::ggbarstats(
  data = data,
  x = Group_E2F6_TOP2A,
  y = Cluster,
  ggtheme = ggplot2::theme_classic(),
  palette = "category10_d3", 
  package = "ggsci", 
  k = 3, 
  perc.k = 1 ,
  results.subtitle = T,xlab ='',
  return='plot',messages=T,
  bf.message = F
) 
ggsave('figure04_E2F6_TOP2A/icgc_e2f6_top2a_subtype_cluster.pdf',p,width = 5,height = 4,units = 'in')
save(meta,file = 'figure04_E2F6_TOP2A/icgc_cluster_clinialinfo.Rdata')
