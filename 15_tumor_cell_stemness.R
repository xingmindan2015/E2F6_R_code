rm(list=ls())
mR=read.csv('figure15_stemness_index/mRNAsi_tcga.csv',row.names = 1)
class(mR$cancer.type)
mR1=mR[mR$cancer.type%in%'LIHC',]
a=rownames(mR1)[str_sub(rownames(mR1),14,15)=='01']
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')

mR1=mR1[str_sub(rownames(mR1),1,12)%in%rownames(meta),]
mR1$ID = str_sub(rownames(mR1),1,12) 
mR1 = mR1[!duplicated(str_sub(rownames(mR1),1,12)),]
rownames(mR1) = mR1$ID
dat = merge(mR1,meta,by = 'ID')

#检测是否符合正态分布
library(lattice)
library(MASS)
histogram(dat$mRNAsi)
histogram(dat$Predict_score)#不符合正态分布
###所以用spearman相关分析进行分析
library(ggplot2)
library(ggpubr)
library(ggpubr)
ggplot(data=dat, aes(x=Predict_score, y=mRNAsi))+
  geom_point(color="#e7857b")+
  stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat, method = 'spearman')+
  theme_bw()
#stat_cor(data=dat, method = "pearson")意为用pearson相关进行相关性分析，可以自行更改方法
ggsave('figure15_stemness_index/tcga_mRNAsi_riskscore.pdf',height = 4,width = 4,units = 'in')
colnames(dat)
ordercolors<-c("#6d9cc5","#e7857b")
library(gghalves)
b1 = ggplot(data = dat,
            aes(x=Group, y=mRNAsi, fill=Group)) +
  geom_half_violin(side = "r", color=NA, alpha=0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width=0.2, linewidth=0.5) +
  geom_half_point_panel(side = "l", shape=21, size=3, color="white") +
  scale_fill_manual(values = ordercolors) +
  stat_compare_means(label = 'p.format')+
  labs(x=NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black"))

# b1 = ggviolin(dat = dat,x='Group',y='mRNAsi',fill = 'Group',add = 'boxplot',add.params = list(fill='white'))+
#   scale_fill_manual(values = c("#6d9cc5","#e7857b"))+
#   stat_compare_means(label.y = 0.65)+
#   theme_bw()+
#   labs(x=' ')+
#   theme(legend.position = "bottom")
ggsave('figure15_stemness_index/tcga_mRNAsi_rigroup.pdf',b1,height = 4,width = 4,units = 'in')
b2 = ggplot(data = dat,
            aes(x=Cluster, y=mRNAsi, fill=Cluster)) +
  geom_half_violin(side = "r", color=NA, alpha=0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width=0.2, linewidth=0.5) +
  geom_half_point_panel(side = "l", shape=21, size=3, color="white") +
  scale_fill_manual(values = ordercolors) +
  stat_compare_means(label = 'p.format')+
  labs(x=NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black"))
# b2 = ggviolin(dat = dat,x='Cluster',y='mRNAsi',fill = 'Cluster',add = 'boxplot',add.params = list(fill='white'))+
#   scale_fill_manual(values = c("#6d9cc5","#e7857b"))+
#   stat_compare_means(label.y = 0.65)+
#   theme_bw()+
#   labs(x=' ')+
#   theme(legend.position = "bottom")
ggsave('figure15_stemness_index/tcga_mRNAsi_cluster.pdf',b2,height = 4,width = 4,units = 'in')




###mDNAsi
rm(list = ls())
mR=read.csv('figure15_stemness_index/mDNAsi_tcga.csv',row.names = 1)
class(mR$cancer.type)
mR1=mR[mR$cancer.type%in%'LIHC',]
a=rownames(mR1)[str_sub(rownames(mR1),14,15)=='01']
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')

mR1=mR1[str_sub(rownames(mR1),1,12)%in%rownames(meta),]
mR1$ID = str_sub(rownames(mR1),1,12) 
mR1 = mR1[!duplicated(str_sub(rownames(mR1),1,12)),]
rownames(mR1) = mR1$ID
dat = merge(mR1,meta,by = 'ID')

#检测是否符合正态分布
library(lattice)
library(MASS)
histogram(dat$mDNAsi)
histogram(dat$Predict_score)#不符合正态分布
###所以用spearman相关分析进行分析
library(ggplot2)
library(ggpubr)
library(ggpubr)
ggplot(data=dat, aes(x=Predict_score, y=mDNAsi))+
  geom_point(color="#e7857b")+
  stat_smooth(method="lm",se=FALSE)+stat_cor(data=dat, method = 'spearman')+
  theme_bw()
#stat_cor(data=dat, method = "pearson")意为用pearson相关进行相关性分析，可以自行更改方法
ggsave('figure15_stemness_index/tcga_mDNAsi_riskscore.pdf',height = 4,width = 4,units = 'in')
colnames(dat)
ordercolors<-c("#6d9cc5","#e7857b")
library(gghalves)
b1 = ggplot(data = dat,
            aes(x=Group, y=mDNAsi, fill=Group)) +
  geom_half_violin(side = "r", color=NA, alpha=0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width=0.2, linewidth=0.5) +
  geom_half_point_panel(side = "l", shape=21, size=3, color="white") +
  scale_fill_manual(values = ordercolors) +
  stat_compare_means(label = 'p.format')+
  labs(x=NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black"))
# b1 = ggviolin(dat = dat,x='Group',y='mDNAsi',fill = 'Group',add = 'boxplot',add.params = list(fill='white'))+
#   scale_fill_manual(values = c("#6d9cc5","#e7857b"))+
#   stat_compare_means(label.y = 0.65)+
#   theme_bw()+
#   labs(x=' ')+
#   theme(legend.position = "bottom")
ggsave('figure15_stemness_index/tcga_mDNAsi_rigroup.pdf',b1,height = 4,width = 4,units = 'in')
# b2 = ggviolin(dat = dat,x='Cluster',y='mDNAsi',fill = 'Cluster',add = 'boxplot',add.params = list(fill='white'))+
#   scale_fill_manual(values = c("#6d9cc5","#e7857b"))+
#   stat_compare_means(label.y = 0.65)+
#   theme_bw()+
#   labs(x=' ')+
#   theme(legend.position = "bottom")
b2 = ggplot(data = dat,
            aes(x=Cluster, y=mDNAsi, fill=Cluster)) +
  geom_half_violin(side = "r", color=NA, alpha=0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width=0.2, linewidth=0.5) +
  geom_half_point_panel(side = "l", shape=21, size=3, color="white") +
  scale_fill_manual(values = ordercolors) +
  stat_compare_means(label = 'p.format')+
  labs(x=NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black"))
ggsave('figure15_stemness_index/tcga_mDNAsi_cluter.pdf',b2,height = 4,width = 4,units = 'in')
