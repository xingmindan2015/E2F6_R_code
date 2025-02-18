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
exp=exp[,str_sub(rownames(data),14,15)!='02']
Group=ifelse(str_sub(rownames(data),14,15)=='01','Tumor','Non-tumor')
Group=factor(Group,levels = c('Non-tumor','Tumor'))
data$Group=Group
colnames(data)
b1 = ggplot(data = data,aes(Group,E2F6,fill=Group ))+
  geom_boxplot()+
  #scale_fill_jco()+
  scale_fill_manual(values = c('#6d9cc5','#e7857b'))+
  stat_compare_means(label = 'p.format')+
  theme_bw()+
  xlab('')
ggsave('figure02_E2F6/E2F6_boxplot_normal_tumor.pdf',b1,width = 3,height = 3,units = 'in')
b2 = ggplot(data = data,aes(Group,TOP2A,fill=Group ))+
  geom_boxplot()+
  scale_fill_manual(values = c('#6d9cc5','#e7857b'))+
  stat_compare_means(label = 'p.format')+
  theme_bw()+
  xlab('')
ggsave('figure03_TOP2A/TOP2A_boxplot_normal_tumor.pdf',b2,width = 3,height = 3,units = 'in')
library(lattice)
library(MASS)
colnames(data)
histogram(data$TOP2A)
histogram(data$E2F6)
f1=ggplot(data=data, aes(x=TOP2A, y=E2F6))+
  geom_point(color="#e7857b")+
  stat_smooth(method="lm",se=T)+
  stat_cor(data=data, method = 'spearman')+
  theme_bw()
ggsave('figure02_E2F6/point_line_E2F6_top2a.pdf',f1,width = 3,height = 3,units = 'in',dpi = 320)
data$Group1=ifelse(data$Group=='Tumor',1,0)
write.csv(data,file = 'figure02_E2F6/tcga_e2f6_top2a_normal_tumor_roc.csv')

rm(list = ls())
load("/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/data_symbols_lihc.Rdata")
load('/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/clinical_lihc_postclear.Rdata')

library(tinyarray)
library(tidyverse)
exp=count1
exp[1:4,1:4]
exp=exp[,str_sub(colnames(exp),14,15)=='01']
exp=log2(cpm(exp)+1)
dim(exp)
colnames(exp)=str_sub(colnames(exp),1,12)
exprSet=exp
# 以病人为中心，表达矩阵按病人ID去重复
exprSet = exprSet[,colnames(exprSet)%in%rownames(meta)]
exprSet=as.data.frame(exprSet)
meta=meta[rownames(meta)%in%colnames(exprSet),]
#调整meta的ID顺序与exprSet列名一致
meta=meta[match(colnames(exprSet),meta$ID),]
identical(meta$ID,colnames(exprSet))
save(meta,exprSet,file = "figure02_E2F6/lihc_sur_data.Rdata")


#E2F6的生存分析
identical(colnames(exprSet),rownames(meta))
a=c('E2F6','TOP2A')
exp=exprSet[a,]

data=as.data.frame(t(exp))

meta$E2F6 = ifelse(data$E2F6> median(data$E2F6),'High','Low')
identical(rownames(meta),rownames(data))
sfit1=survfit(Surv(Time, Event)~E2F6, data=meta)
g1=ggsurvplot(fit = sfit1,palette = c("#e7857b", '#6d9cc5'),
              risk.table =TRUE,
              pval =TRUE,
              conf.int =TRUE,
              ncensor.plot = TRUE,data = meta)
library(patchwork)
p1=g1$plot
p2=g1$table
p_sur1=p1 /p2  + plot_layout(design = 'A
                         A
                         A
                         A
                         B')
ggsave('figure02_E2F6/survival_E2F6.pdf',p_sur1,width = 6,height = 5,units = 'in')

meta$TOP2A = ifelse(data$TOP2A > median(data$TOP2A),'High','Low')
identical(rownames(meta),rownames(data))
sfit1=survfit(Surv(Time, Event)~TOP2A, data=meta)
g1=ggsurvplot(fit = sfit1,palette = c("#e7857b", '#6d9cc5'),
              risk.table =TRUE,
              pval =TRUE,
              conf.int =TRUE,
              ncensor.plot = TRUE,data = meta)
library(patchwork)
p1=g1$plot
p2=g1$table
p_sur1=p1 /p2  + plot_layout(design = 'A
                         A
                         A
                         A
                         B')
ggsave('figure03_TOP2A/survival_TOP2A.pdf',p_sur1,width = 6,height = 5,units = 'in')


#E2F6在肿瘤分期中的差异
identical(rownames(meta),rownames(data))
data$Group=meta$Stage
b1=ggboxplot(data= data,x = "Group",y = "E2F6",color = "Group",add = "jitter")+
  stat_compare_means(comparisons = list(c("I","II"),
                                        c("I","III"),
                                        c("I","IV"),
                                        c('II','III'),
                                        c("II","IV"),
                                        c("III","IV")))+
  scale_color_manual(values = c('#6d9cc5','#a5d3ed','#b5aad4',"#e7857b"))+
  stat_compare_means(label.y = 8)+xlab('')
ggsave('figure02_E2F6/boxplot_e2f6_stage.pdf',b1,width = 4,height = 4,units = 'in')

b2=ggboxplot(data= data,x = "Group",y = "TOP2A",color = "Group",add = "jitter")+
  stat_compare_means(comparisons = list(c("I","II"),
                                        c("I","III"),
                                        c("I","IV"),
                                        c('II','III'),
                                        c("II","IV"),
                                        c("III","IV")))+
  scale_color_manual(values = c('#6d9cc5','#a5d3ed','#b5aad4',"#e7857b"))+
  stat_compare_means(label.y = 15)+xlab('')
ggsave('figure03_TOP2A/boxplot_top2a_stage.pdf',b2,width = 4,height = 4,units = 'in')
#绘制时间依赖性ROC曲线
###绘制timeROC
colnames(meta)
identical(rownames(meta),rownames(data))
dat_time = meta[,c(10,12)]
identical(rownames(dat_time),rownames(data))
dat_time$E2F6=data$E2F6
dat_time$TOP2A=data$TOP2A
head(dat_time)  
table(dat_time$Event)
library(survminer)
library(survival)
library(timeROC)
colnames(dat_time)

result <-with(dat_time, timeROC(T=Time,
                                delta=Event,
                                marker=E2F6,
                                cause=1,
                                times=c(12,24,36),
                                iid = TRUE))



plot_dat = data.frame(fpr = as.numeric(result$FP),
                      tpr = as.numeric(result$TP),
                      time = rep(as.factor(c(12,24,36)),each = nrow(result$TP)))
library(ggplot2)
p_time_roc2=ggplot() + 
  geom_line(data = plot_dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#6d9cc5", "#e7857a","#a5d3ed"),
                     labels = paste0("AUC of ",c(1,2,3),"-y survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
dev.off()
p_time_roc2
ggsave('figure02_E2F6/tcga_e2f6_timeroc.pdf',width = 5,height = 5,dpi = 320,units = 'in')

result <-with(dat_time, timeROC(T=Time,
                                delta=Event,
                                marker=TOP2A,
                                cause=1,
                                times=c(12,24,36),
                                iid = TRUE))



plot_dat = data.frame(fpr = as.numeric(result$FP),
                      tpr = as.numeric(result$TP),
                      time = rep(as.factor(c(12,24,36)),each = nrow(result$TP)))
library(ggplot2)
p_time_roc2=ggplot() + 
  geom_line(data = plot_dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#6d9cc5", "#e7857a","#a5d3ed"),
                     labels = paste0("AUC of ",c(1,2,3),"-y survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
dev.off()
p_time_roc2
ggsave('figure03_TOP2A/tcga_top2a_timeroc.pdf',width = 5,height = 5,dpi = 320,units = 'in')


m_E2F6=coxph(Surv(Time, Event) ~ E2F6, data =  dat_time)
beta <- coef(m_E2F6)
se <- sqrt(diag(vcov(m_E2F6)))
HR <- exp(beta)
HRse <- HR * se
tmp_E2F6 <- round(cbind(coef = beta, 
                       se = se, z = beta/se, 
                       p = 1 - pchisq((beta/se)^2, 1),
                       HR = HR, HRse = HRse,
                       HRz = (HR - 1) / HRse, 
                       HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                       HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                       HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
colnames(dat_time)

m_TOP2A=coxph(Surv(Time, Event) ~ TOP2A, data =  dat_time)
beta <- coef(m_TOP2A)
se <- sqrt(diag(vcov(m_TOP2A)))
HR <- exp(beta)
HRse <- HR * se
tmp_TOP2A <- round(cbind(coef = beta, 
                         se = se, z = beta/se, 
                         p = 1 - pchisq((beta/se)^2, 1),
                         HR = HR, HRse = HRse,
                         HRz = (HR - 1) / HRse, 
                         HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                         HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                         HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)

identical(colnames(tmp_E2F6),colnames(tmp_TOP2A))
tmp_uni_E2F6_TOP2A=rbind(tmp_E2F6,tmp_TOP2A)

rownames(tmp_uni_E2F6_TOP2A)
write.csv(tmp_uni_E2F6_TOP2A,file = 'figure02_E2F6/tcga_unicox_e2f6_top2a.csv')



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
b1 = ggplot(data = data,aes(Group,E2F6,fill=Group ))+
  geom_boxplot()+
  scale_fill_manual(values = c('#6d9cc5','#e7857a'))+
  stat_compare_means(label = 'p.format')+
  theme_bw()+
  xlab('')
ggsave('figure02_E2F6/ICGC_E2F6_boxplot_normal_tumor.pdf',b1,width = 3,height = 3,units = 'in')
b2 = ggplot(data = data,aes(Group,TOP2A,fill=Group ))+
  geom_boxplot()+
  #scale_fill_jco()+
  scale_fill_manual(values = c('#6d9cc5','#e7857a'))+
  stat_compare_means(label = 'p.format')+
  theme_bw()+
  xlab('')
ggsave('figure03_TOP2A/ICGC_TOP2A_boxplot_normal_tumor.pdf',b2,width = 3,height = 3,units = 'in')

library(lattice)
library(MASS)
colnames(data)
histogram(data$TOP2A)
histogram(data$E2F6)
f1=ggplot(data=data, aes(x=TOP2A, y=E2F6))+
  geom_point(color="#e7857a")+
  stat_smooth(method="lm",se=T)+stat_cor(data=data, method = 'spearman')+
  theme_bw()
ggsave('figure02_E2F6/icgc_point_line_E2F6_top2a.pdf',f1,width = 3,height = 3,units = 'in',dpi = 320)

data$Group1=ifelse(data$Group=='Tumor',1,0)
write.csv(data,file = 'figure02_E2F6/icgc_e2f6_top2a_normal_tumor_roc.csv')


colnames(exp_t)=str_split(colnames(exp_t),'-',simplify = T)[,1]
identical(rownames(cl_tumor),colnames(exp_t))
exp=exp_t
exp[1:4,1:4]
exp=log2(cpm(exp)+1)
dim(exp)
exprSet=exp

#调整meta的ID顺序与exprSet列名一致
meta=cl_tumor
identical(rownames(meta),colnames(exprSet))
save(meta,exprSet,file = "figure02_E2F6/icgc_lihc_sur_data.Rdata")


#E2F6的生存分析
identical(colnames(exprSet),rownames(meta))

a=c('E2F6','TOP2A')
exp=exprSet[a,]

data=as.data.frame(t(exp))

meta$E2F6 = ifelse(data$E2F6> median(data$E2F6),'High','Low')
identical(rownames(meta),rownames(data))
colnames(meta)=c("Gender","Age","Event","Stage","Time","E2F6" )
sfit1=survfit(Surv(Time, Event)~E2F6, data=meta)
g1=ggsurvplot(fit = sfit1,palette = c('#e7857a', '#6d9cc5'),
              risk.table =TRUE,
              pval =TRUE,
              conf.int =TRUE,
              ncensor.plot = TRUE,data = meta)
library(patchwork)
p1=g1$plot
p2=g1$table
p_sur1=p1 /p2  + plot_layout(design = 'A
                         A
                         A
                         A
                         B')
ggsave('figure02_E2F6/icgc_survival_E2F6.pdf',p_sur1,width = 6,height = 5,units = 'in')

meta$TOP2A = ifelse(data$TOP2A > median(data$TOP2A),'High','Low')
identical(rownames(meta),rownames(data))
sfit1=survfit(Surv(Time, Event)~TOP2A, data=meta)
g1=ggsurvplot(fit = sfit1,palette = c('#e7857a', '#6d9cc5'),
              risk.table =TRUE,
              pval =TRUE,
              conf.int =TRUE,
              ncensor.plot = TRUE,data = meta)
library(patchwork)
p1=g1$plot
p2=g1$table
p_sur1=p1 /p2  + plot_layout(design = 'A
                         A
                         A
                         A
                         B')
ggsave('figure03_TOP2A/icgc_survival_TOP2A.pdf',p_sur1,width = 6,height = 5,units = 'in')


#E2F6在肿瘤分期中的差异
identical(rownames(meta),rownames(data))
table(meta$Stage)
data$Group=ifelse(meta$Stage==1,'I',
                  ifelse(meta$Stage==2,'II',
                         ifelse(meta$Stage==3,'III','IV')))
b1=ggboxplot(data= data,x = "Group",y = "E2F6",color = "Group",add = "jitter")+
  stat_compare_means(comparisons = list(c("I","II"),
                                        c("I","III"),
                                        c("I","IV"),
                                        c('II','III'),
                                        c("II","IV"),
                                        c("III","IV")))+
  scale_color_manual(values = c('#6d9cc5','#a5d3ed','#b5aad4',"#e7857b"))+
  stat_compare_means(label.y = 8)+xlab('')
ggsave('figure02_E2F6/icgc_boxplot_e2f6_stage.pdf',b1,width = 4,height = 4,units = 'in')

b2=ggboxplot(data= data,x = "Group",y = "TOP2A",color = "Group",add = "jitter")+
  stat_compare_means(comparisons = list(c("I","II"),
                                        c("I","III"),
                                        c("I","IV"),
                                        c('II','III'),
                                        c("II","IV"),
                                        c("III","IV")))+
  scale_color_manual(values = c('#6d9cc5','#a5d3ed','#b5aad4',"#e7857b"))+
  stat_compare_means(label.y = 15)+xlab('')
ggsave('figure03_TOP2A/icgc_boxplot_top2a_stage.pdf',b2,width = 4,height = 4,units = 'in')


#绘制时间依赖性ROC曲线
###绘制timeROC
colnames(meta)
identical(rownames(meta),rownames(data))
dat_time = meta[,c(3,5)]
identical(rownames(dat_time),rownames(data))
dat_time$E2F6=data$E2F6
dat_time$TOP2A=data$TOP2A
head(dat_time)  
table(dat_time$Event)
library(survminer)
library(survival)
library(timeROC)
colnames(dat_time)

result <-with(dat_time, timeROC(T=Time,
                                delta=Event,
                                marker=E2F6,
                                cause=1,
                                times=c(12,24,36),
                                iid = TRUE))



plot_dat = data.frame(fpr = as.numeric(result$FP),
                      tpr = as.numeric(result$TP),
                      time = rep(as.factor(c(12,24,36)),each = nrow(result$TP)))
library(ggplot2)
p_time_roc2=ggplot() + 
  geom_line(data = plot_dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c(c("#6d9cc5", "#e7857a","#a5d3ed")),
                     labels = paste0("AUC of ",c(1,2,3),"-y survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
dev.off()
p_time_roc2
ggsave('figure02_E2F6/icgc_e2f6_timeroc.pdf',width = 6,height = 5,dpi = 320,units = 'in')

result <-with(dat_time, timeROC(T=Time,
                                delta=Event,
                                marker=TOP2A,
                                cause=1,
                                times=c(12,24,36),
                                iid = TRUE))



plot_dat = data.frame(fpr = as.numeric(result$FP),
                      tpr = as.numeric(result$TP),
                      time = rep(as.factor(c(12,24,36)),each = nrow(result$TP)))
library(ggplot2)
p_time_roc2=ggplot() + 
  geom_line(data = plot_dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#6d9cc5", "#e7857a","#a5d3ed"),
                     labels = paste0("AUC of ",c(1,2,3),"-y survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
dev.off()
p_time_roc2
ggsave('figure03_TOP2A/icgc_top2a_timeroc.pdf',width = 6,height = 5,dpi = 320,units = 'in')

m_E2F6=coxph(Surv(Time, Event) ~ E2F6, data =  dat_time)
beta <- coef(m_E2F6)
se <- sqrt(diag(vcov(m_E2F6)))
HR <- exp(beta)
HRse <- HR * se
tmp_E2F6 <- round(cbind(coef = beta, 
                        se = se, z = beta/se, 
                        p = 1 - pchisq((beta/se)^2, 1),
                        HR = HR, HRse = HRse,
                        HRz = (HR - 1) / HRse, 
                        HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                        HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                        HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
colnames(dat_time)

m_TOP2A=coxph(Surv(Time, Event) ~ TOP2A, data =  dat_time)
beta <- coef(m_TOP2A)
se <- sqrt(diag(vcov(m_TOP2A)))
HR <- exp(beta)
HRse <- HR * se
tmp_TOP2A <- round(cbind(coef = beta, 
                         se = se, z = beta/se, 
                         p = 1 - pchisq((beta/se)^2, 1),
                         HR = HR, HRse = HRse,
                         HRz = (HR - 1) / HRse, 
                         HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                         HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                         HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)

identical(colnames(tmp_E2F6),colnames(tmp_TOP2A))
tmp_uni_E2F6_TOP2A=rbind(tmp_E2F6,tmp_TOP2A)

rownames(tmp_uni_E2F6_TOP2A)
write.csv(tmp_uni_E2F6_TOP2A,file = 'figure02_E2F6/icgc_unicox_e2f6_top2a.csv')


# rm(list = ls())
# load('database/GEO/GSE14520/exp1_pd1.Rdata')
# a=c('E2F6','TOP2A')
# exp=exp1
# meta=pd1
# data=as.data.frame(t(exp[a,]))
# identical(rownames(data),rownames(meta))
# 
# data$Group=meta$Group
# colnames(data)
# b1 = ggplot(data = data,aes(Group,E2F6,fill=Group ))+
#   geom_boxplot()+
#   scale_fill_jco()+
#   stat_compare_means(label = 'p.format')+
#   theme_bw()+
#   xlab('')
# ggsave('figure02_E2F6/14520_01_E2F6_boxplot_normal_tumor.pdf',b1,width = 4,height = 3,units = 'in')
# b2 = ggplot(data = data,aes(Group,TOP2A,fill=Group ))+
#   geom_boxplot()+
#   scale_fill_jco()+
#   stat_compare_means(label = 'p.format')+
#   theme_bw()+
#   xlab('')
# ggsave('figure03_TOP2A/14520_01_TOP2A_boxplot_normal_tumor.pdf',b2,width = 4,height = 3,units = 'in')
# 
# library(lattice)
# library(MASS)
# colnames(data)
# histogram(data$TOP2A)
# histogram(data$E2F6)
# f1=ggplot(data=data, aes(x=TOP2A, y=E2F6))+
#   geom_point(color="red")+
#   stat_smooth(method="lm",se=T)+stat_cor(data=data, method = 'spearman')+
#   theme_bw()
# ggsave('figure02_E2F6/14520_01_point_line_E2F6_top2a.pdf',f1,width = 3,height = 3,units = 'in',dpi = 320)
# 
# data$Group1=ifelse(data$Group=='Tumor',1,0)
# write.csv(data,file = 'figure02_E2F6/14520_01_e2f6_top2a_normal_tumor_roc.csv')
# 
# 
# 
# rm(list = ls())
# load('database/GEO/GSE14520/exp2_pd2.Rdata')
# a=c('E2F6','TOP2A')
# exp=exp2
# meta=pd2
# data=as.data.frame(t(exp[a,]))
# identical(rownames(data),rownames(meta))
# 
# data$Group=meta$Group
# colnames(data)
# b1 = ggplot(data = data,aes(Group,E2F6,fill=Group ))+
#   geom_boxplot()+
#   scale_fill_jco()+
#   stat_compare_means(label = 'p.format')+
#   theme_bw()+
#   xlab('')
# ggsave('figure02_E2F6/14520_02_E2F6_boxplot_normal_tumor.pdf',b1,width = 4,height = 3,units = 'in')
# b2 = ggplot(data = data,aes(Group,TOP2A,fill=Group ))+
#   geom_boxplot()+
#   scale_fill_jco()+
#   stat_compare_means(label = 'p.format')+
#   theme_bw()+
#   xlab('')
# ggsave('figure03_TOP2A/14520_02_TOP2A_boxplot_normal_tumor.pdf',b2,width = 4,height = 3,units = 'in')
# f1=ggplot(data=data, aes(x=TOP2A, y=E2F6))+
#   geom_point(color="red")+
#   stat_smooth(method="lm",se=T)+stat_cor(data=data, method = 'spearman')+
#   theme_bw()
# ggsave('figure02_E2F6/14520_02_point_line_E2F6_top2a.pdf',f1,width = 3,height = 3,units = 'in',dpi = 320)
# data$Group1=ifelse(data$Group=='Tumor',1,0)
# write.csv(data,file = 'figure02_E2F6/14520_02_e2f6_top2a_normal_tumor_roc.csv')
# 
# 
# rm(list = ls())
# load('database/GEO/GSE20140/exp1_pd1.Rdata')
# a=c('E2F6','TOP2A')
# exp=exp1
# meta=pd1
# data=as.data.frame(t(exp[a,]))
# identical(rownames(data),rownames(meta))
# 
# data$Group=meta$Group
# colnames(data)
# b1 = ggplot(data = data,aes(Group,E2F6,fill=Group ))+
#   geom_boxplot()+
#   scale_fill_jco()+
#   stat_compare_means(label = 'p.format')+
#   theme_bw()+
#   xlab('')
# ggsave('figure02_E2F6/20140_01_E2F6_boxplot_normal_tumor.pdf',b1,width = 4,height = 3,units = 'in')
# b2 = ggplot(data = data,aes(Group,TOP2A,fill=Group ))+
#   geom_boxplot()+
#   scale_fill_jco()+
#   stat_compare_means(label = 'p.format')+
#   theme_bw()+
#   xlab('')
# ggsave('figure03_TOP2A/20140_01_TOP2A_boxplot_normal_tumor.pdf',b2,width = 4,height = 3,units = 'in')
# f1=ggplot(data=data, aes(x=TOP2A, y=E2F6))+
#   geom_point(color="red")+
#   stat_smooth(method="lm",se=T)+stat_cor(data=data, method = 'spearman')+
#   theme_bw()
# ggsave('figure02_E2F6/20140_01_point_line_E2F6_top2a.pdf',f1,width = 3,height = 3,units = 'in',dpi = 320)
# data$Group1=ifelse(data$Group=='Tumor',1,0)
# write.csv(data,file = 'figure02_E2F6/20140_01_e2f6_top2a_normal_tumor_roc.csv')
# 
# 
# table(meta$Group)
# table(meta$Location)
# meta=meta[meta$Location!='NA',]
# data=data[rownames(data)%in%rownames(meta),]
# identical(rownames(data),rownames(meta))
# data$Location=meta$Location
# data$Location=ifelse(data$Location=='tumor peripherial','Tumor_peripherial','Tumor_central')
# b1 = ggplot(data = data,aes(Location,E2F6,fill=Location ))+
#   geom_boxplot()+
#   scale_fill_jco()+
#   stat_compare_means(label = 'p.format')+
#   theme_bw()+
#   xlab('')
# ggsave('figure02_E2F6/20140_01_E2F6_boxplot_cen_peri.pdf',b1,width = 4,height = 3,units = 'in')
# b2 = ggplot(data = data,aes(Location,TOP2A,fill=Location ))+
#   geom_boxplot()+
#   scale_fill_jco()+
#   stat_compare_means(label = 'p.format')+
#   theme_bw()+
#   xlab('')
# ggsave('figure03_TOP2A/20140_01_TOP2A_boxplot_cen_peri.pdf',b2,width = 4,height = 3,units = 'in')
# 
# rm(list = ls())
# load('database/GEO/GSE20140/exp3_pd3.Rdata')
# a=c('E2F6','TOP2A')
# exp=exp3
# meta=pd3
# table(rownames(exp)%in%'TOP2A')
# class(exp)
# data=as.data.frame(t(exp[a,]))
# 
# identical(rownames(data),rownames(meta))
# colnames(meta)
# data$Vascular_invasion=meta$Vascular_invasion
# colnames(data)
# b1 = ggplot(data = data,aes(Vascular_invasion,E2F6,fill=Vascular_invasion ))+
#   geom_boxplot()+
#   scale_fill_jco()+
#   stat_compare_means(label = 'p.format')+
#   theme_bw()+
#   xlab('')
# ggsave('figure02_E2F6/20140_03_E2F6_boxplot_Vascular_invasion.pdf',b1,width = 4,height = 3,units = 'in')
# b2 = ggplot(data = data,aes(Vascular_invasion,TOP2A,fill=Vascular_invasion ))+
#   geom_boxplot()+
#   scale_fill_jco()+
#   stat_compare_means(label = 'p.format')+
#   theme_bw()+
#   xlab('')
# ggsave('figure03_TOP2A/20140_03_TOP2A_boxplot_Vascular_invasion.pdf',b2,width = 4,height = 3,units = 'in')
# f1=ggplot(data=data, aes(x=TOP2A, y=E2F6))+
#   geom_point(color="red")+
#   stat_smooth(method="lm",se=T)+stat_cor(data=data, method = 'spearman')+
#   theme_bw()
# ggsave('figure02_E2F6/20140_03_point_line_E2F6_top2a.pdf',f1,width = 3,height = 3,units = 'in',dpi = 320)

rm(list = ls())
a=read.csv('figure02_E2F6/tcga_icgc_unicox_e2f6_top2a.csv')
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
pdf('figure02_E2F6/tcga_icgc_e2f6_top2a_uniforest_clinical.pdf',width = 6,height = 4)
forestplot(
  tmp[, c(1, 2, 6,7)],
  mean = tmp[, 3],
  lower = tmp[, 4],
  upper = tmp[, 5],
  zero = 1,
  boxsize = 0.1,
  col = fpColors(box = '#6d9cc5', lines = 'black', zero = 'grey'),
  lty.ci = "solid",
  graph.pos = 3,
  #xticks = F,
  is.summary = c(T, rep(F, 13)),
  align = "l",
  hrzl_lines = list(
    "1" = gpar(lty=1.5),
    "2" = gpar(lty=1.5),
    "8"= gpar(lty=1.5)),
  colgap = unit(5, 'mm')
)
dev.off()
