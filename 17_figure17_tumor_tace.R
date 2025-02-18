
# GSE109211----------------------------------------------------------

#GSE109211无正常对照，有生存信息，准备数据####
rm(list = ls())
options(stringsAsFactors = F)#3.6版本以后的需要用
library(GEOquery)
gse_number = "GSE109211"
eSet <- getGEO(gse_number, 
               filename = 'database/GEO/GSE109211/GSE109211_series_matrix.txt.gz', 
               getGPL = F)
class(eSet)
length(eSet)
exp = eSet@assayData$exprs
#(1)提取表达矩阵exp

exp[1:4,1:4]#注意看一下，注意不要log两次，取过log的表达量的值在0-20之间
exp1=log2(exp+1)
exp1[1:4,1:4]
boxplot(exp1[,1:40])

#(2)提取临床信息
pd<- eSet@phenoData@data

colnames(pd)
pd1=pd[,c(2,33,35)]
colnames(pd1)
colnames(pd1)=c('ID','Outcome','Treatment')

exp1[1:4,1:4]
exp=as.data.frame(exp1)

identical(pd1$ID,colnames(exp))

#(4)提取芯片平台编号
gpl_number <- eSet@annotation

#(5)下载注释文件
a = getGEO(gpl_number,destdir = "database/GEO/GSE109211/")
b = a@dataTable@table
colnames(b)
ids2 = b[,c("ID","Symbol")]
colnames(ids2) = c("probe_id","symbol")
ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),]

table(rownames(exp)%in%ids2$probe_id)
exp=exp[rownames(exp)%in%ids2$probe_id,]
table(duplicated(ids2$symbol))
ids2=ids2[!duplicated(ids2$symbol),]
table(rownames(exp)%in%ids2$probe_id)
exp=exp[rownames(exp)%in%ids2$probe_id,]
class(ids2$probe_id)
ids2=ids2[ids2$probe_id%in%rownames(exp),]
dim(exp)
exp=exp[match(ids2$probe_id,rownames(exp)),]
identical(ids2$probe_id,rownames(exp))
rownames(exp)=ids2$symbol

identical(colnames(exp),rownames(pd1))
save(pd1,exp,file = 'database/GEO/GSE109211/exp_pd1.Rdata')
#模型验证####
rm(list = ls())
load('figure06_two_group/train_e2f6_top2a_logistic_model_genes.Rdata')
load('database/GEO/GSE109211/exp_pd1.Rdata')
a = c('E2F6','TOP2A')
exp = exp [a,]
test1 = as.data.frame(t(exp))
fp <- predict(model2,newdata = test1)
names(fp)
identical(names(fp),pd1$ID)
pd1$Predict_Score = fp

test1$Group_E2F6=ifelse(test1$E2F6>median(test1$E2F6),'High','Low')
test1$Group_TOP2A=ifelse(test1$TOP2A>median(test1$TOP2A),'High','Low')
test1$Group_E2F6_TOP2A=ifelse(test1$Group_E2F6=='High'&test1$Group_TOP2A=='High','E2F6_H&TOP2A_H','Others')
test1$Group_E2F6_TOP2A=factor(test1$Group_E2F6_TOP2A,levels = c('Others','E2F6_H&TOP2A_H'))
identical(rownames(test1),pd1$ID)
test1$Predict_Score = pd1$Predict_Score
pd1$Group = test1$Group_E2F6_TOP2A
write.csv(test1,file='figure06_two_group/GSE109211_roc_data.csv')
library(rms)
fit1 <- lrm(Group_E2F6_TOP2A ~ E2F6+TOP2A, data = test1, x=TRUE,y=TRUE)
cal1 <- calibrate(fit1, method = "boot", B=1000)
pdf('figure17_tace/GSE109211_cal_curv.pdf',width = 5,height = 5)
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
library(ggDCA)
library(rms)
test1$Group_E2F6_TOP2A=ifelse(test1$Group_E2F6_TOP2A=='Others',0,1)
PredicScore <- rms::lrm(Group_E2F6_TOP2A ~ E2F6+TOP2A, data = test1)
E2F6=rms::lrm(Group_E2F6_TOP2A ~ E2F6 , data = test1)
TOP2A=rms::lrm(Group_E2F6_TOP2A ~  TOP2A , data = test1)
#remove.packages('dcurves')
library(ggDCA)
library(rms)
data  <- dca(PredicScore,E2F6,TOP2A)
ggplot(data)
ggsave('figure17_tace/GSE109211_dca.pdf',height = 5,width = 6,dpi=320,units = 'in')
ordercolors<-c("#6d9cc5","#e7857b")
library(gghalves)
colnames(pd1)
b1 = ggplot(data = pd1,
            aes(x=Outcome, y=Predict_Score, fill=Outcome)) +
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

ggsave('figure17_tace/GSE109211_predictscore_boxplot.pdf',b1,width = 4,height = 4,units = 'in')
b2=ggstatsplot::ggbarstats(
  data = pd1,
  x = Outcome,
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
ggsave('figure17_tace/GSE109211_ration.pdf',b2,width = 4,height = 4,units = 'in')

pd3 = pd1[pd1$Treatment == 'Sor',]
b3=ggstatsplot::ggbarstats(
  data = pd3,
  x = Outcome,
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
ggsave('figure17_tace/GSE109211_ration_sor.pdf',b2,width = 4,height = 4,units = 'in')

b4 = ggplot(data = pd3,
            aes(x=Outcome, y=Predict_Score, fill=Outcome)) +
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

ggsave('figure17_tace/GSE109211_predictscore_boxplot_sor.pdf',b4,width = 4,height = 4,units = 'in')


# GSE104580 ---------------------------------------------------------------
#GSE104580无正常对照，有生存信息,准备数据####
rm(list = ls())
options(stringsAsFactors = F)#3.6版本以后的需要用
library(GEOquery)
gse_number = "GSE104580"
eSet <- getGEO(gse_number, 
               filename = 'database/GEO/GSE104580/GSE104580_series_matrix.txt.gz', 
               getGPL = F)
class(eSet)
length(eSet)
exp = eSet@assayData$exprs
#(1)提取表达矩阵exp

exp[1:4,1:4]#注意看一下，注意不要log两次，取过log的表达量的值在0-20之间
boxplot(exp[,1:100])
exp = normalizeBetweenArrays(exp)
boxplot(exp[,1:100])
#(2)提取临床信息
pd<- eSet@phenoData@data
table(pd$`subject subgroup:ch1`)
pd$Group = str_split(pd$`subject subgroup:ch1`,' ',simplify = T)[,2]
colnames(pd)
pd1=pd[,c(2,38)]
colnames(pd1)
colnames(pd1)=c('ID','Outcome')

exp=as.data.frame(exp)

identical(pd1$ID,colnames(exp))

#(4)提取芯片平台编号
gpl_number <- eSet@annotation
if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)#如果是另外的包，可以用查找全部，然后全部用另外一个包替换即可
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)
ids2 = ids
table(rownames(exp)%in%ids2$probe_id)
exp=exp[rownames(exp)%in%ids2$probe_id,]
table(duplicated(ids2$symbol))
ids2=ids2[!duplicated(ids2$symbol),]
table(rownames(exp)%in%ids2$probe_id)
exp=exp[rownames(exp)%in%ids2$probe_id,]
class(ids2$probe_id)
ids2=ids2[ids2$probe_id%in%rownames(exp),]
dim(exp)
exp=exp[match(ids2$probe_id,rownames(exp)),]
identical(ids2$probe_id,rownames(exp))
rownames(exp)=ids2$symbol

identical(colnames(exp),rownames(pd1))
save(pd1,exp,file = 'database/GEO/GSE104580/exp_pd1.Rdata')

#模型验证####
rm(list = ls())
load('figure06_two_group/train_e2f6_top2a_logistic_model_genes.Rdata')
load('database/GEO/GSE104580/exp_pd1.Rdata')
a = c('E2F6','TOP2A')
exp = exp [a,]
test1 = as.data.frame(t(exp))
fp <- predict(model2,newdata = test1)
names(fp)
identical(names(fp),pd1$ID)
pd1$Predict_Score = fp

test1$Group_E2F6=ifelse(test1$E2F6>median(test1$E2F6),'High','Low')
test1$Group_TOP2A=ifelse(test1$TOP2A>median(test1$TOP2A),'High','Low')
test1$Group_E2F6_TOP2A=ifelse(test1$Group_E2F6=='High'&test1$Group_TOP2A=='High','E2F6_H&TOP2A_H','Others')
test1$Group_E2F6_TOP2A=factor(test1$Group_E2F6_TOP2A,levels = c('Others','E2F6_H&TOP2A_H'))
identical(rownames(test1),pd1$ID)
test1$Predict_Score = pd1$Predict_Score
pd1$Group = test1$Group_E2F6_TOP2A
write.csv(test1,file='figure06_two_group/GSE104580_roc_data.csv')
library(rms)
fit1 <- lrm(Group_E2F6_TOP2A ~ E2F6+TOP2A, data = test1, x=TRUE,y=TRUE)
cal1 <- calibrate(fit1, method = "boot", B=1000)
pdf('figure17_tace/GSE104580_cal_curv.pdf',width = 5,height = 5)
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
library(ggDCA)
library(rms)
test1$Group_E2F6_TOP2A=ifelse(test1$Group_E2F6_TOP2A=='Others',0,1)
PredicScore <- rms::lrm(Group_E2F6_TOP2A ~ E2F6+TOP2A, data = test1)
E2F6=rms::lrm(Group_E2F6_TOP2A ~ E2F6 , data = test1)
TOP2A=rms::lrm(Group_E2F6_TOP2A ~  TOP2A , data = test1)
#remove.packages('dcurves')
library(ggDCA)
library(rms)
data  <- dca(PredicScore,E2F6,TOP2A)
ggplot(data)
ggsave('figure17_tace/GSE104580_dca.pdf',height = 5,width = 6,dpi=320,units = 'in')
ordercolors<-c("#6d9cc5","#e7857b")
library(gghalves)
colnames(pd1)
b1 = ggplot(data = pd1,
            aes(x=Outcome, y=Predict_Score, fill=Outcome)) +
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

ggsave('figure17_tace/GSE104580_predictscore_boxplot.pdf',b1,width = 4,height = 4,units = 'in')
b2=ggstatsplot::ggbarstats(
  data = pd1,
  x = Outcome,
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
ggsave('figure17_tace/GSE104580_ration.pdf',b2,width = 4,height = 4,units = 'in')


