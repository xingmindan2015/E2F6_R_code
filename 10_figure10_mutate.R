##正式工作
rm(list = ls())
library(maftools)
load('figure10_mutate/lihc_mutation_dataframe.Rdata')
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
table(meta$Group)
maf$Tumor_Sample_Barcode[1:4]
table(str_sub(maf$Tumor_Sample_Barcode,14,15))
maf1=maf[str_sub(maf$Tumor_Sample_Barcode,14,15)=='01',]
maf1$Tumor_Sample_Barcode=str_sub(maf1$Tumor_Sample_Barcode,1,12)
colnames(meta)
meta1=meta
length(table(maf1$Tumor_Sample_Barcode))
colnames(meta1)=c( "Tumor_Sample_Barcode","Stage","Last_followup","Age","T_stage","N_stage","M_stage","Race",
                   "Gender","Event","Death","Time", "Age_group","Cluster","Group", "Predict_score")
maf1=maf1[maf1$Tumor_Sample_Barcode%in%meta1$Tumor_Sample_Barcode,]
meta1=meta1[meta1$Tumor_Sample_Barcode%in%maf1$Tumor_Sample_Barcode,]
# HR_name=meta1$Tumor_Sample_Barcode[meta1$Riskgroup=='Highrisk']
# LR_name=meta1$Tumor_Sample_Barcode[meta1$Riskgroup=='Lowrisk']
table(meta1$Tumor_Sample_Barcode%in%maf1$Tumor_Sample_Barcode)
length(table(maf1$Tumor_Sample_Barcode))

maf_clinical <- read.maf(maf = maf1,
                         clinicalData = meta1,
                         isTCGA = TRUE
)
save(maf_clinical,file = 'figure10_mutate/lihc_mutation_clinical.Rdata')

# 展示样本层面的信息
getSampleSummary(maf_clinical)
# 展示基因层面的信息
getGeneSummary(maf_clinical)
# 展示临床相关信息
getClinicalData(maf_clinical)
# 获取 MAF 文件的所有列字段
getFields(maf_clinical)
# 保存为文件
write.mafSummary(maf = maf_clinical, basename = 'figure10_mutate/lihc_mutation_clinical_summary.maf')
library(maftools)
class(meta1$Group)
table(meta1$Group)
par('mar')
par(mar=c(1,1,1,1))
pdf('figure10_mutate/oncoplot_diffriskgroup_top20.pdf',width = 8,height = 8)
oncoplot(maf = maf_clinical,top = 20,legend_height = 6,
         clinicalFeatures = "Group",sortByAnnotation = T,
         draw_titv = T)
dev.off()

###计算TMB值
rm(list = ls())
load('figure10_mutate/lihc_mutation_clinical.Rdata')
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
library(maftools)
x = tmb(maf = maf_clinical)
x[1:4,]
class(x)
x=as.data.frame(x)
rownames(x)=x$Tumor_Sample_Barcode
x=x[rownames(x)%in%rownames(meta),]
meta1=meta[rownames(meta)%in%rownames(x),]
meta1=meta1[match(rownames(x),rownames(meta1)),]
identical(rownames(meta1),rownames(x))
x$Group=meta1$Group
colnames(x)=c('ID','total','TMB','total_perMB_log','Group')
ordercolors<-c("#6d9cc5","#e7857b")
library(gghalves)
p1 = ggplot(data = x,
            aes(x=Group, y=TMB, fill=Group)) +
  geom_half_violin(side = "r", color=NA, alpha=0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width=0.2, linewidth=0.5) +
  geom_half_point_panel(side = "l", shape=21, size=3, color="white") +
  scale_fill_manual(values = ordercolors) +
  stat_compare_means()+
  labs(x=NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black"))

ggsave('figure10_mutate/boxplot_diffrisk.pdf',p1,width = 4,height = 3.5,dpi = 320,units = 'in')

identical(rownames(x),rownames(meta1))
x$Predict_score=meta1$Predict_score
p1=ggplot(data=x, aes(x=Predict_score , y=TMB))+
  geom_point(color="#e7857b")+
  stat_smooth(method="lm",se=T)+
  stat_cor(data=x, method = "spearman")+theme_bw()
p1
ggsave('figure10_mutate/tmb_line_point.pdf',width = 4,height = 3.5,dpi = 320,units = 'in')

meta1$TMB=x$TMB

meta1$TMBgroup=ifelse(meta1$TMB<median(meta1$TMB),'LowTMB','HighTMB')
meta1$TMBgroup=factor(meta1$TMBgroup,levels = c('LowTMB','HighTMB'))
sfit <- survfit(Surv(Time, Event)~TMBgroup, data=meta1)
summary(sfit)$table
p2=ggsurvplot(sfit,palette = c("#6d9cc5","#e7857b"),
              risk.table =TRUE,pval =TRUE,
              conf.int =TRUE,xlab ="Time in months", 
              ggtheme =theme_light(), 
              ncensor.plot = TRUE)

library(patchwork)
p1=p2$plot
p2_2=p2$table
p3=p2$ncensor.plot
p_sur=p1 /p2_2  + plot_layout(design = 'A
                         A
                         A
                         A
                         B')
p_sur
ggsave('figure10_mutate/tmb_survival.pdf',width = 6,height = 5,dpi = 320,units = 'in')
save(meta1,file = 'figure10_mutate/TMB_sur_clinical.Rdata')




###比较高低风险组两个队列
# 从临床数据中提取不同风险组对应的"Tumor_Sample_Barcode"
rm(list = ls())
load('figure10_mutate/TMB_sur_clinical.Rdata')
load('figure10_mutate/lihc_mutation_clinical.Rdata')
table(meta1$Group)
clinO=meta1$ID[meta1$Group=='Others']
clinET=meta1$ID[meta1$Group=='E2F6_H&TOP2A_H']
# 使用subsetMaf构建高风险组和低风阻的MAF对象
lihcO=subsetMaf(maf = maf_clinical,tsb = clinO,isTCGA = TRUE)
lihcET=subsetMaf(maf = maf_clinical,tsb = clinET,isTCGA = T)

# 使用mafCompare比较差异突变基因
fvsm <- mafCompare(m1=lihcET, m2=lihcO, m1Name="E2F6_H&TOP2A_H", m2Name="Others", minMut=5)
# 结果保存到文件"highrisk_vs_lowrisk.tsv"
write.csv(fvsm$results, file="figure10_mutate/Others_vs_ET.csv")
#1.绘制森林图
pdf('figure10_mutate/forest_diffriskgroup_degmutation.pdf',width = 10,height = 8)
forestPlot(mafCompareRes=fvsm, pVal=0.05, color=c("#6d9cc5", "#e7857b"), geneFontSize=0.8)
dev.off()
#2.比较两个队列的oncoplot
#使用coOncoplot并排绘制两个队列的oncoplot：
a=fvsm$results
table(a$pval<0.05)
a$Hugo_Symbol[a$pval<0.05]
genes <- a$Hugo_Symbol[a$pval<0.05]
coOncoplot(m1=lihcET, m2=lihcO, m1Name="E2F6_H&TOP2A_H", m2Name="Others", genes=genes)
dev.off()
pdf('figure10_mutate/oncoplot_diffriskgroup_degmutation.pdf',width = 10,height = 12)
oncoplot(maf = maf_clinical,legend_height = 6,genes=genes,
         clinicalFeatures = "Group",sortByAnnotation = T,
         draw_titv = T)
dev.off()

###临床富集分析
clin_enrich <- clinicalEnrichment(maf=maf_clinical, clinicalFeature="Group")
# clinicalEnrichment会进行两两配对（比如female vs male）和分组的Fisher精确检验（female vs rest），结果分别储存于$pairwise_comparision和$groupwise_comparision之中。
# 这里就不贴返回结果了，用writetable保存到文件
write.table(clin_enrich$pairwise_comparision, file="figure10_mutate/clin_enrich_pair.tsv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(clin_enrich$groupwise_comparision, file="figure10_mutate/clin_enrich_group.tsv", quote=FALSE, row.names=FALSE, sep="\t")
write.csv(clin_enrich$pairwise_comparision, file="figure10_mutate/clin_enrich_pair.csv")
write.csv(clin_enrich$groupwise_comparision, file="figure10_mutate/clin_enrich_group.csv")
#绘制富集结果，pVal用来控制输出pvalue的阈值，似乎目前还不能按照fdr进行筛选：
dev.off()
par('mar')
par(mar=c(1,1,1,1))
pdf('figure10_mutate/clin_enrich.pdf',width = 5,height = 5)
plotEnrichmentResults(enrich_res=clin_enrich, pVal=0.05)
dev.off()
#这里的分析做得比较粗糙，因为没有去掉未报道性别的case，
#建议先对临床特征进行筛选再使用，否则会影响groupwise的结果
a=clin_enrich$pairwise_comparision
b=clin_enrich$groupwise_comparision

