rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 200)

#1. segment file数据下载和处理
##1.1 从TCGA下载数据
# library(SummarizedExperiment)
# library(TCGAbiolinks)
# 
# query <- GDCquery(project = "TCGA-LIHC",
#                   data.category = "Copy Number Variation",
#                   data.type = "Masked Copy Number Segment")
# GDCdownload(query, method = "api")
# LIHC_CNV_download <- GDCprepare(query = query,
#                                 save = TRUE, save.filename = "database/data/LIHC_CNV_download.Rdata")
##1.2 数据处理
rm(list = ls())
A = load("figure11_signature/LIHC_CNV_download.Rdata")
tumorCNV <- eval(parse(text = A))

###改名
tumorCNV <- tumorCNV[,2:7]
tumorCNV <- tumorCNV[,c("Sample", "Chromosome",
                        "Start", "End", "Num_Probes", "Segment_Mean")]
##1.3 提取01A肿瘤样本
tum.sam <- substr(tumorCNV$Sample,14,16) == "01A"
table(tum.sam)
tumorCNV <- tumorCNV[tum.sam,]
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
table(meta$Group)
ET_name=rownames(meta)[meta$Group=='E2F6_H&TOP2A_H']
OR_name=rownames(meta)[meta$Group=='Others']
tumorCNV[1:4,]
length(table(str_sub(tumorCNV$Sample,1,12)))
tumorCNV_ET=tumorCNV[str_sub(tumorCNV$Sample,1,12)%in%ET_name,]
tumorCNV_OR=tumorCNV[str_sub(tumorCNV$Sample,1,12)%in%OR_name,]
tumorCNV=tumorCNV[str_sub(tumorCNV$Sample,1,12)%in%rownames(meta),]
write.table(tumorCNV, file = "figure12_cnv/all/segment_file.txt", sep = "\t", row.names = F, quote = F)
write.table(tumorCNV_ET, file = "figure12_cnv/E2F6&TOP2A/segment_file_ET.txt", sep = "\t", row.names = F, quote = F)
write.table(tumorCNV_OR, file = "figure12_cnv/Others/segment_file_OR.txt", sep = "\t", row.names = F, quote = F)


library(data.table)
Marker <- read.delim("figure12_cnv/snp6.na35.remap.hg38.subset.txt.gz")
head(Marker )

fre <- Marker$freqcnv == "FALSE"
table(fre)
Marker <- Marker[fre,1:3]
colnames(Marker) <- c("Marker_Name", "Chromosome", "Marker_Position")
write.table(Marker, file = "figure12_cnv/marker_file.txt", sep = "\t", row.names = F, quote = F)

rm(list = ls())
library(maftools)
lihc.gistic <- readGistic(gisticAllLesionsFile="figure12_cnv/all/all_lesions.conf_90.txt", 
                          gisticAmpGenesFile="figure12_cnv/all/amp_genes.conf_90.txt", 
                          gisticDelGenesFile="figure12_cnv/all/del_genes.conf_90.txt",
                          gisticScoresFile="figure12_cnv/all/scores.gistic", isTCGA=TRUE)
#直接键入GISTIC对象 lihc.gistic可以得到一些简单的统计：
lihc.gistic
#使用getSampleSummary/getGeneSummary/getCytobandSummary进行统计，并用write.GisticSummary保存统计结果：
getSampleSummary(lihc.gistic)
getGeneSummary(lihc.gistic)
getCytobandSummary(lihc.gistic)
write.GisticSummary(gistic = lihc.gistic,basename = 'figure12_cnv/all/lihc_gistic2')

#cnv数据的可视化
#GISTIC的输出结果可以用染色体图、气泡图（点图）、oncoplot进行可视化
#1.染色体图
#使用gisticChromPlot 可以绘制染色各位点G-score
pdf('figure12_cnv/all/boxplot.pdf',width = 10,height = 6)
gisticChromPlot(gistic = lihc.gistic, ref.build = 'hg38',
                color = c("#e7857b","#6d9cc5")
)
dev.off()
mycolor = c("#e7857b","#6d9cc5")
names(mycolor) = c("Amp","Del")
pdf('figure12_cnv/all/heatmap.pdf',width = 10,height = 6)
gisticOncoPlot(gistic = lihc.gistic,
               top = 20,colors = mycolor)
dev.off()

###ET风险组
rm(list = ls())
library(maftools)

lihcHr.gistic <- readGistic(gisticAllLesionsFile="figure12_cnv/E2F6&TOP2A/all_lesions.conf_90.txt", 
                            gisticAmpGenesFile="figure12_cnv/E2F6&TOP2A/amp_genes.conf_90.txt", 
                            gisticDelGenesFile="figure12_cnv/E2F6&TOP2A/del_genes.conf_90.txt",
                            gisticScoresFile="figure12_cnv/E2F6&TOP2A/scores.gistic", isTCGA=TRUE)
#直接键入GISTIC对象 lihc.gistic可以得到一些简单的统计：
lihcHr.gistic
#使用getSampleSummary/getGeneSummary/getCytobandSummary进行统计，并用write.GisticSummary保存统计结果：
getSampleSummary(lihcHr.gistic)
getGeneSummary(lihcHr.gistic)
getCytobandSummary(lihcHr.gistic)
write.GisticSummary(gistic = lihcHr.gistic,basename = 'figure12_cnv/E2F6&TOP2A/lihcHr_gistic2')

#cnv数据的可视化
#GISTIC的输出结果可以用染色体图、气泡图（点图）、oncoplot进行可视化
#1.染色体图
#使用gisticChromPlot 可以绘制染色各位点G-score
pdf('figure12_cnv/E2F6&TOP2A/boxplot_high.pdf',width = 10,height = 6)
gisticChromPlot(gistic = lihcHr.gistic, ref.build = 'hg38',
                color =c("#e7857b","#6d9cc5")
)
dev.off()
pdf('figure12_cnv/E2F6&TOP2A/bubbleplot_high.pdf',width = 6,height = 6)
gisticBubblePlot(gistic = lihcHr.gistic,
                 c("#e7857b","#6d9cc5")
)
dev.off()
mycolor =c("#e7857b","#6d9cc5")
names(mycolor) = c("Amp","Del")
pdf('figure12_cnv/E2F6&TOP2A/heatmap_high.pdf',width = 10,height = 6)
gisticOncoPlot(gistic = lihcHr.gistic,
               top = 20,colors = mycolor)
dev.off()

###OR风险组
rm(list = ls())
library(maftools)
lihcLr.gistic <- readGistic(gisticAllLesionsFile="figure12_cnv/Others/all_lesions.conf_90.txt", 
                            gisticAmpGenesFile="figure12_cnv/Others/amp_genes.conf_90.txt", 
                            gisticDelGenesFile="figure12_cnv/Others/del_genes.conf_90.txt",
                            gisticScoresFile="figure12_cnv/Others/scores.gistic", isTCGA=TRUE)
#直接键入GISTIC对象 lihc.gistic可以得到一些简单的统计：
lihcLr.gistic
#使用getSampleSummary/getGeneSummary/getCytobandSummary进行统计，并用write.GisticSummary保存统计结果：
getSampleSummary(lihcLr.gistic)
getGeneSummary(lihcLr.gistic)
getCytobandSummary(lihcLr.gistic)
write.GisticSummary(gistic = lihcLr.gistic,basename = 'figure12_cnv/Others/lihcLr_gistic2')

#cnv数据的可视化
#GISTIC的输出结果可以用染色体图、气泡图（点图）、oncoplot进行可视化
#1.染色体图
#使用gisticChromPlot 可以绘制染色各位点G-score
pdf('figure12_cnv/Others/boxplot_low.pdf',width = 10,height = 6)

gisticChromPlot(gistic = lihcLr.gistic, ref.build = 'hg38',
                color = c("#e7857b","#6d9cc5")
)
dev.off()
pdf('figure12_cnv/Others/bubbleplot_low.pdf',width = 6,height = 6)
gisticBubblePlot(gistic = lihcLr.gistic,
                 color = c("#e7857b","#6d9cc5")
)
dev.off()
mycolor = c("#e7857b","#6d9cc5")
names(mycolor) = c("Amp","Del")
pdf('figure12_cnv/Others/heatmap_low.pdf',width = 10,height = 6)
gisticOncoPlot(lihcLr.gistic,top = 20,colors = mycolor)
dev.off()


# CNV-柱状图 -----------------------------------------------------------------
rm(list = ls())
library(magrittr)
library(data.table)
library(maftools)
## cnv数据清洗----
dat_cnv <- fread("figure12_cnv/TCGA-LIHC.gistic.tsv.gz", data.table = F)
rownames(dat_cnv) <- dat_cnv$`Gene Symbol`
dat_cnv = dat_cnv[,str_sub(colnames(dat_cnv),14,15)=='01']
colnames(dat_cnv)=str_sub(colnames(dat_cnv),1,12)
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
dat_cnv = dat_cnv[,colnames(dat_cnv)%in%rownames(meta)]
group=meta
library(data.table)
ann <- fread('figure12_cnv/gencode.v22.annotation.gene.probeMap') %>% dplyr::select(1,2)
colnames(ann) <- c('id','symbol')
dat_cnv = rownames_to_column(dat_cnv,var = 'id')
dat <- dplyr::inner_join(dat_cnv,ann,by = 'id') %>% dplyr::select(symbol, everything())
dat <- dat %>% distinct(symbol,.keep_all = TRUE)
rownames(dat) <- dat$symbol
dat <- dat[,-c(1:2)]
gene <- c('E2F6','TOP2A')
dat <- dat[gene,]
input <- dat %>% 
  t() %>% 
  data.frame() %>% 
  melt(variable.name = "gene", value.name = "group") %>% 
  mutate(group = dplyr::case_when(
    group == "1" ~ "amp",
    group == "-1" ~ "del",
    group == "0" ~ NA
  )) %>% 
  arrange(gene) %>% 
  table() %>% 
  data.frame() %>% 
  mutate(f = Freq/337 * 100 * ifelse(group == "del", -1, 1))#337是所有样本的数量
library(ggplot2)
ggplot(data = input, mapping = aes(x = gene, y = f, fill=group))+
  geom_bar(stat = "identity", width = 0.6, alpha = .75, 
           position = position_dodge(0))+
  scale_fill_manual(values = c("#e7857b","#6d9cc5"), )+
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.position = "top", legend.title = element_blank())+
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "",y = "CNV frequency(%)")+
  theme_bw()
ggsave("figure12_cnv/cnv.pdf", width = 3, height = 3,units = 'in')

rm(list = ls())
load('figure06_two_group/e2f6_top2a_two_group_deg_correlation.Rdata')
library(magrittr)
library(data.table)
library(maftools)
## cnv数据清洗----
dat_cnv <- fread("figure12_cnv/TCGA-LIHC.gistic.tsv.gz", data.table = F)
rownames(dat_cnv) <- dat_cnv$`Gene Symbol`
dat_cnv = dat_cnv[,str_sub(colnames(dat_cnv),14,15)=='01']
colnames(dat_cnv)=str_sub(colnames(dat_cnv),1,12)
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
dat_cnv = dat_cnv[,colnames(dat_cnv)%in%rownames(meta)]
group=meta
library(data.table)
ann <- fread('figure12_cnv/gencode.v22.annotation.gene.probeMap') %>% dplyr::select(1,2)
colnames(ann) <- c('id','symbol')
dat_cnv = rownames_to_column(dat_cnv,var = 'id')
dat <- dplyr::inner_join(dat_cnv,ann,by = 'id') %>% dplyr::select(symbol, everything())
dat <- dat %>% distinct(symbol,.keep_all = TRUE)
rownames(dat) <- dat$symbol
dat <- dat[,-c(1:2)]
gene <- type_df$nodes
table(gene%in%rownames(dat))
dat <- dat[rownames(dat)%in%gene,]
input <- dat %>% 
  t() %>% 
  data.frame() %>% 
  melt(variable.name = "gene", value.name = "Group") %>% 
  mutate(Group = dplyr::case_when(
    Group == "1" ~ "amp",
    Group == "-1" ~ "del",
    Group == "0" ~ NA
  )) %>% 
  arrange(gene) %>% 
  table() %>% 
  data.frame() %>% 
  mutate(f = Freq/337 * 100 * ifelse(Group == "del", -1, 1))#337是所有样本的数量
library(ggplot2)
ggplot(data = input, mapping = aes(x = gene, y = f, fill=Group))+
  geom_bar(stat = "identity", width = 0.6, alpha = .75, 
           position = position_dodge(0))+
  scale_fill_manual(values = c("#EE0000","#3B4992"), )+
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.position = "top", legend.title = element_blank())+
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "",y = "CNV frequency(%)")+
  theme_bw()+
  theme(axis.text.x = element_text(vjust = 1,hjust = 1,angle = 80))
ggsave("figure12_cnv/cnv.pdf", width = 6, height = 4,units = 'in')

