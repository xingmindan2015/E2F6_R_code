
# top2a_network -----------------------------------------------------------

rm(list = ls())
a=readxl::read_xlsx('figure01_top2a/LIHC_triplet.xlsx')

top2a_network = filter(a,`#Gene Symbol` =='TOP2A' ) %>% 
  .[,c(2,4,6)]
table(top2a_network$`#TF Symbol`)
colnames(top2a_network)=c('LncRNA_Symbol', "TF_Symbol","Gene_Symbol")
node_lncRNA=data.frame(node = top2a_network$LncRNA_Symbol,
                       type = 'lncRNA') %>% .[!duplicated(.$node),]
node_tf = data.frame(node = top2a_network$TF_Symbol,
                     type = 'TF')%>% .[!duplicated(.$node),]
node_gene = data.frame(node = top2a_network$Gene_Symbol,
                       type = 'Gene')%>% .[!duplicated(.$node),]
node_type = do.call(rbind,list(node_gene,node_lncRNA,node_tf))
node_type = node_type[!duplicated(node_type$node),]
write.csv(node_type,file = 'figure01_top2a/node_type.csv',quote = F,row.names = F)
write.csv(node_lncRNA,file = 'figure01_top2a/node_lncRNA.csv',quote = F,row.names = F)
write.csv(node_tf,file = 'figure01_top2a/node_tf.csv',quote = F,row.names = F)
node_network1 = data.frame(node1 = top2a_network$LncRNA_Symbol,
                           node2 = top2a_network$TF_Symbol)
node_network2 = data.frame(node1 = top2a_network$LncRNA_Symbol,
                           node2 = top2a_network$Gene_Symbol)
node_network3 = data.frame(node1 = top2a_network$TF_Symbol,
                           node2 = top2a_network$Gene_Symbol)
node_network = do.call(rbind,list(node_network1,node_network2,node_network3))
write.csv(node_network,file = 'figure01_top2a/node_network.csv',quote = F,row.names = F)


# TCGA数据库中top2a与TF的相关性分析 ----------------------------------------------------------

rm(list = ls())
load('/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/data_symbols_lihc.Rdata')
load('/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/clinical_lihc_postclear.Rdata')
fpkm_lihc1[1:4,1:4]
boxplot(fpkm_lihc1[,1:10])
exp=log2(fpkm_lihc1+1)
a=c('E2F6','TOP2A','SP1','SMARCC2','SMARCC1','CEBPB','NFYA')
#table(a%in%rownames(exp))
exp=exp[a,]
library(corrplot)
M = cor(t(exp))
testRes = cor.mtest(t(exp), conf.level = 0.95)$p
library(tidyverse)
g = pivot_longer(rownames_to_column(as.data.frame(M),var = "from"),
                 cols = 2:(ncol(M)+1),
                 names_to = "to",
                 values_to = "cor")
gp = pivot_longer(rownames_to_column(as.data.frame(testRes)),
                  cols = 2:(ncol(M)+1),
                  names_to = "gene",
                  values_to = "p")
g$p = gp$p
g = g[g$from!=g$to,]
g$group = case_when(g$cor>0.3 & g$p<0.05 ~ "positive",
                    g$cor< -0.3 & g$p<0.05 ~ "negative",
                    T~"not" )
g = g[g$to=='TOP2A',]
write.csv(g,file = 'figure01_top2a/tcga_tf_top2a_cor.csv',quote = F,row.names = F )
node1 = data.frame(node = g$from,
                   type = 'TF')
node2 = data.frame(node = g$to,
                   type = 'Gene')
node_type = do.call(rbind,list(node1,node2))
node_type = node_type[!duplicated(node_type$node),]
write.csv(node_type,file = 'figure01_top2a/tcga_tf_top2a_type.csv',quote = F,row.names = F)

# ICGC数据库中TOP2A与TF的相关性分析 --------------------------------------------------

rm(list = ls())
load('/home/data/t0202010/R/gongzuo01/ganai/danxibao/livercancersinglecell/database/data/icgc_tumor_normal_rawcount.Rdata')
colnames(exp)
exp[1:4,1:4]
exp=log2(cpm(exp)+1)
a=c('E2F6','TOP2A','SP1','SMARCC2','SMARCC1','CEBPB','NFYA')
table(a%in%rownames(exp))
exp=exp[a,]
library(corrplot)
M = cor(t(exp))
testRes = cor.mtest(t(exp), conf.level = 0.95)$p
library(tidyverse)
g = pivot_longer(rownames_to_column(as.data.frame(M),var = "from"),
                 cols = 2:(ncol(M)+1),
                 names_to = "to",
                 values_to = "cor")
gp = pivot_longer(rownames_to_column(as.data.frame(testRes)),
                  cols = 2:(ncol(M)+1),
                  names_to = "gene",
                  values_to = "p")
g$p = gp$p
g = g[g$from!=g$to,]
g$group = case_when(g$cor>0.3 & g$p<0.05 ~ "positive",
                    g$cor< -0.3 & g$p<0.05 ~ "negative",
                    T~"not" )
g = g[g$to=='TOP2A',]
write.csv(g,file = 'figure01_top2a/icgc_tf_top2a_cor.csv',quote = F,row.names = F )
node1 = data.frame(node = g$from,
                   type = 'TF')
node2 = data.frame(node = g$to,
                   type = 'Gene')
node_type = do.call(rbind,list(node1,node2))
node_type = node_type[!duplicated(node_type$node),]
write.csv(node_type,file = 'figure01_top2a/icgc_tf_top2a_type.csv',quote = F,row.names = F)

# TOP2A与E2F6在正常组与疾病组中的甲基化水平 -----------------------------------------------
#云雨图####
rm(list = ls())
library(tidyverse)
library(data.table)
library(dplyr)
dat <- fread('figure14_methy/TCGA-LIHC.methylation450.tsv.gz') %>% as.data.frame()

colnames(dat)[1] <- 'id'

ann <- fread('figure14_methy/illuminaMethyl450_hg38_GDC') %>% dplyr::select(1,2)

colnames(ann) <- c('id','symbol')
table(duplicated(ann$symbol))
# a <- dat %>% separate(id,c('id','version'),"[.]")

# dat$id <- a$id
merge <- inner_join(dat,ann,by = 'id') %>% dplyr::select(id,symbol,everything())

# merge[,4:ncol(merge)] <- 2^(merge[,4:ncol(merge)]) - 1

merge <- merge %>% distinct(symbol,.keep_all = TRUE)
merge <- merge[,-1]
rownames(merge) = merge$symbol
a = c('TOP2A','E2F6')
table(a %in% rownames(merge))
exp = merge[a,-1]
dat = as.data.frame(t(exp))
table(str_sub(rownames(dat),14,15))
dat = dat[str_sub(rownames(dat),14,15)!='02',]
dat$Group = ifelse(str_sub(rownames(dat),14,15)=='01','Tumor','Non-tumor')
table(dat$Group)
ordercolors<-c("#6d9cc5","#e7857b")
library(gghalves)
p1 = ggplot(data = dat,
           aes(x=Group, y=TOP2A, fill=Group)) +
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
ggsave('figure01_top2a/tcga_methv_top2a.pdf',p1,width = 3,height = 4,units = 'in')
p2 = ggplot(data = dat,
            aes(x=Group, y=E2F6, fill=Group)) +
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
ggsave('figure01_top2a/tcga_methv_e2f6.pdf',p2,width = 3,height = 4,units = 'in')
#分面小提琴图或者豆荚图####
rm(list = ls())
library(tidyverse)
library(data.table)
library(dplyr)
dat <- fread('figure14_methy/TCGA-LIHC.methylation450.tsv.gz') %>% as.data.frame()

colnames(dat)[1] <- 'id'

ann <- fread('figure14_methy/illuminaMethyl450_hg38_GDC') %>% dplyr::select(1,2)

colnames(ann) <- c('id','symbol')
table(duplicated(ann$symbol))
# a <- dat %>% separate(id,c('id','version'),"[.]")

# dat$id <- a$id
merge <- inner_join(dat,ann,by = 'id') %>% dplyr::select(id,symbol,everything())

# merge[,4:ncol(merge)] <- 2^(merge[,4:ncol(merge)]) - 1
a = c('E2F6','TOP2A')
merge <- merge [merge$symbol%in%a,]
rownames(merge) = merge$id
merge_E2F6 = merge[merge$symbol=='E2F6',] %>% 
  .[,-c(1,2)] %>% 
  .[,str_sub(colnames(.),14,15)!='02'] 
rownames(merge_E2F6)
nop_E2F6 = c('cg06239478','cg07890032','cg12320676','cg15066463','cg19388373')
merge_E2F6 = merge_E2F6[!rownames(merge_E2F6)%in%nop_E2F6,]
dat_E2F6 = as.data.frame(t(merge_E2F6))
dat_E2F6$Group = ifelse(str_sub(rownames(dat_E2F6),14,15)=='01','Tumor','Non-tumor')
dat_E2F6$Group = factor(dat_E2F6$Group,levels = c('Non-tumor','Tumor'))
merge_TOP2A = merge[merge$symbol=='TOP2A',]%>% 
  .[,-c(1,2)] %>% 
  .[,str_sub(colnames(.),14,15)!='02']
rownames(merge_TOP2A)
nop_TOP2A = c('cg07388903','cg15769921','cg15797971','cg21499137','cg22912359','cg25581784','cg26164164')
merge_TOP2A = merge_TOP2A[!rownames(merge_TOP2A)%in%nop_TOP2A,]
dat_TOP2A = as.data.frame(t(merge_TOP2A))
dat_TOP2A$Group = ifelse(str_sub(rownames(dat_TOP2A),14,15)=='01','Tumor','Non-tumor')
dat_TOP2A$Group = factor(dat_TOP2A$Group,levels = c('Non-tumor','Tumor'))
library("RColorBrewer")
table(duplicated(colnames(dat_E2F6)))
dat_E2F6 = pivot_longer(dat_E2F6,cols = -Group,names_to = 'Genes',values_to = 'Value')
dat_TOP2A = pivot_longer(dat_TOP2A,cols = -Group,names_to = 'Genes',values_to = 'Value')
colnames(dat_E2F6)
library(devtools)

#install_github("JanCoUnchained/ggunchained")

library(ggunchained)
su = group_by(dat_E2F6,Group,Genes) %>% 
  summarise(s = sd(Value),Value = mean(Value))
p_E2F6 = ggplot(dat_E2F6, aes(x = Genes, y = Value,fill=Group))+
  geom_split_violin(color = NA)+
  guides()+
  stat_compare_means(label = 'p.signif')+
  scale_fill_manual(values = c('#6d9cc5','#e7857b'))+
  geom_point(data = su,aes(x=Genes, y=Value),position=position_dodge(0.1),size=0.1)+
  geom_errorbar(data = su,
                  mapping = aes(ymin = Value -s, ymax = Value+s), #误差条表示95%的置信区间
                  width=0, 
                  position=position_dodge(0.1), 
                  color="black",
                  size=0.1)+  
  theme_bw()+
  xlab('')+
  ylab('E2F6_Beta')+
  theme(axis.text.x = element_text(vjust = 1,hjust = 1,angle = 80))
#p_E2F6 = p_E2F6+facet_wrap(~Genes,scales = 'free',nrow = 1)
su = group_by(dat_TOP2A,Group,Genes) %>% 
  summarise(s = sd(Value),Value = mean(Value))
p_TOP2A = ggplot(dat_TOP2A, aes(x = Genes, y = Value,fill=Group))+
  geom_split_violin(color = NA)+
  guides()+
  stat_compare_means(label = 'p.signif')+
  scale_fill_manual(values = c('#6d9cc5','#e7857b'))+
  geom_point(data = su,aes(x=Genes, y=Value),position=position_dodge(0.1),size=0.1)+
  geom_errorbar(data = su,
                mapping = aes(ymin = Value -s, ymax = Value+s), #误差条表示95%的置信区间
                width=0, 
                position=position_dodge(0.1), 
                color="black",
                size=0.1)+  
  theme_bw()+
  xlab('')+
  ylab('TOP2A_Beta')+
  theme(axis.text.x = element_text(vjust = 1,hjust = 1,angle = 80))
#p_TOP2A = p_TOP2A+facet_wrap(~Genes,scales = 'free',nrow = 2)
ggsave('figure01_top2a/E2F6_geneid_methsurv.pdf',p_E2F6,width = 8,height = 4,units = 'in')
ggsave('figure01_top2a/TOP2A_geneid_methsurv.pdf',p_TOP2A,width = 8,height = 4,units = 'in')
save(merge_E2F6,merge_TOP2A,file = 'figure01_top2a/E2F6_TOP2A_tumor_normal_meth.Rdata')

# E2F6_TOP2A_预后分析 ---------------------------------------------------------

rm(list = ls())
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
load('figure01_top2a/E2F6_TOP2A_tumor_normal_meth.Rdata')
merge_E2F6 = merge_E2F6[,str_sub(colnames(merge_E2F6),14,15)=='01']
colnames(merge_E2F6)=str_sub(colnames(merge_E2F6),1,12)
table(colnames(merge_E2F6)%in%meta$ID)
merge_E2F6=merge_E2F6[,meta$ID]
merge_TOP2A = merge_TOP2A[,str_sub(colnames(merge_TOP2A),14,15)=='01']
colnames(merge_TOP2A)=str_sub(colnames(merge_TOP2A),1,12)
table(colnames(merge_TOP2A)%in%meta$ID)
merge_TOP2A=merge_TOP2A[,meta$ID]
merge_E2F6 = merge_E2F6[,match(meta$ID,colnames(merge_E2F6))]
merge_TOP2A = merge_TOP2A[,match(meta$ID,colnames(merge_TOP2A))]
identical(meta$ID,colnames(merge_E2F6))
identical(meta$ID,colnames(merge_TOP2A))
coxfile = "figure01_top2a/e2f6_uni_cox.Rdata"
exprSet = merge_E2F6
if(!file.exists(coxfile)){
  cox_results <-apply(exprSet , 1 , function(gene){
    meta$gene = gene
    #可直接使用连续型变量
    m = coxph(Surv(Time, Event) ~ gene, data =  meta)
    beta <- coef(m)
    se <- sqrt(diag(vcov(m)))
    HR <- exp(beta)
    HRse <- HR * se
    tmp <- round(cbind(coef = beta, 
                       se = se, z = beta/se, 
                       p = 1 - pchisq((beta/se)^2, 1),
                       HR = HR, HRse = HRse,
                       HRz = (HR - 1) / HRse, 
                       HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                       HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                       HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
    
    return(tmp['gene',]) 
  })
  cox_results=as.data.frame(t(cox_results))
  save(cox_results,file = coxfile)
}
load(coxfile)
table(cox_results$p<0.01)
table(cox_results$p<0.05)
cox = rownames(cox_results)[cox_results$p<0.05];length(cox)
write.csv(cox_results,file = 'figure01_top2a/E2F6_unicox.csv',quote = F)

coxfile = "figure01_top2a/top2a_uni_cox.Rdata"
exprSet = merge_TOP2A
if(!file.exists(coxfile)){
  cox_results <-apply(exprSet , 1 , function(gene){
    meta$gene = gene
    #可直接使用连续型变量
    m = coxph(Surv(Time, Event) ~ gene, data =  meta)
    beta <- coef(m)
    se <- sqrt(diag(vcov(m)))
    HR <- exp(beta)
    HRse <- HR * se
    tmp <- round(cbind(coef = beta, 
                       se = se, z = beta/se, 
                       p = 1 - pchisq((beta/se)^2, 1),
                       HR = HR, HRse = HRse,
                       HRz = (HR - 1) / HRse, 
                       HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                       HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                       HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
    
    return(tmp['gene',]) 
  })
  cox_results=as.data.frame(t(cox_results))
  save(cox_results,file = coxfile)
}
load(coxfile)
table(cox_results$p<0.01)
table(cox_results$p<0.05)
cox = rownames(cox_results)[cox_results$p<0.05];length(cox)
write.csv(cox_results,file = 'figure01_top2a/TOP2A_unicox.csv',quote = F)



# E2F6_TOP2A_CNV ----------------------------------------------------------

rm(list = ls())
library(magrittr)
library(data.table)
library(maftools)
## cnv数据清洗----
dat_cnv <- fread("figure12_cnv/TCGA-LIHC.gistic.tsv.gz", data.table = F)
rownames(dat_cnv) <- dat_cnv$`Gene Symbol`
library(data.table)
ann <- fread('figure12_cnv/gencode.v22.annotation.gene.probeMap') %>% dplyr::select(1,2)
colnames(ann) <- c('id','symbol')
colnames(dat_cnv)[1] = c('id')
dat_cnv = dat_cnv[,!duplicated(colnames(dat_cnv))]
table(duplicated(colnames(dat_cnv)))

dat <- inner_join(dat_cnv,ann,by = 'id') %>% dplyr::select(symbol, everything())
dat <- dat %>% distinct(symbol,.keep_all = TRUE)
rownames(dat) <- dat$symbol
dat <- dat[,-c(1:2)]
gene <- c('TOP2A','E2F6')
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
ggsave("figure01_top2a/top2a_e2f6_cnv.pdf", width = 3, height = 3,units = 'in')


# E2F6_TOP2A_突变瀑布图 --------------------------------------------------------
#从cBioportal数据库中获得



