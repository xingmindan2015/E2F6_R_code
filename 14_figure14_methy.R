rm(list = ls())
library(tidyverse)
library(data.table)

dat <- fread('figure14_methy/TCGA-LIHC.methylation450.tsv.gz') %>% as.data.frame()

colnames(dat)[1] <- 'id'

ann <- fread('figure14_methy/illuminaMethyl450_hg38_GDC') %>% select(1,2)

colnames(ann) <- c('id','symbol')

# a <- dat %>% separate(id,c('id','version'),"[.]")

# dat$id <- a$id
merge <- inner_join(dat,ann,by = 'id') %>% select(id,symbol,everything())

# merge[,4:ncol(merge)] <- 2^(merge[,4:ncol(merge)]) - 1

merge <- merge %>% distinct(symbol,.keep_all = TRUE)
merge <- merge[,-1]
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')

sample <- rownames(meta)

b <- merge[,1,drop = F]

merge_t <- merge[,str_sub(colnames(merge),14,15) == "01"]
colnames(merge_t)=sample
merge_t <- merge_t[,colnames(merge_t) %in% sample]
merge_t1 <- cbind(b, merge_t) %>% select(symbol,everything())
write.csv(merge_t1, "figure14_methy/lihc_dnam_t.csv", row.names = F)
save(merge_t1,file = 'figure14_methy/lihc_dnam_t.Rdata')

library(data.table)
library(tibble)
library(reshape2)
dat <- merge_t1 %>% column_to_rownames("symbol")
identical(colnames(dat),rownames(meta))
Group <- meta$Group
load('figure06_two_group/e2f6_top2a_two_group_deg_correlation.Rdata')
Gene <- type_df$nodes
table(Gene%in%rownames(dat))
dat <- dat[rownames(dat)%in%Gene,]
identical(rownames(meta),colnames(dat))
draw_boxplot(dat,Group,ylab = 'Beta_value',xlab = '',sort = F)
ggsave("figure14_methy/methy_exp_barplot_Group.pdf",width = 8,height = 5,units = 'in')
dev.off()
Group = meta$Cluster
draw_boxplot(dat,Group,ylab = 'Beta_value',xlab = '',sort = F)
ggsave("figure14_methy/methy_exp_barplot_Cluster.pdf",width = 8,height = 5,units = 'in')
dev.off()
# input <- dat %>% 
#   t() %>% 
#   data.frame() %>% 
#   mutate(Group = Group) %>% 
#   melt(id.vars = 'Group', variable.name = "Gene", value.name = "Value") %>% 
#   arrange(Gene)
# class(input$Group)
# 
# library(ggplot2)
# library(corrgram)
# library(ggthemes)
# library(ggpubr)
# input_pvalue <- compare_means(Value ~ Group, data = input, group.by = "Gene")
# input_pvalue <- input_pvalue %>% arrange(Gene)
# write.csv(input_pvalue, "figure14_methy/dnam_exp_pvalue.csv")
# 
# p<-ggplot(input,aes(x=Gene,y=Value,fill=Group))+
#   geom_boxplot(width=0.7,size=0.3,outlier.color = NA,linewidth=0.1,fatten=1,position=position_dodge(0.85))+
#   theme_bw()+scale_fill_manual(values = c("#00468B","#ED0000"))+#颜色顺序同前面group因子化顺序一致。
#   theme(panel.grid = element_blank())+
#   stat_compare_means(symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1),
#                                         symbols = c("***", "**", "*", "")),label = "p.signif",size=6,hjust = 0.5,vjust = 0)+
#   theme(axis.text.x = element_text(angle = 45,hjust = 1))+
#   theme(legend.position = 'top')+xlab('')+ylab('value')+labs(fill='Group')
# 
# p
 ggsave("output/dnam_exp_barplot.pdf",plot=p,width = 12,height = 8)
dev.off()