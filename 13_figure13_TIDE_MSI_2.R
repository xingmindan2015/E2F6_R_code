rm(list = ls())
tide <- read.csv("figure13_TIDE_MSI/Tide_result2.csv",row.names = 1,check.names = F)
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
meta = meta[match(rownames(tide),meta$ID),]
identical(meta$ID,rownames(tide))
tide$Group = meta$Group
table(tide$Group)
a1 = tide[tide$Group=='E2F6_H&TOP2A_H',]
a2 = tide[tide$Group=='Others',]
table(a1$Responder)
table(a2$Responder)
