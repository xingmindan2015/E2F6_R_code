rm(list = ls())
load("figure02_E2F6/lihc_sur_data.Rdata")
load('figure06_two_group/tcga_predict_factor_clinicalinformation.Rdata')
exp = exprSet

exp = exp[c('E2F6','TOP2A'),]
write.table(exp,file = 'figure18_durg/cMAP/exp_e2f6_top2a.txt')
dat = as.data.frame(t(exp))
write.table(dat,file = 'figure18_durg/cMAP/dat_e2f6_top2a.txt')
