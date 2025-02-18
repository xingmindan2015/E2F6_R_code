#GSE10141无正常对照，有生存信息
rm(list = ls())
options(stringsAsFactors = F)#3.6版本以后的需要用
library(GEOquery)
gse_number = "GSE10141"
eSet <- getGEO(gse_number, 
               destdir = 'database/GEO/GSE10141', 
               getGPL = F)
class(eSet)
length(eSet)
eSet = eSet[[1]]
#(1)提取表达矩阵exp
exp <- exprs(eSet)
exp[1:4,1:4]#注意看一下，注意不要log两次，取过log的表达量的值在0-20之间
exp1=log2(exp+1)
exp1[1:4,1:4]
boxplot(exp1)

#(2)提取临床信息
pd<- pData(eSet)

colnames(pd)
pd1=pd[,c(2,45:51)]
colnames(pd1)
colnames(pd1)=c('ID','Alcohol','HBV','HCV','Microvascular_invasion','Satellite_lesions','Event','Time')
pd1$Event=str_split(pd1$Event,':',simplify = T)[,3]
#mfs：无转移生存期；os：总体生存期。
pd2=pd1[,c(1,7:8)]
class(pd2$Event)
pd2$Event=as.numeric(pd2$Event)
class(pd2$Time)
pd2$Time=as.numeric(pd2$Time)
pd2$Time = as.numeric(pd2$Time)/30
# 去掉生存信息不全或者生存时间小于0.1(月)的病人，样本纳排标准不唯一，且差别很大
k1 = pd2$Time>=0.1;table(k1)##去掉1个
k2 = !(is.na(pd2$Time)|is.na(pd2$Event));table(k2)
pd2 = pd2[k1&k2,]
exp1[1:4,1:4]
exp=as.data.frame(exp1)
exp=exp[,colnames(exp)%in%pd2$ID]
exp=exp[,match(pd2$ID,colnames(exp))]
identical(pd2$ID,colnames(exp))

#(4)提取芯片平台编号
gpl_number <- eSet@annotation

#(5)下载注释文件
a = getGEO(gpl_number,destdir = "database/GEO/GSE10141/")
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

identical(colnames(exp),rownames(pd2))
save(pd2,exp,file = 'database/GEO/GSE10141/exp_pd2.Rdata')
colnames(pd2)
a=pd2[,c('ID','Event','Time')]
a$Event=ifelse(a$Event==0,'Alive','Death')
a$Event=factor(a$Event,levels = c('Alive','Death'))
library('table1')
table=table1(~ Time|Event, data=a, overall="Total")
save(table,file='database/GEO/GSE10141/table_GSE10141_HCC_clinicalinformation.Rdata')

#GSE10143
rm(list = ls())
options(stringsAsFactors = F)#3.6版本以后的需要用
library(GEOquery)
gse_number = "GSE10143"
eSet <- getGEO(gse_number, 
               destdir = 'database/GEO/GSE10143', 
               getGPL = F)
class(eSet)
length(eSet)
eSet = eSet[[1]]
#(1)提取表达矩阵exp
exp <- exprs(eSet)
exp[1:4,1:4]#注意看一下，注意不要log两次，取过log的表达量的值在0-20之间
exp1=log2(exp+1)
exp1[1:4,1:4]
boxplot(exp1[,1:162])
exp1=exp1[,1:162]
exp1 = limma::normalizeBetweenArrays(exp1)
#(2)提取临床信息
pd<- pData(eSet)
pd1=pd[pd$geo_accession%in%colnames(exp1),]
colnames(pd1)
pd1=pd1[,c(2,8,47,54,55,57,59,60,61)]
colnames(pd1)
pd1$Group=ifelse(str_detect(pd1$source_name_ch1,'cirrhotic'),'Non-Tumor','Tumor')
pd1=pd1[,c(1,3:10)]
colnames(pd1)=c('ID','Alcohol','HBV','HCV','Microvascular_invasion','Satellite_lesions','Event','Time','Group')
pd1$Event=str_split(pd1$Event,':',simplify = T)[,3]
#mfs：无转移生存期；os：总体生存期。
pd2=pd1[,c(1,7:9)]
class(pd2$Event)
pd2$Event=as.numeric(pd2$Event)
class(pd2$Time)
pd2$Time=as.numeric(pd2$Time)
pd2$Time = as.numeric(pd2$Time)/30
# 去掉生存信息不全或者生存时间小于0.1(月)的病人，样本纳排标准不唯一，且差别很大
# k1 = pd2$Time>=0.1;table(k1)##去掉1个
# k2 = !(is.na(pd2$Time)|is.na(pd2$Event));table(k2)
# pd2 = pd2[k1&k2,]
exp1[1:4,1:4]
exp1=as.data.frame(exp1)
exp1=exp1[,colnames(exp1)%in%pd2$ID]
exp1=exp1[,match(pd2$ID,colnames(exp1))]
identical(pd2$ID,colnames(exp1))

#(4)提取芯片平台编号
gpl_number <- eSet@annotation

#(5)下载注释文件
a = getGEO(gpl_number,destdir = "database/GEO/GSE10143/")
b = a@dataTable@table
colnames(b)
ids2 = b[,c("ID","Symbol")]
colnames(ids2) = c("probe_id","symbol")
ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),]

table(rownames(exp1)%in%ids2$probe_id)
exp1=exp1[rownames(exp1)%in%ids2$probe_id,]
table(duplicated(ids2$symbol))
ids2=ids2[!duplicated(ids2$symbol),]
table(rownames(exp1)%in%ids2$probe_id)
exp1=exp1[rownames(exp1)%in%ids2$probe_id,]
class(ids2$probe_id)
ids2=ids2[ids2$probe_id%in%rownames(exp1),]
dim(exp1)
exp1=exp1[match(ids2$probe_id,rownames(exp1)),]
identical(ids2$probe_id,rownames(exp1))
rownames(exp1)=ids2$symbol

identical(colnames(exp1),rownames(pd2))
save(pd2,exp1,file = 'database/GEO/GSE10143/exp1_pd2.Rdata')
colnames(pd2)
a=pd2[,c('ID','Event','Time','Group')]
a$Event=ifelse(a$Event==0,'Alive','Death')
a$Event=factor(a$Event,levels = c('Alive','Death'))
library('table1')
table=table1(~ Time+Event|Group, data=a, overall="Total")
save(table,file='database/GEO/GSE10143/table_GSE10143_HCCvsLC_clinicalinformation.Rdata')

#GSE14520
rm(list = ls())
options(stringsAsFactors = F)#3.6版本以后的需要用
library(GEOquery)
gse_number = "GSE14520"
eSet <- getGEO(gse_number, 
               destdir = 'database/GEO/GSE14520', 
               getGPL = F)
class(eSet)
length(eSet)
eSet1 = eSet[[1]]
#(1)提取表达矩阵exp
exp1 <- exprs(eSet1)
exp1[1:4,1:4]#注意看一下，注意不要log两次，取过log的表达量的值在0-20之间
exp1[1:4,1:4]
boxplot(exp1[,1:100])
exp1 = limma::normalizeBetweenArrays(exp1)
#(2)提取临床信息
pd1<- pData(eSet1)
colnames(pd1)
pd1[1:4,]

pd1=pd1[,c(1,2,46)]
colnames(pd1)=c('Title','ID','Group')
pd1$Title
pd1$Group=ifelse(str_detect(pd1$Title,'Non-Tumor'),'Non-Tumor','Tumor')
pd1$Group=factor(pd1$Group,levels = c('Non-Tumor','Tumor'))
identical(pd1$ID,colnames(exp1))
table(pd1$Group)
#(4)提取芯片平台编号
gpl_number <- eSet1@annotation

#(5)下载注释文件

if(!require(hthgu133a.db))BiocManager::install("hthgu133a.db")
library(hthgu133a.db)
ls("package:hthgu133a.db")
ids <- toTable(hthgu133aSYMBOL)
head(ids)
ids2=ids
exp1[1:4,]
table(rownames(exp1)%in%ids2$probe_id)
exp1=exp1[rownames(exp1)%in%ids2$probe_id,]
table(duplicated(ids2$symbol))
ids2=ids2[!duplicated(ids2$symbol),]
table(rownames(exp1)%in%ids2$probe_id)
exp1=exp1[rownames(exp1)%in%ids2$probe_id,]
class(ids2$probe_id)
ids2=ids2[ids2$probe_id%in%rownames(exp1),]
dim(exp1)
exp1=exp1[match(ids2$probe_id,rownames(exp1)),]
identical(ids2$probe_id,rownames(exp1))
rownames(exp1)=ids2$symbol

identical(colnames(exp1),rownames(pd1))
save(pd1,exp1,file = 'database/GEO/GSE14520/exp1_pd1.Rdata')

eSet2 = eSet[[2]]
#(1)提取表达矩阵exp
exp2 <- exprs(eSet2)
exp2[1:4,1:4]#注意看一下，注意不要log两次，取过log的表达量的值在0-20之间
exp2[1:4,1:4]
boxplot(exp2)
exp2 = limma::normalizeBetweenArrays(exp2)
#(2)提取临床信息
pd2<- pData(eSet2)
table(pd2$`Tissue:ch1`)
colnames(pd2)
pd2=pd2[,c(1,2)]
colnames(pd2)=c('Title','ID')
pd2$Title
pd2$Group=ifelse(str_detect(pd2$Title,'healthy'),'Non-Tumor',
                 ifelse(str_detect(pd2$Title,'Non-Tumor'),'Non-Tumor','Tumor'))

pd2$Group=factor(pd2$Group,levels = c('Non-Tumor','Tumor'))
identical(pd2$ID,colnames(exp2))
table(pd2$Group)
#(4)提取芯片平台编号
gpl_number <- eSet2@annotation

#(5)下载注释文件

if(!require(hgu133a2.db))BiocManager::install("hgu133a2.db")
library(hgu133a2.db)
ls("package:hgu133a2.db")
ids <- toTable(hgu133a2SYMBOL)
head(ids)
ids2=ids
exp2[1:4,]
table(rownames(exp2)%in%ids2$probe_id)
exp2=exp2[rownames(exp2)%in%ids2$probe_id,]
table(duplicated(ids2$symbol))
ids2=ids2[!duplicated(ids2$symbol),]
table(rownames(exp2)%in%ids2$probe_id)
exp2=exp2[rownames(exp2)%in%ids2$probe_id,]
class(ids2$probe_id)
ids2=ids2[ids2$probe_id%in%rownames(exp2),]
dim(exp2)
exp2=exp2[match(ids2$probe_id,rownames(exp2)),]
identical(ids2$probe_id,rownames(exp2))
rownames(exp2)=ids2$symbol

identical(colnames(exp2),rownames(pd2))
save(pd2,exp2,file = 'database/GEO/GSE14520/exp2_pd2.Rdata')

rm(list = ls())
options(stringsAsFactors = F)#3.6版本以后的需要用
library(GEOquery)
gse_number = "GSE20140"
eSet <- getGEO(gse_number, 
               destdir = 'database/GEO/GSE20140', 
               getGPL = F)
class(eSet)
length(eSet)
eSet1 = eSet[[1]]
#(1)提取表达矩阵exp
exp1 <- exprs(eSet1)
exp1[1:4,1:4]#注意看一下，注意不要log两次，取过log的表达量的值在0-20之间
exp1=log2(exp1+1)
exp1[1:4,1:4]
boxplot(exp1)

#(2)提取临床信息
pd1<- pData(eSet1)
colnames(pd1)
pd1[1:4,]

pd1=pd1[,c(1,2,31,32)]
colnames(pd1)=c('Title','ID','Location','Group')
pd1$Title
table(pd1$Group)
pd1$Group=ifelse(pd1$Group=='cirrhosis','Non-Tumor','Tumor')
pd1$Group=factor(pd1$Group,levels = c('Non-Tumor','Tumor'))
identical(pd1$ID,colnames(exp1))
table(pd1$Group)
#(4)提取芯片平台编号
gpl_number <- eSet1@annotation

#(5)下载注释文件
a = getGEO(gpl_number,destdir = "database/GEO/GSE20140/")
b = a@dataTable@table

colnames(b)
ids2 = b[,c("ID","SYMBOL")]
colnames(ids2) = c("probe_id","symbol")
ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),]

table(rownames(exp1)%in%ids2$probe_id)
exp1=exp1[rownames(exp1)%in%ids2$probe_id,]
table(duplicated(ids2$symbol))
ids2=ids2[!duplicated(ids2$symbol),]
table(rownames(exp1)%in%ids2$probe_id)
exp1=exp1[rownames(exp1)%in%ids2$probe_id,]
class(ids2$probe_id)
ids2=ids2[ids2$probe_id%in%rownames(exp1),]
dim(exp1)
exp1=exp1[match(ids2$probe_id,rownames(exp1)),]
identical(ids2$probe_id,rownames(exp1))
rownames(exp1)=ids2$symbol

identical(colnames(exp1),rownames(pd1))
save(pd1,exp1,file = 'database/GEO/GSE20140/exp1_pd1.Rdata')


eSet2 = eSet[[2]]
#(1)提取表达矩阵exp
exp2 <- exprs(eSet2)
exp2[1:4,1:4]#注意看一下，注意不要log两次，取过log的表达量的值在0-20之间
exp2=log2(exp2+1)
exp2[1:4,1:4]
boxplot(exp2)
exp2 = limma::normalizeBetweenArrays(exp2)
#(2)提取临床信息
pd2<- pData(eSet2)
colnames(pd2)
pd2[1:4,]
a=pd2[,51:64]
pd2=pd2[,c(1,2,8,58,61:63)]
colnames(pd2)=c('Title','ID','Group','Recurrence_status','Survival_time','Survival_status','Recurrence_time')
pd2$Title
table(pd2$Group)
pd2$Group=ifelse(str_detect(pd2$Title,'cirrhotic'),'Non-Tumor','Tumor')
pd2$Group=factor(pd2$Group,levels = c('Non-Tumor','Tumor'))
identical(pd2$ID,colnames(exp2))
table(pd2$Group)
#(4)提取芯片平台编号
gpl_number <- eSet2@annotation

#(5)下载注释文件
a = getGEO(gpl_number,destdir = "database/GEO/GSE20140/")
b = a@dataTable@table

colnames(b)
ids2 = b[,c("ID","Symbol")]
colnames(ids2) = c("probe_id","symbol")
ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),]

table(rownames(exp2)%in%ids2$probe_id)
exp2=exp2[rownames(exp2)%in%ids2$probe_id,]
table(duplicated(ids2$symbol))
ids2=ids2[!duplicated(ids2$symbol),]
table(rownames(exp2)%in%ids2$probe_id)
exp2=exp2[rownames(exp2)%in%ids2$probe_id,]
class(ids2$probe_id)
ids2=ids2[ids2$probe_id%in%rownames(exp2),]
dim(exp2)
exp2=exp2[match(ids2$probe_id,rownames(exp2)),]
identical(ids2$probe_id,rownames(exp2))
rownames(exp2)=ids2$symbol

identical(colnames(exp2),rownames(pd2))
save(pd2,exp2,file = 'database/GEO/GSE20140/exp2_pd2.Rdata')

eSet3 = eSet[[3]]

#(2)提取临床信息
pd3<- pData(eSet3)
colnames(pd3)
table(pd3$source_name_ch1)
pd3=pd3[,c(1,2,36)]

colnames(pd3)=c('Title','ID','Vascular_invasion')
pd3$Title
table(pd3$Vascular_invasion)
pd3=pd3[pd3$Vascular_invasion!='N/A',]
pd3$Vascular_invasion=factor(pd3$Vascular_invasion,levels = c('No','Yes'))

#(1)提取表达矩阵exp
exp3 <- exprs(eSet3)
exp3[1:4,1:4]#注意看一下，注意不要log两次，取过log的表达量的值在0-20之间
dim(exp3)
exp3=exp3[,colnames(exp3)%in%pd3$ID]
dim(exp3)
dim(pd3)
min(exp3)
exp3=log2(exp3+23)
exp3[1:4,1:4]
boxplot(exp3)

identical(pd3$ID,colnames(exp3))

# #(4)提取芯片平台编号
gpl_number <- eSet3@annotation
# 
# #(5)下载注释文件
a = getGEO(gpl_number,destdir = "database/GEO/GSE20140/")
b = a@dataTable@table
# 
colnames(b)
ids2 = b[,c("ID","Symbol")]
colnames(ids2) = c("probe_id","symbol")
ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),]
table(rownames(exp3)%in%ids2$probe_id)
exp3=exp3[rownames(exp3)%in%ids2$probe_id,]
table(duplicated(ids2$symbol))
ids2=ids2[!duplicated(ids2$symbol),]
table(rownames(exp3)%in%ids2$probe_id)
exp3=exp3[rownames(exp3)%in%ids2$probe_id,]
class(ids2$probe_id)
ids2=ids2[ids2$probe_id%in%rownames(exp3),]
dim(exp3)
exp3=exp3[match(ids2$probe_id,rownames(exp3)),]
identical(ids2$probe_id,rownames(exp3))
rownames(exp3)=ids2$symbol
identical(colnames(exp3),rownames(pd3))
save(pd3,exp3,file = 'database/GEO/GSE20140/exp3_pd3.Rdata')

#GSE25097
rm(list = ls())
options(stringsAsFactors = F)#3.6版本以后的需要用
library(GEOquery)
gse_number = "GSE25097"
eSet <- getGEO(gse_number, 
               destdir = 'database/GEO/GSE25097', 
               getGPL = F)
class(eSet)
length(eSet)
eSet1 = eSet[[1]]
#(1)提取表达矩阵exp
exp1 <- exprs(eSet1)
exp1[1:4,1:4]#注意看一下，注意不要log两次，取过log的表达量的值在0-20之间
boxplot(exp1[,1:10])
exp1=log2(exp1+1)
boxplot(exp1[,1:10])
#(2)提取临床信息
pd1<- pData(eSet1)
colnames(pd1)
table(pd1$source_name_ch1)
table(pd1$`tissue:ch1`)
pd1=pd1[,c(1:2,8)]
colnames(pd1)=c('Title','ID','Group1')
table(pd1$Group2)
pd1$Group1=ifelse(pd1$Group1=='cirrhotic','Cirrhotic',
                  ifelse(pd1$Group1=='healthy','Healthy',
                         ifelse(pd1$Group1=='non_tumor','Non-tumor','Tumor')))
pd1$Group2=ifelse(pd1$Group1=='Tumor','Tumor','Non-tumor')
pd1$Group1=factor(pd1$Group1,levels = c('Healthy','Cirrhotic','Non-tumor','Tumor'))
pd1$Group2=factor(pd1$Group2,levels = c('Non-tumor','Tumor'))
table(pd1$Group1)

identical(pd1$ID,colnames(exp1))
table(pd1$Group1)
table(pd1$Group2)
#(4)提取芯片平台编号
gpl_number <- eSet1@annotation

#(5)下载注释文件
a = getGEO(gpl_number,destdir = "database/GEO/GSE25097/")
b = a@dataTable@table

colnames(b)
ids2 = b[,c("ID","GeneSymbol")]
colnames(ids2) = c("probe_id","symbol")
ids2 = ids2[ids2$symbol!="" & !str_detect(ids2$symbol,"///"),]

table(rownames(exp1)%in%ids2$probe_id)
exp1=exp1[rownames(exp1)%in%ids2$probe_id,]
table(duplicated(ids2$symbol))
ids2=ids2[!duplicated(ids2$symbol),]
table(rownames(exp1)%in%ids2$probe_id)
exp1=exp1[rownames(exp1)%in%ids2$probe_id,]
class(ids2$probe_id)
ids2=ids2[ids2$probe_id%in%rownames(exp1),]
dim(exp1)
exp1=exp1[match(ids2$probe_id,rownames(exp1)),]
identical(ids2$probe_id,rownames(exp1))
rownames(exp1)=ids2$symbol

identical(colnames(exp1),rownames(pd1))
save(pd1,exp1,file = 'database/GEO/GSE25097/exp1_pd1.Rdata')

rm(list = ls())
options(stringsAsFactors = F)#3.6版本以后的需要用
library(GEOquery)
gse_number = "GSE104580"
eSet <- getGEO(gse_number, 
               destdir = 'database/GEO/GSE104580', 
               getGPL = F)
class(eSet)
length(eSet)
eSet = eSet[[1]]
#(1)提取表达矩阵exp
exp <- exprs(eSet)
exp[1:4,1:4]#注意看一下，注意不要log两次，取过log的表达量的值在0-20之间
