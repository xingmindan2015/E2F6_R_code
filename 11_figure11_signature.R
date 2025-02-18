#突变风险特征
rm(list = ls())
library(maftools) 
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(stringr)
library(deconstructSigs)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
load('figure10_mutate/lihc_mutation_clinical.Rdata')
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE)
laml=maf_clinical
mut=laml@data
mut[1:4,1:2]
mut=mut[mut$Variant_Type=='SNP',]
a = mut[,c(16,5,6,12,13)]
colnames(a)=c( "Sample","chr", "pos","ref",  "alt")
a$Sample=as.character(a$Sample)

sigs.input <- mut.to.sigs.input(mut.ref = a, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
class(sigs.input)

sigs.input[1:4,1:4]
if(!file.exists('figure11_signature/signature_xiaohua.Rdata')){
  w=lapply(unique(a$Sample), function(i){
    ## signatures.cosmic signatures.nature2013
    sample_1 = whichSignatures(tumor.ref = sigs.input, 
                               signatures.ref = signatures.cosmic, 
                               sample.id =  i, 
                               contexts.needed = TRUE,
                               tri.counts.method = 'default')
    print(c(i,which(unique(a$Sample)==i)))
    return(sample_1$weights)
  })
  w=do.call(rbind,w)
  save(w,file = 'figure11_signature/signature_xiaohua.Rdata')
}
load('figure11_signature/signature_xiaohua.Rdata')
mat = t(w)
pd=laml@clinical.data
head(pd)
colnames(pd)
pd=pd[,c('Tumor_Sample_Barcode',"Group","Cluster","Stage","Gender","Event")]
#调整pd的行顺序
pd = arrange(pd,Group,Cluster,Stage,Gender,Event)
#mat与pd一一对应
mat = mat[,match(pd$Tumor_Sample_Barcode,colnames(mat))]
identical(pd$Tumor_Sample_Barcode,colnames(mat))
#只取一部分signiture来画
s = head(names(sort(rowSums(mat),decreasing = T)),8)
mat = mat[s,]
#顶部直方图
for_bar = table(mut$Tumor_Sample_Barcode)
for_bar = for_bar[match(colnames(mat),names(for_bar))]
for_bar = as.integer(for_bar)

annotation = HeatmapAnnotation(
  mut = anno_barplot(for_bar),
  df = pd[,-1])

ht_list = Heatmap(mat, name = "mat", 
                  cluster_rows = F,
                  cluster_columns = F,
                  show_column_names = F,
                  column_split = pd$Riskgroup,
                  top_annotation = annotation)

pdf('figure11_signature/signature8.pdf',height = 10,width = 10)
draw(ht_list, heatmap_legend_side = "left", annotation_legend_side = "bottom")
dev.off()


