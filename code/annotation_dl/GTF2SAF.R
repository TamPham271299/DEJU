library(Rsubread)

setwd("../../annotation")
SAF <- flattenGTF("gencode.vM32.annotation.gtf", 
                  GTF.featureType="exon", 
                  GTF.attrType="gene_id", 
                  method="merge")

write.table(SAF, "gencode.vM32.annotation.flattened.exon.saf", 
            quote=F, row.names=F, sep="\t")