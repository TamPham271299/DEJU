OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)

GTF_f <- as.character(args[['GTF']])
# GTF_f <- "../../annotation/gencode.vM32.primary_assembly.annotation.gtf.gz"

message('Extracting exon information from GTF ...')
GTF <- import(GTF_f)
exon_dt <- subset(GTF, type == "exon")
exon_info <- data.frame(chr = seqnames(exon_dt),
                        start = start(exon_dt),
                        end = end(exon_dt),
                        strand = strand(exon_dt),
                        geneID = mcols(exon_dt)$gene_id,
                        transcriptID = mcols(exon_dt)$transcript_id)

# write.table(exon_info, "../../annotation/gencode.vM32.primary_assembly.annotation.exonInfo.tsv",
#             quote=F, row.names=F, sep="\t")

message('Converting exon annotation SAF to SJ database ...')
# exon_f <- "/vast/projects/lab_chen/tam/ref_genome/Mus_musculus/Gencode/gencode.vM32.primary_assembly.annotation.exonInfo.tsv"
# exon_info <- read.table(exon_f, header=T)

SJ <- exon_info %>%
  group_by(transcriptID) %>%
  arrange(chr, start, end) %>%
  mutate(start = lead(start)) %>%
  slice(-n()) %>%
  ungroup()

SJ <- SJ[c("geneID", "chr", "end", "start", "strand")]
colnames(SJ)[3:4] <- c("start", "end")
SJ$juncID <- paste(SJ$chr, SJ$start, SJ$end, sep="_")

SJ.1 <- SJ %>%
  group_by(geneID) %>%
  distinct(juncID, .keep_all = TRUE) %>%
  ungroup() %>%
  count(juncID, name = "freq") %>%
  arrange(geneID, chr, start, end)

write.table(SJ.1, "../../annotation/gencode.vM32.primary_assembly.annotation.SJdatabase.tsv",
            quote=F, row.names=F, sep="\t")

