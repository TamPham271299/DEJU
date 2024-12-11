DEU_analysis_limma <- function(r_count, annot, group, design, mode, REF, p) {
  
  message("Constructing DGElist object ...")  
  y <- DGEList(counts=r_count, genes=annot, group=group)
  colnames(y) <- gsub("[.].*$", "", colnames(y))
  
  if (mode == "case_study") y <- filter_gene(y, REF)
  
  message("Filtering exons with low mapping reads ...")
  keep <- filterByExpr(y, group=group)
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  message("Normalizing lib sizes ...")
  y <- normLibSizes(y)
  
  message("Transforming data for LM ...")
  v <- voom(y, design, plot=FALSE)
  
  message("Fitting LM for design matrix ...")
  fit <- lmFit(v, design)
  
  message("Making contrast ...")
  cmd <- paste("makeContrasts(", p, ",levels=design)", sep="")
  contr <- eval(parse(text=cmd))
  cfit <- contrasts.fit(fit, contr)
  
  message("Running diffSplice ..")
  sp <- diffSplice(cfit, geneid="GeneID", robust=TRUE, exonid="Start")
  
  return(list('sp'=sp, 'y'=y, 'v'=v))
  
}