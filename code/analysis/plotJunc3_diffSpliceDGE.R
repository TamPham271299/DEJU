# This code is written originally by Yunshun Chen for diffSplice.
# Edited by Tam Pham for diffSpliceDGE on 06 March 2024.
plotJunc <- function(fit, geneid, genecolname=NULL, FDR=0.05, annotation=NULL)
{
  if(is.null(genecolname)) 
    genecolname <- fit$genecolname
  else
    genecolname <- as.character(genecolname)
  
  geneid <- as.character(geneid)
  i <- fit$genes[, genecolname]==geneid
  i[is.na(i)] <- FALSE
  if(!any(i)) stop(paste0(geneid, " not found."))
  
  # Subsetting
  genes <- fit$genes[i,]
  strand <- genes$Strand[1]
  p.value <- fit$exon.p.value[i]
  coefficients <- fit$coefficients[i]
  
  # Sorting	
  o <- order(genes$Start, genes$End)
  genes <- genes[o,]
  p.value <- p.value[o]
  coefficients <- coefficients[o]
  
  # Exons
  IsExon <- genes$Length > 1
  genes.e <- genes[IsExon, ]
  
  # Missing exons
  if(!is.null(annotation)){
    if(is.null(annotation$Length)) annotation$Length <- annotation$End - annotation$Start + 1L
    sel <- annotation$GeneID == genes$GeneID[1]
    genes.e2 <- annotation[sel, ]
    genes.e2 <- genes.e2[order(genes.e2$Start), ]
    m <- match(genes.e$Start, genes.e2$Start)
    genes.e <- genes.e2
  }
  
  # Introns
  Start.i <- genes.e$End[-nrow(genes.e)]+1
  End.i <- genes.e$Start[-1]-1
  genes.i <- genes.e[-1,]
  genes.i$Start <- Start.i
  genes.i$End <- End.i
  genes.i$Length <- genes.i$End - genes.i$Start + 1
  
  # Extend the gene range
  if(any(!IsExon)){
    genes.j <- genes[!IsExon, ]
    if(min(genes.j$Start) < min(genes.e$Start)){
      intron <- genes.i[1,,drop=FALSE]
      intron$Start <- min(genes.j$Start)
      intron$End <- min(genes.e$Start) - 1
      intron$Length <- intron$End - intron$Start + 1
      genes.i <- rbind(intron, genes.i)
    }
    if(max(genes.j$End) > max(genes.e$End)){
      intron <- genes.i[1,,drop=FALSE]
      intron$Start <- max(genes.e$End) + 1
      intron$End <- max(genes.j$End)
      intron$Length <- intron$End - intron$Start + 1
      genes.i <- rbind(genes.i, intron)
    }
  }
  
  # Combine introns and exons
  genes.ie <- rbind(cbind(genes.e, Flag="Exon"), cbind(genes.i, Flag="Intron"))
  genes.ie <- genes.ie[order(genes.ie$Start), ]
  
  # Pseudo length
  pseudo.length <- (genes.ie$Length)^.5
  pseudo.pos <- cumsum((genes.ie$Length)^.5)
  pseudo.start <- c(0, pseudo.pos[-nrow(genes.ie)])
  pseudo.end <- pseudo.pos
  genes.ie <- cbind(genes.ie, pseudo.start=pseudo.start, pseudo.end=pseudo.end, pseudo.length=pseudo.length)
  
  # Junctions
  if(any(!IsExon)){
    genes.j <- cbind(genes.j, pseudo.start=0, pseudo.end=0)
    for(j in 1:nrow(genes.j)){
      k <- which(genes.j$Start[j] <= genes.ie$End)[1]
      genes.j$pseudo.start[j] <- genes.ie$pseudo.end[k] - (genes.ie$End[k] - genes.j$Start[j]) / genes.ie$Length[k] * genes.ie$pseudo.length[k]
      k <- which(genes.j$End[j] <= genes.ie$End)[1]
      genes.j$pseudo.end[j] <- genes.ie$pseudo.end[k] - (genes.ie$End[k] - genes.j$End[j]) / genes.ie$Length[k] * genes.ie$pseudo.length[k]
    }
  }
  
  # Setup plot
  GeneStart <- min(genes.ie$pseudo.start)
  GeneEnd <- max(genes.ie$pseudo.end)
  gene.length <- GeneEnd - GeneStart
  frame()
  plot.window(xlim=c(GeneStart, GeneEnd), ylim=c(-0.7, 0.7))
  title(main=paste0(geneid, " (", genes$Strand[1], ")"))
  
  # Plot gene range
  rect(xleft=GeneStart, xright=GeneEnd, ybottom=-0.02, ytop=0.02, col="gray", border="gray")
  if(strand=="+"){
    tx.left <- "5'"
    tx.right <- "3'"
  } else {
    tx.left <- "3'"
    tx.right <- "5'"
  }
  text(x=-0.02*gene.length, y=0.1, labels=tx.left)
  text(x=1.02*gene.length, y=0.1, labels=tx.right)
  
  # Significance
  up <- coefficients > 0
  down <- coefficients < 0
  IsSig <- p.adjust(p.value, method="holm") < FDR
  IsSig.j <- IsSig[!IsExon]
  down.j <- down[!IsExon]
  
  # Colouring
  col <- rep("black", sum(i))
  col[up & IsSig] <- "#D62728"
  col[down & IsSig] <- "#1F77B4"
  col.e <- col[IsExon]
  col.j <- col[!IsExon]
  
  # Plot exons
  ex <- genes.ie$Flag=="Exon"
  if(!is.null(annotation)){
    col.e2 <- rep("grey", nrow(genes.e))
    col.e2[m] <- col.e
    col.e <- col.e2
  }
  rect(xleft=genes.ie$pseudo.start[ex], xright=genes.ie$pseudo.end[ex], ybottom=-0.1,ytop=0.1, col=col.e, border=col.e)
  
  # Plot junctions
  if(any(!IsExon)){
    y <- rep(0.11, sum(!IsExon))
    h <- rep(0.25, sum(!IsExon))
    h[IsSig.j] <- 0.4
    y[down.j] <- -y[down.j]
    h[down.j] <- -h[down.j]
    for(i in 1:sum(!IsExon)) {
      x1 <- genes.j$pseudo.start[i]
      x2 <- genes.j$pseudo.end[i]
      curve(4*h[i]/(2*x1*x2 - x1^2 - x2^2)*(x-x1)*(x-x2)+y[i], from=x1, to=x2, add=TRUE, col=col.j[i], lwd=2)
    }
  }
  
  # axis
  if(genes$Strand[1]=="+")
    labels <- paste0("Exon.", 1:length(col.e))
  else
    labels <- paste0("Exon.", length(col.e):1)
  axis(side=1, at=(genes.ie$pseudo.start[ex]+genes.ie$pseudo.end[ex])/2, las=2, labels=labels)
  
  invisible()
}

