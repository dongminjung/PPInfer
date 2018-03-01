# ORA barplot
ORA.barplot <- function (object, category, size, count, pvalue, top = 10, sort = NULL, 
                         decreasing = FALSE, p.adjust.methods = NULL, numChar = NULL, 
                         title = NULL, transparency = 0.5, plot = TRUE) 
{
  object <- data.frame(object)
  if (!is.null(sort))
  {
    index <- order(object[, sort], decreasing = decreasing)
    object <- object[index, ]
  }
  EnrichTab <- object[1:min(top, nrow(object)), ]
  EnrichTab <- EnrichTab[, c(category, size, count, pvalue)]
  colnames(EnrichTab) <- c("Category", "Size", "Count", "Pvalue")
  if (is.null(numChar))
  {
    numChar <- max(nchar(as.character(EnrichTab$Category)))
  }
  else
  {
    if (length(unique(substr(EnrichTab$Category, 1, numChar))) < nrow(EnrichTab))
    {
      numChar <- max(nchar(as.character(EnrichTab$Category)))
      message("Note : numChar is too small.", "\n")
    }
  }
  EnrichTab$Category <- paste(substr(EnrichTab$Category, 1, 
                                     numChar), ifelse(nchar(as.character(EnrichTab$Category)) > 
                                                        numChar, "...", ""), sep = "")
  EnrichTab[, 1] <- factor(EnrichTab[, 1], levels = EnrichTab[nrow(EnrichTab):1,1])
  GeneRatio <- EnrichTab$Count/EnrichTab$Size
  index <- which(colnames(object) == pvalue)
  if (!is.null(p.adjust.methods))
  {
    neg.log10.p.value = "-log10(adj.Pvalue)"
    adj.Pvalue <- p.adjust(object[, index], p.adjust.methods)
    adj.Pvalue <- adj.Pvalue[1:nrow(EnrichTab)]
    EnrichTab <- data.frame(EnrichTab, GeneRatio, adj.Pvalue)
    neg.log10.adj.Pvalue <- -log10(EnrichTab$adj.Pvalue)
    new.EnrichTab <- data.frame(EnrichTab, neg.log10.adj.Pvalue)
  }
  else
  {
    neg.log10.p.value = "-log10(Pvalue)"
    EnrichTab <- data.frame(EnrichTab, GeneRatio)
    neg.log10.Pvalue <- -log10(EnrichTab[,"Pvalue"])
    new.EnrichTab <- data.frame(EnrichTab, neg.log10.Pvalue)
  }
  p <- ggplot(new.EnrichTab, aes_string(x = "Category",
                                        y = ifelse(!is.null(p.adjust.methods) &
                                                     is.numeric(object[, index]),
                                                   "neg.log10.adj.Pvalue", "neg.log10.Pvalue"))) +
    geom_bar(stat = "identity", alpha = transparency) + 
    aes_string(fill = "GeneRatio") +  xlab("") + coord_flip() +
    ylab(neg.log10.p.value) + labs(title = title) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  if (plot == TRUE)
  {
    print(EnrichTab)
    return(p)
  }
  else
  {
    print(p)
    return(EnrichTab)
  }
}