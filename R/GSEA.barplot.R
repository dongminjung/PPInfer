# GSEA barplot
GSEA.barplot <- function (object, category, score, pvalue, top = 10, sort = NULL,
                          decreasing = FALSE, numChar = NULL, title = NULL,
                          transparency = 0.5, plot = TRUE)
{
  object <- data.frame(object)
  if (!is.null(sort))
  {
    index <- order(object[, sort], decreasing = decreasing)
    object <- object[index, ]
  }
  EnrichTab <- object[1:min(top, nrow(object)), ]
  EnrichTab <- EnrichTab[, c(category, score, pvalue)]
  EnrichTab.colnames <- colnames(EnrichTab)
  if (is.null(numChar))
  {
    numChar <- max(nchar(as.character(EnrichTab[, 1])))
  }
  else
  {
    if (length(unique(substr(EnrichTab[, 1], 1, numChar))) < nrow(EnrichTab))
    {
      numChar <- max(nchar(as.character(EnrichTab[, 1])))
      message("Note : numChar is too small.", "\n")
    }
  }
  EnrichTab[, 1] <- paste(substr(EnrichTab[, 1], 1, numChar), 
                          ifelse(nchar(as.character(EnrichTab[, 1])) > numChar, 
                                 "...", ""), sep = "")
  EnrichTab[, 1] <- factor(EnrichTab[, 1], levels = EnrichTab[nrow(EnrichTab):1,1])
  neg.log10.pvalue <- -log10(EnrichTab[, 3])
  p <- ggplot(data.frame(EnrichTab, neg.log10.pvalue),
              aes_string(x = EnrichTab.colnames[1], y = "neg.log10.pvalue")) + 
    geom_bar(stat = "identity", alpha = transparency) + coord_flip() + 
    aes_string(fill = EnrichTab.colnames[2]) + 
    scale_fill_gradient2(low = 'blue', mid = "green", high = 'red',
                         limits = c(min(EnrichTab[, score]), max(EnrichTab[, score]))) +
    xlab("") + ylab(paste("-log10(", EnrichTab.colnames[3], ")", sep = "")) +
    labs(title = title) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  if (plot == TRUE) {
    print(EnrichTab)
    return(p)
  }
  else {
    print(p)
    return(EnrichTab)
  }
}