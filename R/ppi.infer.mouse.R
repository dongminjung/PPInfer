# inference for mouse
ppi.infer.mouse <- function (target, kernel, top = 10, classifier = net.infer,
                             input = "mgi_symbol", output = "mgi_symbol", ...) 
{
  ensembl <- useMart("ensembl")
  mouse.ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
  
  # input
  new.list <- getBM(attributes = c("ensembl_peptide_id", input), 
                    filters = input, values = target, mart = mouse.ensembl)[,1]
  new.list <- na.omit(new.list)
  
  # main
  ppi.pred.10090 <- classifier(new.list, kernel, ...)
  
  # output
  protein_score <- data.frame(ppi.pred.10090$top, ppi.pred.10090$score)
  format.protein <- getBM(attributes = c("ensembl_peptide_id", output),
                          filters = "ensembl_peptide_id", values = ppi.pred.10090$top, 
                          mart = mouse.ensembl)
  index <- match(protein_score[, 1], format.protein[, 1])
  new.protein_score <- cbind(protein_score, format.protein[index, 2])
  new.protein_score[new.protein_score == ""] <- NA
  na.omit.new.protein_score <- na.omit(new.protein_score)
  ppi.pred.10090$top <- as.vector(na.omit.new.protein_score[1:top, 3])
  ppi.pred.10090$score <- as.numeric(na.omit.new.protein_score[1:top, 2])
  ppi.pred.10090
}