# inference for human
ppi.infer.human <- function (target, kernel, top = 10, classifier = net.infer,
                             input = "hgnc_symbol", output = "hgnc_symbol", ...) 
{
  ensembl <- useMart("ensembl")
  human.ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
  
  # input
  new.list <- getBM(attributes = c("ensembl_peptide_id", input), 
                    filters = input, values = target, mart = human.ensembl)[, 1]
  new.list <- na.omit(new.list)
  
  # main
  ppi.pred.9606 <- classifier(new.list, kernel, ...)
  
  # output
  protein_score <- data.frame(ppi.pred.9606$top, ppi.pred.9606$score)
  format.protein <- getBM(attributes = c("ensembl_peptide_id", output),
                          filters = "ensembl_peptide_id", values = ppi.pred.9606$top, 
                          mart = human.ensembl)
  index <- match(protein_score[, 1], format.protein[, 1])
  new.protein_score <- cbind(protein_score, format.protein[index, 2])
  new.protein_score[new.protein_score == ""] <- NA
  na.omit.new.protein_score <- na.omit(new.protein_score)
  ppi.pred.9606$top <- as.vector(na.omit.new.protein_score[1:top, 3])
  ppi.pred.9606$score <- as.numeric(na.omit.new.protein_score[1:top, 2])
  ppi.pred.9606
}