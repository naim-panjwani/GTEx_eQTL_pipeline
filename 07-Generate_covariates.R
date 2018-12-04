#!/usr/bin/R

library(argparser, quietly=TRUE)
library(data.table)

p <- arg_parser("Generate covariates file")
p <- add_argument(p, "famfile", help="Fam file containing of subset of individuals for the tissue of interest; if Omni is needed, ensure a merge with the WGS bplink file has been made")
p <- add_argument(p, "covariates_file", help="The processed covariates file with platform, sex, and age information")
p <- add_argument(p, "pca_file", help="PC-AiR run eigenvectors file")
p <- add_argument(p, "peers_file", help="PEER factors covariates file")
p <- add_argument(p, "outfilename", help="Desired output filename")
argv <- parse_args(p)


#====================================================================
transpose_extra <- function(mat) {
  tmat <- transpose(mat)
  colnames(tmat) <- as.character(tmat[1,])
  tmat <- tmat[-1,]
  tmat <- cbind(ID=colnames(mat)[2:ncol(mat)], tmat)
  return(tmat)
}
#=====================================================================


fam <- fread(argv$famfile, header=F, stringsAsFactors=F)
covariates <- fread(argv$covariates_file, header=T, stringsAsFactors=F)
pca <- fread(argv$pca_file, header=F, stringsAsFactors=F)[,1:4,with=F] # including 3 PC's only
colnames(pca) <- c("IID", "PC1", "PC2", "PC3")
peers <- fread(argv$peers_file, header=T, stringsAsFactors=F)
peers <- transpose_extra(peers)

if((nrow(pca) != nrow(fam)) | (nrow(peers) != nrow(fam))) {
  warning("PCA sample size differs from fam file input")
} else if (nrow(peers) != nrow(fam)) {
  warning("PEER factors sample size differs from fam file input")
}

x1 <- merge(fam[,2,with=F], covariates, by.x="V2", by.y="SUBJID", all.x=T, sort=F)
colnames(x1)[1] <- "IID"
x2 <- merge(x1, pca, by="IID", all.x=T, sort=F)
x3 <- merge(x2, peers, by.x="IID", by.y="ID", all.x=T, sort=F)

final_covars <- transpose_extra(x3)

write.table(final_covars, argv$outfilename, quote=F, row.names=F, col.names=T, sep="\t")

