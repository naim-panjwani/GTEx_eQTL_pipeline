#!/usr/bin/R

library(argparser, quietly=TRUE)
p <- arg_parser("Add standard error column to file")
p <- add_argument(p, "filename", help="File name")
argv <- parse_args(p)

data <- read.table(argv$filename, header=T, stringsAsFactors=F)
standard_error <- abs(data$slope) / qnorm(data$pval_nominal / 2, lower=F)

newdata <- cbind(data, slope_se=standard_error)
fileoutname = gsub(".txt.gz", "_slopeSE.txt", argv$filename)
write.table(newdata, fileoutname, quote=F, row.names=F, col.names=T, sep="\t")

