library(openxlsx)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0){
    stop("No input file given\n", call.=FALSE)
}
if (length(args)==1){
    parts = strsplit(args[1], split="\\.")[[1]]
    args[2] = paste(c(parts[-length(parts)], "tsv"), collapse=".")
}

df <- read.xlsx(args[1])
write.table(df, args[2], sep="\t", quote=FALSE, row.names=FALSE)

