savePathways <- function(pathways, fileName) {
    res <- sapply(names(pathways),function(x) paste(x,paste(pathways[[x]],collapse='\t'), collapse = '\t'))
    lapply(res, write, fileName, append=TRUE)
}
