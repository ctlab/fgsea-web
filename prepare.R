source("./functions.R")
library(org.Mm.eg.db)
org.Mm.eg.anno <- makeOrgAnnotation(org.Mm.eg.db)
save(org.Mm.eg.anno, file="org.Mm.eg.anno.rda")

library(org.Hs.eg.db)
org.Hs.eg.anno <- makeOrgAnnotation(org.Hs.eg.db)
save(org.Hs.eg.anno, file="org.Hs.eg.anno.rda")
