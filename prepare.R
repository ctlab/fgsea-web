source("./functions.R")

dir.create("./data")

library(org.Mm.eg.db)
org.Mm.eg.anno <- makeOrgAnnotation(org.Mm.eg.db)
#save(org.Mm.eg.anno, file="data/org.Mm.eg.anno.rda")

library(org.Hs.eg.db)
org.Hs.eg.anno <- makeOrgAnnotation(org.Hs.eg.db)
#save(org.Hs.eg.anno, file="data/org.Hs.eg.anno.rda")

org.annos <- list("Mus musculus"=org.Mm.eg.anno, "Homo sapiens"=org.Hs.eg.anno)
save(org.annos, file="./data/org.annos.rda")

dir.create("data/msigdb", recursive = T)
download.file("http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p2.rdata", destfile = "./data/msigdb/human_H_v5p2.rdata")

for (x in paste0("c", 1:7)) {
    download.file(paste0("http://bioinf.wehi.edu.au/software/MSigDB/human_", x, "_v5p2.rdata"),
                  destfile = paste0("./data/msigdb/human_", x ,"_v5p2.rdata"))
}

download.file("http://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata", destfile = "./data/msigdb/mouse_H_v5p2.rdata")
for (x in paste0("c", 1:7)) {
    download.file(paste0("http://bioinf.wehi.edu.au/software/MSigDB/mouse_", x, "_v5p2.rdata"),
                  destfile = paste0("./data/msigdb/mouse_", x ,"_v5p2.rdata"))
    
}


for (f in list.files("./data/msigdb", pattern = ".rdata", full.names = T)) {
    load(f)
}

pathways.Hs <- list("hallmark"=Hs.H,
                    "c1"=Hs.c1,
                    "c2.cp.biocarta"=Hs.c2[grep("^BIOCARTA_", names(Hs.c2))],
                    "c2.cp.reactome"=Hs.c2[grep("^REACTOME_", names(Hs.c2))],
                    "c2.cp.kegg"=Hs.c2[grep("^KEGG_", names(Hs.c2))],
                    "c2.cgp"=Hs.c2[grep("^(KEGG|REACTOME|BIOCARTA)_", names(Hs.c2), invert=T)],
                    "c3.mir"=Hs.c3[grep(",MIR-", names(Hs.c3))],
                    "c3.tft"=Hs.c3[grep(",MIR-", names(Hs.c3), invert=T)],
                    "c4.cm"=Hs.c4[grep("^MODULE_", names(Hs.c4))],
                    "c4.cgn"=Hs.c4[grep("^MODULE_", names(Hs.c4), invert=T)],
                    "c5"=Hs.c5,
                    "c6"=Hs.c6,
                    "c7"=Hs.c7
)

pathways.Mm <- list("hallmark"=Mm.H,
                    "c1"=Mm.c1,
                    "c2.cp.biocarta"=Mm.c2[grep("^BIOCARTA_", names(Mm.c2))],
                    "c2.cp.reactome"=Mm.c2[grep("^REACTOME_", names(Mm.c2))],
                    "c2.cp.kegg"=Mm.c2[grep("^KEGG_", names(Mm.c2))],
                    "c2.cgp"=Mm.c2[grep("^(KEGG|REACTOME|BIOCARTA)_", names(Mm.c2), invert=T)],
                    "c3.mir"=Mm.c3[grep(",MIR-", names(Mm.c3))],
                    "c3.tft"=Mm.c3[grep(",MIR-", names(Mm.c3), invert=T)],
                    "c4.cm"=Mm.c4[grep("^MODULE_", names(Mm.c4))],
                    "c4.cgn"=Mm.c4[grep("^MODULE_", names(Mm.c4), invert=T)],
                    "c5"=Mm.c5,
                    "c6"=Mm.c6,
                    "c7"=Mm.c7
)


pathways.msigdb <- list("Mus musculus"=pathways.Mm, "Homo sapiens"=pathways.Hs)
save(pathways.msigdb, file="./data/pathways.msigdb.rda")
