library(tidyverse)
library(data.table)
library(magrittr)

options(scipen = 100, stringsAsFactors = FALSE )

# annotate the Differential expression data with stuff from GFF

# read in the required tables
gff_merge <- read.table("PATH/TO/GFF_ANNOTATION", sep = "\t", quote = "")
dif_expr <- read.table("PATH/TO/TABULAR/DESEQ2", sep = "\t")

# Filter the annotation file for entries that are genes (ie select things in column 3 that are genes)
gff_merge_forDEseq <- gff_merge[gff_merge$V3 == "gene",]

# generate a column containing the gene ID that can be used to map to the DEseq2 table, and call this "map"
gff_merge_forDEseq$map <- gsub(".*gene_id=", "", gff_merge_forDEseq$V9)
gff_merge_forDEseq$map <- gsub(";.*", "", gff_merge_forDEseq$map)

# rename the first column of the DEseq2 output to "map" as it contains the same information as the GFF
# and can be used to join the two tables
colnames(dif_expr) <- c("map", "c2", "c3", "c4", "c5", "c6", "c7")

# join the columns we want from the GFF to the DEseq2 table 
dif_expr_app <- left_join(dif_expr, gff_merge_forDEseq[,c("V1","V4","V5", "V7", "V3", "map", "V9")])

# delete all the information in the attributes column up to the bit we want (which is right at the end)
dif_expr_app$V9 <- gsub(".*Drosophila", "Drosophila", dif_expr_app$V9)

# save the table in tab delimited format under the name specified in the path
write.table(dif_expr_app,
            "PATH/TO/SAVE/FILE.txt",
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)