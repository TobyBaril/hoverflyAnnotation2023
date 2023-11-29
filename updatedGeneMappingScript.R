# load libraries
library(tidyverse)
library(data.table)
library(magrittr)

# set options
options(stringsAsFactors = FALSE, scipen = 100)

# read in top hits table and set column names
tophit <- read.table("top blast hit.txt", sep = "\t", header = T)

# read in gene annotation GFF
gff <- read.table("Episyrphus_balteatus_-GCA_945859705.1-2023_05-genes.gff3", sep = "\t", quote = "")

# column 9 has a match to the tophit "Query" column (minus the trailing .1)
# make a new column with the ID matching in both tophit and gff
# add column containing just gene ID that we can use to map our gene names to
# how do we deal with the multiple transcripts...they must match to the same gene
tophit$map <- gsub("\\.1$", "", tophit$Query)

gff$map <- gsub("ID=CDS:", "", gff$V9)
gff$map <- gsub(";.*", "", gff$map)
gff$map <- gsub("ID=.*", NA, gff$map)
gff$map <- gsub("Parent=.*", NA, gff$map)

gff$tra <- gsub(".*transcript:", "", gff$V9)
gff$tra <- gsub(";.*", "", gff$tra)
gff$tra <- gsub("ID=.*", NA, gff$tra)
gff$tra <- gsub("Parent=.*", NA, gff$tra)

# merge the tables 
gff_merge <- left_join(gff, 
                       tophit,
                       by = "map")

# add information to all rows based on transcript ID
tojoin <- gff_merge[,11:13] %>%
  distinct()
tojoin <- tojoin[! is.na(tojoin$tra),]
tojoin <- tojoin[! is.na(tojoin$Query),]

gff_merge <- gff_merge[,1:11]
gff_merge <- left_join(gff_merge,
                       tojoin)

# add information to all rows based on gene ID
gff_merge$gen <- gsub(".*gene:", "", gff_merge$V9)
gff_merge$gen <- gsub(";.*", "", gff_merge$gen)
gff_merge$gen <- gsub("ID=.*", NA, gff_merge$gen)
gff_merge$gen <- gsub("Parent=.*", NA, gff_merge$gen)

tojoin <- gff_merge[,11:14] %>%
  distinct()
tojoin <- tojoin[! is.na(tojoin$tra),]
tojoin <- tojoin[! is.na(tojoin$gen),]

for (i in 1:length(gff_merge$V1)) {
  if (is.na(gff_merge$tra[i]) & ! is.na(gff_merge$gen[i])) {
    gff_merge$Query[i] <- tojoin$Query[tojoin$gen == gff_merge$gen[i]][[1]]
    gff_merge$hit[i] <- tojoin$hit[tojoin$gen == gff_merge$gen[i]][[1]]
  }
}

# generate new GFF with the extra information in it
gff_merge$V9 <- paste0(gff_merge$V9, ";Drosophila_melanogaster_best_hit_name=", gff_merge$hit, sep = "")
gff_merge <- gff_merge[,1:9]

# write the table to a new gtf file
write.table(gff_merge,
            "Episyrphus_balteatus_annotations-GCA_945859705.1-2023_05-genes_dmelBestHitsAdded.gff",
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)



