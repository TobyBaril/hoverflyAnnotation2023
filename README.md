# hoverflyAnnotation2023
Workflow for annotation of Hoverfly genome (GCA_945859705.1) with D. melanogaster homologs

To match the annotations, BLASTp the D.melanogaster protein sequences against the Hoverfly protein sequences. 

# Step 1. BLASTp the D.melanogaster protein sequences against the Hoverfly protein sequences. 

Here, the query is the E. balteatus proteins, as these are what we want to match up. The subject is the D. melanogaster proteins, as these are what we want to match to. Outfmt is set to 6, which is a tabular format. The std option is to include the standard fields in the output. The qlen and slen options are to include the query and subject lengths in the output.

```
blastp -query /data/toby/hoverflyAnnotations/hoverflyAssembly/epiBalProteins.fa -subject /data/toby/hoverflyAnnotations/droMelAnnotations/dmel-all-translation-r6.46.fasta -outfmt "6 std qlen slen" -out /data/toby/hoverflyAnnotations/analysis/Dmelanogaster_protein_vs_Hoverfly_protein.blastp
```

# Step 2. Generate mapping table from fasta header to work out gene names in D. melanogaster

Take the fasta header sequences, strip out the information we don't need so we have four columns, one with gene ID (e.g DBpp0070000), one with D. melanogaster gene name (e.g. Nep3-PA), one with genbank ID (e.g AAF45370), and one with Uniprot ID (e.g Q9W5Y0)

Explanation below:
- ```grep ">" /data/toby/hoverflyAnnotations/droMelAnnotations/dmel-all-translation-r6.46.fasta```:This command uses grep to search for lines in the file dmel-all-translation-r6.46.fasta that contain the ">" character. These lines typically represent sequence identifiers in a FASTA file format.
- ```sed 's/>//'```: This command uses sed to remove the ">" character from the lines that were identified in the previous step.
- ```awk '{print $1 "\t" $5 "\t" $7 "\t" $7}'```: This command uses awk to print the first, fifth, and seventh fields from the lines that were identified in the previous step. The first field is the sequence identifier, the fifth field is the gene name, and the seventh field contains the genbank and uniprot identifiers, and so is printed twice. 
- ```awk '{OFS="\t"}{gsub(".*GB_protein:", "", $3); gsub(",.*", "", $3); gsub(".*UniProt[^:]*:", "", $4); gsub("[,;].*", "", $4); gsub(".*=", "", $2); gsub(";", "", $2); print $0}'```: This command uses awk to set the output field separator to a tab character, and then uses gsub to remove the text that we don't need from the fields that were identified in the previous step. The gsub commands are as follows:
  - ```gsub(".*GB_protein:", "", $3)```: This command removes the text "GB_protein:" and everything that precedes it from the third field.
  - ```gsub(",.*", "", $3)```: This command removes the text that follows a comma from the third field.
  - ```gsub(".*UniProt[^:]*:", "", $4)```: This command removes the text "UniProt" and everything that precedes it, and the colon that follows it, from the fourth field.
  - ```gsub("[,;].*", "", $4)```: This command removes the text that follows a comma or semicolon from the fourth field.
  - ```gsub(".*=", "", $2)```: This command removes the text that precedes an equals sign from the second field.
  - ```gsub(";", "", $2)```: This command removes the semicolon that follows the text that was removed in the previous step from the second field.
  - ```print $0```: This command prints the entire line.

```
grep ">" /data/toby/hoverflyAnnotations/droMelAnnotations/dmel-all-translation-r6.46.fasta | sed 's/>//' | awk '{print $1 "\t" $5 "\t" $7 "\t" $7}' | awk '{OFS="\t"}{gsub(".*GB_protein:", "", $3); gsub(",.*", "", $3); gsub(".*UniProt[^:]*:", "", $4); gsub("[,;].*", "", $4); gsub(".*=", "", $2); gsub(";", "", $2); print $0}' > /data/toby/hoverflyAnnotations/analysis/dmelanogaster_protein_to_geneName.txt
```

# Step 3. Read the table into R and filter for the best hits to each query sequence. Below is the Rscript.

```
# load libraries
library(tidyverse)
library(data.table)
library(magrittr)

# set options
options(stringsAsFactors = FALSE, scipen = 100)

# read in BLAST table
hits <- read.table("/data/toby/hoverflyAnnotations/analysis/Dmelanogaster_protein_vs_Hoverfly_protein.blastp",
                   sep = "\t")

# There are 1853630 hits in the raw dataset

# add column names
colnames(hits) <- c("eb_gene", "dmel_gene", "perc_id", "alignment_length", "mismatches", "gapOpen", "eb_start", "eb_end",
                    "dmel_start", "dmel_end", "evalue", "bitscore", "eb_full_length", "dmel_full_length")

# sort table by dmel_gene, dmel_start, and bitscore
hits %<>% 
  arrange(dmel_gene, -alignment_length, -bitscore)

# filter hits for those at least 50% of full dmel length
hits$perc_of_dmel_length <- (hits$alignment_length / hits$dmel_full_length) * 100
hits <- hits[hits$perc_of_dmel_length >= 50.0,]

# This reduces the dataset to 214369 potential matches
# We have removed 1639261 hits that were shorter than 50% of total gene length

# filter hits for those with a perc_id of at least 50%
# logically, this equates to a 25% identity between e_bal genes and dmel genes,
# (i.e 50% alignment length * 50% identity = 25% real matches)

hits <- hits[hits$perc_id >=50.0, ]

# This reduces the dataset to 28132 matches
# We have removed a further 186237 hits

# Now filter for the top hit for each eb_gene by seleting the hit with the highest bitscore and, if there are any duplicates,
# lowest evalue

filtered_hits <- hits %>%
  group_by(eb_gene) %>%
  slice_max(bitscore) %>%
  slice_min(evalue) %>%
  arrange(eb_gene)

# There are still some duplicates where eb_gene has matched multiple dmel genes with the same perc_id etc
# This is likely due to massive gene duplication in d_mel
# filter for distinct hits based on eb_gene, bitscore, and evalue

filtered_hits %<>% 
  distinct(eb_gene, bitscore, evalue, .keep_all = TRUE)

# We have (relatively) good matches to 6587 gene annotations

# read in the names of d_mel genes and match them up
gene_map <- read.table("analysis/dmelanogaster_protein_to_geneName.txt", sep = "\t", header = FALSE, quote = "")

# add column names
colnames(gene_map) <- c("dmel_gene", "dmel_name", "genbank", "uniprot")

# replace strings that don't match uniprot with NA
gene_map$uniprot <- gsub(".*=.*", NA, gene_map$uniprot)

# map common names to eb_genes
filtered_hits_labelled <- merge(filtered_hits,
                                gene_map,
                                by = "dmel_gene")

# read in the GTF file for Ebal and append dmel matches to final column
epibal_gtf <- read.table("hoverflyAssembly/Episyrphus_balteatus_annotations-GCA_945859705.1-2023_05-genes.gtf", sep = "\t")

# I think that gene_id matches the gene id in the fasta file...
# add colum containing just gene ID that we can use to map our gene names to
# how do we deal with the multiple transcripts...they must match to the same gene
epibal_gtf$map <- gsub(";.*", "", epibal_gtf$V9)
epibal_gtf$map <- gsub("gene_id ENSEBLG[0]*", "g", epibal_gtf$map)

# add a map column for filtered_hits_labelled
filtered_hits_labelled$map <- gsub(".t[0-9]*", "", filtered_hits_labelled$eb_gene)

# merge the tables 
epibal_gtf_merge <- left_join(epibal_gtf, 
                          filtered_hits_labelled,
                          by = "map")

# generate new GFF with the extra information in it
epibal_gtf_merge$V9 <- paste0(epibal_gtf_merge$V9, 
                              " Drosophila_melanogaster_best_hit_name=", epibal_gtf_merge$dmel_name,
                              "; Drosophila_melanogaster_best_hit_genbankID=", epibal_gtf_merge$genbank,
                              "; Drosophila_melanogaster_best_hit_uniprotID=", epibal_gtf_merge$uniprot,
                              sep = "")
epibal_gtf_new <- epibal_gtf_merge[,1:9]

# write the table to a new gtf file
write.table(epibal_gtf_new,
            "Episyrphus_balteatus_annotations-GCA_945859705.1-2023_05-genes_dmelBestHitsAdded.gtf",
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)
```

This has generated a new GTF with BLAST best hits added where available (ie for matches >50% of full length protein and >50% percent identity). This file is called ```Episyrphus_balteatus_annotations-GCA_945859705.1-2023_05-genes_dmelBestHitsAdded.gtf``` and can be found in toby's mac. All files except the GTF file needed to reproduce this annotation in the ```/src``` directory.

