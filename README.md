# hoverflyAnnotation2023
Workflow for annotation of Hoverfly genome (GCA_945859705.1) with D. melanogaster homologs

To match the annotations, BLASTp the D.melanogaster protein sequences against the Hoverfly protein sequences. 

# Step 1. BLASTp the D.melanogaster protein sequences against the Hoverfly protein sequences. 

Here, the query is the E. balteatus proteins, as these are what we want to match up. The subject is the D. melanogaster proteins, as these are what we want to match to. Outfmt is set to 6, which is a tabular format. The std option is to include the standard fields in the output. The qlen and slen options are to include the query and subject lengths in the output.

```
blastp -query /data/toby/hoverflyAnnotations/hoverflyAssembly/epiBalProteins.fa -subject /data/toby/hoverflyAnnotations/droMelAnnotations/dmel-all-translation-r6.46.fasta -outfmt "6 std qlen slen" -out /data/toby/hoverflyAnnotations/analysis/Dmelanogaster_protein_vs_Hoverfly_protein.blastp
```

# Step 2. Read the table into R and filter for the best hits to each query sequence. Below is the Rscript.



