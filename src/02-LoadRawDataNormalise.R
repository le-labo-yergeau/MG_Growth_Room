###Load raw data, normalize and save as intermediate files for other scripts

#Mapping file
map <- read.table(file = here("data", "raw", "mapping_file.tsv"), row.names = 1, header = T, sep = "\t", comment.char = "") #20 obs in 3 var
map.s <- map[order(row.names(map)),]
map.s$perc_SWHC <- as.character(map.s$perc_SWHC)
map.s$block <- as.character(map.s$block)
saveRDS(map.s, file = here("data", "intermediate", "map.RDS"))

#Plant fresh biomass
plant <- read.table(file = here("data", "raw", "plant_fresh_biomass.txt"),  header = T, sep = "\t", comment.char = "") #160 obs in 3 var
saveRDS(plant, file = here("data", "intermediate", "plant.RDS"))

#Other plant traits
trait <- read.table(file = here("data", "raw", "Traits.txt"),  header = T, sep = "\t", comment.char = "") #160 obs in #10 var
saveRDS(trait, file = here("data", "intermediate", "trait.RDS"))

#Annotations
annot <- fread(file = here("data", "raw", "annotations.tsv"), header = T, sep = "\t", quote = "") # 13108188 obs of 20 variables
saveRDS(annot, file = here("data", "intermediate", "annot.RDS"))

#Gene abundance
gene <- data.frame(read.table(file = here("data", "raw", "merged_gene_abundance.tsv"), row.names = 1, sep="\t", header=T)) #13,108,187 obs of 20 variables
colnames(gene) <- gsub("X", "", colnames(gene))
gene <- gene[,order(colnames(gene))]
saveRDS(gene, here("data", "intermediate", "gene.RDS"))

##Normalize with internal standard
#Retrieve the genes related to the standard from mapping on Thermus genome
#thermus <- data.frame(read.table(here("data", "raw", "qc_mapping_stats.tsv"),row.names = 1, sep="\t",header=T)) #21 obs. of 10 var
#thermus$properlyPaired <- gsub(",", "", thermus$properlyPaired)
#sum.gene.thermus <- data.frame(as.numeric(thermus$properlyPaired))
#rownames(sum.gene.thermus) <-rownames(thermus)
#sum.gene.thermus$toto <- rep("toto",20)
#sum.gene.thermus <- sum.gene.thermus[order(rownames(sum.gene.thermus)),]

#Get extraction masses
#gram <- data.frame(fread(here("data", "raw", "gram.txt"), sep="\t", header=T)) #20 obs of 3 var
#gram <- gram[order(gram$sample),]

#check ordering of three files
#colnames(gene)==gram$sample
#colnames(gene)==rownames(sum.gene.thermus)
#rownames(sum.gene.thermus)==gram$sample

#Normalize all genes
##          1 ng (amount of ng added per sample)*6.022×10^23 (avogadro nr) / 
#           2127482 (length in bp of genome)*10^9 (to get to g)*650 (weight in Dalton of 1 dsDNA bp)
#           which results in 4.35E+06 genomes added per extraction.
#           This amounts to an addition of 9.2626E+12 bp (genome size*copies added)
#           Which would then have been fragmented in 100bp reads 
#           If sequenced to extinction, we would have found 9.2646E+10 fragments (bp added/100bp)
#           Thus, we multiply actual reads of all genes by 9.2646E+10 and divide by number of Thermus reads retrieved
#           Then divide by the number of gram of soil used in extraction.
#norm.vect <- 9.2646e10/sum.gene.thermus[,1] #fraction of Thermus reads retrieved per sample
#norm.vect <- norm.vect/gram[,3] #fraction of Thermus reads retrieved per gram of soil

#gene.norm <- data.frame(apply(gene, 1, "*", norm.vect)) #13,108,187 obs, 20 var
#gene.norm[1:10,1:10]

#Sanity check
#gene.norm[18,344555] #1882758
#gene[344555,18]*norm.vect[18] #1882758

#save intermediate
#saveRDS(gene.norm, file = here("data", "intermediate", "gene.norm.RDS"))

##Normalize with relative
colsum.vect <- colSums(gene)
gene.rel <- data.frame(apply(gene, 1, "/", colsum.vect)) #20 obs, 13,108,187 var
gene.rel[1:10,1:10]

#Sanity check
gene.rel[18,123457] #3.72065e-08
gene[123457,18]/colsum.vect[18] #3.72065e-08

#save intermediate
saveRDS(gene.rel, file = here("data", "intermediate", "gene.rel.RDS"))
saveRDS(colsum.vect, file = here("data", "intermediate", "colsum.vect.RDS"))

##MAGs
#MAG table
mag.table <- read.table(file = here("data", "raw", "MAG.otu_table.tsv"), row.names = 1, header = T, sep = "\t", comment.char = "") #300 obs in 21 var
colnames(mag.table) <- gsub("X", "", colnames(mag.table)) #Remove X in column names
#move taxonomy to other object
mag.tax <- data.frame(mag.table[,21])
colnames(mag.tax) <- "Taxonomy"
row.names(mag.tax) <- row.names(mag.table)
row.names(mag.tax) <- gsub("-", ".", row.names(mag.tax))
mag.tax.sep <- mag.tax %>% separate(Taxonomy, c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")
mag.tax.sep <- mag.tax.sep[,-1] #Fix bug
mag.tax.sep <- data.frame(lapply(mag.tax.sep, function(x) gsub(".__", "", x)))
row.names(mag.tax.sep) <- row.names(mag.table)
mag.table <- mag.table[,-21] #remove taxonomy column
#Relativize
mag.table.s <- mag.table[,order(colnames(mag.table))]#Sort
colnames(mag.table.s) == colnames(gene) #check ordering
mag.table.rel <- data.frame(apply(mag.table.s, 1, "/", colsum.vect)) #20 obs, 300 var
#Sanity check
mag.table.rel[18,254] #4.46478e-06
mag.table.s[254,18]/colsum.vect[18] #4.46478e-06

#MAG link file (MAG-contig-gene)
mag.link <- read.table(file = here("data", "raw", "link.tsv"), header = F, sep = "\t", comment.char = "") #1,183,664 obs in 3 var
colnames(mag.link) <- c("MAG", "contig", "gene")
mag.link$MAG <- gsub("-",".", mag.link$MAG)

#MAG quality from checkM
mag.qual <- read.table(file = here("data", "raw", "out_checkm.txt"), row.names = 1,  header = T, sep = "\t", comment.char = "") #300 obs of 13 variables
#row.names(mag.qual) <- gsub("-", ".", row.names(mag.qual))
mag.qual <- mag.qual[order(row.names(mag.qual)),]

#save intermediate
saveRDS(mag.table.rel, file = here("data", "intermediate", "mag.table.rel.RDS"))
saveRDS(mag.link, file = here("data", "intermediate", "mag.link.RDS"))
saveRDS(mag.tax.sep, file = here("data", "intermediate", "mag.tax.RDS"))
saveRDS(mag.qual, file = here("data", "intermediate", "mag.qual.RDS"))
