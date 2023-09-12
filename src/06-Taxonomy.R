###Produce stack bar chart for taxonomy
#Import data files
gene.rel <- readRDS(file = here("data","intermediate","gene.rel.RDS"))
annot <- readRDS(file = here("data","intermediate","annot.RDS")) 
map.s <- readRDS(file = here("data","intermediate", "map.RDS"))
#Sort
gene.rel.s <- gene.rel[,order(colnames(gene.rel))]
annot.s <- annot[order(annot$gene_id),] 
#Remove extra line in the annot file
annot.s <- annot.s[annot.s$gene_id %in% colnames(gene.rel.s),]
#Check sorting
sum(annot.s$gene_id == colnames(gene.rel.s)) #13108187 
#Create combined data frame
gene.tax <- data.frame(cbind(t(gene.rel.s), annot.s$tax_phylum)) ##3108187 obs of 21 var
#Remove intermediate files
rm(annot)
rm(annot.s)
rm(gene.rel)
rm(gene.rel.s)

#Make numeric
gene.tax.2 <- data.frame(sapply(gene.tax[,1:20],as.numeric), gene.tax$V21)
#Summarize at the phylum level - will make sum of gene relative abundance for each phyla
tax.phylum.summary <- gene.tax.2 %>%
  group_by(gene.tax.V21) %>%
  summarise(across(.cols=everything(), ~ sum(.x, na.rm = TRUE)))
saveRDS(tax.phylum.summary, file = here("data", "intermediate", "tax.phylum.summary.RDS"))

#Fix the headers
colnames(tax.phylum.summary) <- gsub("X", "", colnames(tax.phylum.summary)) 
colSums(tax.phylum.summary[,2:21]) #Check if in rel abundance - should all be 1
#Keep only above 0.5% on average across all samples
tax.phylum.summary.top <- tax.phylum.summary[rowMeans(tax.phylum.summary[,2:21])>0.005, ] 
#Prepare for ggplot
#Transpose, fix lines and headers
tax.phylum.summary.top.t <- as.data.frame(t(tax.phylum.summary.top)) 
rownames(tax.phylum.summary.top.t) <- gsub("X", "", rownames(tax.phylum.summary.top.t)) 
colnames(tax.phylum.summary.top.t) <- tax.phylum.summary.top.t[1,] 
tax.phylum.summary.top.t <- tax.phylum.summary.top.t[-1,]
#Sort
tax.phylum.summary.top.t.s <- tax.phylum.summary.top.t[order(row.names(tax.phylum.summary.top.t)),]#Sort 
#Check sorting
sum(row.names(map.s)==row.names(tax.phylum.summary.top.t.s))#20
#Make numeric
tax.phylum.summary.top.t.s <- sapply(tax.phylum.summary.top.t.s,as.numeric)
row.names(tax.phylum.summary.top.t.s) <- row.names(tax.phylum.summary.top.t)
#Addd mapping file
tax.map <- data.frame(map.s,tax.phylum.summary.top.t.s)
saveRDS(tax.map, file = here("data", "intermediate", "tax.map.RDS"))
#Make long for ggplot
tax.map.long <- gather(tax.map,Phylum,relabund,4:12) #transform in long format for ggplot
#Remove NULL
tax.map.long.nonull <- tax.map.long[tax.map.long$Phylum != "NULL.",]
tax.map.long.nonull$perc_SWHC <- factor(tax.map.long.nonull$perc_SWHC, c("50", "5"))#Reorder manually

#Plot
palette(c(brewer.pal(n = 9, name = "Set1"),"lightgrey", "black", "darkred", "darkblue", "darkgreen", "purple4", "darkgrey", "white"))
stack.phylum <- ggplot(tax.map.long.nonull, aes(fill = Phylum, y = relabund, x = perc_SWHC)) + 
  geom_bar( stat = "summary", fun ="mean", position = "stack") +
  ylab("Fraction of reads") + 
  scale_fill_manual(values = palette(), guide = guide_legend(label.theme = element_text(face = "italic", size = 8))) +
  theme_bw() +
  scale_y_continuous( expand = c(0,0)) +
  scale_x_discrete(name = "% SWHC") +
  facet_wrap(vars(SoilType), labeller = labeller(SoilType = c("IR" = "intermittent", "NI" = "continuous")))
stack.phylum
