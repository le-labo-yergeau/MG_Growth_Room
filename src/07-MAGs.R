#Find MAGs that are influenced by the treatments
#Then look for traits in the genomes

#Load data
mag.table.rel <- readRDS(file = here("data", "intermediate", "mag.table.rel.RDS")) #20 obs of 714 var
mag.link <- readRDS(file = here("data", "intermediate", "mag.link.RDS")) #1833497 obs of 3 var
mag.tax <- readRDS(file = here("data", "intermediate", "mag.tax.RDS")) #1833497 obs of 3 var
map <- readRDS(file = here("data", "intermediate", "map.RDS")) #20 obs of 3 variables
mag.qual <- readRDS(file = here("data", "intermediate", "mag.qual.RDS")) #714 obs of 26 variables

#Keep only MAGs with >50% completion and <10% contamination (medium quality)
sum(row.names(mag.qual) == colnames(mag.table.rel))#check ordering should be = 714
mag.table.rel <- mag.table.rel[,(mag.qual$Completeness>50 & mag.qual$Contamination<10)] #20obs of 35 variables
mag.table.rel <- mag.table.rel[,-5] #Remove spiked-in Thermus thermophilus, 20 obs of 34 variables
mag.tax <- mag.tax[row.names(mag.tax)%in%colnames(mag.table.rel),] #Do same for mag.tax, 34 obs of 7 var.
mag.link <- mag.link[mag.link$MAG%in%colnames(mag.table.rel),] #Do same for mag.link, 114889 obs of 3 var.
mag.list <- colnames(mag.table.rel)
  
#Check sorting
row.names(map) == row.names(mag.table.rel)

#Anova loop for the 34 MAGs
i <- 1
list.MAG.SWHC <- vector()
list.MAG.inter <- vector()
list.MAG.soil <- vector()

while (i<=ncol(mag.table.rel)) {
  anova <- aov(mag.table.rel[,i]~map$SoilType*map$perc_SWHC+map$block)
  if(summary(anova)[[1]][1,5]<(0.05/34)){
    if(mean(mag.table.rel[map$SoilType=="IR",i])< mean(mag.table.rel[map$SoilType=="NI",i])){
      list.MAG.soil <- append(list.MAG.soil, c(colnames(mag.table.rel)[i], "continuous"))
    }else{list.MAG.soil <- append(list.MAG.soil, c(colnames(mag.table.rel)[i], "intermittent"))}  
  }else{list.MAG.soil <-append(list.MAG.soil, c(colnames(mag.table.rel)[i], "unaffected"))}
  if(summary(anova)[[1]][2,5]<(0.05/34)){
    if(mean(mag.table.rel[map$perc_SWHC=="5",i]) < mean(mag.table.rel[map$perc_SWHC=="50",i])){
      list.MAG.SWHC <- append(list.MAG.SWHC, c(colnames(mag.table.rel)[i], "50"))
    }else{list.MAG.SWHC <- append(list.MAG.SWHC, c(colnames(mag.table.rel)[i], "5"))}  
  }else{list.MAG.SWHC <-append(list.MAG.SWHC, c(colnames(mag.table.rel)[i], "unaffected"))}
  if(summary(anova)[[1]][4,5]<(0.05/34)){
    list.MAG.inter <- append(list.MAG.inter, colnames(mag.table.rel)[i])#Only 2
  }
  i <- i+1
}
df.MAG.soil <- data.frame(matrix(list.MAG.soil, nrow = 34, ncol = 2, byrow = TRUE))
colnames(df.MAG.soil) <- c("MAG", "Status")
df.MAG.SWHC <- data.frame(matrix(list.MAG.SWHC, nrow = 34, ncol = 2, byrow = TRUE))
colnames(df.MAG.SWHC) <- c("MAG", "Status")

#Completion and contamination level, contig length for different status - soil
mean(mag.qual[row.names(mag.qual)%in%df.MAG.soil[df.MAG.soil$Status=="intermittent",]$MAG,]$Completeness) #69.585
mean(mag.qual[row.names(mag.qual)%in%df.MAG.soil[df.MAG.soil$Status=="continuous",]$MAG,]$Completeness) #79.18688
mean(mag.qual[row.names(mag.qual)%in%df.MAG.soil[df.MAG.soil$Status=="unaffected",]$MAG,]$Completeness) #69.593
mean(mag.qual[row.names(mag.qual)%in%df.MAG.soil[df.MAG.soil$Status=="intermittent",]$MAG,]$Contamination) #5.9375
mean(mag.qual[row.names(mag.qual)%in%df.MAG.soil[df.MAG.soil$Status=="continuous",]$MAG,]$Contamination) #4.233125
mean(mag.qual[row.names(mag.qual)%in%df.MAG.soil[df.MAG.soil$Status=="unaffected",]$MAG,]$Contamination) #4.107
mean(mag.qual[row.names(mag.qual)%in%df.MAG.soil[df.MAG.soil$Status=="intermittent",]$MAG,]$total_contig_length) #3,847,075
mean(mag.qual[row.names(mag.qual)%in%df.MAG.soil[df.MAG.soil$Status=="continuous",]$MAG,]$total_contig_length) #3,398,789
mean(mag.qual[row.names(mag.qual)%in%df.MAG.soil[df.MAG.soil$Status=="unaffected",]$MAG,]$total_contig_length) #3,116,643


##MAGs affected by soil history
#Taxonomy
sum(row.names(mag.tax) == df.MAG.soil$MAG)#34
df.mag.soil.tax <- cbind(df.MAG.soil,mag.tax)
#Make a summary
df.mag.summary <- df.mag.soil.tax %>%
  group_by(Phylum, Status) %>%
  summarise(count = n())   %>%
  filter(count > 0)
df.mag.summary$Phylum <- gsub("p__","", df.mag.summary$Phylum)
df.mag.summary$Phylum <- gsub("NULL","unassigned", df.mag.summary$Phylum)
others <- data.frame(matrix(c(rep("Others",3), 
                              "continuous", "intermittent", "unaffected", 
                              16-sum(df.mag.summary[df.mag.summary$Status=="continuous",3]),
                              8-sum(df.mag.summary[df.mag.summary$Status=="intermittent",3]),
                              34-24-sum(df.mag.summary[df.mag.summary$Status=="unaffected",3])),
                              nrow = 3, ncol = 3))
colnames(others) <- colnames(df.mag.summary)
others$count <- as.numeric(others$count)
df.mag.summary.others <- rbind(df.mag.summary,others)

#Plot
palette(c(brewer.pal(n = 9, name = "Set1"),"lightgrey", "black", "darkred", "darkblue", "darkgreen", "purple4", "darkgrey", "white"))
stack.MAG.tax <- ggplot(df.mag.summary.others, aes(fill = Phylum, y = count, x = Status)) + 
  geom_bar( stat = "summary", fun ="mean", position = "stack") +
  ylab("Number of MAGs") + 
  scale_fill_manual(values = palette(), guide = guide_legend(title = "Phylum", label.theme = element_text(face = "italic", size = 8))) +
  theme_bw() +
  scale_y_continuous(limits = c(0,20), expand = c(0,0))+
  scale_x_discrete(name = NULL)
stack.MAG.tax 

#Traits
#IAA
gene.rel.IAA <- readRDS(here("data","intermediate","gene.rel.IAA.RDS"))
mag.link.IAA <- mag.link[mag.link$gene %in% colnames(gene.rel.IAA),] #91 obs of 3 variables
#mag.link.IAA$MAG <- gsub("-",".", mag.link.IAA$MAG)
mag.link.IAA.count <- mag.link.IAA %>%
  mutate(MAG = factor(MAG, levels = (mag.list)))%>%
  group_by(MAG, .drop=F) %>%
  summarise(count = n())

mag.link.IAA.count.status <- cbind(mag.link.IAA.count, df.MAG.soil[df.MAG.soil$MAG %in% mag.link.IAA.count$MAG,2])
colnames(mag.link.IAA.count.status) <- c("MAG","Count","Status")
sum.count.mag.IAA <- mag.link.IAA.count.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = sum(Count!=0))
sum.count.mag.IAA$nbMAGTotal <-c(16,8,10)
sum.count.mag.IAA$nbGenePerPositiveMAG <-sum.count.mag.IAA$nbGenes/sum.count.mag.IAA$nbMAGwithGene
sum.count.mag.IAA$PercentPositiveMAG <-sum.count.mag.IAA$nbMAGwithGene/sum.count.mag.IAA$nbMAGTotal*100
sum.count.mag.IAA$Trait <- rep("IAA",3)

#ACC
gene.rel.ACC <- readRDS(here("data","intermediate","gene.rel.ACC.RDS"))
mag.link.ACC <- mag.link[mag.link$gene %in% colnames(gene.rel.ACC),] #7 obs of 3 variables
#mag.link.ACC$MAG <- gsub("-",".", mag.link.ACC$MAG)
mag.link.ACC.count <- mag.link.ACC %>%
  mutate(MAG = factor(MAG, levels = (mag.list)))%>%
  group_by(MAG, .drop=F) %>%
  summarise(count = n())

mag.link.ACC.count.status <- cbind(mag.link.ACC.count, df.MAG.soil[df.MAG.soil$MAG %in% mag.link.ACC.count$MAG,2])
colnames(mag.link.ACC.count.status) <- c("MAG","Count","Status")
sum.count.mag.ACC <- mag.link.ACC.count.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = sum(Count!=0))
sum.count.mag.ACC$nbMAGTotal <-c(16,8,10)
sum.count.mag.ACC$nbGenePerPositiveMAG <-sum.count.mag.ACC$nbGenes/sum.count.mag.ACC$nbMAGwithGene
sum.count.mag.ACC$PercentPositiveMAG <-sum.count.mag.ACC$nbMAGwithGene/sum.count.mag.ACC$nbMAGTotal*100
sum.count.mag.ACC$Trait <- rep("ACC",3)

#Osmolytes
gene.rel.osmo <- readRDS(here("data","intermediate","gene.rel.osmo.RDS"))
mag.link.osmo <- mag.link[mag.link$gene %in% colnames(gene.rel.osmo),] #2882 obs of 3 var.
#mag.link.osmo$MAG <- gsub("-",".", mag.link.osmo$MAG)
mag.link.osmo.count <- mag.link.osmo %>%
  mutate(MAG = factor(MAG, levels = (mag.list)))%>%
  group_by(MAG, .drop=F) %>%
  summarise(count = n())

mag.link.osmo.count.status <- cbind(mag.link.osmo.count, df.MAG.soil[df.MAG.soil$MAG %in% mag.link.osmo.count$MAG,2])
colnames(mag.link.osmo.count.status) <- c("MAG","Count","Status")
sum.count.mag.osmo <- mag.link.osmo.count.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = sum(Count!=0))
sum.count.mag.osmo$nbMAGTotal <-c(16,8,10)
sum.count.mag.osmo$nbGenePerPositiveMAG <-sum.count.mag.osmo$nbGenes/sum.count.mag.osmo$nbMAGwithGene
sum.count.mag.osmo$PercentPositiveMAG <-sum.count.mag.osmo$nbMAGwithGene/sum.count.mag.osmo$nbMAGTotal*100
sum.count.mag.osmo$Trait <- rep("osmo",3)

#EPS
gene.rel.EPS <- readRDS(here("data","intermediate","gene.rel.EPS.RDS"))
mag.link.EPS <- mag.link[mag.link$gene %in% colnames(gene.rel.EPS),] #57 obs of 3 var
#mag.link.EPS$MAG <- gsub("-",".", mag.link.EPS$MAG)
mag.link.EPS.count <- mag.link.EPS %>%
  mutate(MAG = factor(MAG, levels = (mag.list)))%>%
  group_by(MAG, .drop=F) %>%
  summarise(count = n())

mag.link.EPS.count.status <- cbind(mag.link.EPS.count, df.MAG.soil[df.MAG.soil$MAG %in% mag.link.EPS.count$MAG,2])
colnames(mag.link.EPS.count.status) <- c("MAG","Count","Status")
sum.count.mag.EPS <- mag.link.EPS.count.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = sum(Count!=0))
sum.count.mag.EPS$nbMAGTotal <-c(16,8,10)
sum.count.mag.EPS$nbGenePerPositiveMAG <-sum.count.mag.EPS$nbGenes/sum.count.mag.EPS$nbMAGwithGene
sum.count.mag.EPS$PercentPositiveMAG <-sum.count.mag.EPS$nbMAGwithGene/sum.count.mag.EPS$nbMAGTotal*100
sum.count.mag.EPS$Trait <- rep("EPS",3)

#cytokinines
gene.rel.cyto <- readRDS(here("data","intermediate","gene.rel.cyto.RDS"))
mag.link.cyto <- mag.link[mag.link$gene %in% colnames(gene.rel.cyto),] #29 obs of 3 var
#mag.link.cyto$MAG <- gsub("-",".", mag.link.cyto$MAG)
mag.link.cyto.count <- mag.link.cyto %>%
  mutate(MAG = factor(MAG, levels = (mag.list)))%>%
  group_by(MAG, .drop=F) %>%
  summarise(count = n())

mag.link.cyto.count.status <- cbind(mag.link.cyto.count, df.MAG.soil[df.MAG.soil$MAG %in% mag.link.cyto.count$MAG,2])
colnames(mag.link.cyto.count.status) <- c("MAG","Count","Status")
sum.count.mag.cyto <- mag.link.cyto.count.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = sum(Count!=0))
sum.count.mag.cyto$nbMAGTotal <-c(16,8,10)
sum.count.mag.cyto$nbGenePerPositiveMAG <-sum.count.mag.cyto$nbGenes/sum.count.mag.cyto$nbMAGwithGene
sum.count.mag.cyto$PercentPositiveMAG <-sum.count.mag.cyto$nbMAGwithGene/sum.count.mag.cyto$nbMAGTotal*100
sum.count.mag.cyto$Trait <- rep("cyto",3)

#antiox
gene.rel.antiox <- readRDS(here("data","intermediate","gene.rel.antiox.RDS"))
mag.link.antiox <- mag.link[mag.link$gene %in% colnames(gene.rel.antiox),] #349 obs of 3 var
#mag.link.antiox$MAG <- gsub("-",".", mag.link.antiox$MAG)
mag.link.antiox.count <- mag.link.antiox %>%
  mutate(MAG = factor(MAG, levels = (mag.list)))%>%
  group_by(MAG, .drop=F) %>%
  summarise(count = n())

mag.link.antiox.count.status <- cbind(mag.link.antiox.count, df.MAG.soil[df.MAG.soil$MAG %in% mag.link.antiox.count$MAG,2])
colnames(mag.link.antiox.count.status) <- c("MAG","Count","Status")
sum.count.mag.antiox <- mag.link.antiox.count.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = sum(Count!=0))
sum.count.mag.antiox$nbMAGTotal <-c(16,8,10)
sum.count.mag.antiox$nbGenePerPositiveMAG <-sum.count.mag.antiox$nbGenes/sum.count.mag.antiox$nbMAGwithGene
sum.count.mag.antiox$PercentPositiveMAG <-sum.count.mag.antiox$nbMAGwithGene/sum.count.mag.antiox$nbMAGTotal*100
sum.count.mag.antiox$Trait <- rep("antiox",3)

#Combine and plot
sum.count.mag.soil.all <- rbind(sum.count.mag.IAA, sum.count.mag.ACC, sum.count.mag.osmo, 
                                sum.count.mag.cyto, sum.count.mag.EPS, sum.count.mag.antiox)
#Bar charts
sum.count.mag.soil.all$Status <- factor(sum.count.mag.soil.all$Status, c("continuous", "intermittent", "unaffected"))
bar.mag.perc.soil <- ggplot(sum.count.mag.soil.all, aes(fill = Status, x = Trait, y = PercentPositiveMAG)) + 
  geom_bar(stat = "summary", fun ="mean", position = "dodge" ) +
  ylab("Percent positive MAGs") + 
  scale_fill_manual(values = c("darkred", "darkblue", "darkgreen"), guide = guide_legend(title = "Soil"), labels = c("continuous WSH", "intermittent WSH", "unaffected")) +
  theme_bw() +
  scale_y_continuous(limits = c(0,105), expand = c(0,0))+
  scale_x_discrete(labels = c("ACC", "Antioxidants", "Cytokinins", "EPS", "IAA", "Osmolytes"))
bar.mag.perc.soil

bar.gene.mag.soil <- ggplot(sum.count.mag.soil.all, aes(fill = Status, x = Trait, y = nbGenePerPositiveMAG)) + 
  geom_bar(stat = "summary", fun ="mean", position = "dodge" ) +
  ylab("Nb genes per positive MAGs") + 
  scale_fill_manual(values = c("darkred", "darkblue", "darkgreen"), guide = guide_legend(title = "Soil"), labels = c("continuous WSH", "intermittent WSH", "unaffected")) +
  theme_bw() +
  scale_y_continuous(limits = c(0,105), expand = c(0,0) )+
  scale_x_discrete(labels = c("ACC", "Antioxidants", "Cytokinins", "EPS", "IAA", "Osmolytes"))
bar.gene.mag.soil

#How many functional traits per MAG?
df.MAG.soil.overlap <- df.MAG.soil
df.MAG.soil.overlap$IAA <- df.MAG.soil$MAG %in% mag.link.IAA$MAG 
df.MAG.soil.overlap$ACC <- df.MAG.soil$MAG %in% mag.link.ACC$MAG
df.MAG.soil.overlap$EPS <- df.MAG.soil$MAG %in% mag.link.EPS$MAG
df.MAG.soil.overlap$antiox <- df.MAG.soil$MAG %in% mag.link.antiox$MAG
df.MAG.soil.overlap$osmo <- df.MAG.soil$MAG %in% mag.link.osmo$MAG
df.MAG.soil.overlap$cyto <- df.MAG.soil$MAG %in% mag.link.cyto$MAG
df.MAG.soil.overlap$Sum <- rowSums(df.MAG.soil.overlap[,3:8])
MAG.soil.overlap.count <- df.MAG.soil.overlap %>%  
  mutate(Sum = factor(Sum, levels = c(0,1,2,3,4,5,6)), Status = factor(Status, levels = c("continuous", "intermittent", "unaffected"))) %>%
  group_by(Sum,Status, .drop = F) %>%
  summarise(count = n())
MAG.soil.overlap.count$Percent <- MAG.soil.overlap.count$count/rep(c(16,8,10),7)*100
sum(MAG.soil.overlap.count$Percent)
MAG.soil.overlap.count$Sum <- as.character(MAG.soil.overlap.count$Sum)
MAG.soil.overlap.count$Status <- factor(MAG.soil.overlap.count$Status, c("continuous", "intermittent", "unaffected"))

#Plot
bar.soil.overlap <- ggplot(MAG.soil.overlap.count, aes(fill = Status, x = Sum, y = Percent)) + 
  geom_bar(stat = "summary", fun ="mean", position = "dodge" ) +
  ylab("% of MAGs") + 
  scale_fill_manual(values = c("darkred", "darkblue", "darkgreen"), guide = guide_legend(title = "Soil"), labels = c("continuous WSH", "intermittent WSH", "unaffected")) +
  theme_bw() +
  scale_y_continuous(limits = c(0,70), expand = c(0,0))+
  scale_x_discrete(labels = c("0", "1", "2", "3", "4", "5", "6"), name = "Number of functional traits")
bar.soil.overlap

#Supplementary table of MAG characteristics
#Check ordering
row.names(mag.qual[(mag.qual$Completeness>50 & mag.qual$Contamination<10 & mag.qual$Order!="Thermales"),]) == mag.link.ACC.count.status$MAG
row.names(mag.qual[(mag.qual$Completeness>50 & mag.qual$Contamination<10 & mag.qual$Order!="Thermales"),]) == mag.link.IAA.count.status$MAG
row.names(mag.qual[(mag.qual$Completeness>50 & mag.qual$Contamination<10 & mag.qual$Order!="Thermales"),]) == mag.link.cyto.count.status$MAG
row.names(mag.qual[(mag.qual$Completeness>50 & mag.qual$Contamination<10 & mag.qual$Order!="Thermales"),]) == mag.link.antiox.count.status$MAG
row.names(mag.qual[(mag.qual$Completeness>50 & mag.qual$Contamination<10 & mag.qual$Order!="Thermales"),]) == mag.link.osmo.count.status$MAG
row.names(mag.qual[(mag.qual$Completeness>50 & mag.qual$Contamination<10 & mag.qual$Order!="Thermales"),]) == mag.link.EPS.count.status$MAG
#Create data frame
df.suppl <- data.frame (Status=mag.link.ACC.count.status$Status,
                        ACC=mag.link.ACC.count.status$Count,
                        IAA=mag.link.IAA.count.status$Count,
                        Osmolytes=mag.link.osmo.count.status$Count,
                        EPS=mag.link.EPS.count.status$Count,
                        Cytokinines=mag.link.cyto.count.status$Count,
                        Antioxidants=mag.link.antiox.count.status$Count,
                        mag.qual[(mag.qual$Completeness>50 & mag.qual$Contamination<10 & mag.qual$Order!="Thermales"),c(12,13,15,16,17,19:25)]
                        )
#export
write.table(df.suppl, file = here("output", "tables","SupplementaryTable1.txt"), sep = "\t", eol = "\n")
