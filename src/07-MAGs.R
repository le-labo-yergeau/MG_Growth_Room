#Find MAGs that are influenced by the treatments
#Then look for traits in the genomes

#Load data
mag.table.rel <- readRDS(file = here("data", "intermediate", "mag.table.rel.RDS")) #20 obs of 714 var
mag.link <- readRDS(file = here("data", "intermediate", "mag.link.RDS")) #1833497 obs of 3 var
mag.tax <- readRDS(file = here("data", "intermediate", "mag.tax.RDS")) #1833497 obs of 3 var
map <- readRDS(file = here("data", "intermediate", "map.RDS")) #20 obs of 3 variables

#Check sorting
row.names(map) == row.names(mag.table.rel)

#Anova loop for the 714 MAGs
i <- 1
list.MAG.SWHC <- vector()
list.MAG.inter <- vector()
list.MAG.soil <- vector()

while (i<=ncol(mag.table.rel)) {
  anova <- aov(mag.table.rel[,i]~map$SoilType*map$perc_SWHC+map$block)
  if(summary(anova)[[1]][1,5]<(0.05/714)){
    if(mean(mag.table.rel[map$SoilType=="IR",i])< mean(mag.table.rel[map$SoilType=="NI",i])){
      list.MAG.soil <- append(list.MAG.soil, c(colnames(mag.table.rel)[i], "SWSH"))
    }else{list.MAG.soil <- append(list.MAG.soil, c(colnames(mag.table.rel)[i], "no SWSH"))}  
  }else{list.MAG.soil <-append(list.MAG.soil, c(colnames(mag.table.rel)[i], "unaffected"))}
  if(summary(anova)[[1]][2,5]<(0.05/714)){
    if(mean(mag.table.rel[map$perc_SWHC=="5",i]) < mean(mag.table.rel[map$perc_SWHC=="50",i])){
      list.MAG.SWHC <- append(list.MAG.SWHC, c(colnames(mag.table.rel)[i], "50"))
    }else{list.MAG.SWHC <- append(list.MAG.SWHC, c(colnames(mag.table.rel)[i], "5"))}  
  }else{list.MAG.SWHC <-append(list.MAG.SWHC, c(colnames(mag.table.rel)[i], "unaffected"))}
  if(summary(anova)[[1]][4,5]<(0.05/714)){
    list.MAG.inter <- append(list.MAG.inter, colnames(mag.table.rel)[i])#Only 2
  }
  i <- i+1
}
df.MAG.soil <- data.frame(matrix(list.MAG.soil, nrow = 714, ncol = 2, byrow = TRUE))
colnames(df.MAG.soil) <- c("MAG", "Status")
df.MAG.SWHC <- data.frame(matrix(list.MAG.SWHC, nrow = 714, ncol = 2, byrow = TRUE))
colnames(df.MAG.SWHC) <- c("MAG", "Status")

##MAGs affected by soil history
#Taxonomy
sum(row.names(mag.tax) == df.MAG.soil$MAG)#714
df.mag.soil.tax <- cbind(df.MAG.soil,mag.tax)
#Make a summary
df.mag.summary <- df.mag.soil.tax %>%
  group_by(Phylum, Status) %>%
  summarise(count = n())   %>%
  filter(count > 10)
df.mag.summary$Phylum <- gsub("p__","", df.mag.summary$Phylum)
df.mag.summary$Phylum <- gsub("NULL","unassigned", df.mag.summary$Phylum)
others <- data.frame(matrix(c(rep("Others",3), 
                              "SWSH", "no SWSH", "unaffected", 
                              264-sum(df.mag.summary[df.mag.summary$Status=="SWSH",3]),
                              171-sum(df.mag.summary[df.mag.summary$Status=="no SWSH",3]),
                              714-435-sum(df.mag.summary[df.mag.summary$Status=="unaffected",3])),
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
  scale_y_continuous(limits = c(0,300), expand = c(0,0))+
  scale_x_discrete(name = NULL)
stack.MAG.tax 

#Traits
#IAA
gene.rel.IAA <- readRDS(here("data","intermediate","gene.rel.IAA.RDS"))
mag.link.IAA <- mag.link[mag.link$gene %in% colnames(gene.rel.IAA),]
mag.link.IAA$MAG <- gsub("-",".", mag.link.IAA$MAG)
mag.link.IAA.count <- mag.link.IAA %>%
  group_by(MAG) %>%
  summarise(count = n())

mag.link.IAA.count.status <- cbind(mag.link.IAA.count, df.MAG.soil[df.MAG.soil$MAG %in% mag.link.IAA.count$MAG,2])
colnames(mag.link.IAA.count.status) <- c("MAG","Count","Status")
sum.count.mag.IAA <- mag.link.IAA.count.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = n())
sum.count.mag.IAA$nbMAGTotal <-c(264,171,279)
sum.count.mag.IAA$nbGenePerPositiveMAG <-sum.count.mag.IAA$nbGenes/sum.count.mag.IAA$nbMAGwithGene
sum.count.mag.IAA$PercentPositiveMAG <-sum.count.mag.IAA$nbMAGwithGene/sum.count.mag.IAA$nbMAGTotal*100
sum.count.mag.IAA$Trait <- rep("IAA",3)

#ACC
gene.rel.ACC <- readRDS(here("data","intermediate","gene.rel.ACC.RDS"))
mag.link.ACC <- mag.link[mag.link$gene %in% colnames(gene.rel.ACC),]
mag.link.ACC$MAG <- gsub("-",".", mag.link.ACC$MAG)
mag.link.ACC.count <- mag.link.ACC %>%
  group_by(MAG) %>%
  summarise(count = n())

mag.link.ACC.count.status <- cbind(mag.link.ACC.count, df.MAG.soil[df.MAG.soil$MAG %in% mag.link.ACC.count$MAG,2])
colnames(mag.link.ACC.count.status) <- c("MAG","Count","Status")
sum.count.mag.ACC <- mag.link.ACC.count.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = n())
sum.count.mag.ACC$nbMAGTotal <-c(264,171,279)
sum.count.mag.ACC$nbGenePerPositiveMAG <-sum.count.mag.ACC$nbGenes/sum.count.mag.ACC$nbMAGwithGene
sum.count.mag.ACC$PercentPositiveMAG <-sum.count.mag.ACC$nbMAGwithGene/sum.count.mag.ACC$nbMAGTotal*100
sum.count.mag.ACC$Trait <- rep("ACC",3)

#Osmolytes
gene.rel.osmo <- readRDS(here("data","intermediate","gene.rel.osmo.RDS"))
mag.link.osmo <- mag.link[mag.link$gene %in% colnames(gene.rel.osmo),]
mag.link.osmo$MAG <- gsub("-",".", mag.link.osmo$MAG)
mag.link.osmo.count <- mag.link.osmo %>%
  group_by(MAG) %>%
  summarise(count = n())

mag.link.osmo.count.status <- cbind(mag.link.osmo.count, df.MAG.soil[df.MAG.soil$MAG %in% mag.link.osmo.count$MAG,2])
colnames(mag.link.osmo.count.status) <- c("MAG","Count","Status")
sum.count.mag.osmo <- mag.link.osmo.count.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = n())
sum.count.mag.osmo$nbMAGTotal <-c(264,171,279)
sum.count.mag.osmo$nbGenePerPositiveMAG <-sum.count.mag.osmo$nbGenes/sum.count.mag.osmo$nbMAGwithGene
sum.count.mag.osmo$PercentPositiveMAG <-sum.count.mag.osmo$nbMAGwithGene/sum.count.mag.osmo$nbMAGTotal*100
sum.count.mag.osmo$Trait <- rep("osmo",3)

#EPS
gene.rel.EPS <- readRDS(here("data","intermediate","gene.rel.EPS.RDS"))
mag.link.EPS <- mag.link[mag.link$gene %in% colnames(gene.rel.EPS),]
mag.link.EPS$MAG <- gsub("-",".", mag.link.EPS$MAG)
mag.link.EPS.count <- mag.link.EPS %>%
  group_by(MAG) %>%
  summarise(count = n())

mag.link.EPS.count.status <- cbind(mag.link.EPS.count, df.MAG.soil[df.MAG.soil$MAG %in% mag.link.EPS.count$MAG,2])
colnames(mag.link.EPS.count.status) <- c("MAG","Count","Status")
sum.count.mag.EPS <- mag.link.EPS.count.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = n())
sum.count.mag.EPS$nbMAGTotal <-c(264,171,279)
sum.count.mag.EPS$nbGenePerPositiveMAG <-sum.count.mag.EPS$nbGenes/sum.count.mag.EPS$nbMAGwithGene
sum.count.mag.EPS$PercentPositiveMAG <-sum.count.mag.EPS$nbMAGwithGene/sum.count.mag.EPS$nbMAGTotal*100
sum.count.mag.EPS$Trait <- rep("EPS",3)

#cytokinines
gene.rel.cyto <- readRDS(here("data","intermediate","gene.rel.cyto.RDS"))
mag.link.cyto <- mag.link[mag.link$gene %in% colnames(gene.rel.cyto),]
mag.link.cyto$MAG <- gsub("-",".", mag.link.cyto$MAG)
mag.link.cyto.count <- mag.link.cyto %>%
  group_by(MAG) %>%
  summarise(count = n())

mag.link.cyto.count.status <- cbind(mag.link.cyto.count, df.MAG.soil[df.MAG.soil$MAG %in% mag.link.cyto.count$MAG,2])
colnames(mag.link.cyto.count.status) <- c("MAG","Count","Status")
sum.count.mag.cyto <- mag.link.cyto.count.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = n())
sum.count.mag.cyto$nbMAGTotal <-c(264,171,279)
sum.count.mag.cyto$nbGenePerPositiveMAG <-sum.count.mag.cyto$nbGenes/sum.count.mag.cyto$nbMAGwithGene
sum.count.mag.cyto$PercentPositiveMAG <-sum.count.mag.cyto$nbMAGwithGene/sum.count.mag.cyto$nbMAGTotal*100
sum.count.mag.cyto$Trait <- rep("cyto",3)

#antiox
gene.rel.antiox <- readRDS(here("data","intermediate","gene.rel.antiox.RDS"))
mag.link.antiox <- mag.link[mag.link$gene %in% colnames(gene.rel.antiox),]
mag.link.antiox$MAG <- gsub("-",".", mag.link.antiox$MAG)
mag.link.antiox.count <- mag.link.antiox %>%
  group_by(MAG) %>%
  summarise(count = n())

mag.link.antiox.count.status <- cbind(mag.link.antiox.count, df.MAG.soil[df.MAG.soil$MAG %in% mag.link.antiox.count$MAG,2])
colnames(mag.link.antiox.count.status) <- c("MAG","Count","Status")
sum.count.mag.antiox <- mag.link.antiox.count.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = n())
sum.count.mag.antiox$nbMAGTotal <-c(264,171,279)
sum.count.mag.antiox$nbGenePerPositiveMAG <-sum.count.mag.antiox$nbGenes/sum.count.mag.antiox$nbMAGwithGene
sum.count.mag.antiox$PercentPositiveMAG <-sum.count.mag.antiox$nbMAGwithGene/sum.count.mag.antiox$nbMAGTotal*100
sum.count.mag.antiox$Trait <- rep("antiox",3)

#Combine and plot
sum.count.mag.soil.all <- rbind(sum.count.mag.IAA, sum.count.mag.ACC, sum.count.mag.osmo, 
                                sum.count.mag.cyto, sum.count.mag.EPS, sum.count.mag.antiox)
#Bar charts
sum.count.mag.soil.all$Status <- factor(sum.count.mag.soil.all$Status, c("SWSH", "no SWSH", "unaffected"))
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
  group_by(Sum,Status) %>%
  summarise(count = n())
MAG.soil.overlap.count$Percent <- MAG.soil.overlap.count$count/rep(c(264,171,279),7)*100
sum(MAG.soil.overlap.count$Percent)
MAG.soil.overlap.count$Sum <- as.character(MAG.soil.overlap.count$Sum)
MAG.soil.overlap.count$Status <- factor(MAG.soil.overlap.count$Status, c("SWSH", "no SWSH", "unaffected"))

#Plot
bar.soil.overlap <- ggplot(MAG.soil.overlap.count, aes(fill = Status, x = Sum, y = Percent)) + 
  geom_bar(stat = "summary", fun ="mean", position = "dodge" ) +
  ylab("% of MAGs") + 
  scale_fill_manual(values = c("darkred", "darkblue", "darkgreen"), guide = guide_legend(title = "Soil"), labels = c("continuous WSH", "intermittent WSH", "unaffected")) +
  theme_bw() +
  scale_y_continuous(limits = c(0,27), expand = c(0,0))+
  scale_x_discrete(labels = c("0", "1", "2", "3", "4", "5", "6"), name = "Number of functional traits")
bar.soil.overlap

#How many have both osmolytes and IAA?
sum(df.MAG.soil.overlap[df.MAG.soil.overlap$Status=="no SWSH",]$IAA & df.MAG.soil.overlap[df.MAG.soil.overlap$Status=="no SWSH",]$osmo) #84
sum(df.MAG.soil.overlap[df.MAG.soil.overlap$Status=="SWSH",]$IAA & df.MAG.soil.overlap[df.MAG.soil.overlap$Status=="SWSH",]$osmo) #107
sum(df.MAG.soil.overlap[df.MAG.soil.overlap$Status=="unaffected",]$IAA & df.MAG.soil.overlap[df.MAG.soil.overlap$Status=="unaffected",]$osmo) #97

##MAGs affected by soil water content (actual)
#Taxonomy
sum(row.names(mag.tax) == df.MAG.SWHC$MAG)#714
df.mag.SWHC.tax <- cbind(df.MAG.SWHC,mag.tax)
#Make a summary
df.mag.SWHC.summary <- df.mag.SWHC.tax %>%
  group_by(Phylum, Status) %>%
  summarise(count = n())   %>%
  filter(Phylum %in% c("p__Acidobacteria", "p__Actinobacteria", "p__NULL", "p__Proteobacteria", "p__Verrucomicrobia"))
df.mag.SWHC.summary$Phylum <- gsub("p__","", df.mag.SWHC.summary$Phylum)
df.mag.SWHC.summary$Phylum <- gsub("NULL","unassigned", df.mag.SWHC.summary$Phylum)
others <- data.frame(matrix(c(rep("Others",3), 
                              "5", "50", "unaffected", 
                              3-sum(df.mag.SWHC.summary[df.mag.SWHC.summary$Status=="5",3]),
                              42-sum(df.mag.SWHC.summary[df.mag.SWHC.summary$Status=="50",3]),
                              669-sum(df.mag.SWHC.summary[df.mag.SWHC.summary$Status=="unaffected",3])),
                            nrow = 3, ncol = 3))
colnames(others) <- colnames(df.mag.SWHC.summary)
others$count <- as.numeric(others$count)
df.mag.SWHC.summary.others <- rbind(df.mag.SWHC.summary,others)

#Plot
palette(c(brewer.pal(n = 9, name = "Set1"),"lightgrey", "black", "darkred", "darkblue", "darkgreen", "purple4", "darkgrey", "white"))
stack.MAG.SWHC.tax <- ggplot(df.mag.SWHC.summary.others, aes(fill = Phylum, y = count, x = Status)) + 
  geom_bar( stat = "summary", fun ="mean", position = "stack") +
  ylab("Number of MAGs") + 
  scale_fill_manual(values = palette(), guide = guide_legend(title = "Phylum", label.theme = element_text(face = "italic", size = 8))) +
  theme_bw() +
  scale_y_continuous(limits = c(0,700), expand = c(0,0))+
  scale_x_discrete(name = NULL)
stack.MAG.SWHC.tax 

#Traits
#IAA
mag.link.IAA.count.SWHC.status <- cbind(mag.link.IAA.count, df.MAG.SWHC[df.MAG.SWHC$MAG %in% mag.link.IAA.count$MAG,2])
colnames(mag.link.IAA.count.SWHC.status) <- c("MAG","Count","Status")
sum.count.mag.SWHC.IAA <- mag.link.IAA.count.SWHC.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = n())
sum.count.mag.SWHC.IAA$nbMAGTotal <-c(3,42,669)
sum.count.mag.SWHC.IAA$nbGenePerPositiveMAG <-sum.count.mag.SWHC.IAA$nbGenes/sum.count.mag.SWHC.IAA$nbMAGwithGene
sum.count.mag.SWHC.IAA$PercentPositiveMAG <-sum.count.mag.SWHC.IAA$nbMAGwithGene/sum.count.mag.SWHC.IAA$nbMAGTotal*100
sum.count.mag.SWHC.IAA$Trait <- rep("IAA",3)

#ACC
mag.link.ACC.count.SWHC.status <- cbind(mag.link.ACC.count, df.MAG.SWHC[df.MAG.SWHC$MAG %in% mag.link.ACC.count$MAG,2])
colnames(mag.link.ACC.count.SWHC.status) <- c("MAG","Count","Status")
sum.count.mag.SWHC.ACC <- mag.link.ACC.count.SWHC.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = n())
sum.count.mag.SWHC.ACC[3,] <- list("50",0,0)
sum.count.mag.SWHC.ACC$nbMAGTotal <-c(42,669,3)
sum.count.mag.SWHC.ACC$nbGenePerPositiveMAG <-sum.count.mag.SWHC.ACC$nbGenes/sum.count.mag.SWHC.ACC$nbMAGwithGene
sum.count.mag.SWHC.ACC$PercentPositiveMAG <-sum.count.mag.SWHC.ACC$nbMAGwithGene/sum.count.mag.SWHC.ACC$nbMAGTotal*100
sum.count.mag.SWHC.ACC$Trait <- rep("ACC",3)
sum.count.mag.SWHC.ACC[3,5] <- 0

#Osmolytes
mag.link.osmo.count.SWHC.status <- cbind(mag.link.osmo.count, df.MAG.SWHC[df.MAG.SWHC$MAG %in% mag.link.osmo.count$MAG,2])
colnames(mag.link.osmo.count.SWHC.status) <- c("MAG","Count","Status")
sum.count.mag.SWHC.osmo <- mag.link.osmo.count.SWHC.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = n())
sum.count.mag.SWHC.osmo$nbMAGTotal <-c(3,42,669)
sum.count.mag.SWHC.osmo$nbGenePerPositiveMAG <-sum.count.mag.SWHC.osmo$nbGenes/sum.count.mag.SWHC.osmo$nbMAGwithGene
sum.count.mag.SWHC.osmo$PercentPositiveMAG <-sum.count.mag.SWHC.osmo$nbMAGwithGene/sum.count.mag.SWHC.osmo$nbMAGTotal*100
sum.count.mag.SWHC.osmo$Trait <- rep("osmo",3)

#EPS
mag.link.EPS.count.SWHC.status <- cbind(mag.link.EPS.count, df.MAG.SWHC[df.MAG.SWHC$MAG %in% mag.link.EPS.count$MAG,2])
colnames(mag.link.EPS.count.SWHC.status) <- c("MAG","Count","Status")
sum.count.mag.SWHC.EPS <- mag.link.EPS.count.SWHC.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = n())
sum.count.mag.SWHC.EPS$nbMAGTotal <-c(3,42,669)
sum.count.mag.SWHC.EPS$nbGenePerPositiveMAG <-sum.count.mag.SWHC.EPS$nbGenes/sum.count.mag.SWHC.EPS$nbMAGwithGene
sum.count.mag.SWHC.EPS$PercentPositiveMAG <-sum.count.mag.SWHC.EPS$nbMAGwithGene/sum.count.mag.SWHC.EPS$nbMAGTotal*100
sum.count.mag.SWHC.EPS$Trait <- rep("EPS",3)

#cytokinines
mag.link.cyto.count.SWHC.status <- cbind(mag.link.cyto.count, df.MAG.SWHC[df.MAG.SWHC$MAG %in% mag.link.cyto.count$MAG,2])
colnames(mag.link.cyto.count.SWHC.status) <- c("MAG","Count","Status")
sum.count.mag.SWHC.cyto <- mag.link.cyto.count.SWHC.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = n())
sum.count.mag.SWHC.cyto$nbMAGTotal <-c(3,42,669)
sum.count.mag.SWHC.cyto$nbGenePerPositiveMAG <-sum.count.mag.SWHC.cyto$nbGenes/sum.count.mag.SWHC.cyto$nbMAGwithGene
sum.count.mag.SWHC.cyto$PercentPositiveMAG <-sum.count.mag.SWHC.cyto$nbMAGwithGene/sum.count.mag.SWHC.cyto$nbMAGTotal*100
sum.count.mag.SWHC.cyto$Trait <- rep("cyto",3)

#antiox
mag.link.antiox.count.SWHC.status <- cbind(mag.link.antiox.count, df.MAG.SWHC[df.MAG.SWHC$MAG %in% mag.link.antiox.count$MAG,2])
colnames(mag.link.antiox.count.SWHC.status) <- c("MAG","Count","Status")
sum.count.mag.SWHC.antiox <- mag.link.antiox.count.SWHC.status %>%
  group_by(Status) %>%
  summarise(nbGenes = sum(Count), nbMAGwithGene = n())
sum.count.mag.SWHC.antiox$nbMAGTotal <-c(3,42,669)
sum.count.mag.SWHC.antiox$nbGenePerPositiveMAG <-sum.count.mag.SWHC.antiox$nbGenes/sum.count.mag.SWHC.antiox$nbMAGwithGene
sum.count.mag.SWHC.antiox$PercentPositiveMAG <-sum.count.mag.SWHC.antiox$nbMAGwithGene/sum.count.mag.SWHC.antiox$nbMAGTotal*100
sum.count.mag.SWHC.antiox$Trait <- rep("antiox",3)

#Combine and plot
sum.count.mag.SWHC.all <- rbind(sum.count.mag.SWHC.IAA, sum.count.mag.SWHC.ACC, sum.count.mag.SWHC.osmo, 
                                sum.count.mag.SWHC.cyto, sum.count.mag.SWHC.EPS, sum.count.mag.SWHC.antiox)
#Bar charts
bar.mag.perc.SWHC <- ggplot(sum.count.mag.SWHC.all, aes(fill = Status, x = Trait, y = PercentPositiveMAG)) + 
  geom_bar(stat = "summary", fun ="mean", position = "dodge" ) +
  ylab("Percent positive MAGs") + 
  scale_fill_manual(values = c("red", "blue", "darkgreen"), guide = guide_legend(title = "%SWHC")) +
  theme_bw() +
  scale_y_continuous(limits = c(0,105), expand = c(0,0))+
  scale_x_discrete(labels = c("ACC", "Antioxidants", "Cytokinins", "EPS", "IAA", "Osmolytes"))
bar.mag.perc.SWHC

bar.gene.mag.SWHC <- ggplot(sum.count.mag.SWHC.all, aes(fill = Status, x = Trait, y = nbGenePerPositiveMAG)) + 
  geom_bar(stat = "summary", fun ="mean", position = "dodge" ) +
  ylab("Nb genes per positive MAGs") + 
  scale_fill_manual(values = c("red", "blue", "darkgreen"), guide = guide_legend(title = "%SWHC")) +
  theme_bw() +
  scale_y_continuous(limits = c(0,125), expand = c(0,0))+
  scale_x_discrete(labels = c("ACC", "Antioxidants", "Cytokinins", "EPS", "IAA", "Osmolytes"))
bar.gene.mag.SWHC

#How many functional traits per MAG?
df.MAG.SWHC.overlap <- df.MAG.SWHC
df.MAG.SWHC.overlap$IAA <- df.MAG.SWHC$MAG %in% mag.link.IAA$MAG 
df.MAG.SWHC.overlap$ACC <- df.MAG.SWHC$MAG %in% mag.link.ACC$MAG
df.MAG.SWHC.overlap$EPS <- df.MAG.SWHC$MAG %in% mag.link.EPS$MAG
df.MAG.SWHC.overlap$antiox <- df.MAG.SWHC$MAG %in% mag.link.antiox$MAG
df.MAG.SWHC.overlap$osmo <- df.MAG.SWHC$MAG %in% mag.link.osmo$MAG
df.MAG.SWHC.overlap$cyto <- df.MAG.SWHC$MAG %in% mag.link.cyto$MAG
df.MAG.SWHC.overlap$Sum <- factor(rowSums(df.MAG.SWHC.overlap[,3:8]),levels = c("0", "1", "2", "3", "4", "5", "6"))
df.MAG.SWHC.overlap$Status <- factor(df.MAG.SWHC.overlap$Status,levels = c("5", "50", "unaffected"))

MAG.SWHC.overlap.count <- df.MAG.SWHC.overlap %>%
  group_by(Sum,Status, .drop = FALSE) %>%
  summarise(count = n())
MAG.SWHC.overlap.count$Percent <- MAG.SWHC.overlap.count$count/rep(c(3, 42, 669),7)*100
sum(MAG.SWHC.overlap.count$Percent)
#MAG.SWHC.overlap.count$Sum <- as.character(MAG.SWHC.overlap.count$Sum)
#MAG.SWHC.overlap.count$Status <- factor(MAG.SWHC.overlap.count$Status, c("5", "50", "unaffected"))

#Plot
bar.SWHC.overlap <- ggplot(MAG.SWHC.overlap.count, aes(fill = Status, x = Sum, y = Percent)) + 
  geom_bar(stat = "summary", fun ="mean", position = "dodge" ) +
  ylab("% of MAGs") + 
  scale_fill_manual(values = c("red", "blue", "darkgreen"), guide = guide_legend(title = "%SWHC")) +
  theme_bw() +
  scale_y_continuous(limits = c(0,35), expand = c(0,0))+
  scale_x_discrete(labels = c("0", "1", "2", "3", "4", "5", "6"), name = "Number of functional traits")
bar.SWHC.overlap
