###PCOA and Permanova

##Normalized by relative abundance
#Load data
gene.rel <- readRDS(here("data", "intermediate", "gene.rel.RDS"))
map.s <- readRDS(here("data","intermediate","map.RDS"))

#Calculate Bray-Curtis distance
bray.rel <- vegdist(gene.rel, method="bray")
saveRDS(bray.rel, here("data","intermediate", "bray.rel.RDS"))
rm(gene.rel)

#Principal coordinate analysis
pcoa.rel <- cmdscale(sqrt(bray.rel), k=19, eig=T)

#Sort
pcoa.rel.s <- pcoa.rel$points[order(row.names(pcoa.rel$points)),]
sum(row.names(pcoa.rel.s) == row.names(map.s)) # 20

#create a data frame
pcoa.rel.map <- data.frame(pcoa.rel.s, map.s)
pcoa.rel.map$perc_SWHC <- as.character(pcoa.rel.map$perc_SWHC)

#PCT explained
pcoa.rel$eig[1]/sum(pcoa.rel$eig)*100 #17.1994%
pcoa.rel$eig[2]/sum(pcoa.rel$eig)*100 #5.914303%

#Reorder manually
pcoa.rel.map$perc_SWHC <- factor(pcoa.rel.map$perc_SWHC, c("50", "5"))

#Plot
pcoa.rel.plot <- ggplot(data=pcoa.rel.map, aes(x=X1, y=X2, shape=SoilType, colour=perc_SWHC)) + 
  geom_point() +
  xlab("PCoA axis 1 = 17.2%") + 
  ylab("PCoA axis 2 = 5.9%") + 
  theme_bw()+
  scale_color_discrete(name = "% SWHC", labels = c("50%", "5%"), type = c("blue", "red"))+
  scale_shape_discrete(name = "Soil water stress history", labels = c("intermittent", "continuous"))
pcoa.rel.plot

#Permanova
#Check for betadisp
row.names(data.frame(as.matrix(bray.rel)))==row.names(map.s) #Check ordering
anova(betadisper(bray.rel, group = map.s$SoilType)) #F=0.3893, P=0.5405
anova(betadisper(bray.rel, group = map.s$perc_SWHC)) #F=0.059, P=0.8108

#Specify blocks
perm <-how(nperm = 9999)
setBlocks (perm) <- with(map.s, block)

set.seed(678)
adonis2(bray.rel~map.s$perc_SWHC*map.s$SoilType, permutations = perm)

#Df SumOfSqs      R2      F Pr(>F)    
#map.s$perc_SWHC                 1  0.09299 0.05669 1.4500 0.1602    
#map.s$SoilType                  1  0.44790 0.27307 6.9845 0.0001 ***
#  map.s$perc_SWHC:map.s$SoilType  1  0.07332 0.04470 1.1434 0.2673    
#Residual                       16  1.02604 0.62554                  
#Total                          19  1.64025 1.00000       