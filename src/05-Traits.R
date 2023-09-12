###Look for traits linked to plant resistance to water stress

#Load data
annot <- fread(here("data", "raw", "annotations.tsv"))#13,108,188 obs of 33 variables
gene.norm <- readRDS(here("data", "intermediate", "gene.norm.RDS"))
gene.rel <- readRDS(here("data", "intermediate", "gene.rel.RDS"))
map.s <- readRDS(here("data", "intermediate", "map.RDS"))

##IAA - root elongation
##Searching for KEGG entry of the last step leading to IAA
##K01501: nitrilase - spme other pathways
##K01426: amidase - some other pathways
##K21801: indoleacetamide hydrolase
##K11816: indole-3-pyruvate monooxygenase
##K11817: indole-3-acetaldehyde oxidase
##K00128: aldehyde dehydrogenase - many other pathways

#Create KEGG entry list for IAA
iaa.KO.list <- c("K01501", "K01426", "K21801", "K11816", "K11817", "K00128")
#Get lines that match our IAA list
annot.IAA <- annot[annot$kegg_entry %in% iaa.KO.list,] #11,497 obs. of 33 variables
##Normalized by relative abundance
#Get genes abundance that match the IAA annotations
gene.rel.IAA <- gene.rel[,colnames(gene.rel) %in% annot.IAA$gene_id] #20 obs of 11,497 variables
#Save intermediate files
saveRDS(gene.rel.IAA, here("data","intermediate","gene.rel.IAA.RDS"))
#Look at sum of IAA-related genes
rownames(map.s)==rownames(gene.rel.IAA)#check sorting
sum.map.rel.IAA <- data.frame(map.s, "sum.IAA" = rowSums(gene.rel.IAA), "log.sum.IAA" = log10(rowSums(gene.rel.IAA)))
#Tukey letters - from ANOVA below
tukey.IAA.rel <- data.frame(letters = c("ac", "a", "b", "c"), x = c(0.8, 1.2, 1.8, 2.2), SoilType=c("IR","IR", "NI", "NI") , perc_SWHC= c("50","5","50","5"), y=c(0.00116, 0.00116, 0.00112, 0.00114))
#Plot
sum.map.rel.IAA$perc_SWHC <- factor(sum.map.rel.IAA$perc_SWHC, c("50", "5")) #Reorder manually
IAA.rel.box <- ggplot(data = sum.map.rel.IAA, aes(x=SoilType, y=sum.IAA, fill=perc_SWHC))+
  geom_boxplot()+
  theme_bw()+
  ylab("Relative abundance")+
  scale_fill_manual(name = "% SWHC", values = c("blue", "red"))+
  scale_x_discrete(name = "Soil water stress history", labels = c("intermittent", "continuous"))+
  geom_text(data = tukey.IAA.rel, aes(x=x, y=y, label = letters)) #Add letters from tukey
IAA.rel.box
#Anova
shapiro.test(sum.map.rel.IAA$sum.IAA) #P=0.08849
bartlett.test(sum.IAA~perc_SWHC, data = sum.map.rel.IAA) #P=0.168
bartlett.test(sum.IAA~SoilType, data = sum.map.rel.IAA) #P=0.0.01482
summary(aov(sum.IAA~SoilType*perc_SWHC+block, data = sum.map.rel.IAA))

#Df    Sum Sq   Mean Sq F value   Pr(>F)    
#SoilType            1 7.281e-09 7.281e-09  49.065 1.42e-05 ***
#  perc_SWHC           1 1.697e-09 1.697e-09  11.436  0.00545 ** 
#  block               4 9.720e-10 2.430e-10   1.638  0.22848    
#SoilType:perc_SWHC  1 5.320e-10 5.320e-10   3.582  0.08276 .  
#Residuals          12 1.781e-09 1.480e-10  

TukeyHSD(aov(sum.IAA~SoilType*perc_SWHC+block, data = sum.map.rel.IAA))

#diff           lwr           upr     p adj
#NI:50-IR:50 -4.847270e-05 -7.134689e-05 -2.559851e-05 0.0002028
#IR:5-IR:50   8.112223e-06 -1.476197e-05  3.098641e-05 0.7229628
#NI:5-IR:50  -1.973763e-05 -4.261182e-05  3.136559e-06 0.0997398
#IR:5-NI:50   5.658492e-05  3.371073e-05  7.945911e-05 0.0000459
#NI:5-NI:50   2.873507e-05  5.860876e-06  5.160926e-05 0.0132386
#NI:5-IR:5   -2.784986e-05 -5.072405e-05 -4.975664e-06 0.0161891
#Letters for plot
#IR50 IR5 NI50  NI5
# ac  a   b     c

#Difference - %
#Soil type
(mean(sum.map.rel.IAA$sum.IAA[sum.map.rel.IAA$SoilType=="IR"]) -
  mean(sum.map.rel.IAA$sum.IAA[sum.map.rel.IAA$SoilType=="NI"])) /
  mean(sum.map.rel.IAA$sum.IAA[sum.map.rel.IAA$SoilType=="IR"])*100
#3.336975
#SWHC
(mean(sum.map.rel.IAA$sum.IAA[sum.map.rel.IAA$perc_SWHC=="5"]) -
  mean(sum.map.rel.IAA$sum.IAA[sum.map.rel.IAA$perc_SWHC=="50"])) /
  mean(sum.map.rel.IAA$sum.IAA[sum.map.rel.IAA$perc_SWHC=="5"])*100
#1.625061

#Permanova on the gene table
rownames(map.s)==rownames(gene.rel.IAA)#Check ordering
#Specify blocks
perm <-how(nperm = 9999)
setBlocks (perm) <- with(map.s, block)
set.seed(123)
adonis2(gene.rel.IAA~perc_SWHC*SoilType,data=map.s, permutations = perm, method = "bray")

#Df SumOfSqs      R2      F Pr(>F)    
#perc_SWHC           1  0.08397 0.05960 1.5428 0.1297    
#SoilType            1  0.39036 0.27709 7.1722 0.0001 ***
#  perc_SWHC:SoilType  1  0.06364 0.04517 1.1692 0.2592    
#Residual           16  0.87083 0.61814                  
#Total              19  1.40880 1.00000 


##ACC deaminase - ethylene biosynthesis
##K01505: 1-aminocyclopropane-1-carboxylate deaminase
#Get lines that match K01505
annot.ACC <- annot[annot$kegg_entry=="K01505",] #535 obs. of 33 variables
##Normalized by relative abundance
#Get genes abundance that match the ACC annotations
gene.rel.ACC <- gene.rel[,colnames(gene.rel) %in% annot.ACC$gene_id] #20 obs of 535 variables
#Save intermediate files
saveRDS(gene.rel.ACC, here("data","intermediate","gene.rel.ACC.RDS"))
#Look at sum of ACC-related genes
rownames(map.s)==rownames(gene.rel.ACC)#check sorting
sum.map.rel.ACC <- data.frame(map.s, "sum.ACC" = rowSums(gene.rel.ACC), "log.sum.ACC" = log10(rowSums(gene.rel.ACC)))
#Tukey letters - from ANOVA below
tukey.ACC.rel <- data.frame(letters = c("a", "a", "ab", "b"), x = c(0.8, 1.2, 1.8, 2.2), SoilType=c("IR","IR", "NI", "NI") , perc_SWHC= c("50","5","50","5"), y=c(5.2e-5, 5.3e-5, 5.2e-5, 5.7e-5))
#Plot
sum.map.rel.ACC$perc_SWHC <- factor(sum.map.rel.ACC$perc_SWHC, c("50", "5")) #Reorder manually
ACC.rel.box <- ggplot(data = sum.map.rel.ACC, aes(x=SoilType, y=sum.ACC, fill=perc_SWHC))+
  geom_boxplot()+
  theme_bw()+
  ylab("Relative abundance")+
  scale_fill_manual(name = "% SWHC", values = c("blue", "red"))+
  scale_x_discrete(name = "Soil water stress history", labels = c("intermittent", "continuous"))+
  geom_text(data = tukey.ACC.rel, aes(x=x, y=y, label = letters)) #Add letters from tukey
ACC.rel.box
#Anova
shapiro.test(sum.map.rel.ACC$sum.ACC) #P=0.6948
bartlett.test(sum.ACC~perc_SWHC, data = sum.map.rel.ACC) #P=0.2167
bartlett.test(sum.ACC~SoilType, data = sum.map.rel.ACC) #P=0.9243
summary(aov(sum.ACC~SoilType*perc_SWHC+block, data = sum.map.rel.ACC))

#Df    Sum Sq   Mean Sq F value  Pr(>F)   
#SoilType            1 7.383e-11 7.383e-11  14.064 0.00277 **
#  perc_SWHC           1 2.627e-11 2.627e-11   5.004 0.04505 * 
#  block               4 4.917e-11 1.229e-11   2.342 0.11394   
#SoilType:perc_SWHC  1 5.010e-12 5.010e-12   0.954 0.34796   
#Residuals          12 6.300e-11 5.250e-12              

TukeyHSD(aov(sum.ACC~SoilType*perc_SWHC+block, data = sum.map.rel.ACC))

#diff           lwr          upr     p adj
#NI:50-IR:50  2.841847e-06 -1.460360e-06 7.144054e-06 0.2551676
#IR:5-IR:50   1.291163e-06 -3.011044e-06 5.593370e-06 0.8096234
#NI:5-IR:50   6.134731e-06  1.832524e-06 1.043694e-05 0.0055179
#IR:5-NI:50  -1.550684e-06 -5.852891e-06 2.751523e-06 0.7132567
#NI:5-NI:50   3.292885e-06 -1.009323e-06 7.595092e-06 0.1593403
#NI:5-IR:5    4.843569e-06  5.413616e-07 9.145776e-06 0.0260882
#Letters for plot
#IR50 IR5 NI50  NI5
# a   a   ab     b

#Difference - %
#Soil type
(mean(sum.map.rel.ACC$sum.ACC[sum.map.rel.ACC$SoilType=="NI"]) -
    mean(sum.map.rel.ACC$sum.ACC[sum.map.rel.ACC$SoilType=="IR"])) /
  mean(sum.map.rel.ACC$sum.ACC[sum.map.rel.ACC$SoilType=="NI"])*100
#7.435555
#SWHC
(mean(sum.map.rel.ACC$sum.ACC[sum.map.rel.ACC$perc_SWHC=="5"]) -
    mean(sum.map.rel.ACC$sum.ACC[sum.map.rel.ACC$perc_SWHC=="50"])) /
  mean(sum.map.rel.ACC$sum.ACC[sum.map.rel.ACC$perc_SWHC=="5"])*100
#4.502566

#Permanova on the gene table
rownames(map.s)==rownames(gene.rel.ACC)#Check ordering
#Specify blocks
perm <-how(nperm = 9999)
setBlocks (perm) <- with(map.s, block)
set.seed(123)
adonis2(gene.rel.ACC~perc_SWHC*SoilType,data=map.s, permutations = perm, method = "bray")

#Df SumOfSqs      R2      F Pr(>F)    
#perc_SWHC           1  0.09020 0.05835 1.5745 0.1396    
#SoilType            1  0.46466 0.30060 8.1109 0.0001 ***
#  perc_SWHC:SoilType  1  0.07430 0.04806 1.2969 0.2209    
#Residual           16  0.91662 0.59298                  
#Total              19  1.54578 1.00000   

##Osmolytes
##Use headers of supplementary material 2.xlsx of 
##McParland EL, Alexander H & Johnson WM (2021)
##The Osmolyte Ties That Bind: Genomic Insights Into 
##Synthesis and Breakdown of Organic Osmolytes in Marine Microbes.
## Frontiers in Marine Science 8.
osmo.list <- read.table(file = here("data", "raw", "osmo.list.txt"))
#Get annotation lines that match our IAA list
annot.osmo <- annot[annot$kegg_entry %in% osmo.list$V1,] #326,208 obs. of 33 variables
##Normalized by relative abundance
#Get genes abundance that match the osmo annotations
gene.rel.osmo <- gene.rel[,colnames(gene.rel) %in% annot.osmo$gene_id] #20 obs of 326,208 variables
#Save intermediate files
saveRDS(gene.rel.osmo, here("data","intermediate","gene.rel.osmo.RDS"))
#Look at sum of osmo-related genes
rownames(map.s)==rownames(gene.rel.osmo)#check sorting
sum.map.rel.osmo <- data.frame(map.s, "sum.osmo" = rowSums(gene.rel.osmo), "log.sum.osmo" = log10(rowSums(gene.rel.osmo)))
#Tukey letters - from ANOVA below
tukey.osmo.rel <- data.frame(letters = c("a", "a", "b", "b"), x = c(0.8, 1.2, 1.8, 2.2), SoilType=c("IR","IR", "NI", "NI") , perc_SWHC= c("50","5","50","5"), y=c(0.034, 0.034, 0.033, 0.0335))
#Plot
sum.map.rel.osmo$perc_SWHC <- factor(sum.map.rel.osmo$perc_SWHC, c("50", "5")) #Reorder manually
osmo.rel.box <- ggplot(data = sum.map.rel.osmo, aes(x=SoilType, y=sum.osmo, fill=perc_SWHC))+
  geom_boxplot()+
  theme_bw()+
  ylab("Relative abundance")+
  scale_fill_manual(name = "% SWHC", values = c("blue", "red"))+
  scale_x_discrete(name = "Soil water stress history", labels = c("intermittent", "continuous"))+
  geom_text(data = tukey.osmo.rel, aes(x=x, y=y, label = letters)) #Add letters from tukey
osmo.rel.box
#Anova
shapiro.test(sum.map.rel.osmo$sum.osmo) #P=0.817
bartlett.test(sum.osmo~perc_SWHC, data = sum.map.rel.osmo) #P=0.7955
bartlett.test(sum.osmo~SoilType, data = sum.map.rel.osmo) #P=0.7038
summary(aov(sum.osmo~SoilType*perc_SWHC+block, data = sum.map.rel.osmo))

#Df    Sum Sq   Mean Sq F value   Pr(>F)    
#SoilType            1 3.447e-06 3.447e-06  35.707 6.45e-05 ***
#  perc_SWHC           1 7.400e-08 7.400e-08   0.762    0.400    
#block               4 6.680e-07 1.670e-07   1.730    0.208    
#SoilType:perc_SWHC  1 0.000e+00 0.000e+00   0.001    0.982    
#Residuals          12 1.158e-06 9.700e-08                     

TukeyHSD(aov(sum.osmo~SoilType*perc_SWHC+block, data = sum.map.rel.osmo))

#diff           lwr           upr     p adj
#NI:50-IR:50 -0.0008335263 -0.0014169472 -0.0002501054 0.0054413
#IR:5-IR:50   0.0001180992 -0.0004653218  0.0007015201 0.9297974
#NI:5-IR:50  -0.0007090348 -0.0012924558 -0.0001256139 0.0163769
#IR:5-NI:50   0.0009516255  0.0003682045  0.0015350464 0.0019702
#NI:5-NI:50   0.0001244915 -0.0004589294  0.0007079124 0.9192061
#NI:5-IR:5   -0.0008271340 -0.0014105549 -0.0002437131 0.0057545
#Letters for plot
#IR50 IR5 NI50  NI5
# a   a   b     b

#Difference - %
#Soil type
(mean(sum.map.rel.osmo$sum.osmo[sum.map.rel.osmo$SoilType=="IR"]) -
    mean(sum.map.rel.osmo$sum.osmo[sum.map.rel.osmo$SoilType=="NI"])) /
  mean(sum.map.rel.osmo$sum.osmo[sum.map.rel.osmo$SoilType=="IR"])*100
#2.484511


#Permanova on the gene table
rownames(map.s)==rownames(gene.rel.osmo)#Check ordering
#Specify blocks
perm <-how(nperm = 9999)
setBlocks (perm) <- with(map.s, block)
set.seed(123)
adonis2(gene.rel.osmo~perc_SWHC*SoilType,data=map.s, permutations = perm, method = "bray")

#Df SumOfSqs      R2      F Pr(>F)    
#perc_SWHC           1  0.07562 0.05482 1.4437 0.1682    
#SoilType            1  0.40633 0.29456 7.7578 0.0001 ***
#  perc_SWHC:SoilType  1  0.05946 0.04310 1.1352 0.2779    
#Residual           16  0.83804 0.60752                  
#Total              19  1.37945 1.00000 

##EPS - soil aggregation
#genes in KEGG
#Mode: Single Entry to Database
#From: KEGG PATHWAY map00543
#To: KEGG ORTHOLOGY
#Hits: 73 from 1 database
#
#ID                   Definition
#----------------------------------------------------------------------------------------------------
#K00640               cysE; serine O-acetyltransferase [EC:2.3.1.30] 
#K02851               wecA, tagO, rfe; UDP-GlcNAc:undecaprenyl-phosphate/decaprenyl-phosphate GlcNAc-1-phosphate transfera 
#K02852               wecG, rffM; UDP-N-acetyl-D-mannosaminouronate:lipid I N-acetyl-D-mannosaminouronosyltransferase [EC: 
#K03208               wcaI; putative colanic acid biosynthesis glycosyltransferase WcaI 
#K03606               wcaJ; undecaprenyl-phosphate glucose phosphotransferase [EC:2.7.8.31] 
#K03818               wcaF; putative colanic acid biosynthesis acetyltransferase WcaF [EC:2.3.1.-] 
#K03819               wcaB; putative colanic acid biosynthesis acetyltransferase WcaB [EC:2.3.1.-] 
#K11936               pgaC, icaA; poly-beta-1,6-N-acetyl-D-glucosamine synthase [EC:2.4.1.-] 
#K11937               pgaD; biofilm PGA synthesis protein PgaD 
#K12582               wecF, rffT; dTDP-N-acetylfucosamine:lipid II N-acetylfucosaminyltransferase [EC:2.4.1.325] 
#K13656               gumD; undecaprenyl-phosphate glucose phosphotransferase [EC:2.7.8.31] 
#K13657               gumH, aceA; alpha-1,3-mannosyltransferase [EC:2.4.1.252] 
#K13658               gumI; beta-1,4-mannosyltransferase [EC:2.4.1.251] 
#K13659               gumK; 2-beta-glucuronyltransferase [EC:2.4.1.264] 
#K13660               gumM; beta-1,4-glucosyltransferase [EC:2.4.1.-] 
#K13663               gumF; acyltransferase [EC:2.3.1.-] 
#K13664               gumG; acyltransferase [EC:2.3.1.-] 
#K13665               gumL; pyruvyltransferase 
#K13683               wcaE; putative colanic acid biosynthesis glycosyltransferase WcaE [EC:2.4.-.-] 
#K13684               wcaC; putative colanic acid biosynthesis glycosyltransferase WcaC [EC:2.4.-.-] 
#K16555               exoO; succinoglycan biosynthesis protein ExoO [EC:2.4.-.-] 
#K16556               exoM; succinoglycan biosynthesis protein ExoM [EC:2.4.-.-] 
#K16557               exoA; succinoglycan biosynthesis protein ExoA [EC:2.4.-.-] 
#K16558               exoL; succinoglycan biosynthesis protein ExoL [EC:2.-.-.-] 
#K16560               exoH; succinoglycan biosynthesis protein ExoH 
#K16562               exoW; succinoglycan biosynthesis protein ExoW [EC:2.4.-.-] 
#K16563               exoV; succinoglycan biosynthesis protein ExoV 
#K16564               exoU; succinoglycan biosynthesis protein ExoU [EC:2.4.-.-] 
#K16566               exoY; exopolysaccharide production protein ExoY 
#K16568               exoZ; exopolysaccharide production protein ExoZ 
#K16700               amsB, cpsE; amylovoran/stewartan biosynthesis glycosyltransferase AmsB/CpsE [EC:2.4.-.-] 
#K16701               amsD; amylovoran biosynthesis glycosyltransferase AmsD [EC:2.4.-.-] 
#K16702               amsE; amylovoran biosynthesis glycosyltransferase AmsE [EC:2.4.-.-] 
#K16703               wcaL, amsK, cpsK; colanic acid/amylovoran/stewartan biosynthesis glycosyltransferase WcaL/AmsK/CpsK  
#K16707               amsG, cpsA; UDP-galactose-lipid carrier transferase 
#K16710               wcaK, amsJ; colanic acid/amylovoran biosynthesis protein WcaK/AmsJ 
#K19290               alg8; mannuronan synthase [EC:2.4.1.33] 
#K19291               alg44; mannuronan synthase [EC:2.4.1.33] 
#K19293               algX; alginate biosynthesis protein AlgX 
#K19294               algI; alginate O-acetyltransferase complex protein AlgI 
#K19295               algJ; alginate O-acetyltransferase complex protein AlgJ 
#K19296               algF; alginate O-acetyltransferase complex protein AlgF 
#K20921               vpsD, epsF; polysaccharide biosynthesis protein VpsD 
#K20922               vpsI; polysaccharide biosynthesis protein VpsI 
#K20997               pslA; polysaccharide biosynthesis protein PslA 
#K20999               pslF; polysaccharide biosynthesis protein PslF 
#K21001               pslH; polysaccharide biosynthesis protein PslH 
#K21002               pslI; polysaccharide biosynthesis protein PslI 
#K21154               pssM; exopolysaccharide glucosyl ketal-pyruvate-transferase [EC:2.5.1.98] 
#K21461               icaD; poly-beta-1,6-N-acetyl-D-glucosamine synthesis protein 
#K25205               pslC; polysaccharide biosynthesis protein PslC 
#K25875               wcaA; putative colanic acid biosynthesis glycosyltransferase WcaA 
#K25886               vpsK; polysaccharide biosynthesis protein VpsK 
#K25887               aceA; undecaprenyl-phosphate glucose phosphotransferase [EC:2.7.8.31] 
#K25888               aceB; beta-1,4-glucosyltransferase [EC:2.4.1.-] 
#K25889               aceQ; glucosyltransferase 
#K25890               aceP; beta-1,6-glucosyltransferase 
#K25891               aceR; rhamnosyltransferase 
#K25892               aceI; acyltransferase 
#K25902               pssA; acidic exopolysaccharide biosynthesis protein PssA 
#K25903               pssD; exopolysaccharide biosynthesis glucuronosyltransferase PssD 
#K25904               pssE; exopolysaccharide biosynthesis glucuronosyltransferase PssE 
#K25905               pssC; exopolysaccharide biosynthesis glycosyltransferase PssC 
#K25906               pssS; exopolysaccharide biosynthesis glucosyltransferase PssS [EC:2.4.1.-] 
#K25907               pssF; exopolysaccharide biosynthesis glycosyltransferase PssF 
#K25908               pssI; exopolysaccharide biosynthesis glycosyltransferase PssI 
#K25909               pssG; exopolysaccharide biosynthesis glycosyltransferase PssG 
#K25910               pssH; exopolysaccharide biosynthesis glycosyltransferase PssH 
#K25911               pssJ; exopolysaccharide biosynthesis galactosyltransferase PssJ 
#K25912               pssR; exopolysaccharide acyltransferase PssR 
#K25913               pssK; exopolysaccharide biosynthesis protein PssK 
#K25957               cpsG; stewartan biosynthesis glycosyltransferase CpsG 
#K25958               cpsF; stewartan biosynthesis galactosyltransferase CpsF                               

EPS.list <- read.table(file = here("data", "raw", "EPS.list.txt"))
#Get annotation lines that match our IAA list
annot.EPS <- annot[annot$kegg_entry %in% EPS.list$V1,] #5351 obs. of 33 variables
##Normalized by relative abundance
#Get genes abundance that match the EPS annotations
gene.rel.EPS <- gene.rel[,colnames(gene.rel) %in% annot.EPS$gene_id] #20 obs of 5351 variables
#Save intermediate files
saveRDS(gene.rel.EPS, here("data","intermediate","gene.rel.EPS.RDS"))
#Look at sum of EPS-related genes
rownames(map.s)==rownames(gene.rel.EPS)#check sorting
sum.map.rel.EPS <- data.frame(map.s, "sum.EPS" = rowSums(gene.rel.EPS), "log.sum.EPS" = log10(rowSums(gene.rel.EPS)))
#Tukey letters - from ANOVA below
tukey.EPS.rel <- data.frame(letters = c("a", "a", "b", "a"), x = c(0.8, 1.2, 1.8, 2.2), SoilType=c("IR","IR", "NI", "NI") , perc_SWHC= c("50","5","50","5"), y=c(0.00049, 0.00049, 0.000505, 0.00049))
#Plot
sum.map.rel.EPS$perc_SWHC <- factor(sum.map.rel.EPS$perc_SWHC, c("50", "5")) #Reorder manually
EPS.rel.box <- ggplot(data = sum.map.rel.EPS, aes(x=SoilType, y=sum.EPS, fill=perc_SWHC))+
  geom_boxplot()+
  theme_bw()+
  ylab("Relative abundance")+
  scale_fill_manual(name = "% SWHC", values = c("blue", "red"))+
  scale_x_discrete(name = "Soil water stress history", labels = c("intermittent", "continuous"))+
  geom_text(data = tukey.EPS.rel, aes(x=x, y=y, label = letters)) #Add letters from tukey
EPS.rel.box
#Anova
shapiro.test(sum.map.rel.EPS$sum.EPS) #P=0.2782
bartlett.test(sum.EPS~perc_SWHC, data = sum.map.rel.EPS) #P=0.3793
bartlett.test(sum.EPS~SoilType, data = sum.map.rel.EPS) #P=0.4493
summary(aov(sum.EPS~SoilType*perc_SWHC+block, data = sum.map.rel.EPS))

#Df    Sum Sq   Mean Sq F value   Pr(>F)    
#SoilType            1 3.306e-10 3.306e-10  19.655 0.000816 ***
#  perc_SWHC           1 2.863e-10 2.863e-10  17.019 0.001407 ** 
#  block               4 2.361e-10 5.900e-11   3.510 0.040575 *  
#  SoilType:perc_SWHC  1 3.740e-11 3.740e-11   2.221 0.161965    
#Residuals          12 2.018e-10 1.680e-11  

TukeyHSD(aov(sum.EPS~SoilType*perc_SWHC+block, data = sum.map.rel.EPS))

#diff           lwr           upr     p adj
#NI:50-IR:50  1.086478e-05  3.163869e-06  1.856570e-05 0.0059608
#IR:5-IR:50  -4.833135e-06 -1.253405e-05  2.867779e-06 0.2931644
#NI:5-IR:50   5.649744e-07 -7.135940e-06  8.265888e-06 0.9961482
#IR:5-NI:50  -1.569792e-05 -2.339883e-05 -7.997004e-06 0.0002899
#NI:5-NI:50  -1.029981e-05 -1.800072e-05 -2.598895e-06 0.0086909
#NI:5-IR:5    5.398109e-06 -2.302805e-06  1.309902e-05 0.2138470
#Letters for plot
#IR50 IR5 NI50  NI5
# a   a   b     a

#Difference - %
#Soil type
(mean(sum.map.rel.EPS$sum.EPS[sum.map.rel.EPS$SoilType=="NI"]) -
    mean(sum.map.rel.EPS$sum.EPS[sum.map.rel.EPS$SoilType=="IR"])) /
  mean(sum.map.rel.EPS$sum.EPS[sum.map.rel.EPS$SoilType=="NI"])*100
#1.673764
#SWHC
(mean(sum.map.rel.EPS$sum.EPS[sum.map.rel.EPS$perc_SWHC=="50"]) -
    mean(sum.map.rel.EPS$sum.EPS[sum.map.rel.EPS$perc_SWHC=="5"])) /
  mean(sum.map.rel.EPS$sum.EPS[sum.map.rel.EPS$perc_SWHC=="50"])*100
#1.558377

#Permanova on the gene table
rownames(map.s)==rownames(gene.rel.EPS)#Check ordering
#Specify blocks
perm <-how(nperm = 9999)
setBlocks (perm) <- with(map.s, block)
set.seed(123)
adonis2(gene.rel.EPS~perc_SWHC*SoilType,data=map.s, permutations = perm, method = "bray")

#Df SumOfSqs      R2      F Pr(>F)    
#perc_SWHC           1  0.08281 0.05484 1.4108 0.1733    
#SoilType            1  0.42273 0.27997 7.2022 0.0001 ***
#  perc_SWHC:SoilType  1  0.06526 0.04322 1.1119 0.2825    
#Residual           16  0.93911 0.62197                  
#Total              19  1.50991 1.00000 

##cytokinines - stomatal closure
#genes in KEGG
#Mode: Single Entry to Database
#From: KEGG PATHWAY map00908
#To: KEGG ORTHOLOGY
#Hits: 10 from 1 database
#
#ID                   Definition
#----------------------------------------------------------------------------------------------------
#K00279               CKX; cytokinin dehydrogenase [EC:1.5.99.12] 
#K00791               miaA, TRIT1; tRNA dimethylallyltransferase [EC:2.5.1.75] 
#K10717               CYP735A; cytokinin trans-hydroxylase 
#K10760               IPT; adenylate dimethylallyltransferase (cytokinin synthase) [EC:2.5.1.27 2.5.1.112] 
#K13492               ZOG1; zeatin O-glucosyltransferase [EC:2.4.1.203] 
#K13493               UGT76C1_2; cytokinin-N-glucosyltransferase [EC:2.4.1.-] 
#K13494               ZOX1; zeatin O-xylosyltransferase [EC:2.4.2.40] 
#K13495               CISZOG; cis-zeatin O-glucosyltransferase [EC:2.4.1.215] 
#K13496               UGT73C; UDP-glucosyltransferase 73C [EC:2.4.1.-] 
#K23452               UGT85A; UDP-glucosyltransferase 85A [EC:2.4.1.-] 

cyto.list <- c("K00279", "K00791", "K10717","K10760","K13492","K13493","K13494", "K13495", "K13496","K23452")
#Get annotation lines that match our IAA list
annot.cyto <- annot[annot$kegg_entry %in% cyto.list,] #2,406 obs. of 33 variables
##Normalized by relative abundance
#Get genes abundance that match the cyto annotations
gene.rel.cyto <- gene.rel[,colnames(gene.rel) %in% annot.cyto$gene_id] #20 obs of 2406 variables
#Save intermediate files
saveRDS(gene.rel.cyto, here("data","intermediate","gene.rel.cyto.RDS"))
#Look at sum of cyto-related genes
rownames(map.s)==rownames(gene.rel.cyto)#check sorting
sum.map.rel.cyto <- data.frame(map.s, "sum.cyto" = rowSums(gene.rel.cyto), "log.sum.cyto" = log10(rowSums(gene.rel.cyto)))
#Tukey letters - from ANOVA below
#tukey.cyto.rel <- data.frame(letters = c("a", "a", "b", "b"), x = c(0.8, 1.2, 1.8, 2.2), SoilType=c("IR","IR", "NI", "NI") , perc_SWHC= c("50","5","50","5"), y=c(0.034, 0.034, 0.033, 0.0335))
#Plot
sum.map.rel.cyto$perc_SWHC <- factor(sum.map.rel.cyto$perc_SWHC, c("50", "5")) #Reorder manually
cyto.rel.box <- ggplot(data = sum.map.rel.cyto, aes(x=SoilType, y=sum.cyto, fill=perc_SWHC))+
  geom_boxplot()+
  theme_bw()+
  ylab("Relative abundance")+
  scale_x_discrete(name = "Soil water stress history", labels = c("intermittent", "continuous"))+
  scale_fill_manual(name = "% SWHC", values = c("blue", "red"))#+
  #geom_text(data = tukey.cyto.rel, aes(x=x, y=y, label = letters)) #Add letters from tukey
cyto.rel.box
#Anova
shapiro.test(sum.map.rel.cyto$sum.cyto) #P=0.7469
bartlett.test(sum.cyto~perc_SWHC, data = sum.map.rel.cyto) #P=0.5273
bartlett.test(sum.cyto~SoilType, data = sum.map.rel.cyto) #P=0.9516
summary(aov(sum.cyto~SoilType*perc_SWHC+block, data = sum.map.rel.cyto))

#Df    Sum Sq   Mean Sq F value Pr(>F)
#SoilType            1 1.500e-13 1.520e-13   0.025  0.877
#perc_SWHC           1 5.100e-13 5.060e-13   0.084  0.777
#block               4 3.097e-11 7.743e-12   1.280  0.332
#SoilType:perc_SWHC  1 1.830e-11 1.830e-11   3.024  0.108
#Residuals          12 7.261e-11 6.051e-12                        

TukeyHSD(aov(sum.cyto~SoilType*perc_SWHC+block, data = sum.map.rel.cyto))

#diff           lwr          upr     p adj
#NI:5-IR:5    2.087539e-06 -2.531279e-06 6.706356e-06 0.5558530
#IR:50-IR:5   1.594948e-06 -3.023869e-06 6.213766e-06 0.7384327
#NI:50-IR:5  -1.434749e-07 -4.762292e-06 4.475342e-06 0.9997023
#IR:50-NI:5  -4.925900e-07 -5.111407e-06 4.126227e-06 0.9884535
#NI:50-NI:5  -2.231013e-06 -6.849831e-06 2.387804e-06 0.5034740
#NI:50-IR:50 -1.738423e-06 -6.357240e-06 2.880394e-06 0.6862337
#Letters for plot
#IR50 IR5 NI50  NI5
# a   a   a     a

#Permanova on the gene table
rownames(map.s)==rownames(gene.rel.cyto)#Check ordering
#Specify blocks
perm <-how(nperm = 9999)
setBlocks (perm) <- with(map.s, block)
set.seed(123)
adonis2(gene.rel.cyto~perc_SWHC*SoilType,data=map.s, permutations = perm, method = "bray")

#Df SumOfSqs      R2      F Pr(>F)    
#perc_SWHC           1  0.08727 0.05645 1.4424 0.1593    
#SoilType            1  0.42668 0.27601 7.0522 0.0002 ***
#  perc_SWHC:SoilType  1  0.06389 0.04133 1.0560 0.3344    
#Residual           16  0.96806 0.62621                  
#Total              19  1.54591 1.00000  

##antioxydants - ROS suppression
#genes in KEGG - superoxide dismutase, gluthatione peroxidase and cytochrome c oxidase
#
#K00518
#sodN; nickel superoxide dismutase [EC:1.15.1.1]
#K04564
#SOD2; superoxide dismutase, Fe-Mn family [EC:1.15.1.1]
#K04565
#SOD1; superoxide dismutase, Cu-Zn family [EC:1.15.1.1]
#K04569
#CCS; copper chaperone for superoxide dismutase
#K16627
#SOD3; superoxide dismutase, Cu-Zn family [EC:1.15.1.1]
#K03781
#katE, CAT, catB, srpA; catalase [EC:1.11.1.6]
#K03782
#katG; catalase-peroxidase [EC:1.11.1.21]
#K07217
#ydbD; manganese catalase [EC:1.11.1.6]
#K19885
#rebD; dichlorochromopyrrolate synthase / catalase [EC:1.21.98.2 1.11.1.6]
#K00432
#gpx, btuE, bsaA; glutathione peroxidase [EC:1.11.1.9]
#K05361
#GPX4; phospholipid-hydroperoxide glutathione peroxidase [EC:1.11.1.12]
#K11207
#TDPX; glutathione peroxidase-type tryparedoxin peroxidase [EC:1.11.1.-]
#K23075
#garA; glutathione amide-dependent peroxidase [EC:1.11.1.17]
#K00404
#ccoN; cytochrome c oxidase cbb3-type subunit I [EC:7.1.1.9]
#K00405
#ccoO; cytochrome c oxidase cbb3-type subunit II
#K00406
#ccoP; cytochrome c oxidase cbb3-type subunit III
#K00407
#ccoQ; cytochrome c oxidase cbb3-type subunit IV
#K00424
#cydX; cytochrome bd-I ubiquinol oxidase subunit X [EC:7.1.1.7]
#K00425
#cydA; cytochrome bd ubiquinol oxidase subunit I [EC:7.1.1.7]
#K00426
#cydB; cytochrome bd ubiquinol oxidase subunit II [EC:7.1.1.7]
#K00428
#E1.11.1.5; cytochrome c peroxidase [EC:1.11.1.5]
#K02256
#COX1; cytochrome c oxidase subunit 1 [EC:7.1.1.9]
#K02258
#COX11, ctaG; cytochrome c oxidase assembly protein subunit 11
#K02260
#COX17; cytochrome c oxidase assembly protein subunit 17
#K02261
#COX2; cytochrome c oxidase subunit 2
#K02262
#COX3; cytochrome c oxidase subunit 3
#K02263
#COX4; cytochrome c oxidase subunit 4
#K02264
#COX5A; cytochrome c oxidase subunit 5a
#K02265
#COX5B; cytochrome c oxidase subunit 5b
#K02266
#COX6A; cytochrome c oxidase subunit 6a
#K02267
#COX6B; cytochrome c oxidase subunit 6b
#K02268
#COX6C; cytochrome c oxidase subunit 6c
#K02269
#COX7; cytochrome c oxidase subunit 7
#K02270
#COX7A; cytochrome c oxidase subunit 7a
#K02271
#COX7B; cytochrome c oxidase subunit 7b
#K02272
#COX7C; cytochrome c oxidase subunit 7c
#K02273
#COX8; cytochrome c oxidase subunit 8
#K02274
#coxA, ctaD; cytochrome c oxidase subunit I [EC:7.1.1.9]
#K02275
#coxB, ctaC; cytochrome c oxidase subunit II [EC:7.1.1.9]
#K02276
#coxC, ctaE; cytochrome c oxidase subunit III [EC:7.1.1.9]
#K02277
#coxD, ctaF; cytochrome c oxidase subunit IV [EC:7.1.1.9]
#K02297
#cyoA; cytochrome o ubiquinol oxidase subunit II [EC:7.1.1.3]
#K02298
#cyoB; cytochrome o ubiquinol oxidase subunit I [EC:7.1.1.3]
#K02299
#cyoC; cytochrome o ubiquinol oxidase subunit III
#K02300
#cyoD; cytochrome o ubiquinol oxidase subunit IV
#K02826
#qoxA; cytochrome aa3-600 menaquinol oxidase subunit II [EC:7.1.1.5]
#K02827
#qoxB; cytochrome aa3-600 menaquinol oxidase subunit I [EC:7.1.1.5]
#K02828
#qoxC; cytochrome aa3-600 menaquinol oxidase subunit III [EC:7.1.1.5]
#K02829
#qoxD; cytochrome aa3-600 menaquinol oxidase subunit IV [EC:7.1.1.5]
#K15408
#coxAC; cytochrome c oxidase subunit I+III [EC:7.1.1.9]
#K15862
#ccoNO; cytochrome c oxidase cbb3-type subunit I/II [EC:7.1.1.9]
#K18173
#COA1; cytochrome c oxidase assembly factor 1
#K18174
#COA2; cytochrome c oxidase assembly factor 2
#K18175
#CCDC56, COA3; cytochrome c oxidase assembly factor 3, animal type
#K18176
#COA3; cytochrome c oxidase assembly factor 3, fungi type
#K18177
#COA4; cytochrome c oxidase assembly factor 4
#K18178
#COA5, PET191; cytochrome c oxidase assembly factor 5
#K18179
#COA6; cytochrome c oxidase assembly factor 6
#K18180
#COA7, SELRC1, RESA1; cytochrome c oxidase assembly factor 7
#K18181
#COX14; cytochrome c oxidase assembly factor 14
#K18182
#COX16; cytochrome c oxidase assembly protein subunit 16
#K18183
#COX19; cytochrome c oxidase assembly protein subunit 19
#K18184
#COX20; cytochrome c oxidase assembly protein subunit 20
#K18185
#COX23; cytochrome c oxidase assembly protein subunit 23
#K18189
#TACO1; translational activator of cytochrome c oxidase 1
#K22501
#appX; cytochrome bd-II ubiquinol oxidase subunit AppX [EC:7.1.1.7]
#K24007
#soxD; cytochrome aa3-type oxidase subunit SoxD
#K24008
#soxC; cytochrome aa3-type oxidase subunit III
#K24009
#soxB; cytochrome aa3-type oxidase subunit I [EC:7.1.1.4]
#K24010
#soxA; cytochrome aa3-type oxidase subunit II [EC:7.1.1.4]
#K24011
#soxM; cytochrome aa3-type oxidase subunit I/III [EC:7.1.1.4]

antiox.list <- read.table(file = here("data", "raw", "antiox.list.txt"))
#Get annotation lines that match our IAA list
annot.antiox <- annot[annot$kegg_entry %in% antiox.list$V1,] #35,833 obs. of 33 variables
##Normalized by relative abundance
#Get genes abundance that match the antiox annotations
gene.rel.antiox <- gene.rel[,colnames(gene.rel) %in% annot.antiox$gene_id] #20 obs of 35833 variables
#Save intermediate files
saveRDS(gene.rel.antiox, here("data","intermediate","gene.rel.antiox.RDS"))
#Look at sum of antiox-related genes
rownames(map.s)==rownames(gene.rel.antiox)#check sorting
sum.map.rel.antiox <- data.frame(map.s, "sum.antiox" = rowSums(gene.rel.antiox), "log.sum.antiox" = log10(rowSums(gene.rel.antiox)))
#Tukey letters - from ANOVA below
tukey.antiox.rel <- data.frame(letters = c("ab", "a", "b", "a"), x = c(0.8, 1.2, 1.8, 2.2), SoilType=c("IR","IR", "NI", "NI") , perc_SWHC= c("50","5","50","5"), y=c(0.00345, 0.0034, 0.00352, 0.0034))
#Plot
sum.map.rel.antiox$perc_SWHC <- factor(sum.map.rel.antiox$perc_SWHC, c("50", "5")) #Reorder manually
antiox.rel.box <- ggplot(data = sum.map.rel.antiox, aes(x=SoilType, y=sum.antiox, fill=perc_SWHC))+
  geom_boxplot()+
  theme_bw()+
  ylab("Relative abundance")+
  scale_fill_manual(name = "% SWHC", values = c("blue", "red"))+
  scale_x_discrete(name = "Soil water stress history", labels = c("intermittent", "continuous"))+
  geom_text(data = tukey.antiox.rel, aes(x=x, y=y, label = letters)) #Add letters from tukey
antiox.rel.box
#Anova
shapiro.test(sum.map.rel.antiox$sum.antiox) #P=0.5552
bartlett.test(sum.antiox~perc_SWHC, data = sum.map.rel.antiox) #P=0.5959
bartlett.test(sum.antiox~SoilType, data = sum.map.rel.antiox) #P=0.0259
summary(aov(sum.antiox~SoilType*perc_SWHC+block, data = sum.map.rel.antiox))

#Df    Sum Sq   Mean Sq F value  Pr(>F)   
#SoilType            1 1.025e-08 1.025e-08   3.484 0.08659 . 
#perc_SWHC           1 4.330e-08 4.330e-08  14.709 0.00237 **
#  block               4 6.670e-09 1.670e-09   0.567 0.69168   
#SoilType:perc_SWHC  1 1.042e-08 1.042e-08   3.541 0.08435 . 
#Residuals          12 3.532e-08 2.940e-09                                    

TukeyHSD(aov(sum.antiox~SoilType*perc_SWHC+block, data = sum.map.rel.antiox))

#diff           lwr          upr     p adj
#NI:5-IR:5   -3.697230e-07 -1.022427e-04 0.0001015032 0.9999995
#IR:50-IR:5   4.739905e-05 -5.447390e-05 0.0001492720 0.5332156
#NI:50-IR:5   1.383430e-04  3.647003e-05 0.0002402159 0.0078186
#IR:50-NI:5   4.776878e-05 -5.410418e-05 0.0001496417 0.5270898
#NI:50-NI:5   1.387127e-04  3.683975e-05 0.0002405857 0.0076738
#NI:50-IR:50  9.094393e-05 -1.092903e-05 0.0001928169 0.0860407
#Letters for plot
#IR50 IR5 NI50  NI5
# ab   a   b     a

#Difference - %
#Soil type
(mean(sum.map.rel.antiox$sum.antiox[sum.map.rel.antiox$SoilType=="NI"]) -
    mean(sum.map.rel.antiox$sum.antiox[sum.map.rel.antiox$SoilType=="IR"])) /
  mean(sum.map.rel.antiox$sum.antiox[sum.map.rel.antiox$SoilType=="NI"])*100
#1.342902
#SWHC
(mean(sum.map.rel.antiox$sum.antiox[sum.map.rel.antiox$perc_SWHC=="50"]) -
    mean(sum.map.rel.antiox$sum.antiox[sum.map.rel.antiox$perc_SWHC=="5"])) /
  mean(sum.map.rel.antiox$sum.antiox[sum.map.rel.antiox$perc_SWHC=="50"])*100
#2.739988

#Permanova on the gene table
rownames(map.s)==rownames(gene.rel.antiox)#Check ordering
#Specify blocks
perm <-how(nperm = 9999)
setBlocks (perm) <- with(map.s, block)
set.seed(123)
adonis2(gene.rel.antiox~perc_SWHC*SoilType,data=map.s, permutations = perm, method = "bray")

#Df SumOfSqs      R2      F Pr(>F)    
#perc_SWHC           1  0.07562 0.05482 1.4437 0.1682    
#SoilType            1  0.40633 0.29456 7.7578 0.0001 ***
#  perc_SWHC:SoilType  1  0.05946 0.04310 1.1352 0.2779    
#Residual           16  0.83804 0.60752                  
#Total              19  1.37945 1.00000 