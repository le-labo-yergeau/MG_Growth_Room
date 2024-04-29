##Anova of plant data

#Load data
plant <- readRDS(here("data", "intermediate", "plant.RDS"))
trait <- readRDS(here("data", "intermediate", "trait.RDS"))

#Merge
sum((plant$SWHC==trait$Treatment) & 
      (plant$Cultivar==trait$cultivar) & 
      (plant$Soil_type==trait$Soil_type) &
      (plant$Replicates==trait$Replicates)
      )#160
plant[,7:12] <- trait[,2:7]

#keep only relevant for this study
plant.MG <- plant[plant$SWHC %in% c(5,50) & plant$Cultivar=="ACNass",]
plant.MG$SWHC <- as.character(plant.MG$SWHC)
plant.MG$Replicates <- as.character(plant.MG$Replicates)

#Replace weird value
plant.MG[5,10] <- 100

#Test normality and heteroscedascticity
shapiro.test(plant.MG$Shoot_biomass) #P=0.002085
shapiro.test(plant.MG$Root_biomass) #P=0.00407
shapiro.test(plant.MG$Shoot_dry_biomass) #P=0.001823
shapiro.test(plant.MG$Root_dry_biomass) #P=0.0001028
shapiro.test(plant.MG$Root_length_cm) #P=0.5033
shapiro.test(plant.MG$Leaf_relative_water_content_percentage) #P=0.0097
shapiro.test(plant.MG$Leaf_moisture) #P=0.02015
shapiro.test(plant.MG$leaf_dry_mater_content_mg.g) #P=0.2111

bartlett.test(plant.MG$Shoot_biomass~plant.MG$Soil_type) #P=0.8317
bartlett.test(plant.MG$Root_biomass~plant.MG$Soil_type) #P=0.01504
bartlett.test(plant.MG$Shoot_dry_biomass~plant.MG$Soil_type) #P=0.809
bartlett.test(plant.MG$Root_dry_biomass~plant.MG$Soil_type) #P=0.4458
bartlett.test(plant.MG$Root_length_cm~plant.MG$Soil_type) #P=0.3498
bartlett.test(plant.MG$Leaf_relative_water_content_percentage~plant.MG$Soil_type) #P=0.4736
bartlett.test(plant.MG$Leaf_moisture~plant.MG$Soil_type) #P=0.737
bartlett.test(plant.MG$leaf_dry_mater_content_mg.g~plant.MG$Soil_type) #P=0.4255

bartlett.test(plant.MG$Shoot_biomass~plant.MG$SWHC) #P=0.0001883
bartlett.test(plant.MG$Root_biomass~plant.MG$SWHC) #P=4.641e-06
bartlett.test(plant.MG$Shoot_dry_biomass~plant.MG$SWHC) #P=0.0005088
bartlett.test(plant.MG$Root_dry_biomass~plant.MG$SWHC) #P=3.82e-09
bartlett.test(plant.MG$Root_length_cm~plant.MG$SWHC) #P=0.0009281
bartlett.test(plant.MG$Leaf_relative_water_content_percentage~plant.MG$SWHC) #P=0.001997
bartlett.test(plant.MG$Leaf_moisture~plant.MG$SWHC) #P=0.200
bartlett.test(plant.MG$leaf_dry_mater_content_mg.g~plant.MG$SWHC) #P=0.3127

#Log transform
shapiro.test(log(plant.MG$Shoot_biomass)) #P=0.003372
shapiro.test(log(plant.MG$Root_biomass)) #P=0.06382
shapiro.test(log(plant.MG$Shoot_dry_biomass)) #P=0.003163
shapiro.test(log(plant.MG$Root_dry_biomass)) #P=0.3643
shapiro.test(log(plant.MG$Root_length_cm)) #P=0.2.702e-05
shapiro.test(log(plant.MG$Leaf_relative_water_content_percentage)) #P=0.002189
shapiro.test(log(plant.MG$Leaf_moisture)) #P=0.008252
shapiro.test(log(plant.MG$leaf_dry_mater_content_mg.g)) #P=0.3752

bartlett.test(log(plant.MG$Shoot_biomass)~plant.MG$Soil_type) #P=0.5672
bartlett.test(log(plant.MG$Root_biomass)~plant.MG$Soil_type) #P=0.0864
bartlett.test(log(plant.MG$Shoot_dry_biomass)~plant.MG$Soil_type) #P=0.6215
bartlett.test(log(plant.MG$Root_dry_biomass)~plant.MG$Soil_type) #P=0.5374
bartlett.test(log(plant.MG$Root_length_cm)~plant.MG$Soil_type) #P=0.00693
bartlett.test(log(plant.MG$Leaf_relative_water_content_percentage)~plant.MG$Soil_type) #P=0.2251
bartlett.test(log(plant.MG$Leaf_moisture)~plant.MG$Soil_type) #P=0.3014
bartlett.test(log(plant.MG$leaf_dry_mater_content_mg.g)~plant.MG$Soil_type) #P=0.6127

bartlett.test(log(plant.MG$Shoot_biomass)~plant.MG$SWHC) #P=0.09125
bartlett.test(log(plant.MG$Root_biomass)~plant.MG$SWHC) #P=0.9958
bartlett.test(log(plant.MG$Shoot_dry_biomass)~plant.MG$SWHC) #P=0.8011
bartlett.test(log(plant.MG$Root_dry_biomass)~plant.MG$SWHC) #P=0.1041
bartlett.test(log(plant.MG$Root_length_cm)~plant.MG$SWHC) #P=3.487e-05
bartlett.test(log(plant.MG$Leaf_relative_water_content_percentage)~plant.MG$SWHC) #P=2.455e-05
bartlett.test(log(plant.MG$Leaf_moisture)~plant.MG$SWHC) #P=0.9.235e-05
bartlett.test(log(plant.MG$leaf_dry_mater_content_mg.g)~plant.MG$SWHC) #P=0.9466

#So, 
#Shoot_biomass: log
#Root_biomass: log
#Shoot_dry_biomass: log
#Root_dry_biomass: log
#Root_length: not transformed
#Leaf moisture %: not transformed
#LWC: not transformed
#Leaf dry matter content: not transformed

#Anova
#Shoot fresh biomass
summary(aov(log(plant.MG$Shoot_biomass)~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#Df Sum Sq Mean Sq F value   Pr(>F)    
#plant.MG$SWHC                     1 17.993  17.993 397.465 1.45e-10 ***
#  plant.MG$Soil_type                1  0.282   0.282   6.221   0.0282 *  
#  plant.MG$Replicates               4  0.234   0.058   1.292   0.3274    
#plant.MG$SWHC:plant.MG$Soil_type  1  0.180   0.180   3.984   0.0691 .  
#Residuals                        12  0.543   0.045        
TukeyHSD(aov(log(plant.MG$Shoot_biomass)~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#diff        lwr        upr     p adj
#5:IR           -50:IR             -1.7070488 -2.1065553 -1.3075424 0.0000001
#50:NI           -50:IR            -0.0473959 -0.4469023  0.3521105 0.9842788
#5:NI           -50:IR             -2.1343035 -2.5338099 -1.7347970 0.0000000
#50:NI           -5:IR              1.6596529  1.2601465  2.0591594 0.0000002
#5:NI           -5:IR              -0.4272546 -0.8267610 -0.0277482 0.0349550
#5:NI           -50:NI             -2.0869075 -2.4864140 -1.6874011 0.0000000
#Letters for plot
#50:IR  5:IR  50:NI 5:NI
#  a    b     a     c
group_by(plant.MG, interaction(SWHC,Soil_type)) %>% summarize(mean(Shoot_biomass))
#5.IR                                           0.242
#50.IR                                          1.36 
#5.NI                                           0.162
#50.NI                                          1.26 


#Root fresh biomass
summary(aov(log(plant.MG$Root_biomass)~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#Df Sum Sq Mean Sq F value   Pr(>F)    
#plant.MG$SWHC                     1 16.525  16.525 229.121  3.5e-09 ***
#  plant.MG$Soil_type                1  0.082   0.082   1.132 0.308339    
#plant.MG$Replicates               4  0.251   0.063   0.869 0.510404    
#plant.MG$SWHC:plant.MG$Soil_type  1  1.566   1.566  21.715 0.000551 ***
#  Residuals                        12  0.865   0.072
TukeyHSD(aov(log(plant.MG$Root_biomass)~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#diff        lwr         upr     p adj
#5:IR           -50:IR             -1.2582850 -1.7625522 -0.75401771 0.0000422
#50:NI           -50:IR             0.6874362  0.1831689  1.19170346 0.0076104
#5:NI           -50:IR             -1.6901818 -2.1944491 -1.18591456 0.0000020
#50:NI           -5:IR              1.9457212  1.4414539  2.44998842 0.0000004
#5:NI           -5:IR              -0.4318968 -0.9361641  0.07237041 0.1029240
#5:NI           -50:NI             -2.3776180 -2.8818853 -1.87335075 0.0000000
#Letters for plot
#50:IR  5:IR  50:NI 5:NI
#  a    b     c     b
group_by(plant.MG, interaction(SWHC,Soil_type)) %>% summarize(mean(Root_biomass))
#5.IR                                           0.106
#50.IR                                          0.37
#5.NI                                           0.072
#50.NI                                          0.73 

#Shoot dry biomass
summary(aov(log(plant.MG$Shoot_dry_biomass)~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#Df Sum Sq Mean Sq F value   Pr(>F)    
#plant.MG$SWHC                     1  7.911   7.911 228.162 3.59e-09 ***
#plant.MG$Soil_type                1  0.017   0.017   0.485    0.500    
#plant.MG$Replicates               4  0.085   0.021   0.612    0.662    
#plant.MG$SWHC:plant.MG$Soil_type  1  0.073   0.073   2.120    0.171    
#Residuals                        12  0.416   0.035            
TukeyHSD(aov(log(plant.MG$Shoot_dry_biomass)~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#diff        lwr        upr     p adj
#50:IR-5:IR   1.13663264  0.7869900  1.4862753 0.0000027
#5:NI-5:IR   -0.17921760 -0.5288602  0.1704250 0.4554329
#50:NI-5:IR   1.19988654  0.8502439  1.5495292 0.0000015
#5:NI-50:IR  -1.31585024 -1.6654929 -0.9662076 0.0000006
#50:NI-50:IR  0.06325391 -0.2863887  0.4128965 0.9482471
#50:NI-5:NI   1.37910415  1.0294615  1.7287468 0.0000003
#Letters for plot
#50:IR  5:IR  50:NI 5:NI
#  a    b     a     b
group_by(plant.MG, interaction(SWHC,Soil_type)) %>% summarize(mean(Shoot_dry_biomass))
#5.IR                                           0.062
#50.IR                                          0.196
#5.NI                                           0.052
#50.NI                                          0.206

#Shoot dry biomass
summary(aov(log(plant.MG$Root_dry_biomass)~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#Df Sum Sq Mean Sq F value   Pr(>F)    
#plant.MG$SWHC                     1 18.072  18.072  40.394 3.63e-05 ***
#  plant.MG$Soil_type                1  0.048   0.048   0.108    0.748    
#plant.MG$Replicates               4  0.350   0.087   0.196    0.936    
#plant.MG$SWHC:plant.MG$Soil_type  1  0.493   0.493   1.102    0.315    
#Residuals                        12  5.369   0.447            
TukeyHSD(aov(log(plant.MG$Root_dry_biomass)~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#diff        lwr        upr     p adj
#50:IR-5:IR   1.5871547  0.3312058  2.8431035 0.0127337
#5:NI-5:IR   -0.4122846 -1.6682334  0.8436642 0.7661002
#50:NI-5:IR   1.8029086  0.5469598  3.0588575 0.0052556
#5:NI-50:IR  -1.9994393 -3.2553881 -0.7434905 0.0023908
#50:NI-50:IR  0.2157539 -1.0401949  1.4717028 0.9551242
#50:NI-5:NI   2.2151932  0.9592444  3.4711421 0.0010346
#Letters for plot
#50:IR  5:IR  50:NI 5:NI
#  a    b     a     b
group_by(plant.MG, interaction(SWHC,Soil_type)) %>% summarize(mean(Root_dry_biomass))
#5.IR                                           0.008
#50.IR                                          0.518
#5.NI                                           0.054
#50.NI                                          0.54


#Root length
summary(aov(plant.MG$Root_length_cm~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#Df Sum Sq Mean Sq F value Pr(>F)  
#plant.MG$SWHC                     1  123.8  123.75   3.594 0.0823 .
#plant.MG$Soil_type                1    1.1    1.13   0.033 0.8594  
#plant.MG$Replicates               4  106.4   26.61   0.773 0.5634  
#plant.MG$SWHC:plant.MG$Soil_type  1   16.7   16.65   0.484 0.5000  
#Residuals                        12  413.2   34.43      
TukeyHSD(aov(plant.MG$Root_length_cm~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#diff        lwr       upr     p adj
#50:IR-5:IR   3.15  -7.868317 14.168317 0.8303768
#5:NI-5:IR   -2.30 -13.318317  8.718317 0.9237875
#50:NI-5:IR   4.50  -6.518317 15.518317 0.6310285
#5:NI-50:IR  -5.45 -16.468317  5.568317 0.4843660
#50:NI-50:IR  1.35  -9.668317 12.368317 0.9827477
#50:NI-5:NI   6.80  -4.218317 17.818317 0.3060347
#50:IR  5:IR  50:NI 5:NI
#  a    a     a     a
group_by(plant.MG, interaction(SWHC,Soil_type)) %>% summarize(mean(Root_length_cm))
#1 5.IR                                             13.9
#2 50.IR                                            17.0
#3 5.NI                                             11.6
#4 50.NI                                            18.4

#Leaf relative water content %
summary(aov(plant.MG$Leaf_relative_water_content_percentage~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#Df Sum Sq Mean Sq F value   Pr(>F)    
#plant.MG$SWHC                     1   6377    6377  62.670 4.19e-06 ***
#  plant.MG$Soil_type                1     20      20   0.197    0.665    
#plant.MG$Replicates               4    942     236   2.315    0.117    
#plant.MG$SWHC:plant.MG$Soil_type  1     19      19   0.191    0.670    
#Residuals                        12   1221     102    
TukeyHSD(aov(plant.MG$Leaf_relative_water_content_percentage~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#diff       lwr       upr     p adj
#50:IR-5:IR   33.744  14.80283  52.68517 0.0009505
#5:NI-5:IR    -3.972 -22.91317  14.96917 0.9228524
#50:NI-5:IR   33.710  14.76883  52.65117 0.0009586
#5:NI-50:IR  -37.716 -56.65717 -18.77483 0.0003586
#50:NI-50:IR  -0.034 -18.97517  18.90717 0.9999999
#50:NI-5:NI   37.682  18.74083  56.62317 0.0003615
#50:IR  5:IR  50:NI 5:NI
#  a    b     a     b
group_by(plant.MG, interaction(SWHC,Soil_type)) %>% summarize(mean(Leaf_relative_water_content_percentage))
#1 5.IR                                                                     63.4
#2 50.IR                                                                    97.1
#3 5.NI                                                                     59.4
#4 50.NI                                                                    97.1

#Leaf water content
summary(aov(plant.MG$Leaf_moisture~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#Df Sum Sq Mean Sq F value   Pr(>F)    
#plant.MG$SWHC                     1 116.89  116.89 153.110 3.43e-08 ***
#  plant.MG$Soil_type                1   1.54    1.54   2.017    0.181    
#plant.MG$Replicates               4   4.45    1.11   1.457    0.275    
#plant.MG$SWHC:plant.MG$Soil_type  1   0.32    0.32   0.419    0.530    
#Residuals           
TukeyHSD(aov(plant.MG$Leaf_moisture~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#diff        lwr       upr     p adj
#50:IR-5:IR   33.744   6.326688  61.16131 0.0151122
#5:NI-5:IR    -3.972 -31.389312  23.44531 0.9721571
#50:NI-5:IR   44.462  17.044688  71.87931 0.0020640
#5:NI-50:IR  -37.716 -65.133312 -10.29869 0.0071406
#50:NI-50:IR  10.718 -16.699312  38.13531 0.6612711
#50:NI-5:NI   48.434  21.016688  75.85131 0.0010209
#50:IR  5:IR  50:NI 5:NI
#  a    b     a     b
group_by(plant.MG, interaction(SWHC,Soil_type)) %>% summarize(mean(Leaf_moisture))
#1 5.IR                                            3.67
#2 50.IR                                           8.25
#3 5.NI                                            2.86
#4 50.NI                                           7.95

#Leaf dry matter content
summary(aov(plant.MG$leaf_dry_mater_content_mg.g~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#Df  Sum Sq  Mean Sq F value   Pr(>F)    
#plant.MG$SWHC                     1 0.01352 0.013520  34.228 7.83e-05 ***
#  plant.MG$Soil_type                1 0.00242 0.002420   6.127   0.0292 *  
#  plant.MG$Replicates               4 0.00098 0.000245   0.620   0.6567    
#plant.MG$SWHC:plant.MG$Soil_type  1 0.00032 0.000320   0.810   0.3858    
#Residuals                        12 0.00474 0.000395                                    
TukeyHSD(aov(plant.MG$leaf_dry_mater_content_mg.g~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#diff          lwr          upr     p adj
#50:IR-5:IR  -0.044 -0.081318508 -0.006681492 0.0197792
#5:NI-5:IR    0.030 -0.007318508  0.067318508 0.1328160
#50:NI-5:IR  -0.030 -0.067318508  0.007318508 0.1328160
#5:NI-50:IR   0.074  0.036681492  0.111318508 0.0003723
#50:NI-50:IR  0.014 -0.023318508  0.051318508 0.6883326
#50:NI-5:NI  -0.060 -0.097318508 -0.022681492 0.0022106
#50:IR  5:IR  50:NI 5:NI
#  a    bc     ac     b
group_by(plant.MG, interaction(SWHC,Soil_type)) %>% summarize(mean(leaf_dry_mater_content_mg.g))
#1 5.IR                                                         0.15 
#2 50.IR                                                        0.106
#3 5.NI                                                         0.18 
#4 50.NI                                                        0.12 

#Create object with Tukey letters for plots
tukey.plant <- data.frame(letters = c("a", "b", "c", "b", "a", "b", "a", "c"), Part=c(rep("Root_biomass",4), rep("Shoot_biomass",4)), x = c(rep(c(0.85, 1.25, 1.85, 2.25),2)), Soil_type=c("IR","IR", "NI", "NI","IR","IR", "NI", "NI") , SWHC= c("50","5","50","5","50","5","50","5"), y=c(0.5, 0.25, 1.0, 0.25, 1.80, 0.40, 1.5, 0.35))

#Plot
plant.MG$SWHC <- factor(plant.MG$SWHC, c("50", "5")) #Reorder manually
plant.MG.long <- gather(plant.MG,Part,Biomass,5:6) #make long
plant.box <- ggplot(data = plant.MG.long, aes(x=Soil_type, y=Biomass, fill=SWHC))+
                    geom_boxplot()+
                    facet_wrap(vars(Part), scales = "free_y", labeller = labeller(Part = c("Shoot_biomass" = "Shoot", "Root_biomass" = "Root") ))+
                    theme_bw()+
                    scale_x_discrete(name="Soil water stress history", labels = c("intermittent", "continuous"))+
                    #xlab() +
                    ylab("Fresh biomass")+
                    scale_fill_manual(name = "% SWHC", values = c("blue", "red"))+
                    geom_text(data = tukey.plant, aes(x=x, y=y, label = letters, hjust=1.5)) #Add letters from tukey
plant.box

#Percentage decrease - only significant terms in ANOVA
#Shoot fresh
(mean(plant.MG$Shoot_biomass[plant.MG$SWHC==50])-mean(plant.MG$Shoot_biomass[plant.MG$SWHC==5]))/mean(plant.MG$Shoot_biomass[plant.MG$SWHC==50])*100
#84.58015
(mean(plant.MG$Shoot_biomass[plant.MG$Soil_type=="IR"])-mean(plant.MG$Shoot_biomass[plant.MG$Soil_type=="NI"]))/mean(plant.MG$Shoot_biomass[plant.MG$Soil_type=="IR"])*100
#10.76345
#Root fresh
(mean(plant.MG$Root_biomass[plant.MG$SWHC==50])-mean(plant.MG$Root_biomass[plant.MG$SWHC==5]))/mean(plant.MG$Root_biomass[plant.MG$SWHC==50])*100
#83.81818
#Interaction
(mean(plant.MG$Root_biomass[1:5])-mean(plant.MG$Root_biomass[11:15]))/mean(plant.MG$Root_biomass[1:5])*100
#49.31507
#Shoot dry
(mean(plant.MG$Shoot_dry_biomass[plant.MG$SWHC==50])-mean(plant.MG$Shoot_dry_biomass[plant.MG$SWHC==5]))/mean(plant.MG$Shoot_dry_biomass[plant.MG$SWHC==50])*100
#71.64179
#Leaf moisture
(mean(plant.MG$Leaf_moisture[plant.MG$SWHC==50])-mean(plant.MG$Leaf_moisture[plant.MG$SWHC==5]))/mean(plant.MG$Leaf_moisture[plant.MG$SWHC==50])*100
#59.71347
#Leaf relative water content
(mean(plant.MG$Leaf_relative_water_content_percentage[plant.MG$SWHC==50])-mean(plant.MG$Leaf_relative_water_content_percentage[plant.MG$SWHC==5]))/mean(plant.MG$Leaf_relative_water_content_percentage[plant.MG$SWHC==50])*100
#36.77393
#Leaf dry matter content
(mean(plant.MG$leaf_dry_mater_content_mg.g[plant.MG$SWHC==5])-mean(plant.MG$leaf_dry_mater_content_mg.g[plant.MG$SWHC==50]))/mean(plant.MG$leaf_dry_mater_content_mg.g[plant.MG$SWHC==5])*100
#31.51515 (increased under dry soil)
(mean(plant.MG$leaf_dry_mater_content_mg.g[plant.MG$Soil_type=="NI"])-mean(plant.MG$leaf_dry_mater_content_mg.g[plant.MG$Soil_type=="IR"]))/mean(plant.MG$leaf_dry_mater_content_mg.g[plant.MG$Soil_type=="NI"])*100
#14.66667 (increased under continuous)