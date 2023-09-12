##Anova of plant data

#Load data
plant <- readRDS(here("data", "intermediate", "plant.RDS"))
lwc <- readRDS(here("data", "intermediate", "lwc.RDS"))

#Merge
sum((plant$SWHC==lwc$Treatment) & 
      (plant$Cultivar==lwc$cultivar) & 
      (plant$Soil_type==lwc$Soil_type) &
      (plant$Replicates==lwc$Replicates)
      )#160
plant$LWC <- lwc$Leaf_moisture

#keep only relevant for this study
plant.MG <- plant[plant$SWHC %in% c(5,50) & plant$Cultivar=="ACNass",]
plant.MG$SWHC <- as.character(plant.MG$SWHC)
plant.MG$Replicates <- as.character(plant.MG$Replicates)

#Test normality and heteroscedascticity
shapiro.test(plant.MG$Shoot_biomass) #P=0.002085
shapiro.test(plant.MG$Root_biomass) #P=0.00407
shapiro.test(plant.MG$LWC) #P=0.02017
bartlett.test(plant.MG$Shoot_biomass~plant.MG$Soil_type) #P=0.8317
bartlett.test(plant.MG$Root_biomass~plant.MG$Soil_type) #P=0.01504
bartlett.test(plant.MG$LWC~plant.MG$Soil_type) #P=0.7367
bartlett.test(plant.MG$Shoot_biomass~plant.MG$SWHC) #P=0.0001883
bartlett.test(plant.MG$Root_biomass~plant.MG$SWHC) #P=4.641e-06
bartlett.test(plant.MG$LWC~plant.MG$SWHC) #P=0.2002

#Log transform
shapiro.test(log(plant.MG$Shoot_biomass)) #P=0.003372
shapiro.test(log(plant.MG$Root_biomass)) #P=0.06382
shapiro.test(log(plant.MG$LWC)) #P=0.008259
bartlett.test(log(plant.MG$Shoot_biomass)~plant.MG$Soil_type) #P=0.5672
bartlett.test(log(plant.MG$Root_biomass)~plant.MG$Soil_type) #P=0.08064
bartlett.test(log(plant.MG$LWC)~plant.MG$Soil_type)#P=0.3013
bartlett.test(log(plant.MG$Shoot_biomass)~plant.MG$SWHC) #P=0.09125
bartlett.test(log(plant.MG$Root_biomass)~plant.MG$SWHC) #P=0.9958
bartlett.test(log(plant.MG$LWC)~plant.MG$SWHC) #P=9.252e-05

#Anova
#Shoot
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

#Root
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

#LWC
summary(aov(plant.MG$LWC~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#Df Sum Sq Mean Sq F value   Pr(>F)    
#plant.MG$SWHC                     1 116.90  116.90 153.173 3.43e-08 ***
#  plant.MG$Soil_type                1   1.54    1.54   2.016    0.181    
#plant.MG$Replicates               4   4.46    1.11   1.460    0.275    
#plant.MG$SWHC:plant.MG$Soil_type  1   0.32    0.32   0.420    0.529    
#Residuals                        12   9.16    0.76                     
TukeyHSD(aov(plant.MG$LWC~plant.MG$SWHC*plant.MG$Soil_type+plant.MG$Replicates))
#diff       lwr        upr     p adj
#50:IR-5:IR   4.582000  2.941617  6.2223835 0.0000135
#5:NI-5:IR   -0.808000 -2.448383  0.8323835 0.4877420
#50:NI-5:IR   4.280638  2.640255  5.9210215 0.0000269
#5:NI-50:IR  -5.390000 -7.030383 -3.7496165 0.0000025
#50:NI-50:IR -0.301362 -1.941745  1.3390215 0.9460201
#50:NI-5:NI   5.088638  3.448255  6.7290215 0.0000045
#Letters for plot
#50:IR  5:IR  50:NI 5:NI
#  a    b     a     b

#Create object with Tukey letters
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
#Shoot
(mean(plant.MG$Shoot_biomass[plant.MG$SWHC==50])-mean(plant.MG$Shoot_biomass[plant.MG$SWHC==5]))/mean(plant.MG$Shoot_biomass[plant.MG$SWHC==50])*100
#84.58015
(mean(plant.MG$Shoot_biomass[plant.MG$Soil_type=="IR"])-mean(plant.MG$Shoot_biomass[plant.MG$Soil_type=="NI"]))/mean(plant.MG$Shoot_biomass[plant.MG$Soil_type=="IR"])*100
#10.76345
#Root
(mean(plant.MG$Root_biomass[plant.MG$SWHC==50])-mean(plant.MG$Root_biomass[plant.MG$SWHC==5]))/mean(plant.MG$Root_biomass[plant.MG$SWHC==50])*100
#83.81818
#Interaction
(mean(plant.MG$Root_biomass[1:5])-mean(plant.MG$Root_biomass[11:15]))/mean(plant.MG$Root_biomass[1:5])*100
#49.31507
#LWC
(mean(plant.MG$LWC[plant.MG$SWHC==50])-mean(plant.MG$LWC[plant.MG$SWHC==5]))/mean(plant.MG$LWC[plant.MG$SWHC==50])*100
#59.71506