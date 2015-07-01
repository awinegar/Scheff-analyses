##################################################################################################
####LINEAR MIXED EFFECT MODELS FOR SCHEFFERVILLE CLADOCERANS                                     #
##################################################################################################

##In this script: Prep for LME models and use of lmer package 
#Beta-diversity focused analyses 
#Previous script: 
#R version: 3.1.2 (Pumpkin Helmet)

##Last update: June 15, 2015 
##Associated workspace: workspace_scheff_LMER.RData
##Associated markdown: 
##Github: 

##################################################################################################

##################################################################################################
####PACKAGES                                                                                     #
##################################################################################################

#Data prep
library(permute)
library(vegan)
library(boot)
library(rich)
library(Iso)
library(vegan)
library(plyr)
library(reshape)
library(reshape2)
library(arm)
library(AICcmodavg)

#Linear mixed effect models
library(lme4)

#Plotting
library(ggplot2)
library(gridExtra)

##################################################################################################

##################################################################################################
####RAREFYING CLADOCERAN SPECIES RICHNESS                                                        #
##################################################################################################

clad<- read.csv(file.choose()) #schefferville_lmerdata_MAR2015_clad.csv
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\Schefferville general\Analysis data\Mixed effect model data
#Count matrix of Cladocerans from schefferville_lmerdata_Mar2015.xls
#Counts rounded to nearest whole integer. 

row.names(clad) <- clad[,1]
clad <- clad[,-1]
Srar <- rarefy(clad, min(rowSums(clad)))
Srar

#Export single column for addition to: schefferville_lmerdata_MAR2015.xls

##################################################################################################

##################################################################################################

##################################################################################################
####PCA ON CLADOCERAN ASSEMBLAGE                                                                 #
##################################################################################################

clad.relabund<- read.csv(file.choose()) #schefferville_lmerdata_MAR2015_cladrelabun.csv
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\Schefferville general\Analysis data\Mixed effect model data
#Matrix of relative abundance of Cladocerans from schefferville_lmerdata_Mar2015.xls

row.names(clad.relabund)<- clad.relabund[,1]

clad.pca<- rda(decostand(clad.relabund[,4:41], method="hell")) #2nd column is lake, 3rd estimated year. 
clad.pc1.scr<- as.data.frame(scores(clad.pca, dis="sites", choices=1))

#Create nice looking plot of this. 
#Re-extract scores for PC1 and 2. 

clad.scr<- as.data.frame(scores(clad.pca, dis="sites", choices=1:2))
clad.scr<- as.data.frame(cbind(clad.relabund$Lake, clad.relabund$Est_year, clad.scr))
colnames(clad.scr) [1]<- 'Lake'
colnames(clad.scr) [2]<- 'Year'
  
clad.pca.plot<-ggplot()
clad.pca.plot<- clad.pca.plot + geom_vline(x=0,colour="grey50") 
clad.pca.plot<- clad.pca.plot+ geom_hline(y=0,colour="grey50") 
#clad.pca.plot<- clad.pca.plot + geom_point(data = clad.scr, aes(x = PC1, y = PC2, label=rownames(clad.scr), colour = Lake),pch = 15, size=2) 
clad.pca.plot<- clad.pca.plot + labs(x= "PC1 (32% var explained)", y= "PC2 (23% var exp)") + theme_bw()
clad.pca.plot<- clad.pca.plot + geom_text(data = clad.scr, aes(x = PC1, y = PC2, label=clad.scr$Year, colour = Lake), size=6) 
clad.pca.plot<- clad.pca.plot + scale_colour_brewer(type="qual", palette="Dark2")
clad.pca.plot<- clad.pca.plot + geom_path(data = clad.scr, aes(x = PC1, y= PC2, colour= clad.scr$Lake), size =0.2, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))
clad.pca.plot<- clad.pca.plot + theme(axis.text.x = element_text(colour="black", size=16))
clad.pca.plot<- clad.pca.plot + theme(axis.text.y = element_text(colour="black", size=16))
clad.pca.plot<- clad.pca.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
clad.pca.plot<- clad.pca.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))

#Export single column for addition to: schefferville_lmerdata_MAR2015.xls
write.csv(clad.pc1.scr, "clad.pc1.scr.csv")

##Use this data to make some quick stratigraphies. 

clad.long<- melt(clad.relabund, id.vars = c("Sample_ID_Master", "Lake", "Est_year"))

#Dauriat
clad.long.dar<- as.data.frame(subset(clad.long, Lake == "Dauriat", drop=T))

dar.plot<- ggplot(clad.long.dar, aes(x=Est_year, y=value)) + geom_point() + geom_path() + facet_wrap(~variable)
dar.plot2<- ggplot(clad.long.dar, aes(x=value, y=Est_year)) + geom_point() + geom_path() + facet_wrap(~variable)

#Knob
clad.long.kb<- as.data.frame(subset(clad.long, Lake == "Knob", drop=T))

kb.plot<- ggplot(clad.long.kb, aes(x=value, y=Est_year)) + geom_point() + geom_path() + facet_wrap(~variable)


##More polished stratigraphies using lumped Bosmina data (lumped due to 2 different taxonomists- Bosmina spp. only)

#Bosmina lumped data ("boslp" = Bosmina lumped)
clad.relabund.boslp<- read.csv(file.choose()) #schefferville_lmerdata_MAR2015_cladrelabun_BOSMINALP
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\Schefferville general\Analysis data\Mixed effect model data

clad.long.boslp<- melt(clad.relabund.boslp, id.vars = c("Sample_ID_Master", "Lake", "Est_year"))
#Re-order "variable" column so that species are alphabetical. 
clad.long.boslp<- clad.long.boslp[order(clad.long.boslp$variable),]#Not working, but want to have the 
#species in alphabetical order. 


#Dauriat
clad.longboslp.dar<- as.data.frame(subset(clad.long.boslp, Lake == "Dauriat", drop=T))

dummyYear<- as.factor(clad.longboslp.dar$Est_year) #was trying this for histograms, not working. 

#http://sape.inf.usi.ch/quick-reference/ggplot2/linetype (line types)

dar.plot<- ggplot(clad.longboslp.dar, aes(x=value, y=Est_year)) + geom_point() + geom_path(linetype="solid") + facet_wrap(~variable)
dar.plot<- dar.plot + labs(x="Relative abundance", y="Estimated year") +theme_bw()
dar.plot<- dar.plot + theme(axis.text.x = element_text(colour="black", size=14))
dar.plot<- dar.plot + theme(axis.text.y = element_text(colour="black", size=14))
dar.plot<- dar.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
dar.plot<- dar.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
dar.plot<- dar.plot + theme(strip.text = element_text(size=12))#facet titles
dar.plot<- dar.plot + annotate("rect", xmin=0.0, xmax=1.0, ymin=1954, ymax=1982, alpha=0.2, fill="darkorange3") #highlight mining time period
#Note: plot still doesn't have species in alphabetical order. 



#Knob
clad.longboslp.kb<- as.data.frame(subset(clad.long.boslp, Lake == "Knob", drop=T))

#Cut scale down?
kb.plot<- ggplot(clad.longboslp.kb, aes(x=value, y=Est_year)) + geom_point() + geom_path() + facet_wrap(~variable)
kb.plot<- kb.plot + labs(x="Relative abundance", y="Estimated year") +theme_bw()
kb.plot<- kb.plot + theme(axis.text.x = element_text(colour="black", size=14))
kb.plot<- kb.plot + theme(axis.text.y = element_text(colour="black", size=14))
kb.plot<- kb.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
kb.plot<- kb.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
kb.plot<- kb.plot + theme(strip.text = element_text(size=12))#facet titles
kb.plot<- kb.plot + annotate("rect", xmin=0.0, xmax=0.6, ymin=1954, ymax=1982, alpha=0.2, fill="darkorange3") #highlight mining time period
#Note: plot still doesn't have species in alphabetical order. 


##################################################################################################

##################################################################################################

##################################################################################################
####PCA ON METAL DATA                                                                            #
##################################################################################################

metals.ti<- read.csv(file.choose()) #schefferville_lmerdata_MAR2015_metalTi.csv
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\Schefferville general\Analysis data\Mixed effect model data
#Matrix of raw metal data: Titanium ratios

row.names(metals.ti)<- metals.ti[,1]

metals.ti2<- na.omit(metals.ti) #Remove row 28, "DL4_1" with NA values.

metals.pca<- rda(metals.ti2[,2:33])
metals.pc1.scr<- as.data.frame(scores(metals.pca, dis="sites", choices=1))

#Export single column for addition to: schefferville_lmerdata_MAR2015.xls
#Note: does not include a PC1 value for DL4_1
write.csv(metals.pc1.scr, "metals.pc1.scr.csv")

##################################################################################################


##################################################################################################

##################################################################################################
####LINEAR MIXED EFFECT MODELS                                                                   #
##################################################################################################

####Data for models

scheff.data<- read.csv(file.choose()) #schefferville_lmerdata_MAR2015_formodels.csv
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\Schefferville general\Analysis data\Mixed effect model data
#.csv of data from scheffervilleLmerdata_MAR2015.xls

scheff.data<- na.omit(scheff.data) #Remove DL4_1 observation

#Make Lake and Time_period factors
as.factor(scheff.data$Lake)
as.factor(scheff.data$Time_period)

####Data exploration/transformation
hist(scheff.data$Clad_S)
hist(scheff.data$Metal_EF)
hist(scheff.data$Clad_PC1)
hist(scheff.data$Metal_PC1)
hist(scheff.data$Cd_EF)
hist(scheff.data$Co_EF)
hist(scheff.data$Hg_EF)
#May need to consider transformations. 

####Clad_S through time
s.plot<- ggplot(scheff.data, aes(x=Est_year_numeric, y=Clad_S, colour = Lake)) + geom_point(size=4) + geom_path(size=0.75)
s.plot<- s.plot + labs(x= "Estimated year", y= "Rarefied taxa S") + theme_bw()
s.plot<- s.plot + scale_colour_brewer(type="qual", palette="Dark2")
s.plot<- s.plot + theme(axis.text.x = element_text(colour="black", size=16))
s.plot<- s.plot + theme(axis.text.y = element_text(colour="black", size=16))
s.plot<- s.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
s.plot<- s.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
s.plot<- s.plot + annotate("rect", xmin=1954, xmax=1982, ymin=5, ymax=11.4, alpha=0.2)


####Scaling of the data 
#Z-correct the variables so all on same scale. 

#Clad_S
Clad_S.Z<- (scheff.data$Clad_S-mean(scheff.data$Clad_S))/sd(scheff.data$Clad_S)

#Metal_EF
Metal_EF.Z<- (scheff.data$Metal_EF-mean(scheff.data$Metal_EF))/sd(scheff.data$Metal_EF)

#Clad_PC1
Clad_PC1.Z<- (scheff.data$Clad_PC1-mean(scheff.data$Clad_PC1))/sd(scheff.data$Clad_PC1)

#Metal_PC1
Metal_PC1.Z<- (scheff.data$Metal_PC1-mean(scheff.data$Metal_PC1))/sd(scheff.data$Metal_PC1)

#Cd_EF
Cd_EF.Z<- (scheff.data$Cd_EF-mean(scheff.data$Cd_EF))/sd(scheff.data$Cd_EF)

#Co_EF
Co_EF.Z<- (scheff.data$Co_EF-mean(scheff.data$Co_EF))/sd(scheff.data$Co_EF)

#Hg_EF
Hg_EF.Z<- (scheff.data$Hg_EF-mean(scheff.data$Hg_EF))/sd(scheff.data$Hg_EF)


####Look at linear models to determine if mixed models are necessary. 

##Clad_S and Clas_S.Z by Metal_EF and Metal_EF.Z##
lm.test1<- lm(Clad_S~Metal_EF, data=scheff.data) #Adj R2 = 0.4886
lm.test1.resid<- rstandard(lm.test1)
#Lake effect
plot(lm.test1.resid~ scheff.data$Lake, xlab = "Lake", ylab="Standardized residuals")
abline(0,0, lty=2) #show variation across lakes 
#Time period effect
plot(lm.test1.resid~ scheff.data$Time_period, xlab = "Time_period", ylab="Standardized residuals")
abline(0,0, lty=2) #show some variation across time periods 
  
lm.test2<- lm(Clad_S.Z~Metal_EF.Z, data=scheff.data) #Adj R2 = 0.4886
lm.test2.resid<- rstandard(lm.test2)
#Lake effect
plot(lm.test2.resid~ scheff.data$Lake, xlab = "Lake", ylab="Standardized residuals")
abline(0,0, lty=2) #show variation across lakes 
#Time period effect
plot(lm.test2.resid~ scheff.data$Time_period, xlab = "Time_period", ylab="Standardized residuals")
abline(0,0, lty=2) #show some variation across time periods 

##Clad_S and Clad_S.Z by Metal_PC1 and Metal_PC1.Z##
lm.test3<- lm(Clad_S~Metal_PC1, data=scheff.data) #Adj R2 = 0.139
lm.test3.resid<- rstandard(lm.test3)
#Lake effect
plot(lm.test3.resid~ scheff.data$Lake, xlab = "Lake", ylab="Standardized residuals")
abline(0,0, lty=2) #show variation across lakes 
#Time period effect
plot(lm.test3.resid~ scheff.data$Time_period, xlab = "Time_period", ylab="Standardized residuals")
abline(0,0, lty=2) #show some variation across time periods 
  
lm.test4<- lm(Clad_S.Z~Metal_PC1.Z, data=scheff.data) #Adj R2 = 0.139 
lm.test4.resid<- rstandard(lm.test4)
#Lake effect
plot(lm.test4.resid~ scheff.data$Lake, xlab = "Lake", ylab="Standardized residuals")
abline(0,0, lty=2) #show variation across lakes 
#Time period effect
plot(lm.test4.resid~ scheff.data$Time_period, xlab = "Time_period", ylab="Standardized residuals")
abline(0,0, lty=2) #show some variation across time periods 

##Clad_PC1 and Clad_PC1.Z by Metal_EF and Metal_EF.Z##
lm.test5<- lm(Clad_PC1~Metal_EF, data=scheff.data) #Adj R2 = 0.1033
lm.test5.resid<- rstandard(lm.test5)
#Lake effect
plot(lm.test5.resid~ scheff.data$Lake, xlab = "Lake", ylab="Standardized residuals")
abline(0,0, lty=2) #show variation across lakes 
#Time period effect
plot(lm.test5.resid~ scheff.data$Time_period, xlab = "Time_period", ylab="Standardized residuals")
abline(0,0, lty=2) #show some variation across time periods 
  
lm.test6<- lm(Clad_PC1.Z~Metal_EF.Z, data=scheff.data) #Adj R2 = 0.1033 
lm.test6.resid<- rstandard(lm.test6)
#Lake effect
plot(lm.test6.resid~ scheff.data$Lake, xlab = "Lake", ylab="Standardized residuals")
abline(0,0, lty=2) #show variation across lakes 
#Time period effect
plot(lm.test6.resid~ scheff.data$Time_period, xlab = "Time_period", ylab="Standardized residuals")
abline(0,0, lty=2) #show some variation across time periods 

##Clad_PC1 and Clad_PC1.Z by Metal_PC1 and Metal_PC1.Z## 
lm.test7<- lm(Clad_PC1~Metal_PC1, data=scheff.data) #Adj R2 = -0.03537
lm.test7.resid<- rstandard(lm.test7)
#Lake effect
plot(lm.test7.resid~ scheff.data$Lake, xlab = "Lake", ylab="Standardized residuals")
abline(0,0, lty=2) #show variation across lakes 
#Time period effect
plot(lm.test7.resid~ scheff.data$Time_period, xlab = "Time_period", ylab="Standardized residuals")
abline(0,0, lty=2) #show some variation across time periods 
  
lm.test8<- lm(Clad_PC1.Z~Metal_PC1.Z, data=scheff.data) #Adj R2 = -0.03537
lm.test8.resid<- rstandard(lm.test8)
#Lake effect
plot(lm.test8.resid~ scheff.data$Lake, xlab = "Lake", ylab="Standardized residuals")
abline(0,0, lty=2) #show variation across lakes 
#Time period effect
plot(lm.test8.resid~ scheff.data$Time_period, xlab = "Time_period", ylab="Standardized residuals")
abline(0,0, lty=2) #show some variation across time periods 

##Mixed models with random factors are appropriate for all the models. 


####Base Model 1: Clad_S~Metal_EF + Lake + Time
#Note REML = TRUE for these models 
#Note: did not need to use Z-standardized, all the same (see above)

#Base 1 model 1
#Description: Lake and time period are random factors, varying intercepts for both
b1m1<- lmer(Clad_S ~ Metal_EF + (1|Lake) + (1|Time_period), data=scheff.data, REML=TRUE)
b1m1
summary(b1m1)

#Base 1 model 2
#Description: Lake and time period as random factors, intercepts and slopes vary with respect to metal
b1m2<- lmer(Clad_S~Metal_EF + (1+Metal_EF|Lake) + (1+Metal_EF|Time_period), data=scheff.data, REML=TRUE)
b1m2
summary(b1m2)

#Base 1 model 3
#Description: Lake as a random factor, varying intercept 
b1m3<- lmer(Clad_S~Metal_EF + (1|Lake), data=scheff.data, REML=TRUE)
b1m3
summary(b1m3)

#Base 1 model 4
#Description: Time period as random factor, varying intercept 
b1m4<- lmer(Clad_S~Metal_EF + (1|Time_period), data=scheff.data, REML=TRUE)
b1m4
summary(b1m4)

#Base 1 model 5
#Description: Lake as a random factor, varying slope and intercept with respect to metal 
b1m5<- lmer(Clad_S~Metal_EF + (1+Metal_EF|Lake), data=scheff.data, REML=TRUE)
b1m5
summary(b1m5)

#Base 1 model 6
#Description: Time period as a random factor, varying slope and intercept with respect to metal 
b1m6<- lmer(Clad_S~Metal_EF + (1+Metal_EF|Time_period), data=scheff.data, REML=TRUE)
b1m6
summary(b1m6)

#Base 1 model 7
#Description: Lake and time period as random factors, varying intercept plus slope for lake, varying intercept only for time period 
b1m7<- lmer(Clad_S~Metal_EF + (1+Metal_EF|Lake) + (1|Time_period), data=scheff.data, REML=TRUE)
b1m7
summary(b1m7)

#Base 1 model 8
#Description: Lake and time period as random factors, varying intercept plus slope for time period, varying intercept only for lake
b1m8<- lmer(Clad_S~Metal_EF + (1|Lake) + (1+Metal_EF|Time_period), data=scheff.data, REML=TRUE)
b1m8
summary(b1m8)

#Base 1 null model
#Description: Linear model, cladoceran richness by metal EF, no random factors 
b1null<- lm(Clad_S~Metal_EF, data=scheff.data)
b1null
summary(b1null)

####Base Model 2: Clad_S~Metal_PCA + Lake + Time

#Base 2 model 1
#Description: Lake and time period as random factors, varying intercepts for both 
b2m1<- lmer(Clad_S ~ Metal_PC1 + (1|Lake) + (1|Time_period), data=scheff.data, REML=TRUE)
b2m1
summary(b2m1)

#Base 2 model 2
#Description: Lake and time period as random factors, varying intercepts and slopes for both with respect to metal
b2m2<- lmer(Clad_S~Metal_PC1 + (1+Metal_EF|Lake) + (1+Metal_EF|Time_period), data=scheff.data, REML=TRUE)
b2m2
summary(b2m2)

#Base 2 model 3
#Description: Lake as random factor, varying intercept 
b2m3<- lmer(Clad_S~Metal_PC1 + (1|Lake), data=scheff.data, REML=TRUE)
b2m3
summary(b2m3)

#Base 2 model 4
#Description: Time period as random factor, varying intercept 
b2m4<- lmer(Clad_S~Metal_PC1 + (1|Time_period), data=scheff.data, REML=TRUE)
b2m4
summary(b2m4)

#Base 2 model 5
#Description: Lake as random factor, varying intercept and slope with respect to metal 
b2m5<- lmer(Clad_S~Metal_PC1 + (1+Metal_PC1|Lake), data=scheff.data,REML=TRUE)
b2m5
summary(b2m5)

#Base 2 model 6
#Description: Time period as random factor, varying intercept and slope with respect to metal 
b2m6<- lmer(Clad_S~Metal_PC1 + (1+Metal_PC1|Time_period), data=scheff.data, REML=TRUE)
b2m6
summary(b2m6)

#Base 2 model 7
#Description: Lake and and time period as random factors, varying intercept and slope for lake, varying intercept for time period 
b2m7<- lmer(Clad_S~Metal_PC1 + (1+Metal_PC1|Lake) + (1|Time_period), data=scheff.data, REML=TRUE)
b2m7
summary(b2m7)

#Base 2 model 8
#Description: Lake and time period as random factors, varying intercept and slope for time period, varying intercept for lake
b2m8<- lmer(Clad_S~Metal_PC1 + (1|Lake) + (1+Metal_PC1|Time_period), data=scheff.data, REML=TRUE)
b2m8
summary(b2m8)

#Base 2 null
#Description: Linear model, cladoceran richness by metal PC1, no random factors 
b2null<- lm(Clad_S~Metal_PC1, data=scheff.data)
b2null
summary(b2null)

####Base Model 3: Clad_PCA~Metal_EF + Lake + Time

#Base 3 model 1
#Description: Lake and time period as random factors, varying intercepts 
b3m1<- lmer(Clad_PC1 ~ Metal_EF + (1|Lake) + (1|Time_period), data=scheff.data, REML=TRUE)
b3m1
summary(b3m1)

#Base 3 model 2
#Description: Lake and time period as random factors, varying intercepts and slopes with respect to metal 
b3m2<- lmer(Clad_PC1~Metal_EF + (1+Metal_EF|Lake) + (1+Metal_EF|Time_period), data=scheff.data, REML=TRUE)
b3m2
summary(b3m2)

#Base 3 model 3
#Description: Lake as random factor, varying intercept 
b3m3<- lmer(Clad_PC1~Metal_EF + (1|Lake), data=scheff.data, REML=TRUE)
b3m3
summary(b3m3)

#Base 3 model 4
#Description: Time period as random factor, varying intercept 
b3m4<- lmer(Clad_PC1~Metal_EF + (1|Time_period), data=scheff.data, REML=TRUE)
b3m4
summary(b3m4)

#Base 3 model 5
#Description: Lake as random factor, varying intercept and slope with respect to metal 
b3m5<- lmer(Clad_PC1~Metal_EF + (1+Metal_EF|Lake), data=scheff.data,REML=TRUE)
b3m5
summary(b3m5)

#Base 3 model 6
#Description: Time period as random factor, varying intercept and slope with respect to metal 
b3m6<- lmer(Clad_PC1~Metal_EF + (1+Metal_EF|Time_period), data=scheff.data, REML=TRUE)
b3m6
summary(b3m6)

#Base 3 model 7
#Description: Lake and time period as random factors, varying intercept and slope for lake, varying intercept for time period 
b3m7<- lmer(Clad_PC1~Metal_EF + (1+Metal_EF|Lake) + (1|Time_period), data=scheff.data, REML=TRUE)
b3m7
summary(b3m7)

#Base 3 model 8
#Description: Lake and time period as random factors, varying intercept for lake, varying intercept and slope for time period 
b3m8<- lmer(Clad_PC1~Metal_EF + (1|Lake) + (1+Metal_EF|Time_period), data=scheff.data, REML=TRUE)
b3m8
summary(b3m8)

#Base 3 null
#Description: Linear model, cladoceran PC1 by metal EF, no random factors 
b3null<- lm(Clad_PC1~Metal_EF, data=scheff.data)
b3null
summary(b3null)


####Base Model 4: Clad_PCA~Metal_PCA + Lake + Time

#Base 4 model 1
#Description: Lake and time period as random factors, varying intercepts 
b4m1<- lmer(Clad_PC1 ~ Metal_PC1 + (1|Lake) + (1|Time_period), data=scheff.data, REML=TRUE)
b4m1
summary(b4m1)

#Base 4 model 2
#Description: Lake and time period as random factors, varying intercepts and slopes with respect to metal 
b4m2<- lmer(Clad_PC1~Metal_PC1 + (1+Metal_PC1|Lake) + (1+Metal_PC1|Time_period), data=scheff.data, REML=TRUE)
b4m2
summary(b4m2)

#Base 4 model 4
#Description: Lake as random factor, varying intercept 
b4m3<- lmer(Clad_PC1~Metal_PC1 + (1|Lake), data=scheff.data, REML=TRUE)
b4m3
summary(b4m3)

#Base 4 model 4
#Description: Time period as random factor, varying intercept 
b4m4<- lmer(Clad_PC1~Metal_PC1 + (1|Time_period), data=scheff.data, REML=TRUE)
b4m4
summary(b4m4)

#Base 4 model 5
#Description: Lake as random factor, varying intercept and slope with respect to metal 
b4m5<- lmer(Clad_PC1~Metal_PC1 + (1+Metal_PC1|Lake), data=scheff.data,REML=TRUE)
b4m5
summary(b4m5)

#Base 4 model 6
#Description: Time period as random factor, varying intercept and slope with respect to metal 
b4m6<- lmer(Clad_PC1~Metal_PC1 + (1+Metal_PC1|Time_period), data=scheff.data, REML=TRUE)
b4m6
summary(b4m6)

#Base 4 model 7
#Description: Lake and time period as random factor, varying intercept and slope for lake, varying intercept for time period
b4m7<- lmer(Clad_PC1~Metal_EF + (1+Metal_PC1|Lake) + (1|Time_period), data=scheff.data, REML=TRUE)
b4m7
summary(b4m7)

#Base 4 model 8
#Description: Lake and time period as random factor, varying intercept for lake, varying intercept and slope for time period
b4m8<- lmer(Clad_PC1~Metal_PC1 + (1|Lake) + (1+Metal_PC1|Time_period), data=scheff.data, REML=TRUE)
b4m8
summary(b4m8)

#Base 4 null
#Description: Linear model, cladoceran PC1 by metal PC1, no random factors 
b4null<- lm(Clad_PC1~Metal_PC1, data=scheff.data)
b4null
summary(b4null)


####Model selection 
#Note REML = FALSE for model selection/comparison  

b1m1b<- lmer(Clad_S ~ Metal_EF + (1|Lake) + (1|Time_period), data=scheff.data, REML=FALSE)

b1m2b<- lmer(Clad_S~Metal_EF + (1+Metal_EF|Lake) + (1+Metal_EF|Time_period), data=scheff.data, REML=FALSE)

b1m3b<- lmer(Clad_S~Metal_EF + (1|Lake), data=scheff.data, REML=FALSE)

b1m4b<- lmer(Clad_S~Metal_EF + (1|Time_period), data=scheff.data, REML=FALSE)
 
b1m5b<- lmer(Clad_S~Metal_EF + (1+Metal_EF|Lake), data=scheff.data, REML=FALSE)

b1m6b<- lmer(Clad_S~Metal_EF + (1+Metal_EF|Time_period), data=scheff.data, REML=FALSE)

b1m7b<- lmer(Clad_S~Metal_EF + (1+Metal_EF|Lake) + (1|Time_period), data=scheff.data, REML=FALSE)

b1m8b<- lmer(Clad_S~Metal_EF + (1|Lake) + (1+Metal_EF|Time_period), data=scheff.data, REML=FALSE)

b1nullb<- lm(Clad_S~Metal_EF, data=scheff.data)

b2m1b<- lmer(Clad_S ~ Metal_PC1 + (1|Lake) + (1|Time_period), data=scheff.data, REML=FALSE)

b2m2b<- lmer(Clad_S~Metal_PC1 + (1+Metal_EF|Lake) + (1+Metal_EF|Time_period), data=scheff.data, REML=FALSE)

b2m3b<- lmer(Clad_S~Metal_PC1 + (1|Lake), data=scheff.data, REML=FALSE)

b2m4b<- lmer(Clad_S~Metal_PC1 + (1|Time_period), data=scheff.data, REML=FALSE)

b2m5b<- lmer(Clad_S~Metal_PC1 + (1+Metal_PC1|Lake), data=scheff.data,REML=FALSE)

b2m6b<- lmer(Clad_S~Metal_PC1 + (1+Metal_PC1|Time_period), data=scheff.data, REML=FALSE)

b2m7b<- lmer(Clad_S~Metal_PC1 + (1+Metal_PC1|Lake) + (1|Time_period), data=scheff.data, REML=FALSE)

b2m8b<- lmer(Clad_S~Metal_PC1 + (1|Lake) + (1+Metal_PC1|Time_period), data=scheff.data, REML=FALSE)

b2nullb<- lm(Clad_S~Metal_PC1, data=scheff.data)

b3m1b<- lmer(Clad_PC1 ~ Metal_EF + (1|Lake) + (1|Time_period), data=scheff.data, REML=FALSE)

b3m2b<- lmer(Clad_PC1~Metal_EF + (1+Metal_EF|Lake) + (1+Metal_EF|Time_period), data=scheff.data, REML=FALSE)

b3m3b<- lmer(Clad_PC1~Metal_EF + (1|Lake), data=scheff.data, REML=FALSE)

b3m4b<- lmer(Clad_PC1~Metal_EF + (1|Time_period), data=scheff.data, REML=FALSE)

b3m5b<- lmer(Clad_PC1~Metal_EF + (1+Metal_EF|Lake), data=scheff.data,REML=FALSE)

b3m6b<- lmer(Clad_PC1~Metal_EF + (1+Metal_EF|Time_period), data=scheff.data, REML=FALSE)

b3m7b<- lmer(Clad_PC1~Metal_EF + (1+Metal_EF|Lake) + (1|Time_period), data=scheff.data, REML=FALSE)

b3m8b<- lmer(Clad_PC1~Metal_EF + (1|Lake) + (1+Metal_EF|Time_period), data=scheff.data, REML=FALSE)

b3nullb<- lm(Clad_PC1~Metal_EF, data=scheff.data)

b4m1b<- lmer(Clad_PC1 ~ Metal_PC1 + (1|Lake) + (1|Time_period), data=scheff.data, REML=FALSE)

b4m2b<- lmer(Clad_PC1~Metal_PC1 + (1+Metal_PC1|Lake) + (1+Metal_PC1|Time_period), data=scheff.data, REML=FALSE)

b4m3b<- lmer(Clad_PC1~Metal_PC1 + (1|Lake), data=scheff.data, REML=FALSE)

b4m4b<- lmer(Clad_PC1~Metal_PC1 + (1|Time_period), data=scheff.data, REML=FALSE)

b4m5b<- lmer(Clad_PC1~Metal_PC1 + (1+Metal_PC1|Lake), data=scheff.data,REML=FALSE)

b4m6b<- lmer(Clad_PC1~Metal_PC1 + (1+Metal_PC1|Time_period), data=scheff.data, REML=FALSE)

b4m7b<- lmer(Clad_PC1~Metal_EF + (1+Metal_PC1|Lake) + (1|Time_period), data=scheff.data, REML=FALSE)

b4m8b<- lmer(Clad_PC1~Metal_PC1 + (1|Lake) + (1+Metal_PC1|Time_period), data=scheff.data, REML=FALSE)

b4nullb<- lm(Clad_PC1~Metal_PC1, data=scheff.data)

##Model evaluation
AICc<-c(AICc(b1m1b), AICc(b1m2b), AICc(b1m3b), AICc(b1m4b), AICc(b1m5b), AICc(b1m6b), AICc(b1m7b), AICc(b1m8b), AICc(b1nullb), AICc(b2m1b), AICc(b2m2b), AICc(b2m3b), AICc(b2m4b), AICc(b2m5b), AICc(b2m6b), AICc(b2m7b), AICc(b2m8b), AICc(b2nullb), AICc(b3m1b), AICc(b3m2b), AICc(b3m3b), AICc(b3m4b), AICc(b3m5b), AICc(b3m6b), AICc(b3m7b), AICc(b3m8b), AICc(b3nullb), AICc(b4m1b), AICc(b4m2b), AICc(b4m3b), AICc(b4m4b), AICc(b4m5b), AICc(b4m6b), AICc(b4m7b), AICc(b4m8b), AICc(b4nullb))
# Put values into one table for easy comparision
Model<-c("b1m1", "b1m2", "b1m3", "b1m4", "b1m5", "b1m6", "b1m7", "b1m8", "b1null", "b2m1", "b2m2", "b2m3", "b2m4", "b2m5", "b2m6", "b2m7", "b2m8", "b2null", "b3m1", "b3m2", "b3m3", "b3m4", "b3m5", "b3m6", "b3m7", "b3m8", "b3null", "b4m1", "b4m2", "b4m3", "b4m4", "b4m5", "b4m6", "b4m7", "b4m8", "b4null")
AICtable<-data.frame(Model=Model, AICc=AICc)
AICtable
#Doesn't seem like these are correct- I think accurate for b1 and b2 models, issues with b3 and b4 (but not as predictive)

##
##Use b1m3 (b1m3b)- Clad_S ~ Metal_EF with lake as random factor and eventually compare it to Cd, Co, and Hg models below
##

##Model evaluation of b1m3
#Look at independence, plot fitted values versus residuals 
E1<- resid(b1m3b)
F1<- fitted(b1m3b)
plot(x = F1, 
     y = E1, 
     xlab = "Fitted Values",
     ylab = "Normalized residuals")
abline(h = 0, lty = 2)
#Even spaced residuals- model a good fit for the data.

#Check assumption of homogeneity- residuals vs. covariate 
#Metal_EF
plot(x= scheff.data$Metal_EF, 
     y = E1, 
     xlab = "Metal_EF", 
     ylab = "Normalized residuals")
abline(h = 0, lty = 2)

#Lake
boxplot(E1 ~ Lake,   
        ylab = "Normalized residuals",
        data = scheff.data, xlab = "Lake")
abline(h = 0, lty = 2)
##But didn't get homogeneity. 

#Check normality of residuals 
hist(E1)
#good 

#Determine if slope and therefore effect of Metal_EF on Clad_S is significantly different from zero- calculate CI
#for the slope parameter. If CI overlaps with zero, slope is not significantly different from 0 at 0.05 level. 
#CI = Standard Error of the estimate x 1.96 plus or minus the parameter estimate.
#Slope from b1m3b = -0.055

upperCI<- -0.055 + 0.01148*1.96 #-0.032
lowerCI<- -0.055 - 0.01148*1.96 #-0.0775
#Doesn't overlap zero so significantly different from 0. 


#Visualization of b1m3b.
#Put lake coefficents in the dataframe. 
lake.coef<- as.data.frame(coef(b1m3b)$Lake)
colnames(lake.coef)<- c("Intercept", "Slope")

# Plot 1 - All Data
# Make a plot that includes all the data
b1m3.plot1 <- ggplot(aes(x=Metal_EF, y=Clad_S), data=scheff.data) + geom_point() + xlab("Metal EF") + ylab("Cladoceran rar S")
b1m3.plot1<- b1m3.plot1 + geom_abline(intercept=9.3548, slope=-0.055) #Overall
b1m3.plot1<- b1m3.plot1 + geom_abline(intercept = lake.coef[1,1], slope=lake.coef[1,2], colour="coral2")  #Dauriat
b1m3.plot1<- b1m3.plot1 + geom_abline(intercept = lake.coef[2,1], slope = lake.coef[2,2], colour="red")  #Denault
b1m3.plot1<- b1m3.plot1 + geom_abline(intercept = lake.coef[3,1], slope = lake.coef[3,2], colour="blue")  #Dolly
b1m3.plot1<- b1m3.plot1 + geom_abline(intercept = lake.coef[4,1], slope = lake.coef[4,2], colour="darkgoldenrod") #Knob
b1m3.plot1<- b1m3.plot1 + theme_bw()
b1m3.plot1<- b1m3.plot1 + theme(axis.text.x = element_text(colour="black", size=16))
b1m3.plot1<- b1m3.plot1 + theme(axis.text.y = element_text(colour="black", size=16))
b1m3.plot1<- b1m3.plot1 + theme(axis.title.x = element_text(size = rel(2), angle=00))
b1m3.plot1<- b1m3.plot1 + theme(axis.title.y = element_text(size = rel(2), angle=90))




####Linear models looking at Cd, Co and Hg. 
cd.null<-lm(Clad_S~Cd_EF, data=scheff.data)
cd.null
summary(cd.null)

co.null<-lm(Clad_S~Co_EF, data=scheff.data)
co.null
summary(co.null)
  
hg.null<-lm(Clad_S~Hg_EF, data=scheff.data) 
hg.null
summary(hg.null)

##################################################################################################