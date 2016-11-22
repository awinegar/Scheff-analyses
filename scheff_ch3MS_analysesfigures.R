##################################################################################################
####SCHEFFERVILLE: Ch. 3 manuscript figures                                                      #
##################################################################################################

##In this script: Final analyses and figures for 3rd chapter (general Schefferville manuscript)
#Previous script: Various R scripts from R studio scripts/Schefferville summer 2013
#R version: 3.1.2 (Pumpkin Helmet)

##Last update: Jan.14 , 2016
##Associated workspace: 
##Associated markdown: 
##Associated .txt of R script: 
##Github: 

##################################################################################################

##################################################################################################
####PACKAGES                                                                                     #
##################################################################################################

library(permute)
library(vegan)
library(boot)
library(rich)
library(Iso)
library(vegan)
library(plyr)
library(reshape)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(arm)
library(AICcmodavg)
library(lme4)
library(ggmap)
library(mapproj)
library(rioja) #using C2 and stratigraphic plots in R

##################################################################################################

##################################################################################################
####WORKING DIRECTORY FOR FIGURES                                                                #
##################################################################################################
##Move any necessary data to this directory (remember to update if change any files in alternative paths)

setwd("C:/Users/Winegardner/Documents/MCGILL/PhD chapters and projects/Schefferville general/Ch 3 manuscript figures")


##################################################################################################
####FIG 1ab - MAP OF SCHEFFERVILLE AREA                                                          #
##################################################################################################
##Make a map showing where Schefferville is in Quebec and the location of the 4 lakes in the Schefferville
#region

#Fig 1A - Location of Schefferville 
scheffcoor<- read.csv("scheff_coor.csv") #scheff_coor.csv in set working directory 

reg<- get_map(location = "Schefferville Quebec", zoom = 5, source = 'google')

regMAP1<- ggmap(reg) + geom_point(aes(x=Lon, y=Lat), data = scheffcoor, colour = "black", size=7,alpha=0.5) + labs(x="Longitude (DD)",y="Latitude (DD)")
regMAP1<- regMAP1 + annotate("text", x=-66.8333, y=55.5, label="Schefferville")
regMAP1<- regMAP1 + annotate("text", x=-72, y=48, label="Qu�bec City")
regMAP1<- regMAP1 + theme(axis.text.x = element_text(colour="black",size=10))
regMAP1<- regMAP1 + theme(axis.text.y = element_text(colour="black",size=10))
regMAP1<- regMAP1 + theme(axis.title.x = element_text(size = rel(1.5), angle=00))
regMAP1<- regMAP1 + theme(axis.title.y = element_text(size = rel(1.5), angle=90))
regMAP1<- regMAP1 + annotate("text", x=-78, y=61.2, label="(a)", size=10)

#Fig 1B - Four lakes across the Schefferville landscape 
coor <- read.csv("lakes_coor.csv") #lakes_coor.csv in set working directory 
coor<- coor[-4,]#remove Dolly and Denault lakes 
coor<- coor[-3,]

coor1<-get_map(location = "Schefferville Quebec", zoom = 12, source = 'google')

coorMAP1<-ggmap(coor1)+
  geom_point(aes(x=Lon,y=Lat), data = coor,colour="red",size=7,alpha=0.5) + labs(x="Longitude (DD)",y="Latitude (DD)")
coorMAP1<- coorMAP1 + annotate("text", x=-66.82, y=54.810, label="Dauriat")
coorMAP1<- coorMAP1 + annotate("text", x=-66.80, y=54.786, label="Knob")
#coorMAP1<- coorMAP1 + annotate("text", x=-66.75, y=54.810, label="Dolly")
#coorMAP1<- coorMAP1 + annotate("text", x=-66.87, y=54.84, label="Denault")
coorMAP1<- coorMAP1 + theme(axis.text.x = element_text(colour="black",size=10))
coorMAP1<- coorMAP1 + theme(axis.text.y = element_text(colour="black",size=10))
coorMAP1<- coorMAP1 + theme(axis.title.x = element_text(size = rel(1.5), angle=00))
coorMAP1<- coorMAP1 + theme(axis.title.y = element_text(size = rel(1.5), angle=90))
coorMAP1<- coorMAP1 + annotate("text", x=-66.9, y=54.88, label="(b)", size=10)

#Putting Figure 1a and 1b together
fig1<- grid.arrange(regMAP1, coorMAP1, nrow=1)

##################################################################################################



##################################################################################################
####FIG 2 - METAL EF, DITP AND CORE PICTURES FIGURE                                              #
##################################################################################################

ef<- read.csv("euk_metalEF_full.csv") #includes EFs for a suite of metals as opposed to just the composite. 
#Subset so only including DAR and KB. 
ef<- as.data.frame(subset(ef, Lake == "Dauriat " | Lake == "Knob", drop=T))


error.ef<- aes(xmax = ef$Est_yearCIC + ef$Est_year_error, xmin = ef$Est_yearCIC - ef$Est_year_error)

#Metal EF plot (a) 
ef.plot<- ggplot(ef, aes(x=Est_yearCIC, y=Total_EF, colour = Lake)) + geom_path() + coord_flip() 
ef.plot<- ef.plot + geom_point(aes(x=Est_yearCIC, y=Total_EF, shape=factor(Year_type)), size=4)
ef.plot<- ef.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3"))  
ef.plot<- ef.plot + scale_shape_manual(values = c(16,1))
ef.plot<- ef.plot + labs(x = "Estimated year", y = "Metal EF") + theme_bw() + theme(legend.position="none")
ef.plot<- ef.plot + geom_errorbarh(error.ef, width=1.75, colour="black") 
ef.plot<- ef.plot + theme(axis.text.x = element_text(colour="black", size=16))
ef.plot<- ef.plot + theme(axis.text.y = element_text(colour="black", size=16))
ef.plot<- ef.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
ef.plot<- ef.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 
ef.plot<- ef.plot + annotate("text", x=2000, y=0, label="(a)", size=10)
ef.plot<- ef.plot + annotate("text", x=1880, y=40, label="Dauriat", size=8)
ef.plot<- ef.plot + annotate("text", x=1750, y=12.5, label="Knob", size=8)

#PCA plot (b) 
pcmet.plot<- ggplot(ef, aes(x=Est_yearCIC, y=Metal_PC1, colour = Lake)) + geom_path() + coord_flip()
pcmet.plot<- pcmet.plot + geom_point(aes(x=Est_yearCIC, y=Metal_PC1, shape=factor(Year_type)), size=4)
pcmet.plot<- pcmet.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3"))  
pcmet.plot<- pcmet.plot + scale_shape_manual(values = c(16,1))
pcmet.plot<- pcmet.plot + labs(x = "Estimated year", y = "Metal PC1") + theme_bw() + theme(legend.position="none")
pcmet.plot<- pcmet.plot + geom_errorbarh(error.ef, width=1.75, colour="black") 
pcmet.plot<- pcmet.plot + theme(axis.text.x = element_text(colour="black", size=16))
pcmet.plot<- pcmet.plot + theme(axis.text.y = element_text(colour="black", size=16))
pcmet.plot<- pcmet.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
pcmet.plot<- pcmet.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
pcmet.plot<- pcmet.plot + annotate("text", x=2000, y=-3.5, label="(b)", size=10)

#DITP plot (c) 
#DI-TP for Dauriat (from Laperriere et al. 2008)
ditp<- read.csv("DITP_DAR_Lap.csv")

ditp.plot<- ggplot(ditp, aes(x=Year, y=DITP.ugL)) + geom_point(size=4, colour="black") + geom_path(colour="black") + coord_flip()
ditp.plot<- ditp.plot + theme_bw() 
ditp.plot<- ditp.plot + xlab("Estimated year")
ditp.plot<- ditp.plot + ylab(expression(paste(DITP (mu*L))))
ditp.plot<- ditp.plot + theme(axis.text.x = element_text(colour="black", size=16))
ditp.plot<- ditp.plot + theme(axis.text.y = element_text(colour="black", size=16))
ditp.plot<- ditp.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
ditp.plot<- ditp.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
ditp.plot<- ditp.plot + annotate("text", x=2010, y=20, label="(c)", size=10)

#Lithology plot (d) 
#Use Lithology column that is now in euk_metalEF_full.csv 
as.factor(ef$Lithology)
lith.plot<- ggplot(ef, aes(x=Lake, y=Est_yearCIC, colour=Lithology)) + geom_point(size=4)
#can add shape for black and white plots, shape=Lithology
lith.plot<- lith.plot + theme_bw()
lith.plot<- lith.plot + scale_colour_manual(values = c("Black" = "black", "Dark_brown" = "saddlebrown", "Light_brown" = "peru", "Orange_brown" = "orangered1"))
lith.plot<- lith.plot + labs(x="Lake", y="Estimated year")
lith.plot<- lith.plot + theme(axis.text.x = element_text(colour="black", size=12, angle= 45))
lith.plot<- lith.plot + theme(axis.text.y = element_text(colour="black", size=16))
lith.plot<- lith.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
lith.plot<- lith.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
lith.plot<- lith.plot + annotate("text", x=0.6, y=2010, label="(d)", size=10)


#Putting panel together
fig2<- grid.arrange(ef.plot, pcmet.plot, ditp.plot, lith.plot, nrow=2) #Need to export once (c) and (d) done. 

#Experimenting with lithology plot- want to do a stacked histogram. 
#lith2.plot<- ggplot(ef, aes(x=Est_year, fill=Lithology), stat="identity") + geom_bar() + facet_wrap(~Lake) #getting closer 


#lithology<- as.factor(ef$Lithology)
#lith3.plot<- ggplot(ef, aes(x=Lake, y=Depth_cm, fill=lithology)) + geom_bar(stat="identity") 
#+ coord_flip()
#title.grob<- textGrob(
 # label = "(d)",
  #x = unit(0, "lines"),
  #y = unit(0, "lines"), 
  #hjust = -1, vjust = 0, 
  #gp = gpar(fontsize = 28))

#l.p<- arrangeGrob(lith.plot, main = title.grob)




##################################################################################################

##################################################################################################
####FIG 3/4 - CLADOCERAN STRATIGRAPHIES FOR DAR AND KB (full profiles)   + PCAs                  #
##################################################################################################
#######Use rioja package for C2 plots. Only for DAR and KB since have cladoceran profiles##### 

#Count and relative abundance data
#clad<- read.csv("clad_DARKB.csv")
#clad.rel<- read.csv("clad_DARKB_relabund.csv")

#Data with only species that had at least 5% relative abundance in at least one interval. 
#Species richness stil takes account ALL species. 
#clad.dom_relabund<- read.csv("clad_DARKB_relabund_greater5per.csv") 

#####Stratigraphies produced in R - refer to C2 plots instead#####
#Data to use for this figure
clad.dom_percent<- read.csv("clad_DARKB_relabund_percent.csv") #relative abundance from _greater5per multiplied by 100. 
#Includes unsplit Bosmina as well
#Species are in order of littoral then both littoral and pelagic/littoral (if plot in order, can then add a label in illustrator)

#Littoral species
#Acroperus.harpae
#Alona.affinis
#Alona.circumfimbria_guttata
#Alona.quandrangularis
#Alona.spp
#Alonella.nana
#Chydorus.cf.sphaericus
#Chydorus.gibbus
#Eurycercus.spp
#Paralona.piger

#Both pelagic/littoral species
#Bosmina.spp.total
#Bosmina.spp.und (undefined)
#Bosmina.longirostris
#Eubosmina.longispina
#Daphnia.longispina 


#Subset into DAR and KB
clad.dom.DAR<- as.data.frame(subset(clad.dom_percent, Lake == "Dauriat ", drop=T))
clad.dom.KB<- as.data.frame(subset(clad.dom_percent, Lake == "Knob", drop=T))

#Dar plots - FIGURE 3
names.sp<- c("A. harpae", "A. affinis", "A. circumfimbria", "A. qaudrangularis", "Alona spp.", 
          "A. nana", "C. sphaericus", "C. gibbus", "Eurycercus spp.", "P. piger",
          "Bosmina.spp.total", "Bosmina spp.und", "B. longirostris", "E. longispina", "D. longispina") #short forms of names--> write out in legend

names.S<- c("Cladoceran richness", "Cladoceran PC1")

Years.DAR<- clad.dom.DAR$Est_year
Depths.DAR<- clad.dom.DAR$Depth_cm

#Plot for richness so that uses different symbols and extended axis. 
dar.strat.S<-  strat.plot(clad.dom.DAR[,7:8], yvar=Years.DAR, plot.bar=FALSE, plot.line=TRUE,
                          col.bar="darkgrey", lwd.bar=8, ylabel = "Estimated year", cex.ylabel=1, cex.xlabel=1,
                          srt.xlabel=0, x.names=names.S)
#Export as .eps with dimensions w = 950, h=546 for Illustrator 
#Clean up all fonts- change to Arial/size = regular 14
#Add Depth axis 
#Need depths for years 2000, 1950, 1900 and 1850- find by isolating x from cross multiplication of nearest years
#Year 2000 (x = 2000 *0.75/2006.907) = ~0.75
#Year 1950 (x = 1950*17.75/1960.458) = ~17.7
#Year 1900 (x = 1900*23.5/1929.691) = ~23
#Year 1850 (x = 1850*26.5/1896.23) = ~26
#Add (a) - size 16 
                  
#Plot for species profiles. 
dar.strat.sp<- strat.plot(clad.dom.DAR[,9:23], yvar=Years.DAR, scale.percent=TRUE, plot.bar=TRUE, plot.line=FALSE,
                       col.bar="darkgrey", lwd.bar=8, ylabel = "Estimated year", cex.ylabel=1, cex.xlabel=1,
                       srt.xlabel=80, x.names=names.sp, xLeft=0.15) #xLeft=0.15 creates room for a 2ndary y-axis. 
#species are 9:23
#Export as .eps with dimensions w = 1050, h=850 for Illustrator 
#Clean up all fonts- change to Arial/size = regular 14
#Add Depth axis- same as (a) 
#Add (b) - size 16

######
#Trying to figure out how to add a secondary axis for depths (cm). _ For now, add in Illustrator 
#https://craticula.wordpress.com/2013/03/07/test-post/
#dYears<- c(2006, 1998, 1986, 1977, 1969, 1966, 1960, 1941, 1929, 1921, 1921, 1896, 1850, 1840, 1830, 1820) #years rounded down- verify double 1921
#dDepths<- approx(Years, Depths, xout=dYears, rule=2)$y #not exactly sure what is going on here? 

#Add littorial/pelagic designations - this can be done in illustrator now that plotted in correct order. 

#Split Bosmina - not sure how to do this in C2- can make it happen in ggplot- for now these are all represented separately in the plot. 
######


#Knob plots - FIGURE 4
Years.KB<- clad.dom.KB$Est_year
Depths.KB<- clad.dom.KB$Depth_cm

kb.strat.S<- strat.plot(clad.dom.KB[,7:8], yvar=Years.KB, plot.bar=FALSE, plot.line=TRUE,
                        col.bar="darkgrey", lwd.bar=8, ylabel = "Estimated year", cex.ylabel=1, cex.xlabel=1,
                        srt.xlabel=0, x.names=names.S)
#Export as .eps with dimensions w = 950, h=546 for Illustrator 
#Clean up all fonts- change to Arial/size = regular 14
#Add Depth axis 
#Need depths for years 2000, 1950, 1900, 1850, 1800, 1750, 1700- find by isolating x from cross multiplication of nearest years
#Year 2000 (x = 2000 *0.875/2003.007) = ~0.87
#Year 1950 (x = 1950*5.625/1954.2) = ~5.6
#Year 1900 (x = 1900*10.5/1904.09) = ~10.5
#Year 1850 (x = 1850*12.5/1883.531) = ~12.3
#Year 1800 (x = 1800*21.125/1794.9) = ~21
#Year 1750 (x = 1750*27.375/1730.65) = ~28
#Year 1700 (x = 1700*33.675/1665.9) = ~34
#Add (a) - size 16 



kb.strat.sp<- strat.plot(clad.dom.KB[,9:23], yvar=Years.KB, scale.percent=TRUE, plot.bar=TRUE, plot.line=FALSE,
                         col.bar="darkgrey", lwd.bar=8, ylabel = "Estimated year", cex.ylabel=1, cex.xlabel=1,
                         srt.xlabel=80, x.names=names.sp, xLeft=0.15) #xLeft=0.15 creates room for a 2ndary y-axis.
#Export as .eps with dimensions w = 1050, h=850 for Illustrator 
#Clean up all fonts- change to Arial/size = regular 14
#Add Depth axis- same as (a) 
#Add (b) - size 16


##################################################################################################
####FIG 5 - RICHNESS THROUGH TIME                                                                #
##################################################################################################

mem<- read.csv("schefferville_lmerdata_MAR2015_formodels.csv") #schefferville_lmerdata_MAR2015_formodels.csv
mem<- na.omit(mem) #Remove DL4_1 observation
#Remove rows 29-31 (other Dolly and Denault observations- so just Dauriat and Knob)
mem<- mem[-30,]
mem<- mem[-29,]
mem<- mem[-28,] #now just 27 observations and 2 lakes. 

mem.DAR<- as.data.frame(subset(mem, Lake == "Dauriat ", drop=T))


mem.KB<- as.data.frame(subset(mem, Lake == "Knob", drop=T))

#DAR richness plot
DAR.S.plot<- ggplot(aes(x=Est_yearCIC, y=Clad_S), data=mem.DAR) + geom_point(size=4) + geom_path()
DAR.S.plot<- DAR.S.plot + xlab("Estimated year") + ylab("Cladoceran rar S")
DAR.S.plot<- DAR.S.plot + theme_bw()
DAR.S.plot<- DAR.S.plot + theme(axis.text.x = element_text(colour="black", size=16))
DAR.S.plot<- DAR.S.plot + theme(axis.text.y = element_text(colour="black", size=16))
DAR.S.plot<- DAR.S.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
DAR.S.plot<- DAR.S.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
DAR.S.plot<- DAR.S.plot + annotate("rect", xmin=1939, xmax=1977, ymin=4, ymax=12, alpha=0.2)
DAR.S.plot<- DAR.S.plot + annotate("text", x=1850, y=12, label = "(a)", size=12)

#KB richness plot 
KB.S.plot<- ggplot(aes(x=Est_yearCIC, y=Clad_S), data=mem.KB) + geom_point(size=4) + geom_path()
KB.S.plot<- KB.S.plot + xlab("Estimated year") + ylab("Cladoceran rar S")
KB.S.plot<- KB.S.plot + theme_bw()
KB.S.plot<- KB.S.plot + theme(axis.text.x = element_text(colour="black", size=16))
KB.S.plot<- KB.S.plot + theme(axis.text.y = element_text(colour="black", size=16))
KB.S.plot<- KB.S.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
KB.S.plot<- KB.S.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
KB.S.plot<- KB.S.plot + annotate("rect", xmin=1939, xmax=1977, ymin=7, ymax=12, alpha=0.2)
KB.S.plot<- KB.S.plot + annotate("text", x=1660, y=12, label = "(b)", size=12)

#combo plot
S.plot.combo<- grid.arrange(DAR.S.plot, KB.S.plot, nrow=2) #could work, but on different scale. 

#Try with facet
S.plot.facet<- ggplot(aes(x=Est_yearCIC, y=Clad_S), data=mem) + geom_point(size=4) + geom_path() + facet_wrap(~Lake, nrow=2) 
S.plot.facet<- S.plot.facet + theme(strip.text.x = element_text(size=25, face="bold"), panel.background=element_rect(fill="white")) 
S.plot.facet<- S.plot.facet + xlab("Estimated year") + ylab("Cladoceran rar S")
S.plot.facet<- S.plot.facet + theme(axis.text.x = element_text(colour="black", size=16))
S.plot.facet<- S.plot.facet + theme(axis.text.y = element_text(colour="black", size=16))
S.plot.facet<- S.plot.facet + theme(axis.title.x = element_text(size = rel(2), angle=00))
S.plot.facet<- S.plot.facet + theme(axis.title.y = element_text(size = rel(2), angle=90))
S.plot.facet<- S.plot.facet + annotate("rect", xmin=1939, xmax=1977, ymin=4, ymax=12, alpha=0.2)


#####
##################################################################################################
####FIG 6 - PCA OF ASSEMBLAGES                                                                   #
##################################################################################################

####Creating a PCA of the assemblages. 
clad.dom_relabund<- read.csv("clad_DARKB_relabund_greater5per_Boscorr.csv") 
#relative abundance for taxa with greater than 5 percent. 
#columns 7:21 are the species rows. 
#this file has had Bosmina counts corrected. 

#Remove columns of B.long and Eubosmina that represented proportion of total that were assigned by interpolation. 
#Removing columns: B.longisostris.assigned and Eubosmina.longispina.assigned
clad.dom_relabund<- clad.dom_relabund[,-15]
clad.dom_relabund<- clad.dom_relabund[,-16]
#Now species columns 7:19

#Remove row 11 --> duplicate time point in DAR? Keep out til sort that out. 
clad.dom_relabund<- clad.dom_relabund[-11,]

#Rename the taxa columns so that they are shorter- will show up better in PCA
colnames(clad.dom_relabund)[8]<- 'A.harpae' #Acroperus harpae
colnames(clad.dom_relabund)[9]<- 'A.affinis' #Alona affinis
colnames(clad.dom_relabund)[10]<- 'A.circumfimbria' #Alona circumfimbria_guttata
colnames(clad.dom_relabund)[11]<- 'A.quadrangularis' #Alona quadrangularis 
colnames(clad.dom_relabund)[12]<- 'Alona spp.' #Alona.spp
colnames(clad.dom_relabund)[13]<- 'A.nana' #Alonella nana
colnames(clad.dom_relabund)[14]<- 'B.longirostris' #Bosmina longirostris total 
colnames(clad.dom_relabund)[15]<- 'E.longispina' #Eubosmina longispina total
colnames(clad.dom_relabund)[16]<- 'C.cf.sphaericus' #Chydorus cf. sphaericus 
colnames(clad.dom_relabund)[17]<- 'C.gibbus' #Chydorus gibbus 
colnames(clad.dom_relabund)[18]<- 'D.longispina' #Daphnia longispina 
colnames(clad.dom_relabund)[19]<- 'Eurycercus spp.' #Eurycercus.spp
colnames(clad.dom_relabund)[20]<- 'P.piger' #Paralona piger

#Run PCA and extract scores
clad.pca<- rda(decostand(clad.dom_relabund[,8:20], method="hell")) #Hellinger transformation 
clad.scr.sit<- as.data.frame(scores(clad.pca, dis="sites", choices=1:2))
clad.scr.sp<- as.data.frame(scores(clad.pca, dis="sp", choices=1:2))

#Combine site score matrix with other descriptive columns from cladoceran data set
clad.scr.sit<- as.data.frame(cbind(clad.dom_relabund$Lake, clad.dom_relabund$Est_year, clad.scr.sit))
colnames(clad.scr.sit)[1]<- 'Lake' 
colnames(clad.scr.sit)[2]<- 'Est_year' #NOTE: REMEMBER TO FIX EXTRAPOLATIONS IF NEEDED. 


#Figure
pca.plot<-ggplot()
pca.plot<- pca.plot + geom_vline(x=0,colour="grey50") 
pca.plot<- pca.plot + geom_hline(y=0,colour="grey50") 
pca.plot<- pca.plot + geom_point(data = clad.scr.sit, aes(x = PC1, y = PC2, shape = Lake), size=4, colour = "dimgrey") + theme_bw()
pca.plot<- pca.plot + geom_path(data = subset(clad.scr.sit, Lake == "Dauriat ", drop=T), aes(x = PC1, y= PC2), linetype = 'solid', size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))
pca.plot<- pca.plot + geom_path(data = subset(clad.scr.sit, Lake == "Knob", drop=T), aes(x = PC1, y= PC2), linetype = 'longdash', size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))
#trajectories mapped- but add year labels for start and end years. 
#start at oldest time point, end at most recent --> double check though. 
pca.plot<- pca.plot + geom_text(data = clad.scr.sp, aes(x = PC1, y = PC2, label = rownames(clad.scr.sp)), size = 5, angle = 0,  vjust = 1, colour = 'mediumblue')
pca.plot<- pca.plot + labs(x= "PC1 = 0.46", y= "PC2 = 0.22") #adjust if data changes 
#NB: total inertia = 0.26
pca.plot<- pca.plot + theme(axis.text.x = element_text(colour="black", size=16))
pca.plot<- pca.plot + theme(axis.text.y = element_text(colour="black", size=16))
pca.plot<- pca.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
pca.plot<- pca.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 
pca.plot<- pca.plot + annotate("text", x=0.36, y=0.24, label="1920") 
pca.plot<- pca.plot + annotate("text", x=0.20, y=0.54, label="2011")
pca.plot<- pca.plot + annotate("text", x=0.32, y=-0.05, label="2003")
pca.plot<- pca.plot + annotate("text", x=0.44, y=-0.088, label="1665")
#ready for export- pending any data changes 


##################################################################################################

##################################################################################################
####FIG 8 - MIXED EFFECT MODEL FIGURE                                                            #
##################################################################################################
##Mixed effect models already computed and tested in LMER script. 
#Just import data in order to make MEM figure. 

mem<- read.csv("schefferville_lmerdata_MAR2015_formodels.csv") #schefferville_lmerdata_MAR2015_formodels.csv
mem<- na.omit(mem) #Remove DL4_1 observation
#Remove rows 29-31 (other Dolly and Denault observations- so just Dauriat and Knob)
mem<- mem[-30,]
mem<- mem[-29,]
mem<- mem[-28,] #now just 27 observations and 2 lakes. 

#MEM model (lake)
mem.mod<- lmer(Clad_S~Metal_EF + (1|Lake), data=mem, REML=FALSE)

#Determine if slope and therefore effect of Metal_EF on Clad_S is significantly different from zero- calculate CI
#for the slope parameter. If CI overlaps with zero, slope is not significantly different from 0 at 0.05 level. 
#CI = Standard Error of the estimate x 1.96 plus or minus the parameter estimate.
#Slope from b1m3b = -0.070766

upperCI<- -0.070766 + 0.008868*1.96 #-0.05338472
lowerCI<- -0.070766 - 0.008868*1.96 #-0.08814728
#Doesn't overlap zero so significantly different from 0. 

lake.coef<- as.data.frame(coef(mem.mod)$Lake)
colnames(lake.coef)<- c("Intercept", "Slope")


#Make a plot just showing the linear (null model)
lm.plot <- ggplot(aes(x=Metal_EF, y=Clad_S, colour=Lake), data=mem) + geom_point() + xlab("Metal EF") + ylab("Cladoceran rar S")
lm.plot<- lm.plot + stat_smooth(method="lm", se=FALSE, size=1)
lm.plot<- lm.plot + scale_colour_manual(values = c("Dauriat " = "darkblue", "Knob" = "firebrick3"))
lm.plot<- lm.plot + theme_bw()
lm.plot<- lm.plot + theme(axis.text.x = element_text(colour="black", size=16))
lm.plot<- lm.plot + theme(axis.text.y = element_text(colour="black", size=16))
lm.plot<- lm.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
lm.plot<- lm.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
lm.plot<- lm.plot + annotate("text", x=0, y=12, label = "(a)")
lm.plot<- lm.plot +theme(legend.position="none")


# Make a plot that includes all the data (mem with lake)
mem.plot <- ggplot(aes(x=Metal_EF, y=Clad_S), data=mem) + geom_point() + xlab("Metal EF") + ylab("Cladoceran rar S")
mem.plot<- mem.plot + geom_abline(intercept=10.832548, slope=-0.070766, linetype="dotdash", size=2, colour="dimgrey") #Overall
#mem.plot<- mem.plot + geom_abline(intercept = lake.coef[1,1], slope=lake.coef[1,2], colour="darkblue", size=2)  #Dauriat
#mem.plot<- mem.plot + geom_abline(intercept = lake.coef[2,1], slope = lake.coef[2,2], colour="firebrick3", size=2) #Knob
mem.plot<- mem.plot + theme_bw()
mem.plot<- mem.plot + theme(axis.text.x = element_text(colour="black", size=16))
mem.plot<- mem.plot + theme(axis.text.y = element_text(colour="black", size=16))
mem.plot<- mem.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
mem.plot<- mem.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
mem.plot<- mem.plot + annotate("text", x=0, y=12, label = "(b)")
#From when had all 4 lakes in
#mem.plot<- mem.plot + geom_abline(intercept = lake.coef[2,1], slope = lake.coef[2,2], colour="darkgoldenrod4", size=2)  #Denault
#mem.plot<- mem.plot + geom_abline(intercept = lake.coef[3,1], slope = lake.coef[3,2], colour="black", size=2)  #Dolly

mod.plots.combo<- grid.arrange(lm.plot, mem.plot, nrow=2)

##################################################################################################

##################################################################################################
####FIG 7 - TEMPORAL BETA-DIVERSITY FIGURES                                                      #
##################################################################################################

######Source functions######
decompose.D2 <- function(Y1, Y2, den.type=2)
  # Compare two surveys:
  # Decompose the Ruzicka and percentage difference dissimilarities into A, B and C.
  #
  # Parameters --
  # 
  # Y1 : First survey data with sites in rows and species in columns. 
  # Y2 : Second survey data with sites in rows and species in columns. 
  # The sites and species must be the same in the two tables, and in the same order.
  # The files may contain species presence-absence or quantitative abundance data.
  # 
  # den.type -- Denominator type for the indices
#     1 : (A+B+C) as in the Ruzicka dissimilarity
#     2 : (2A+B+C) as in the Percentage difference dissimilarity (alias Bray-Curtis)
#
# Value (output of the function) --
#
# mat1 : A, B and C results, with sites in rows and A, B and C in columns.
# mat2 : A,B,C,D divided by a denominator [either (A+B+C) or (2A+B+C)]; D = (B+C).
#
# Details: Numerical results in output matrices --
# A -- aj is the part of the abundance of species j that is common to the two survey vectors: aj = min(y1j, y2j). A is the sum of the aj values for all species in the functional group under study.
# B -- bj is the part of the abundance of species j that is higher in survey 1 than in survey 2: bj = y1j – y2j. B is the sum of the bj values for all species in the functional group under study.
# C -- cj is the part of the abundance of species j that is higher in survey 2 than in survey 1: cj = y2j – y1j. C is the sum of the cj values for all species in the functional group under study.
#
# Example --
# test1 = matrix(runif(50,0,100),10,5)
# test2 = matrix(runif(50,0,100),10,5)
# (res = decompose.D2(test1, test2, den.type=1))# License: GPL-2
#
# License: GPL-2
# Author:: Pierre Legendre
{
  ### Internal function
  den <- function(A,B,C,den.type) if(den.type==1) den=(A+B+C) else den=(2*A+B+C) #specifying the demoninator to use
  ### End internal function
  #
  n = nrow(Y1)
  p = ncol(Y1)
  if(nrow(Y2)!=n) stop("The data tables do not have the same number of rows") #warning messages in case matrices don't match
  if(ncol(Y2)!=p) stop("The data tables do not have the same number of columns")
  #
  ABC = c("A","B","C")        # A = similarity #column headings for output matrix 1
  ABCD = c("A","B","C","D")   # D = dissimilarity #column headings for output matrix 2
  #
  mat1 = matrix(NA,n,3) #output matrix 1
  colnames(mat1) = ABC
  mat2 = matrix(NA,n,4) #output matrix 2 
  colnames(mat2) = ABCD
  if(!is.null(rownames(Y1))) { 
    rownames(mat1)=rownames(mat2)=rownames(Y1) 
  } else {
    rownames(mat1)=rownames(mat2)=paste("Site",1:n,sep=".") }
  #
  for(i in 1:n) {
    YY = rbind(Y1[i,], Y2[i,])
    A = sum(apply(YY,2,min))
    tmp = YY[1,] - YY[2,]
    B = sum(tmp[tmp>0])
    C = -sum(tmp[tmp<0])
    D = B+C
    mat1[i,] = c(A,B,C)
    mat2[i,] = c(A,B,C,D)/den(A,B,C,den.type) #uses the internal function from above
  }
  #
  list(mat1=mat1, mat2=mat2)
}
########End of decompose.D2 function 

#######See BetaDiv script##### 



#Raw data
#bd<- read.csv("schefferville_clad_bd_JULY2015.csv") #schefferville_clad_bd_JULY2015.csv

#Data by intervals
bdtemp.interval.res<- read.csv("schefferville_bdtemp_byinterval_res.csv")

interval.long<- melt(bdtemp.interval.res, id.vars=c("Lake", "Interval_comparison", "Comparison_order", "Year_comparisons", "Year_midpoint_comparison", "Year_midpoint_comparisonCIC"))
colnames(interval.long)[7]<- 'Beta_component'
colnames(interval.long)[8]<- 'Value'

#Subset
interval.long.DAR<- as.data.frame(subset(interval.long, Lake == "DAR", drop=T))
interval.long.KB<- as.data.frame(subset(interval.long, Lake == "KB", drop=T))

#Horizontal figures 
#DAR
bd.DAR<- ggplot(interval.long.DAR, aes(x=Comparison_order, y=Value, colour=Beta_component)) + geom_point(size=4) + geom_path(size=2) 
bd.DAR<- bd.DAR + labs(x= "Comparison", y= "Beta value")
bd.DAR<- bd.DAR + theme_bw()
bd.DAR<- bd.DAR + scale_color_manual(values = c("Total_beta" = "black", "Species_loss" = "darkred", "Species_gain" = "dimgrey"))
bd.DAR<- bd.DAR + theme(axis.text.x = element_text(colour="black", size=16, angle=45))
bd.DAR<- bd.DAR + theme(axis.text.y = element_text(colour="black", size=16))
bd.DAR<- bd.DAR + theme(axis.title.x = element_text(size = rel(2), angle=00))
bd.DAR<- bd.DAR + theme(axis.title.y = element_text(size = rel(2), angle=90))
#bd.DAR<- bd.DAR + scale_x_discrete(labels=c("Pre 1850 - Pre 1850","Pre 1850 - Pre 1850","Pre 1850 - Pre 1850","Pre 1850 - 1896", "1896 - 1921", "1921 - 1925", "1925 - 1929", "1929 - 1941", "1941 - 1960", "1960 - 1966", "1966 - 1969", "1969 - 1977", "1977 - 1986", "1986 - 1998", "1998 - 2006"))
#changed [1921-1921 and 1921-1929] to [1921-1925 and 1925-1929] --> sort this out.  
bd.DAR<- bd.DAR + annotate("text", x=1, y=0.7, label="(a)", size=12)
bd.DAR<- bd.DAR + theme(legend.position="none")
bd.DAR<- bd.DAR + annotate("text", x=15, y=0.45, label="Total", size=11, colour="black")
bd.DAR<- bd.DAR + annotate("text", x=15, y=0.18, label="Loss", size=11, colour="darkred")
bd.DAR<- bd.DAR + annotate("text", x=15, y=0.23, label="Gain", size=11, colour="dimgrey")

#KB
bd.KB<- ggplot(interval.long.KB, aes(x=Comparison_order, y=Value, colour=Beta_component)) + geom_point(size=4) + geom_path(size=2) 
bd.KB<- bd.KB + labs(x= "Comparison", y= "Beta value")
bd.KB<- bd.KB + theme_bw()
bd.KB<- bd.KB + scale_color_manual(values = c("Total_beta" = "black", "Species_loss" = "darkred", "Species_gain" = "dimgrey"))
bd.KB<- bd.KB + theme(axis.text.x = element_text(colour="black", size=16, angle=45))
bd.KB<- bd.KB + theme(axis.text.y = element_text(colour="black", size=16))
bd.KB<- bd.KB + theme(axis.title.x = element_text(size = rel(2), angle=00))
bd.KB<- bd.KB + theme(axis.title.y = element_text(size = rel(2), angle=90))
bd.KB<- bd.KB + scale_x_discrete(labels=c("1665 - 1730","1730 - 1794","1794 - 1883","1883 - 1904", "1904 - 1945", "1945 - 1954", "1954 - 1967", "1967 - 1974", "1974 - 1991", "1991 - 2003"))
bd.KB<- bd.KB + annotate("text", x=1, y=0.4, label="(b)", size=12)
bd.KB<- bd.KB + theme(legend.position="none")
bd.KB<- bd.KB + annotate("text", x=10, y=0.25, label="Total", size=11, colour="black")
bd.KB<- bd.KB + annotate("text", x=10, y=0.14, label="Loss", size=11, colour="darkred")
bd.KB<- bd.KB + annotate("text", x=10, y=0.08, label="Gain", size=11, colour="dimgrey")

fig6combo<- grid.arrange(bd.DAR, bd.KB, nrow=2)


#Vertical figures + faceted so on the same scale. 
#Need to standardize the age points if plotting both lakes on same scale. 
bd.plot<- ggplot(interval.long, aes(x=Year_midpoint_comparison, y=Value, colour=Beta_component)) + geom_point(size=4) + geom_path(size=1) + coord_flip() + facet_wrap(~Lake)
bd.plot<- bd.plot + labs(x= "Midpoint of comparison", y= "Beta value")
bd.plot<- bd.plot + theme_bw()
bd.plot<- bd.plot + scale_color_manual(values = c("Total_beta" = "black", "Species_loss" = "darkred", "Species_gain" = "dimgrey"))
bd.plot<- bd.plot + theme(axis.text.x = element_text(colour="black", size=16))
bd.plot<- bd.plot + theme(axis.text.y = element_text(colour="black", size=16))
bd.plot<- bd.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
bd.plot<- bd.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))

#Re-make vertical figures with only total beta and species loss. 
#Remove species gain component. 
interval.long.red<- as.data.frame(subset(interval.long, Beta_component == "Total_beta" | Beta_component == "Species_loss"))


bd.plot2<- ggplot(interval.long.red, aes(x=Year_midpoint_comparisonCIC, y=Value, colour=Beta_component)) + geom_point(size=4) + geom_path(size=1) + coord_flip() + facet_wrap(~Lake)
bd.plot2<- bd.plot2 + labs(x= "Midpoint of comparison", y= "Beta value")
bd.plot2<- bd.plot2 + theme_bw()
bd.plot2<- bd.plot2 + scale_color_manual(values = c("Total_beta" = "black", "Species_loss" = "darkred"))
bd.plot2<- bd.plot2 + theme(axis.text.x = element_text(colour="black", size=16))
bd.plot2<- bd.plot2 + theme(axis.text.y = element_text(colour="black", size=16))
bd.plot2<- bd.plot2 + theme(axis.title.x = element_text(size = rel(2), angle=00))
bd.plot2<- bd.plot2 + theme(axis.title.y = element_text(size = rel(2), angle=90))
bd.plot2<- bd.plot2 + annotate("rect", xmin=1939, xmax=1977, ymin=0, ymax=0.8, alpha=.2)
bd.plot2<- bd.plot2 + theme(strip.text.x = element_text(size=25, face="bold"), panel.background=element_rect(fill="white"))

#Use legend instead of annotation 
#bd.plot<- bd.plot + annotate("text", x=15, y=0.45, label="Total", size=11, colour="black")
#bd.plot<- bd.plot + annotate("text", x=15, y=0.18, label="Loss", size=11, colour="darkred")
#bd.plot<- bd.plot + annotate("text", x=15, y=0.23, label="Gain", size=11, colour="dimgrey")
#bd.plot<- bd.plot + theme(legend.position="none")
##################################################################################################

##################################################################################################

##################################################################################################
####FIG 9 - CALCULATION OF HYPOTHETICAL TEMPORAL BD                                              #
##################################################################################################

fig9dat<- read.csv("fig9_hypotheticaldata.csv")

scenarioA<- as.data.frame(subset(fig9dat, Scenario == "A", drop=T))

scenarioB<- as.data.frame(subset(fig9dat, Scenario == "B", drop=T))

scenarioC<- as.data.frame(subset(fig9dat, Scenario == "C", drop=T))

#Temporal BD for Scenario A

#Scenario A - T0
A.T0<- as.data.frame(subset(scenarioA, Time_point == "1865", drop=T))

#Scenario A - T1
A.T1<- as.data.frame(subset(scenarioA, Time_point == "1965", drop=T))

#Scenario A - T2
A.T2<- as.data.frame(subset(scenarioA, Time_point == "2016", drop=T))

#Calculate temporal BD
#A T0 to T1
A.bd.T0T1<- decompose.D2(A.T0[,3:10], A.T1[,3:10], den.type=2)  

#A T1 to T2
A.bd.T1T2<- decompose.D2(A.T1[,3:10], A.T2[,3:10], den.type=2)  

#Temporal BD for Scenario B

#Scenario B - T0
B.T0<- as.data.frame(subset(scenarioB, Time_point == "1865", drop=T))

#Scenario B - T1
B.T1<- as.data.frame(subset(scenarioB, Time_point == "1965", drop=T))

#Scenario B - T2
B.T2<- as.data.frame(subset(scenarioB, Time_point == "2016", drop=T))

#Calculate temporal BD
#B T0 to T1
B.bd.T0T1<- decompose.D2(B.T0[,3:10], B.T1[,3:10], den.type=2)  

#B T1 to T2
B.bd.T1T2<- decompose.D2(B.T1[,3:10], B.T2[,3:10], den.type=2) 


#Temporal BD for Scenario C

#Scenario C - T0
C.T0<- as.data.frame(subset(scenarioC, Time_point == "1865", drop=T))

#Scenario C - T1
C.T1<- as.data.frame(subset(scenarioC, Time_point == "1965", drop=T))

#Scenario C - T2
C.T2<- as.data.frame(subset(scenarioC, Time_point == "2016", drop=T))

#C T0 to T1
C.bd.T0T1<- decompose.D2(C.T0[,3:10], C.T1[,3:10], den.type=2)  

#C T1 to T2
C.bd.T1T2<- decompose.D2(C.T1[,3:10], C.T2[,3:10], den.type=2) 


##################################################################################################
####S3 - FULL METAL PROFILES                                                                     #
##################################################################################################

metsp<- read.csv("metalsp_full.csv") #metalsp_full.csv #with NAs for age errors as 0s (then changed back to NAs)
#metsp<- na.omit(metsp)#remove NAs
##Top panel 
metsp #use 
##BUT NEED TO ADD  error for age estimates.
#PLUS confirm that these are the metals I want to showcase. 

#Subset so only uaing DAR and KB. 
metsp<- as.data.frame(subset(metsp, Lake == "Dauriat " | Lake == "Knob", drop=T))

error<- aes(xmax = metsp$Est_yearCIC + metsp$Est_year_error, xmin = metsp$Est_yearCIC - metsp$Est_year_error)

##May be easiest to do this using facet_wrap 
metsp.long<- melt(metsp, id.vars=c("Sample_ID_Master", "Lake", "Est_year", "Est_yearCIC", "Year_type", "Est_year_error"))
colnames(metsp.long)[7]<- 'Metal'
colnames(metsp.long)[8]<- 'Conc'

##Cut out metals that I am not interested in. 
metsp.long.high<- as.data.frame(subset(metsp.long, Metal == "Al.Ti" | Metal == "Fe.Ti" | Metal == "Mn.Ti" | Metal == "Zn.Ti", drop=T))

metsp.long.low<- as.data.frame(subset(metsp.long, Metal == "As.Ti" | Metal == "Cd.Ti" | Metal == "Co.Ti" | Metal == "Cr.Ti" | Metal == "Cu.Ti" | Metal == "Hg.Ti" |
                                         Metal == "Ni.Ti" | Metal == "Sb.Ti", drop=T))

#Panel 1- high conc metals
met.plot<- ggplot(metsp.long.high, aes(x=Est_yearCIC, y=Conc, colour=Lake)) + geom_path() + facet_wrap(~Metal)
met.plot<- met.plot + geom_point(aes(x=Est_yearCIC, y=Conc, shape=factor(Year_type)), size=2)
met.plot<- met.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3"))  
met.plot<- met.plot + scale_shape_manual(values = c(16,1))
met.plot<- met.plot + labs(x = "Estimated year", y = "[Metal] (ppm)") + theme_bw() + theme(legend.position="bottom") + coord_flip()
met.plot<- met.plot + theme(axis.text.x = element_text(colour="black", size=16))
met.plot<- met.plot + theme(axis.text.y = element_text(colour="black", size=16))
met.plot<- met.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
met.plot<- met.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

#Panel 2- low conc metals 
met.plot<- ggplot(metsp.long.low, aes(x=Est_yearCIC, y=Conc, colour=Lake)) + geom_path() + facet_wrap(~Metal)
met.plot<- met.plot + geom_point(aes(x=Est_yearCIC, y=Conc, shape=factor(Year_type)), size=2)
met.plot<- met.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3"))  
met.plot<- met.plot + scale_shape_manual(values = c(16,1))
met.plot<- met.plot + labs(x = "Estimated year", y = "[Metal] (ppm)") + theme_bw() + theme(legend.position="bottom") + coord_flip()
met.plot<- met.plot + theme(axis.text.x = element_text(colour="black", size=16))
met.plot<- met.plot + theme(axis.text.y = element_text(colour="black", size=16))
met.plot<- met.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
met.plot<- met.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 


#########
##Top panel with all 4 lakes - no age error bars. 
cu.plot<- ggplot(metsp, aes(x=Est_year, y=Cu.Ti, colour=Lake)) + geom_path() 
#+ facet_wrap(~Lake)
cu.plot<- cu.plot + geom_point(aes(x=Est_year, y=Cu.Ti, shape=factor(Year_type)), size=4)
cu.plot<- cu.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3", "Dolly" = "black", "Denault" = "darkgoldenrod4"))  
cu.plot<- cu.plot + scale_shape_manual(values = c(16,1))
cu.plot<- cu.plot + labs(x = "Estimated year", y = "[Cu]") + theme_bw() + theme(legend.position="none")
cu.plot<- cu.plot + geom_errorbarh(error, width=0.25, colour="black") + coord_flip()
cu.plot<- cu.plot + theme(axis.text.x = element_text(colour="black", size=16))
cu.plot<- cu.plot + theme(axis.text.y = element_text(colour="black", size=16))
cu.plot<- cu.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
cu.plot<- cu.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

fe.plot<- ggplot(metsp, aes(x=Est_year, y=Fe.Ti, colour=Lake)) + geom_path() 
#+ facet_wrap(~Lake)
fe.plot<- fe.plot + geom_point(aes(x=Est_year, y=Fe.Ti, shape=factor(Year_type)), size=4)
fe.plot<- fe.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3", "Dolly" = "black", "Denault" = "darkgoldenrod4"))   
fe.plot<- fe.plot + scale_shape_manual(values = c(16,1))
fe.plot<- fe.plot + labs(x = "Estimated year", y = "[Fe]") + theme_bw() + theme(legend.position="none")
fe.plot<- fe.plot + geom_errorbarh(error, width=0.25, colour="black") + coord_flip()
fe.plot<- fe.plot + theme(axis.text.x = element_text(colour="black", size=16))
fe.plot<- fe.plot + theme(axis.text.y = element_text(colour="black", size=16))
fe.plot<- fe.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
fe.plot<- fe.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

al.plot<- ggplot(metsp, aes(x=Est_year, y=Al.Ti, colour=Lake)) + geom_path() 
#+ facet_wrap(~Lake)
al.plot<- al.plot + geom_point(aes(x=Est_year, y=Al.Ti, shape=factor(Year_type)), size=4)
al.plot<- al.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3", "Dolly" = "black", "Denault" = "darkgoldenrod4"))   
al.plot<- al.plot + scale_shape_manual(values = c(16,1))
al.plot<- al.plot + labs(x = "Estimated year", y = "[Al]") + theme_bw() + theme(legend.position="none")
al.plot<- al.plot + geom_errorbarh(error, width=0.25, colour="black") + coord_flip()
al.plot<- al.plot + theme(axis.text.x = element_text(colour="black", size=16))
al.plot<- al.plot + theme(axis.text.y = element_text(colour="black", size=16))
al.plot<- al.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
al.plot<- al.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

cd.plot<- ggplot(metsp, aes(x=Est_year, y=Cd.Ti, colour=Lake)) + geom_path() 
#+ facet_wrap(~Lake)
cd.plot<- cd.plot + geom_point(aes(x=Est_year, y=Cd.Ti, shape=factor(Year_type)), size=4)
cd.plot<- cd.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3", "Dolly" = "black", "Denault" = "darkgoldenrod4"))   
cd.plot<- cd.plot + scale_shape_manual(values = c(16,1))
cd.plot<- cd.plot + labs(x = "Estimated year", y = "[Cd]") + theme_bw() + theme(legend.position="none")
cd.plot<- cd.plot + geom_errorbarh(error, width=0.25, colour="black") + coord_flip()
cd.plot<- cd.plot + theme(axis.text.x = element_text(colour="black", size=16))
cd.plot<- cd.plot + theme(axis.text.y = element_text(colour="black", size=16))
cd.plot<- cd.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
cd.plot<- cd.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

#Putting top panel together
fig2top<- grid.arrange(cu.plot, fe.plot, al.plot, cd.plot, nrow=1)


##################################################################################################

##################################################################################################
####RESTING STAGE PLOTS                                                                          #
##################################################################################################

dar.rstage<- read.csv("dauriat_restingstages_Mar2016_CICmodels.csv")

kb.rstage<- read.csv("knob_restingstages_Mar2016_CICmodels.csv")

#DAR
dar.core2<- ggplot(dar.rstage, aes(x=Est_ageCIC, y=Count_per_g, fill=Specimen_type)) + geom_bar(stat = 'identity', width = 4.00) + coord_flip()
dar.core2<- dar.core2 + labs(y= "Counts per gram", x= "Estimated year") + theme_bw()
dar.core2<- dar.core2 + theme(legend.position = "top")
dar.core2<- dar.core2 + scale_fill_manual(values = c("CAS_1" = "coral3", "CAS_2" = "darkgoldenrod1", "EPH_1" = "chocolate4", "EPH_2" = "burlywood4", "EPH_X" = "aquamarine4"))
dar.core2<- dar.core2 + theme(axis.text.x = element_text(colour="black", size=16))
dar.core2<- dar.core2 + theme(axis.text.y = element_text(colour="black", size=16))
dar.core2<- dar.core2 + theme(axis.title.x = element_text(size = rel(2), angle=00))
dar.core2<- dar.core2 + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

#KB
kb.core2<- ggplot(kb.rstage, aes(x=Est_ageCIC, y=Count_per_g, fill=Specimen_type)) + geom_bar(stat = 'identity', width = 4.00) + coord_flip()
kb.core2<- kb.core2 + labs(y= "Counts per gram", x= "Estimated year") + theme_bw()
kb.core2<- kb.core2 + scale_fill_manual(values = c("CAS_1" = "coral3", "CAS_2" = "darkgoldenrod1", "CAS_und" = "chartreuse4", "CAS_X" = "bisque4", "EPH_1" = "chocolate4", "EPH_2" = "burlywood4", "EPH_X" = "aquamarine4"))
kb.core2<- kb.core2 + theme(legend.position = "top")
kb.core2<- kb.core2 + theme(axis.text.x = element_text(colour="black", size=16))
kb.core2<- kb.core2 + theme(axis.text.y = element_text(colour="black", size=16))
kb.core2<- kb.core2 + theme(axis.title.x = element_text(size = rel(2), angle=00))
kb.core2<- kb.core2 + theme(axis.title.y = element_text(size = rel(2), angle=90)) 















###These were replaced in their own script. 

##################################################################################################

##################################################################################################
####AGE PLOTS FOR TECHNICAL APPENDIX (S2) - DAR AND KB ONLY                                      #
##################################################################################################

##From R script "schefferville_age_models.R" in the Resting stages R folder. 

##For now, use those created for technical appendix- may need error bars afterwards. 

#####EXTRA########################################################################################

##################################################################################################
####PCA SCORES FOR DOLLY AND DENAULT                                                             #
##################################################################################################

#metsp<- read.csv("metalsp_full.csv") #metalsp_full.csv #with NAs for age errors as 0s (then changed back to NAs)
#metsp<- na.omit(metsp)#remove NAs

#metsp.pca<- rda(metsp[,6:37])
#pc.scores<- as.data.frame(scores(metsp.pca, dis="sites", choices=1))
#write.csv(pc.scores, "pc.scores.csv")

##################################################################################################

###################################################################################################

##Making beta diversity figures by time period as opposed to interval. 
#Temporal beta diversity results 
#temp.BD<- read.csv("schefferville_bdtemp_res.csv") #schefferville_bdtemp_res.csv

#Subset into DAR and KB
#temp.BD.dar<- as.data.frame(subset(temp.BD, Lake == "DAR", drop=T))
#temp.BD.kb<- as.data.frame(subset(temp.BD, Lake == "KB", drop=T))

#Make long
#temp.BD.dar.long<- melt(temp.BD.dar, id.vars=c("Lake", "Comparison"))
#colnames(temp.BD.dar.long) [3] <- 'Component'
#colnames(temp.BD.dar.long) [4] <- 'Value'

#temp.BD.kb.long<- melt(temp.BD.kb, id.vars=c("Lake", "Comparison"))
#colnames(temp.BD.kb.long) [3] <- 'Component'
#colnames(temp.BD.kb.long) [4] <- 'Value'

#Make time comparisons into factors and re-order them 
#Dar
#temp.BD.dar.long$Comparison<- factor(temp.BD.dar.long$Comparison, levels=c("T0_T1", "T1_T2", "T0_T2"))

#Kb

#DAR plot
#bd.dar.plot<- ggplot(temp.BD.dar.long, aes(x=Comparison, y=Value, colour=Component)) + geom_point(size=4)
#bd.dar.plot<- bd.dar.plot + labs(x= "Time comparison", y= "Beta diversity") + theme_bw()
#bd.dar.plot<- bd.dar.plot + theme(axis.text.x = element_text(colour="black", size=16))
#bd.dar.plot<- bd.dar.plot + theme(axis.text.y = element_text(colour="black", size=16))
#bd.dar.plot<- bd.dar.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
#bd.dar.plot<- bd.dar.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
#bd.dar.plot<- bd.dar.plot + ggtitle("DAR")

#KB plot
#bd.kb.plot<- ggplot(temp.BD.kb.long, aes(x=Comparison, y=Value, colour=Component)) + geom_point(size=4)
#bd.kb.plot<- bd.kb.plot + labs(x= "Time comparison", y= "Beta diversity") + theme_bw()
#bd.kb.plot<- bd.kb.plot + theme(axis.text.x = element_text(colour="black", size=16))
#bd.kb.plot<- bd.kb.plot + theme(axis.text.y = element_text(colour="black", size=16))
#bd.kb.plot<- bd.kb.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
#bd.kb.plot<- bd.kb.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
#bd.kb.plot<- bd.kb.plot + ggtitle("KNOB")

#Combo
#bd.plot<- grid.arrange(bd.dar.plot, bd.kb.plot, nrow=2)
