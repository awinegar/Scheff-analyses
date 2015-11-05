##################################################################################################
####SCHEFFERVILLE: Ch. 3 manuscript figures                                                      #
##################################################################################################

##In this script: Final analyses and figures for 3rd chapter (general Schefferville manuscript)
#Previous script: Various R scripts from R studio scripts/Schefferville summer 2013
#R version: 3.1.2 (Pumpkin Helmet)

##Last update: Nov. 5, 2015
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
regMAP1<- regMAP1 + theme(axis.text.x = element_text(colour="black",size=10))
regMAP1<- regMAP1 + theme(axis.text.y = element_text(colour="black",size=10))
regMAP1<- regMAP1 + theme(axis.title.x = element_text(size = rel(1.5), angle=00))
regMAP1<- regMAP1 + theme(axis.title.y = element_text(size = rel(1.5), angle=90))
regMAP1<- regMAP1 + annotate("text", x=-78, y=61.2, label="(a)", size=10)

#Fig 1B - Four lakes across the Schefferville landscape 
coor <- read.csv("lakes_coor.csv") #lakes_coor.csv in set working directory 

coor1<-get_map(location = "Schefferville Quebec", zoom = 12, source = 'google')

coorMAP1<-ggmap(coor1)+
  geom_point(aes(x=Lon,y=Lat), data = coor,colour="red",size=7,alpha=0.5) + labs(x="Longitude (DD)",y="Latitude (DD)")
coorMAP1<- coorMAP1 + annotate("text", x=-66.82, y=54.810, label="Dauriat")
coorMAP1<- coorMAP1 + annotate("text", x=-66.80, y=54.786, label="Knob")
coorMAP1<- coorMAP1 + annotate("text", x=-66.75, y=54.810, label="Dolly")
coorMAP1<- coorMAP1 + annotate("text", x=-66.87, y=54.84, label="Denault")
coorMAP1<- coorMAP1 + theme(axis.text.x = element_text(colour="black",size=10))
coorMAP1<- coorMAP1 + theme(axis.text.y = element_text(colour="black",size=10))
coorMAP1<- coorMAP1 + theme(axis.title.x = element_text(size = rel(1.5), angle=00))
coorMAP1<- coorMAP1 + theme(axis.title.y = element_text(size = rel(1.5), angle=90))
coorMAP1<- coorMAP1 + annotate("text", x=-66.9, y=54.88, label="(b)", size=10)

#Putting Figure 1a and 1b together
fig1<- grid.arrange(regMAP1, coorMAP1, nrow=1)

##################################################################################################

##################################################################################################
####FIG 2 - PCA SCORES FOR DOLLY AND DENAULT                                                     #
##################################################################################################

metsp<- read.csv("metalsp_full.csv") #metalsp_full.csv #with NAs for age errors as 0s (then changed back to NAs)
#metsp<- na.omit(metsp)#remove NAs

#metsp.pca<- rda(metsp[,6:37])
#pc.scores<- as.data.frame(scores(metsp.pca, dis="sites", choices=1))
#write.csv(pc.scores, "pc.scores.csv")

##################################################################################################

##################################################################################################
####FIG 2 - METAL LOADING FIGURE                                                                 #
##################################################################################################
##Plots of metal concentrations and metal loading. 

##Top panel 
metsp #use 
##BUT NEED TO ADD  error for age estimates.
#PLUS confirm that these are the metals I want to showcase. 

error<- aes(xmax = metsp$Est_year + metsp$Est_year_error, xmin = metsp$Est_year - metsp$Est_year_error)

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


##Bottom panel 
ef<- read.csv("euk_metalEF_full.csv") #includes EFs for a suite of metals as opposed to just the composite. 

error.ef<- aes(xmax = ef$Est_year + ef$Est_year_error, xmin = ef$Est_year - ef$Est_year_error)

#Plot- works if use just total metal_EF - just 2 lakes
ef.plot<- ggplot(ef, aes(x=Est_year, y=Total_EF, colour = Lake)) + geom_path() + coord_flip() #+ facet_wrap(~Lake)
ef.plot<- ef.plot + geom_point(aes(x=Est_year, y=Total_EF, shape=factor(Year_type)), size=4)
ef.plot<- ef.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3", "Dolly" = "black", "Denault" = "darkgoldenrod4"))  
ef.plot<- ef.plot + scale_shape_manual(values = c(16,1))
ef.plot<- ef.plot + labs(x = "Estimated year", y = "Metal EF") + theme_bw() + theme(legend.position="none")
ef.plot<- ef.plot + geom_errorbarh(error.ef, width=0.25, colour="black") 
ef.plot<- ef.plot + theme(axis.text.x = element_text(colour="black", size=16))
ef.plot<- ef.plot + theme(axis.text.y = element_text(colour="black", size=16))
ef.plot<- ef.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
ef.plot<- ef.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

#PC plots
pcmet.plot<- ggplot(ef, aes(x=Est_year, y=Metal_PC1, colour = Lake)) + geom_path() + coord_flip() #+ facet_wrap(~Lake)
pcmet.plot<- pcmet.plot + geom_point(aes(x=Est_year, y=Metal_PC1, shape=factor(Year_type)), size=4)
pcmet.plot<- pcmet.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3", "Dolly" = "black", "Denault" = "darkgoldenrod4"))  
pcmet.plot<- pcmet.plot + scale_shape_manual(values = c(16,1))
pcmet.plot<- pcmet.plot + labs(x = "Estimated year", y = "Metal PC1") + theme_bw() + theme(legend.position="none")
pcmet.plot<- pcmet.plot + geom_errorbarh(error.ef, width=0.25, colour="black") 
pcmet.plot<- pcmet.plot + theme(axis.text.x = element_text(colour="black", size=16))
pcmet.plot<- pcmet.plot + theme(axis.text.y = element_text(colour="black", size=16))
pcmet.plot<- pcmet.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
pcmet.plot<- pcmet.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))

#Putting bottom panel together
fig2bottom<- grid.arrange(ef.plot, pcmet.plot, nrow=1)

##################################################################################################

##################################################################################################
####FIG 3/4 - CLADOCERAN STRATIGRAPHIES FOR DAR AND KB (full profiles)                             #
##################################################################################################
##Use rioja package for C2 plots. Only for DAR and KB since have cladoceran profiles. 

#Count and relative abundance data
clad<- read.csv("clad_DARKB.csv")
clad.rel<- read.csv("clad_DARKB_relabund.csv")

#Subset into DAR and KB
clad.DAR<- as.data.frame(subset(clad, Lake == "Dauriat ", drop=T))
clad.KB<- as.data.frame(subset(clad, Lake == "Knob", drop=T))

clad.rel.DAR<- as.data.frame(subset(clad.rel, Lake == "Dauriat ", drop=T))
clad.rel.KB<- as.data.frame(subset(clad.rel, Lake == "Knbo", drop=T))

#Data with only species that had at least 5% relative abundance in at least one interval. 
#Species richness stil takes account ALL species. 
clad.dom<- read.csv("clad_DARKB_greater5per.csv")
clad.dom_relabund<- read.csv("clad_DARKB_relabund_greater5per.csv") 
clad.dom_percent<- read.csv("clad_DARKB_relabund_percent.csv") #relative abundance from _greater5per multiplied by 100. 

#Subset into DAR and KB
clad.dom.DAR<- as.data.frame(subset(clad.dom_percent, Lake == "Dauriat ", drop=T))
clad.dom.KB<- as.data.frame(subset(clad.dom_percent, Lake == "Knob", drop=T))

#Dar plots
names<- c("Cladoceran richness", "A. harpae", "A. affinis", "A. circumfimbria", "A. qaudrangularis", "Alona spp.", "A. nana", "Bosmina spp.", "C. sphaericus", "C. gibbus", "D. longispina", "Eurycercus spp.", "P. piger") #short forms of names--> write out in legend
dar.strat<- strat.plot(clad.dom.DAR[,6:18], yvar=clad.dom.DAR$Est_year, scale.percent=TRUE, plot.bar=TRUE, plot.line=FALSE,
                       col.bar="darkgrey", lwd.bar=8, ylabel = "Estimated year", cex.ylabel=1, cex.xlabel=1,
                       srt.xlabel=80, x.names=names)

#Knob plots 
kb.strat<- strat.plot(clad.dom.KB[,6:18], yvar=clad.dom.DAR$Est_year, scale.percent=TRUE, plot.bar=TRUE, plot.line=FALSE,
                       col.bar="darkgrey", lwd.bar=8, ylabel = "Estimated year", cex.ylabel=1, cex.xlabel=1,
                       srt.xlabel=80, x.names=names)

##################################################################################################

##################################################################################################
####FIG 5 - MIXED EFFECT MODEL FIGURE                                                            #
##################################################################################################
##Mixed effect models already computed and tested in LMER script. 
#Just import data in order to make MEM figure. 

mem<- read.csv("schefferville_lmerdata_MAR2015_formodels.csv") #schefferville_lmerdata_MAR2015_formodels.csv
mem<- na.omit(mem) #Remove DL4_1 observation

mem.mod<- lmer(Clad_S~Metal_EF + (1|Lake), data=mem, REML=FALSE)

#Determine if slope and therefore effect of Metal_EF on Clad_S is significantly different from zero- calculate CI
#for the slope parameter. If CI overlaps with zero, slope is not significantly different from 0 at 0.05 level. 
#CI = Standard Error of the estimate x 1.96 plus or minus the parameter estimate.
#Slope from b1m3b = -0.055

upperCI<- -0.055 + 0.01148*1.96 #-0.032
lowerCI<- -0.055 - 0.01148*1.96 #-0.0775
#Doesn't overlap zero so significantly different from 0. 

lake.coef<- as.data.frame(coef(mem.mod)$Lake)
colnames(lake.coef)<- c("Intercept", "Slope")

# Make a plot that includes all the data
mem.plot <- ggplot(aes(x=Metal_EF, y=Clad_S), data=mem) + geom_point() + xlab("Metal EF") + ylab("Cladoceran rar S")
mem.plot<- mem.plot + geom_abline(intercept=9.3548, slope=-0.055, linetype="dotdash", size=2, colour="dimgrey") #Overall
mem.plot<- mem.plot + geom_abline(intercept = lake.coef[1,1], slope=lake.coef[1,2], colour="darkblue", size=2)  #Dauriat
mem.plot<- mem.plot + geom_abline(intercept = lake.coef[2,1], slope = lake.coef[2,2], colour="darkgoldenrod4", size=2)  #Denault
mem.plot<- mem.plot + geom_abline(intercept = lake.coef[3,1], slope = lake.coef[3,2], colour="black", size=2)  #Dolly
mem.plot<- mem.plot + geom_abline(intercept = lake.coef[4,1], slope = lake.coef[4,2], colour="firebrick3", size=2) #Knob
mem.plot<- mem.plot + theme_bw()
mem.plot<- mem.plot + theme(axis.text.x = element_text(colour="black", size=16))
mem.plot<- mem.plot + theme(axis.text.y = element_text(colour="black", size=16))
mem.plot<- mem.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
mem.plot<- mem.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))

##################################################################################################

##################################################################################################
####FIG 6 - TEMPORAL BETA-DIVERSITY FIGURES                                                      #
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
bd<- read.csv("schefferville_clad_bd_JULY2015.csv") #schefferville_clad_bd_JULY2015.csv

#Data by intervals
bdtemp.interval.res<- read.csv("schefferville_bdtemp_byinterval_res.csv")

interval.long<- melt(bdtemp.interval.res, id.vars=c("Lake", "Interval_comparison", "Comparison_order", "Year_comparisons"))
colnames(interval.long)[5]<- 'Beta_component'
colnames(interval.long)[6]<- 'Value'

#Subset
interval.long.DAR<- as.data.frame(subset(interval.long, Lake == "DAR", drop=T))
interval.long.KB<- as.data.frame(subset(interval.long, Lake == "KB", drop=T))

#DAR
bd.DAR<- ggplot(interval.long.DAR, aes(x=Comparison_order, y=Value, colour=Beta_component)) + geom_point(size=4) + geom_path(size=2) 
bd.DAR<- bd.DAR + labs(x= "Comparison", y= "Beta value")
bd.DAR<- bd.DAR + theme_bw()
bd.DAR<- bd.DAR + scale_color_manual(values = c("Total_beta" = "black", "Species_loss" = "darkred", "Species_gain" = "dimgrey"))
bd.DAR<- bd.DAR + theme(axis.text.x = element_text(colour="black", size=16, angle=45))
bd.DAR<- bd.DAR + theme(axis.text.y = element_text(colour="black", size=16))
bd.DAR<- bd.DAR + theme(axis.title.x = element_text(size = rel(2), angle=00))
bd.DAR<- bd.DAR + theme(axis.title.y = element_text(size = rel(2), angle=90))
bd.DAR<- bd.DAR + scale_x_discrete(labels=c("Pre 1850 - Pre 1850","Pre 1850 - Pre 1850","Pre 1850 - Pre 1850","Pre 1850 - 1896", "1896 - 1921", "1921 - 1921", "1921 - 1929", "1929 - 1941", "1941 - 1960", "1960 - 1966", "1966 - 1969", "1969 - 1977", "1977 - 1986", "1986 - 1998", "1998 - 2006"))
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

##################################################################################################



#####EXTRA#####
#Panel with only 2 lakes - probably will not use. 
#cu.plot<- ggplot(metsp, aes(x=Est_year, y=Cu.Ti, colour=Lake)) + geom_path() + coord_flip() 
#+ facet_wrap(~Lake)
#cu.plot<- cu.plot + geom_point(aes(x=Est_year, y=Cu.Ti, shape=factor(Year_type)), size=4)
#cu.plot<- cu.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3"))  
#cu.plot<- cu.plot + scale_shape_manual(values = c(16,1))
#cu.plot<- cu.plot + labs(x = "Estimated year", y = "[Cu]") + theme_bw() + theme(legend.position="none")
#cu.plot<- cu.plot + theme(axis.text.x = element_text(colour="black", size=16))
#cu.plot<- cu.plot + theme(axis.text.y = element_text(colour="black", size=16))
#cu.plot<- cu.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
#cu.plot<- cu.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

#fe.plot<- ggplot(metsp, aes(x=Est_year, y=Fe.Ti, colour=Lake)) + geom_path() + coord_flip() 
#+ facet_wrap(~Lake)
#fe.plot<- fe.plot + geom_point(aes(x=Est_year, y=Fe.Ti, shape=factor(Year_type)), size=4)
#fe.plot<- fe.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3"))  
#fe.plot<- fe.plot + scale_shape_manual(values = c(16,1))
#fe.plot<- fe.plot + labs(x = "Estimated year", y = "[Fe]") + theme_bw() + theme(legend.position="none")
#fe.plot<- fe.plot + theme(axis.text.x = element_text(colour="black", size=16))
#fe.plot<- fe.plot + theme(axis.text.y = element_text(colour="black", size=16))
#fe.plot<- fe.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
#fe.plot<- fe.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

#al.plot<- ggplot(metsp, aes(x=Est_year, y=Al.Ti, colour=Lake)) + geom_path() + coord_flip() 
#+ facet_wrap(~Lake)
#al.plot<- al.plot + geom_point(aes(x=Est_year, y=Al.Ti, shape=factor(Year_type)), size=4)
#al.plot<- al.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3"))  
#al.plot<- al.plot + scale_shape_manual(values = c(16,1))
#al.plot<- al.plot + labs(x = "Estimated year", y = "[Al]") + theme_bw() + theme(legend.position="none")
#al.plot<- al.plot + theme(axis.text.x = element_text(colour="black", size=16))
#al.plot<- al.plot + theme(axis.text.y = element_text(colour="black", size=16))
#al.plot<- al.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
#al.plot<- al.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

#cd.plot<- ggplot(metsp, aes(x=Est_year, y=Cd.Ti, colour=Lake)) + geom_path() + coord_flip() 
#+ facet_wrap(~Lake)
#cd.plot<- cd.plot + geom_point(aes(x=Est_year, y=Cd.Ti, shape=factor(Year_type)), size=4)
#cd.plot<- cd.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3"))  
#cd.plot<- cd.plot + scale_shape_manual(values = c(16,1))
#cd.plot<- cd.plot + labs(x = "Estimated year", y = "[Cd]") + theme_bw() + theme(legend.position="none")
#cd.plot<- cd.plot + theme(axis.text.x = element_text(colour="black", size=16))
#cd.plot<- cd.plot + theme(axis.text.y = element_text(colour="black", size=16))
#cd.plot<- cd.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
#cd.plot<- cd.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

#Putting top panel together
#fig2top<- grid.arrange(cu.plot, fe.plot, al.plot, cd.plot, nrow=1)

#Plot- works if use just total metal_EF - just 2 lakes
#ef.plot<- ggplot(ef, aes(x=Est_year, y=Total_EF, colour = Lake)) + geom_path() + coord_flip() + facet_wrap(~Lake)
#ef.plot<- ef.plot + geom_point(aes(x=Est_year, y=Total_EF, shape=factor(Year_type)), size=4)
#ef.plot<- ef.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3"))  
#ef.plot<- ef.plot + scale_shape_manual(values = c(16,1))
#ef.plot<- ef.plot + labs(x = "Estimated year", y = "Metal EF") + theme_bw() + theme(legend.position="none")
#ef.plot<- ef.plot + theme(axis.text.x = element_text(colour="black", size=16))
#ef.plot<- ef.plot + theme(axis.text.y = element_text(colour="black", size=16))
#ef.plot<- ef.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
#ef.plot<- ef.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

#PC plots
#pcmet.plot<- ggplot(ef, aes(x=Est_year, y=Metal_PC1, colour = Lake)) + geom_path() + coord_flip() + facet_wrap(~Lake)
#pcmet.plot<- pcmet.plot + geom_point(aes(x=Est_year, y=Metal_PC1, shape=factor(Year_type)), size=4)
#pcmet.plot<- pcmet.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3"))  
#pcmet.plot<- pcmet.plot + scale_shape_manual(values = c(16,1))
#pcmet.plot<- pcmet.plot + labs(x = "Estimated year", y = "Metal PC1") + theme_bw() + theme(legend.position="none")
#pcmet.plot<- pcmet.plot + theme(axis.text.x = element_text(colour="black", size=16))
#pcmet.plot<- pcmet.plot + theme(axis.text.y = element_text(colour="black", size=16))
#pcmet.plot<- pcmet.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
#pcmet.plot<- pcmet.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
###################################################################################################

##Making beta diversity figures by time period as opposed to interval. 
#Temporal beta diversity results 
temp.BD<- read.csv("schefferville_bdtemp_res.csv") #schefferville_bdtemp_res.csv

#Subset into DAR and KB
temp.BD.dar<- as.data.frame(subset(temp.BD, Lake == "DAR", drop=T))
temp.BD.kb<- as.data.frame(subset(temp.BD, Lake == "KB", drop=T))

#Make long
temp.BD.dar.long<- melt(temp.BD.dar, id.vars=c("Lake", "Comparison"))
colnames(temp.BD.dar.long) [3] <- 'Component'
colnames(temp.BD.dar.long) [4] <- 'Value'

temp.BD.kb.long<- melt(temp.BD.kb, id.vars=c("Lake", "Comparison"))
colnames(temp.BD.kb.long) [3] <- 'Component'
colnames(temp.BD.kb.long) [4] <- 'Value'

#Make time comparisons into factors and re-order them 
#Dar
temp.BD.dar.long$Comparison<- factor(temp.BD.dar.long$Comparison, levels=c("T0_T1", "T1_T2", "T0_T2"))

#Kb

#DAR plot
bd.dar.plot<- ggplot(temp.BD.dar.long, aes(x=Comparison, y=Value, colour=Component)) + geom_point(size=4)
bd.dar.plot<- bd.dar.plot + labs(x= "Time comparison", y= "Beta diversity") + theme_bw()
bd.dar.plot<- bd.dar.plot + theme(axis.text.x = element_text(colour="black", size=16))
bd.dar.plot<- bd.dar.plot + theme(axis.text.y = element_text(colour="black", size=16))
bd.dar.plot<- bd.dar.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
bd.dar.plot<- bd.dar.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
bd.dar.plot<- bd.dar.plot + ggtitle("DAR")

#KB plot
bd.kb.plot<- ggplot(temp.BD.kb.long, aes(x=Comparison, y=Value, colour=Component)) + geom_point(size=4)
bd.kb.plot<- bd.kb.plot + labs(x= "Time comparison", y= "Beta diversity") + theme_bw()
bd.kb.plot<- bd.kb.plot + theme(axis.text.x = element_text(colour="black", size=16))
bd.kb.plot<- bd.kb.plot + theme(axis.text.y = element_text(colour="black", size=16))
bd.kb.plot<- bd.kb.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
bd.kb.plot<- bd.kb.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
bd.kb.plot<- bd.kb.plot + ggtitle("KNOB")

#Combo
bd.plot<- grid.arrange(bd.dar.plot, bd.kb.plot, nrow=2)
