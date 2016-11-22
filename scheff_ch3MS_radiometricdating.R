##################################################################################################
####SCHEFFERVILLE: Final radiometric dating plots for technical appendix                         #
##################################################################################################

##In this script: Re-do of radiometric dating plots for Knob and Dauriat- based on feedback
#from Jean-Philippe Jenny (March 2016)
#Previous script: 
#R version: 3.1.2 (Pumpkin Helmet)

##Last update: Mar. 25 , 2016
##Associated workspace: 
##Associated markdown: 
##Associated .txt of R script: 
##Github: 

##################################################################################################

##################################################################################################
####PACKAGES                                                                                     #
##################################################################################################

library(plyr)
library(reshape)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(scales)

##################################################################################################

##################################################################################################
####LAKE KNOB - WORKING DIRECTORY                                                                #
##################################################################################################

setwd("C:/Users/Winegardner/Documents/MCGILL/Quebec regional data and field data/Pb210 dating/Knob2_GEOTOP/Final - Oct 2014")

#Data
kbdat<- read.csv("KB2_AGEMODELS_FINAL_COMPARE_Oct2014.csv")
kbdat.CIC<- as.data.frame(subset(kbdat, Model == "CIC", drop=T))
  
#Plot (a) - 210Pb activity vs. depth (log-scale)

plot.a<- ggplot(kbdat.CIC, aes(x=Mean_depth, y=Pb210_dpm_g)) + geom_point(size=2) + geom_path() + coord_flip() + scale_x_reverse() 
plot.a<- plot.a + scale_y_continuous(trans=log10_trans()) #log-scale for y
plot.a<- plot.a + labs(x = "Depth in core (cm)", y = "210Pb activity (dpm/g)") + theme_bw() 
plot.a<- plot.a + theme(axis.text.x = element_text(colour="black", size=16))
plot.a<- plot.a + theme(axis.text.y = element_text(colour="black", size=16))
plot.a<- plot.a + theme(axis.title.x = element_text(size = rel(2), angle=00))
plot.a<- plot.a + theme(axis.title.y = element_text(size = rel(2), angle=90)) 
plot.a<- plot.a + geom_errorbar(aes(ymax = Pb210_dpm_g + Pb210_error, ymin = Pb210_dpm_g - Pb210_error))
plot.a<- plot.a + ggtitle("(a)")

#Plot (b) - 210Pb excess vs. depth (log-scale)
plot.b<- ggplot(kbdat.CIC, aes(x=Mean_depth, y=Pb210_excess)) + geom_point(size=2) + geom_path() + coord_flip() + scale_x_reverse() 
plot.b<- plot.b + scale_y_continuous(trans=log10_trans()) #log-scale for y
plot.b<- plot.b + labs(x = "Depth in core (cm)", y = "210Pb excess (dpm/g)") + theme_bw() 
plot.b<- plot.b + theme(axis.text.x = element_text(colour="black", size=16))
plot.b<- plot.b + theme(axis.text.y = element_text(colour="black", size=16))
plot.b<- plot.b + theme(axis.title.x = element_text(size = rel(2), angle=00))
plot.b<- plot.b + theme(axis.title.y = element_text(size = rel(2), angle=90)) 
plot.b<- plot.b + ggtitle("(b)")
  
#Plot (c) - 137Cs activity vs. depth (log-scale)
plot.c<- ggplot(kbdat.CIC, aes(x=Mean_depth, y=Cs137_dpm_g)) + geom_point(size=2) + geom_path() + coord_flip() + scale_x_reverse() 
plot.c<- plot.c + scale_y_continuous(trans=log10_trans()) #log-scale for y
plot.c<- plot.c + labs(x = "Depth in core (cm)", y = "137Cs activity (dpm/g)") + theme_bw() 
plot.c<- plot.c + theme(axis.text.x = element_text(colour="black", size=16))
plot.c<- plot.c + theme(axis.text.y = element_text(colour="black", size=16))
plot.c<- plot.c + theme(axis.title.x = element_text(size = rel(2), angle=00))
plot.c<- plot.c + theme(axis.title.y = element_text(size = rel(2), angle=90)) 
plot.c<- plot.c + geom_errorbar(aes(ymax = Cs137_dpm_g + Cs137_error, ymin = Cs137_dpm_g - Cs137_error))
plot.c<- plot.c + ggtitle("(c)")
  
#Plot (d) - 210Pb excess vs. age from CIC model 
plot.d<- ggplot(kbdat.CIC, aes(x=Age, y=Pb210_excess)) + geom_point(size=2) + geom_path() + coord_flip()  
plot.d<- plot.d + labs(x = "Age from CIC model", y = "210Pb excess (dpm/g)") + theme_bw() 
plot.d<- plot.d + theme(axis.text.x = element_text(colour="black", size=16))
plot.d<- plot.d + theme(axis.text.y = element_text(colour="black", size=16))
plot.d<- plot.d + theme(axis.title.x = element_text(size = rel(2), angle=00))
plot.d<- plot.d + theme(axis.title.y = element_text(size = rel(2), angle=90)) 
plot.d<- plot.d + ggtitle("(a)")
  
#Plot (e) - 137Cs excess vs. age from CIC model 
plot.e<- ggplot(kbdat.CIC, aes(x=Age, y=Cs137_dpm_g)) + geom_point(size=2) + geom_path() + coord_flip()  
plot.e<- plot.e + labs(x = "Age from CIC model", y = "137Cs activity (dpm/g)") + theme_bw() 
plot.e<- plot.e + theme(axis.text.x = element_text(colour="black", size=16))
plot.e<- plot.e + theme(axis.text.y = element_text(colour="black", size=16))
plot.e<- plot.e + theme(axis.title.x = element_text(size = rel(2), angle=00))
plot.e<- plot.e + theme(axis.title.y = element_text(size = rel(2), angle=90)) 
plot.e<- plot.e + ggtitle("(b)")

#Put plots together
plot.kb<- grid.arrange(plot.a, plot.b, plot.c, nrow=1)

plot.kb2<- grid.arrange(plot.d, plot.e, nrow=1)

##################################################################################################

##################################################################################################
####LAKE DAURIAT - WORKING DIRECTORY                                                                #
##################################################################################################

setwd("C:/Users/Winegardner/Documents/MCGILL/Quebec regional data and field data/Pb210 dating/Dauriat2_GEOTOP/FINAL- Oct 2014")

#Data
dardat<- read.csv("DAR2_AGEMODELS_FINAL_COMPARE_Oct2014.csv")
dardat.CIC<- as.data.frame(subset(dardat, Model == "CIC", drop=T))

#Plot (a) - 210Pb activity vs. depth (log-scale)

plot.a<- ggplot(dardat.CIC, aes(x=Mean_depth, y=Pb210_dpm_g)) + geom_point(size=2) + geom_path() + coord_flip() + scale_x_reverse() 
plot.a<- plot.a + scale_y_continuous(trans=log10_trans()) #log-scale for y
plot.a<- plot.a + labs(x = "Depth in core (cm)", y = "210Pb activity (dpm/g)") + theme_bw() 
plot.a<- plot.a + theme(axis.text.x = element_text(colour="black", size=16))
plot.a<- plot.a + theme(axis.text.y = element_text(colour="black", size=16))
plot.a<- plot.a + theme(axis.title.x = element_text(size = rel(2), angle=00))
plot.a<- plot.a + theme(axis.title.y = element_text(size = rel(2), angle=90)) 
plot.a<- plot.a + geom_errorbar(aes(ymax = Pb210_dpm_g + Pb210_error, ymin = Pb210_dpm_g - Pb210_error))
plot.a<- plot.a + ggtitle("(a)")

#Plot (b) - 210Pb excess vs. depth (log-scale)
plot.b<- ggplot(dardat.CIC, aes(x=Mean_depth, y=Pb210_excess)) + geom_point(size=2) + geom_path() + coord_flip() + scale_x_reverse() 
plot.b<- plot.b + scale_y_continuous(trans=log10_trans()) #log-scale for y
plot.b<- plot.b + labs(x = "Depth in core (cm)", y = "210Pb excess (dpm/g)") + theme_bw() 
plot.b<- plot.b + theme(axis.text.x = element_text(colour="black", size=16))
plot.b<- plot.b + theme(axis.text.y = element_text(colour="black", size=16))
plot.b<- plot.b + theme(axis.title.x = element_text(size = rel(2), angle=00))
plot.b<- plot.b + theme(axis.title.y = element_text(size = rel(2), angle=90)) 
plot.b<- plot.b + ggtitle("(b)")

#Plot (c) - 137Cs activity vs. depth (log-scale)
plot.c<- ggplot(dardat.CIC, aes(x=Mean_depth, y=Cs137_dpm_g)) + geom_point(size=2) + geom_path() + coord_flip() + scale_x_reverse() 
plot.c<- plot.c + scale_y_continuous(trans=log10_trans()) #log-scale for y
plot.c<- plot.c + labs(x = "Depth in core (cm)", y = "137Cs activity (dpm/g)") + theme_bw() 
plot.c<- plot.c + theme(axis.text.x = element_text(colour="black", size=16))
plot.c<- plot.c + theme(axis.text.y = element_text(colour="black", size=16))
plot.c<- plot.c + theme(axis.title.x = element_text(size = rel(2), angle=00))
plot.c<- plot.c + theme(axis.title.y = element_text(size = rel(2), angle=90)) 
plot.c<- plot.c + geom_errorbar(aes(ymax = Cs137_dpm_g + Cs137_error, ymin = Cs137_dpm_g - Cs137_error))
plot.c<- plot.c + ggtitle("(c)")

#Plot (d) - 210Pb excess vs. age from CIC model 
plot.d<- ggplot(dardat.CIC, aes(x=Age, y=Pb210_excess)) + geom_point(size=2) + geom_path() + coord_flip()  
plot.d<- plot.d + labs(x = "Age from CIC model", y = "210Pb excess (dpm/g)") + theme_bw() 
plot.d<- plot.d + theme(axis.text.x = element_text(colour="black", size=16))
plot.d<- plot.d + theme(axis.text.y = element_text(colour="black", size=16))
plot.d<- plot.d + theme(axis.title.x = element_text(size = rel(2), angle=00))
plot.d<- plot.d + theme(axis.title.y = element_text(size = rel(2), angle=90)) 
plot.d<- plot.d + ggtitle("(a)")

#Plot (e) - 137Cs excess vs. age from CIC model 
plot.e<- ggplot(dardat.CIC, aes(x=Age, y=Cs137_dpm_g)) + geom_point(size=2) + geom_path() + coord_flip()  
plot.e<- plot.e + labs(x = "Age from CIC model", y = "137Cs activity (dpm/g)") + theme_bw() 
plot.e<- plot.e + theme(axis.text.x = element_text(colour="black", size=16))
plot.e<- plot.e + theme(axis.text.y = element_text(colour="black", size=16))
plot.e<- plot.e + theme(axis.title.x = element_text(size = rel(2), angle=00))
plot.e<- plot.e + theme(axis.title.y = element_text(size = rel(2), angle=90)) 
plot.e<- plot.e + ggtitle("(b)")

#Put plots together
plot.dar<- grid.arrange(plot.a, plot.b, plot.c, nrow=1)

plot.dar2<- grid.arrange(plot.d, plot.e, nrow=1)
