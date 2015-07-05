##################################################################################################
####BETA-DIVERSITY ANALYSES FOR SCHEFFERVILLE CLADOCERANS                                        #
##################################################################################################

##In this script: Beta diversity analyses (mostly temporal), alpha and gamma diversity 
#Previous script: scheff_LMER_March2015.R
#R version: 3.1.2 (Pumpkin Helmet)

##Last update: July 5, 2015 
##Associated workspace: workspace_scheff_BetaDiv.RData 
##Associated markdown: 
##Github: Scheff-analyses repository 

##################################################################################################

##################################################################################################
####PACKAGES                                                                                     #
##################################################################################################

#Data prep and diversity 
library(permute)
library(vegan)
library(boot)
library(rich)
library(Iso)
library(vegan)
library(plyr)
library(reshape)
library(reshape2)

#Beta diversity 
#No specific packages needed

##################################################################################################

##################################################################################################
####SOURCE FUNCTIONS                                                                             #
##################################################################################################

####BETA.DIV FUNCTION- for computing beta diversity (also computes LCBD and SCBD) (spatial)
##################################################################################################
beta.div <- function(Y, method="hellinger", sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
  #
  # Compute estimates of total beta diversity as the total variance in Y, 
  # for 20 dissimilarity coefficients or analysis of raw data (not recommended). 
  # LCBD indices are tested by permutation within columns of Y.
  # This version includes direct calculation of the Jaccard, Sorensen and Ochiai 
  # coefficients for presence-absence data.
  #
  # Arguments --
  # 
  # Y : community composition data matrix.
  # method : name of one of the 20 dissimilarity coefficients, or "none" for
#          direct calculation on Y (also the case with method="euclidean").
# sqrt.D : If sqrt.D=TRUE, the distances in matrix D are square-rooted before 
#          computation of SStotal, BDtotal and LCBD. 
# samp : If samp=TRUE, the abundance-based distances (ab.jaccard, ab.sorensen,
#        ab.ochiai, ab.simpson) are computed for sample data. If samp=FALSE, 
#        they are computed for true population data.
# nperm : Number of permutations for test of LCBD.
# save.D : If save.D=TRUE, the distance matrix will appear in the output list.
# clock : If clock=TRUE, the computation time is printed in the R console.
#
# Reference --
#
# Legendre, P. and M. De Cáceres. 2013. Beta diversity as the variance of 
# community data: dissimilarity coefficients and partitioning. 
# Ecology Letters 16: 951-963. 
#
# License: GPL-2 
# Author:: Pierre Legendre, December 2012, April-May 2013
{
  ### Internal functions
  centre <- function(D,n)
    # Centre a square matrix D by matrix algebra
    # mat.cen = (I - 11'/n) D (I - 11'/n)
  {  One <- matrix(1,n,n)
  mat <- diag(n) - One/n
  mat.cen <- mat %*% D %*% mat
  }
  ###
  BD.group1 <- function(Y, method, save.D, per, n)
  {
    if(method=="profiles") Y = decostand(Y, "total")
    if(method=="hellinger") Y = decostand(Y, "hellinger")
    if(method=="chord") Y = decostand(Y, "norm")
    if(method=="chisquare") Y = decostand(Y, "chi.square")
    #
    s <- scale(Y, center=TRUE, scale=FALSE)^2   # eq. 1
    SStotal <- sum(s)          # eq. 2
    BDtotal <- SStotal/(n-1)   # eq. 3
    if(!per) { SCBD<-apply(s,2,sum)/SStotal }else{ SCBD<-NA }  # eqs. 4a and 4b
    LCBD <- apply(s, 1, sum)/SStotal  # eqs. 5a and 5b
    #
    D <- NA
    if(!per & save.D)   D <- dist(Y)
    #
    out <- list(SStotal_BDtotal=c(SStotal,BDtotal), SCBD=SCBD, LCBD=LCBD, 
                method=method, D=D)
  }
  ###
  BD.group2 <- function(Y, method, sqrt.D, n)
  {
    if(method == "divergence") {
      D = D11(Y)		
      
    } else if(any(method == 
                  c("jaccard","sorensen","ochiai"))) 
    {
      if(method=="jaccard") D = dist.binary(Y, method=1) # ade4 takes sqrt(D)
      if(method=="sorensen")  D = dist.binary(Y, method=5) #ade4 takes sqrt(D)
      if(method=="ochiai") D = dist.binary(Y, method=7) # ade4 takes sqrt(D)
      
    } else if(any(method == 
                  c("manhattan","canberra","whittaker","percentagedifference","wishart"))) 
    {
      if(method=="manhattan") D = vegdist(Y, "manhattan")
      if(method=="canberra")  D = vegdist(Y, "canberra")
      if(method=="whittaker") D = vegdist(decostand(Y,"total"),"manhattan")/2
      if(method=="percentagedifference") D = vegdist(Y, "bray")
      if(method=="wishart")   D = WishartD(Y)
    } else {
      if(method=="modmeanchardiff") D = D19(Y)
      if(method=="kulczynski")  D = vegdist(Y, "kulczynski")
      if(method=="ab.jaccard")  D = chao(Y, coeff="Jaccard", samp=samp)
      if(method=="ab.sorensen") D = chao(Y, coeff="Sorensen", samp=samp)
      if(method=="ab.ochiai")   D = chao(Y, coeff="Ochiai", samp=samp)
      if(method=="ab.simpson")  D = chao(Y, coeff="Simpson", samp=samp)
    }
    #
    if(sqrt.D) D = sqrt(D)
    SStotal <- sum(D^2)/n      # eq. 8
    BDtotal <- SStotal/(n-1)   # eq. 3
    delta1 <- centre(as.matrix(-0.5*D^2), n)   # eq. 9
    LCBD <- diag(delta1)/SStotal               # eq. 10b
    #
    out <- list(SStotal_BDtotal=c(SStotal,BDtotal), LCBD=LCBD, 
                method=method, D=D)
  }
  ###
  ###
  epsilon <- sqrt(.Machine$double.eps)
  method <- match.arg(method, c("euclidean", "manhattan", "modmeanchardiff", "profiles", "hellinger", "chord", "chisquare", "divergence", "canberra", "whittaker", "percentagedifference", "wishart", "kulczynski", "ab.jaccard", "ab.sorensen","ab.ochiai","ab.simpson","jaccard","sorensen","ochiai","none"))
  #
  if(any(method == c("profiles", "hellinger", "chord", "chisquare", "manhattan", "modmeanchardiff", "divergence", "canberra", "whittaker", "percentagedifference", "kulczynski"))) require(vegan)
  if(any(method == c("jaccard","sorensen","ochiai"))) require(ade4)
  #
  if(is.table(Y)) Y <- Y[1:nrow(Y),1:ncol(Y)]    # In case class(Y) is "table"
  n <- nrow(Y)
  if((n==2)&(dist(Y)[1]<epsilon)) stop("Y contains two identical rows, hence BDtotal = 0")
  #
  aa <- system.time({
    if(any(method == 
           c("euclidean", "profiles", "hellinger", "chord", "chisquare","none"))) {
      note <- "Info -- This coefficient is Euclidean"
      res <- BD.group1(Y, method, save.D, per=FALSE, n)
      #
      # Permutation test for LCBD indices, distances group 1
      if(nperm>0) {
        p <- ncol(Y)
        nGE.L = rep(1,n)
        for(iperm in 1:nperm) {
          Y.perm = apply(Y,2,sample)
          res.p <- BD.group1(Y.perm, method, save.D, per=TRUE, n)
          ge <- which(res.p$LCBD+epsilon >= res$LCBD)
          nGE.L[ge] <- nGE.L[ge] + 1
        }
        p.LCBD <- nGE.L/(nperm+1)
      } else { p.LCBD <- NA }
      #
      if(save.D) { D <- res$D } else { D <- NA }
      #
      out <- list(SStotal_BDtotal=res$SStotal_BDtotal, SCBD=res$SCBD, 
                  LCBD=res$LCBD, p.LCBD=p.LCBD, method=method, note=note, D=D)
      
    } else {
      #
      if(method == "divergence") {
        note = "Info -- This coefficient is Euclidean"
      } else if(any(method == c("jaccard","sorensen","ochiai"))) {
        note = c("Info -- This coefficient is Euclidean because dist.binary ",
                 "of ade4 computes it as sqrt(D). Use beta.div with option sqrt.D=FALSE")
      } else if(any(method == 
                    c("manhattan","canberra","whittaker","percentagedifference","wishart"))) {
        if(sqrt.D) {
          note = "Info -- This coefficient, in the form sqrt(D), is Euclidean"
        } else {
          note = c("Info -- For this coefficient, sqrt(D) would be Euclidean", 
                   "Use is.euclid(D) of ade4 to check Euclideanarity of this D matrix")
        }
      } else {
        note = c("Info -- This coefficient is not Euclidean", 
                 "Use is.euclid(D) of ade4 to check Euclideanarity of this D matrix")
      }
      #
      res <- BD.group2(Y, method, sqrt.D, n)
      #
      # Permutation test for LCBD indices, distances group 2
      if(nperm>0) {
        nGE.L = rep(1,n)
        for(iperm in 1:nperm) {
          Y.perm = apply(Y,2,sample)
          res.p <- BD.group2(Y.perm, method, sqrt.D, n)
          ge <- which(res.p$LCBD+epsilon >= res$LCBD)
          nGE.L[ge] <- nGE.L[ge] + 1
        }
        p.LCBD <- nGE.L/(nperm+1)
      } else { p.LCBD <- NA }
      #
      if(sqrt.D) note.sqrt.D<-"sqrt.D=TRUE"  else  note.sqrt.D<-"sqrt.D=FALSE"
      if(save.D) { D <- res$D } else { D <- NA }
      #
      out <- list(SStotal_BDtotal=res$SStotal_BDtotal, LCBD=res$LCBD,  
                  p.LCBD=p.LCBD, method=c(method,note.sqrt.D), note=note, D=D)
    }
    #
  })
  aa[3] <- sprintf("%2f",aa[3])
  if(clock) cat("Time for computation =",aa[3]," sec\n")
  #
  class(out) <- "beta.div"
  out
}

D11 <- function(Y, algo=1)
  #
  # Compute Clark's coefficient of divergence. 
  # Coefficient D11 in Legendre and Legendre (2012, eq. 7.51).
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, April 2011
{
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  # Prepare to divide by pp = (p-d) = no. species present at both sites
  Y.ap <- 1 - decostand(Y, "pa")
  d <- Y.ap %*% t(Y.ap)
  pp <- p-d   # n. species present at the two compared sites
  #
  if(algo==1) {   # Faster algorithm
    D <- matrix(0, n, n)
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        num <- (Y[i,]-Y[j,])
        den <- (Y[i,]+Y[j,])
        sel <- which(den > 0)
        D[i,j] = sqrt(sum((num[sel]/den[sel])^2)/pp[i,j])
      }
    }
    #
  } else {   # Slower algorithm 
    D <- matrix(0, n, n)
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        temp = 0
        for(p2 in 1:p) {
          den = Y[i,p2] + Y[j,p2]
          if(den > 0) {
            temp = temp + ((Y[i,p2] - Y[j,p2])/den)^2
          }
        }
        D[i,j] = sqrt(temp/pp[i,j])
      }
    }
    #
  }	
  DD <- as.dist(D)
}

D19 <- function(Y)
  #
  # Compute the Modified mean character difference.
  # Coefficient D19 in Legendre and Legendre (2012, eq. 7.46).
  # Division is by pp = number of species present at the two compared sites
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, April 2011
{
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  # Prepare to divide by pp = (p-d) = n. species present at both sites
  Y.ap <- 1 - decostand(Y, "pa")
  d <- Y.ap %*% t(Y.ap)
  pp <- p-d   # n. species present at the two compared sites
  #
  D <- vegdist(Y, "manhattan")
  DD <- as.dist(as.matrix(D)/pp)
}

WishartD <- function(Y)
  #
  # Compute dissimilarity - 1 - Wishart similarity ratio (Wishart 1969).
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, August 2012
{
  CP = crossprod(t(Y))
  SS = apply(Y^2,1,sum)
  n = nrow(Y)
  mat.sq = matrix(0, n, n)
  for(i in 2:n) {
    for(j in 1:(n-1)) { mat.sq[i,j] = CP[i,j]/(SS[i] + SS[j] - CP[i,j]) }
  }
  mat = 1 - as.dist(mat.sq)
}

chao <- function(mat, coeff="Jaccard", samp=TRUE)
  #
  # Compute Chao et al. (2006) abundance-based indices.
  #
  # Arguments -
  # mat = data matrix, species abundances
  # coef = "Jaccard" : modified abundance-based Jaccard index
  #        "Sorensen": modified abundance-based Sørensen index
  #        "Ochiai"  : modified abundance-based Ochiai index
  #        "Simpson" : modified abundance-based Simpson index
  # samp=TRUE : Compute dissimilarities for sample data
  #     =FALSE: Compute dissimilarities for true population data
#
# Details -
# For coeff="Jaccard", the output values are identical to those
# produced by vegan's function vegdist(mat, "chao").
#
# Help received from A. Chao and T. C. Hsieh in July 2012 for the computation  
# of dissimilarities for true population data is gratefully acknowledged.
#
# Reference --
# Chao, A., R. L. Chazdon, R. K. Colwell and T. J. Shen. 2006. 
# Abundance-based similarity indices and their estimation when there 
# are unseen species in samples. Biometrics 62: 361–371.
#
# License: GPL-2 
# Author:: Pierre Legendre, July 2012
{
  require(vegan)
  nn = nrow(mat)
  res = matrix(0,nn,nn)
  if(samp) {   # First for sample data
    for(k in 2:nn) {
      for(j in 1:(k-1)) {
        #cat("k =",k,"  j =",j,"\n")
        v1 = mat[j,]   # Vector 1
        v2 = mat[k,]   # Vector 2
        v1.pa = decostand(v1,"pa")   # Vector 1 in presence-absence form
        v2.pa = decostand(v2,"pa")   # Vector 2 in presence-absence form
        N.j = sum(v1)   # Sum of abundances in vector 1
        N.k = sum(v2)   # Sum of abundances in vector 2
        shared.sp = v1.pa * v2.pa   # Vector of shared species ("pa")
        if(sum(shared.sp) == 0) { 
          res[k,j] = 1
        } else {
          C.j = sum(shared.sp * v1)   # Sum of shared sp. abundances in v1
          C.k = sum(shared.sp * v2)   # Sum of shared sp. abundances in v2
          # a1.j = sum(shared.sp * v1.pa)
          # a1.k = sum(shared.sp * v2.pa)
          a1.j = length(which((shared.sp * v2) == 1)) # Singletons in v2
          a1.k = length(which((shared.sp * v1) == 1)) # Singletons in v1
          a2.j = length(which((shared.sp * v2) == 2)) # Doubletons in v2
          if(a2.j == 0) a2.j <- 1
          a2.k = length(which((shared.sp * v1) == 2)) # Doubletons in v1
          if(a2.k == 0) a2.k <- 1
          # S.j = sum(v1[which(v2 == 1)]) # Sum abund. in v1 for singletons in v2
          # S.k = sum(v2[which(v1 == 1)]) # Sum abund. in v2 for singletons in v1
          sel2 = which(v2 == 1)
          sel1 = which(v1 == 1)
          if(length(sel2)>0) S.j = sum(v1[sel2]) else S.j = 0
          if(length(sel1)>0) S.k = sum(v2[sel1]) else S.k = 0
          
          U.j = (C.j/N.j) + ((N.k-1)/N.k) * (a1.j/(2*a2.j)) * (S.j/N.j) # Eq. 11
          if(U.j > 1) U.j <- 1
          U.k = (C.k/N.k) + ((N.j-1)/N.j) * (a1.k/(2*a2.k)) * (S.k/N.k) # Eq. 12
          if(U.k > 1) U.k <- 1
          
          if(coeff == "Jaccard") {                     # "Jaccard"
            res[k,j] = 1 - (U.j*U.k/(U.j + U.k - U.j*U.k))
          } else if(coeff == "Sorensen") {         # "Sorensen"
            res[k,j] = 1 - (2*U.j*U.k/(U.j + U.k))
          } else if(coeff == "Ochiai") {           # "Ochiai"
            res[k,j] = 1 - (sqrt(U.j*U.k))
          } else if(coeff == "Simpson") { 
            # Simpson (1943), or Lennon et al. (2001) in Chao et al. (2006)
            res[k,j] = 1 -
              (U.j*U.k/(U.j*U.k+min((U.j-U.j*U.k),(U.k-U.j*U.k))))
          } else { # 
            stop("Incorrect coefficient name")
          }
        }
      }
    }
    
  } else {   # Now for complete population data
    
    for(k in 2:nn) {
      for(j in 1:(k-1)) {
        v1 = mat[j,]   # Vector 1
        v2 = mat[k,]   # Vector 2
        v1.pa = decostand(v1,"pa")   # Vector 1 in presence-absence form
        v2.pa = decostand(v2,"pa")   # Vector 2 in presence-absence form
        shared.sp = v1.pa * v2.pa    # Vector of shared species ("pa")
        if(sum(shared.sp) == 0) { 
          res[k,j] = 1
        } else {
          N1 = sum(v1)   # Sum of abundances in vector 1
          N2 = sum(v2)   # Sum of abundances in vector 2
          U = sum(shared.sp * v1)/N1   # Sum of shared sp. abundances in v1
          V = sum(shared.sp * v2)/N2   # Sum of shared sp. abundances in v2
          
          if(coeff == "Jaccard") {                     # "Jaccard"
            res[k,j] = 1 - (U*V/(U + V - U*V))
          } else if(coeff == "Sorensen") {         # "Sorensen"
            res[k,j] = 1 - (2*U*V/(U + V))
          } else if(coeff == "Ochiai") {           # "Ochiai"
            res[k,j] = 1 - (sqrt(U*V))
          } else if(coeff == "Simpson") { # "Simpson"
            res[k,j] = 1 - (U*V/(U*V+min((U-U*V),(V-U*V)))) # Eq. ?
          } else { # 
            stop("Incorrect coefficient name")
          }
        }
      }
    }
  }
  res <- as.dist(res)
}

######## End of beta.div function
##################################################################################################

####BETA.DIV.COMP SOURCE FUNCTION- for partitioning beta diversity into components (spatial)
##################################################################################################
beta.div.comp <- function(mat, coef="J", quant=FALSE, save.abc=FALSE)
  #
  # Description --
  # 
  # Podani-family and Baselga-family decompositions of the Jaccard and Sørensen 
  # dissimilarity coefficients into replacement and richness difference 
  # components, for species presence-absence or abundance data, as described  
  # in Legendre (2014).
  #
  # Usage --
  #
  # beta.div.comp(mat, coef="J", quant=FALSE, save.abc=FALSE)
#
# Arguments --
#
# mat : Data in matrix or data.frame form.
# coef : Family of coefficients to be computed --
#        "S" or "Sorensen": Podani family, Sørensen-based indices
#        "J" or "Jaccard" : Podani family, Jaccard-based indices
#        "BS" : Baselga family, Sørensen-based indices
#        "BJ" : Baselga family, Sørensen-based indices
#        "N" : Podani & Schmera (2011) relativized nestedness index.
#        The quantitative form in Sørensen family is the percentage difference.
#        The quantitative form in the Jaccard family is the Ruzicka index.
#
# quant=TRUE : Compute the quantitative form of replacement, nestedness and D.
#      =FALSE: Compute the presence-absence form of the coefficients.
# save.abc=TRUE : Save the matrices of parameters a, b and c used in the
#      presence-absence calculations.
#
# Details --
#
#    For species presence-absence data, the distance coefficients are 
# Jaccard=(b+c)/(a+b+c) and Sørensen=(b+c)/(2*a+b+c) with usual abc notation.
#
#    For species abundance data, the distance coefficients are 
# the Ruzicka index = (B+C)/(A+B+C) and Odum's percentage difference 
# (incorrectly called Bray-Curtis) = (B+C)/(2A+B+C), where  
# A = sum of the intersections (or minima) of species abundances at two sites,
# B = sum at site 1 minus A, C = sum at site 2 minus A.
#
#    The binary (quant=FALSE) and quantitative (quant=TRUE) forms of the S and  
# J indices return the same values when computed for presence-absence data.
#
# Value --
#
# repl : Replacement matrix, class = 'dist'.
# rich : Richness/abundance difference or nestedness matrix, class = 'dist'.
#        With options "BJ", "BS" and "N", 'rich' contains nestedness indices.
#        With option "N", the 'repl' and 'rich' values do not add up to 'D'.
# D    : Dissimilarity matrix, class = 'dist'.
# part : Beta diversity partitioning -- 
#        1. Total beta div. = sum(D.ij)/(n*(n-1)) (Legendre & De Cáceres 2013)
#        2. Total replacement diversity 
#        3. Total richness difference diversity (or nestedness)
#        4. Total replacement div./Total beta div.
#        5. Total richness difference div. (or nestedness)/Total beta div.
# Note : Name of the dissimilarity coefficient.
#
# References --
#
# Baselga, A. (2010) Partitioning the turnover and nestedness components of beta 
# diversity. Global Ecology and Biogeography, 19, 134–143.
#
# Baselga, A. (2012) The relationship between species replacement, dissimilarity 
# derived from nestedness, and nestedness. Global Ecology and Biogeography, 21, 
# 1223–1232. 
#
# Baselga, A. (2013) Separating the two components of abundance-based 
# dissimilarity: balanced changes in abundance vs. abundance gradients. Methods 
# in Ecology and Evolution, 4, 552–557.
#
# Carvalho, J.C., Cardoso, P., Borges, P.A.V., Schmera, D. & Podani, J. (2013)
# Measuring fractions of beta diversity and their relationships to nestedness: 
# a theoretical and empirical comparison of novel approaches. Oikos, 122, 
# 825–834.
#
# Legendre, P. 2014. Interpreting the replacement and richness difference   
# components of beta diversity. Global Ecology and Biogeography 23: (in press).
#
# Podani, J., Ricotta, C. & Schmera, D. (2013) A general framework for analyzing 
# beta diversity, nestedness and related community-level phenomena based on 
# abundance data. Ecological Complexity, 15, 52-61.
#
# Podani, J. & Schmera, D. 2011. A new conceptual and methodological framework 
# for exploring and explaining pattern in presence-absence data. Oikos, 120, 
# 1625–1638.
#
# License: GPL-2 
# Author:: Pierre Legendre
{
  coef <- pmatch(coef, c("S", "J", "BS", "BJ", "N"))
  if(coef==5 & quant) stop("coef='N' and quant=TRUE: combination not programmed")
  mat <- as.matrix(mat)
  n <- nrow(mat)
  if(is.null(rownames(mat))) noms <- paste("Site",1:n,sep="")
  else noms <- rownames(mat)
  #
  if(!quant) {      # Binary data provided, or make the data binary
    if(coef==1) form="Podani family, Sorensen" 
    if(coef==2) form="Podani family, Jaccard"
    if(coef==3) form="Baselga family, Sorensen" 
    if(coef==4) form="Baselga family, Jaccard"
    if(coef==5) form="Podani & Schmera (2011) relativized nestedness"
    mat.b <- ifelse(mat>0, 1, 0)
    a <- mat.b %*% t(mat.b)
    b <- mat.b %*% (1 - t(mat.b))
    c <- (1 - mat.b) %*% t(mat.b)
    min.bc <- pmin(b,c)
    #
    if(coef==1 || coef==2) {
      repl <- 2*min.bc   # replacement, turnover, beta-3
      rich <- abs(b-c)   # nestedness, richness diff., beta-rich
      #
      # Add the denominators
      if(coef==1) {                # Sørensen-based components
        repl <- repl/(2*a+b+c)
        rich <- rich/(2*a+b+c)
        D <- (b+c)/(2*a+b+c)
      } else if(coef==2) {     # Jaccard-based components
        repl <- repl/(a+b+c)
        rich <- rich/(a+b+c)
        D <- (b+c)/(a+b+c)
      }
    } else if(coef==3) {     # Baselga 2010 components based on Sørensen
      D <- (b+c)/(2*a+b+c)             # Sørensen dissimilarity
      repl <- min.bc/(a+min.bc)        # replacement, turnover
      rich <- D-repl                   # richness difference
      
    } else if(coef==4) {      # Baselga 2012 components based on Jaccard
      D <- (b+c)/(a+b+c)               # Jaccard dissimilarity
      repl <- 2*min.bc/(a+2*min.bc)    # replacement, turnover
      rich <- D-repl                   # richness difference
    } else if(coef==5) {      # rich = Podani N = nestdness based on Jaccard
      repl <- 2*min.bc/(a+b+c)
      D <- (b+c)/(a+b+c)
      rich <- matrix(0,n,n)
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          aa = a[i,j]; bb = b[i,j]; cc = c[i,j]
          if(a[i,j] == 0)  rich[i,j] <- 0  
          else  rich[i,j] <- (aa + abs(bb-cc))/(aa+bb+cc) 
        }
      }
    }
    
    rownames(repl) <- rownames(rich) <- rownames(D) <- noms
    D <- as.dist(D)
    repl <- as.dist(repl)
    rich <- as.dist(rich)
    total.div <- sum(D)/(n*(n-1))
    repl.div <- sum(repl)/(n*(n-1))
    rich.div <- sum(rich)/(n*(n-1))
    part <- c(total.div,repl.div,rich.div,repl.div/total.div,rich.div/total.div)
    #
    if(save.abc) {
      res <- list(repl=repl, rich=rich, D=D, part=part, Note=form, 
                  a=as.dist(a), b=as.dist(b), c=as.dist(c))
    } else { 
      res <- list(repl=repl, rich=rich, D=D, part=part, Note=form)
    }
    #
  } else {      # Quantitative data
    # Calculations based on individuals.within.species
    if(coef==1) form<-"Podani family, percentage difference" 
    if(coef==2) form<-"Podani family, Ruzicka"
    if(coef==3) form<-"Baselga family, percentage difference"
    if(coef==4) form<-"Baselga family, Ruzicka"
    # Baselga (2013) notation:
    # A = W = sum of minima in among-site comparisons
    # B = site.1 sum - W = K.1 - W
    # C = site.2 sum - W = K.2 - W
    K <- vector("numeric", n)   # site (row) sums
    W <- matrix(0,n,n)
    repl <- matrix(0,n,n)
    rich <- matrix(0,n,n)
    D <- matrix(0,n,n)
    rownames(repl) <- rownames(rich) <- rownames(D) <- noms
    K <- apply(mat,1,sum)         # Row sums
    for(i in 2:n) for(j in 1:(i-1)) W[i,j] <- sum(pmin(mat[i,], mat[j,]))
    #
    # Quantitative extensions of the S and J decompositions
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        repl[i,j] <- 2*(min(K[i],K[j])-W[i,j]) # 2*min(B,C)
        rich[i,j] <- abs(K[i]-K[j])            # abs(B-C)
      }
    }
    #
    # Add the denominators
    if(coef==1) {         # Sørensen-based (% difference) components
      for(i in 2:n) {
        for(j in 1:(i-1)) {	                        # Baselga 2013 notation:
          repl[i,j] <- repl[i,j]/(K[i]+K[j])          # 2min(B,C)/(2A+B+C)
          rich[i,j] <- rich[i,j]/(K[i]+K[j])          # abs(B-C)/(2A+B+C)
          # cat(K[i], K[j], W[i,j],"\n")
          D[i,j] <- (K[i]+K[j]-2*W[i,j])/(K[i]+K[j])  # (B+C)/(2A+B+C)
        }
      }
    } else if(coef==2) {    # Jaccard-based (Ruzicka) components
      for(i in 2:n) {
        for(j in 1:(i-1)) {                         # Baselga 2013 notation:
          repl[i,j] <- repl[i,j]/(K[i]+K[j]-W[i,j])   # 2min(B,C)/(A+B+C)
          rich[i,j] <- rich[i,j]/(K[i]+K[j]-W[i,j])   # abs(B-C)/(A+B+C)
          # cat(K[i], K[j], W[i,j],"\n")
          D[i,j]<-(K[i]+K[j]-2*W[i,j])/(K[i]+K[j]-W[i,j]) # (B+C)/(A+B+C)
        }
      }
    }
    #
    # Baselga (2013): quantitative extensions of the Baselga (2010) indices
    if(coef==3) {   # Baselga (2013) indices decomposing percentage difference
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          repl[i,j] <- (min(K[i],K[j])-W[i,j])/min(K[i],K[j])
          rich[i,j] <- abs(K[i]-K[j])*W[i,j]/((K[i]+K[j])*min(K[i],K[j]))
          # cat(K[i], K[j], W[i,j],"\n")
          D[i,j] <- (K[i]+K[j]-2*W[i,j])/(K[i]+K[j])
        }
      }
    }	
    if(coef==4) {   # Decomposing Ruzicka in the spirit of Baselga 2013
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          repl[i,j] <- 
            2*(min(K[i],K[j])-W[i,j])/(2*min(K[i],K[j])-W[i,j])
          rich[i,j] <- abs(K[i]-K[j])*W[i,j]/
            ((K[i]+K[j]-W[i,j])*(2*min(K[i],K[j])-W[i,j]))
          # cat(K[i], K[j], W[i,j],"\n")
          D[i,j] <- (K[i]+K[j]-2*W[i,j])/(K[i]+K[j]-W[i,j])
        }
      }
    }	
    #
    repl <- as.dist(repl)
    rich <- as.dist(rich)
    D <- as.dist(D)
    repl.div <- sum(repl)/(n*(n-1))
    rich.div <- sum(rich)/(n*(n-1))
    total.div <- sum(D)/(n*(n-1))
    part <- c(total.div,repl.div,rich.div,repl.div/total.div,rich.div/total.div)
    #
    res <- list(repl=repl, rich=rich, D=D, part=part, Note=form)
  }
  res
}
########End of beta.div.comp function
##################################################################################################

####LCBD COMP SOURCE FUNCTION- partition LCBD into components (spatial)
##################################################################################################
##LCBD.comp() source function
LCBD.comp <- function(x, sqrt.x=TRUE)
{
  ### Internal function
  centre <- function(D,n)
    # Centre a square matrix D by matrix algebra
    # mat.cen = (I - 11'/n) D (I - 11'/n)
  {
    One <- matrix(1,n,n)
    mat <- diag(n) - One/n
    mat.cen <- mat %*% D %*% mat
  }
  ###
  n <- nrow(as.matrix(x))
  if(sqrt.x) {
    # x = sqrt(x)
    SStotal <- sum(x)/n # eq. 8
    BDtotal <- SStotal/(n-1) # eq. 3
    G <- centre(as.matrix(-0.5*x), n) # Gower-centred matrix
  } else {
    SStotal <- sum(x^2)/n # eq. 8
    BDtotal <- SStotal/(n-1) # eq. 3
    G <- centre(as.matrix(-0.5*x^2), n) # Gower-centred matrix
  }
  LCBD <- diag(G)/SStotal # Legendre & De Caceres (2013), eq. 10b
  out <- list(SStotal_BDtotal=c(SStotal,BDtotal), LCBD=LCBD, D=x)
} 

# Arguments --
#
# x : D or beta diversity component matrix, class=dist.
# sqrt.x : Take sqrt() of components before computing LCBD.comp. Use
# sqrt.x=TRUE for the replacement and richness/abundance difference indices
# computed by beta.div.comp(), as well as for the corresponding D matrices.

########End of LCBD.comp function 
##################################################################################################

####DECOMPOSE.D2 TEMPORAL BD SOURCE FUNCTION- for computing temporal beta diversity and gain and loss components
##################################################################################################
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
##################################################################################################

####PAIRED.DIFF2_2.R TEMPORAL EXCEPTIONAL SITES FUNCTION (with the BDC matrix)- calculating significant temporal change
##################################################################################################
paired.diff2 <- function(mat1,mat2,method="hellinger", pa.tr=FALSE, nperm=99, permute.sp=1, BCD=TRUE, replace=FALSE, clock=FALSE)
  #
  # Compute and test differences between pairs of data vectors of observations at T1 and T2.
  #
  # Test hypothesis (H0) that an object is not exceptionally different between T1 and T2.
  # Example in palaeoecology: ancient and modern diatom communities in sediment cores.
  # Example in sequence data: "regions" before and after treatment.
  # 
  # Arguments --
  # mat1, mat2: two matrices (or data.frames) with the same number of rows and columns. 
  # The rows must correspond to the same objects (e.g. sites) and the colums to the same 
  # variables (e.g. species).
# method={"hellinger", "chord", "ruzicka", "%difference", "euclidean"}. 
#   Methods {"hellinger", "chord"} are obtained by transformation of the  
#     species data followed by calculation of the Euclidean distance. These distances 
#     have the Euclidean property. 
#     If pa.tr=TRUE, sqrt(2)*sqrt(1-Ochiai) is computed.
#   Methods {"ruzicka", "%difference"} are obtained by computing a dissimilarity function. 
#     It is recommended to take the square root of these dissimilarities before computing 
#     ordinations by principal coordinate analysis. However, that precaution is not 
#     important here; the results of the permutation tests will be the same for these
#     dissimilarities square-rooted or not.
#     If pa.tr=TRUE, either the Jaccard or the Sørensen coefficient is computed.
# pa.tr=FALSE: do NOT transform the data to presence-absence.
#      =TRUE : transform the data to binary (i.e. presence-absence) form.
# nperm = number of permutations for the permutation test.
##
# This version of the function contains three permutation methods --
# permute.sp=1 : permute data separately in each column, both matrices in the same way.
#           =2 : permute data separately in each column. Do not force the permutations to 
#  			 start at the same point in the two matrices.
#           =3 : permute entire rows in each matrix separately (suggestion D. Borcard).
##
# BCD=TRUE  : Compute and save the B and C components of the %difference and Ruzicka D.
#             For the %difference, they are expressed as B/(2A+B+C) and C/(2A+B+C).
#             For the Ruzicka D, they are expressed as B/(A+B+C) and C/(A+B+C).
# BCD=FALSE : Do not compute the components. BCD=FALSE for D other than %diff and Ruzicka.
# replace=FALSE : sampling without replacement for regular permutation test.
#        =TRUE  : sampling with replacement. The testing method is then bootstrapping.
# clock=FALSE : Do not print the computation time.
#      =TRUE  : Print the time (in sec) used for computation.
#
# Details --
# H0: in each matrix (e.g. each time), the sites do not differ in species composition. 
#     They only differ by random sampling of each species' statistical population.
# H1: Some sites are exceptionally different between T1 and T2.
# 
# The randomization procedures are the following:
# 1. In each matrix, the original values (e.g. species abundances) are permuted at random, 
# independently in each column. Permutation of the two matrices is started with the same 
# random seed, so that the values in each column (e.g. species) are permuted in the same 
# way in mat1.perm and mat2.perm. 
# 2. The transformation, if any, is recomputed on the permuted data matrices. This is 
# necessary to make sure that the permuted data are transformed in the same way as the 
# initial data, with row sums or row lengths of 1. In this way, the D of the permuted data 
# will be comparable to the reference D.
# 3. The distances between T1 and T2 are recomputed, for each site separately.
#
# For presence-absence data, this function computes the binary forms of the quantitative 
# coefficients listed under the 'method' parameter. The "hellinger" and "chord" 
# transformations produce the Ochiai distance, or more precisely: 
# D.Hellinger = D.chord = sqrt(2) * sqrt(1 - S.Ochiai) 
# where "S.Ochiai" designates the Ochiai similarity coefficient.  
# The "%difference" dissimilarity produces (1 – S.Sørensen) 
# whereas the "ruzicka" dissimilarity produces (1 – S.Jaccard).
#
# Community composition data could be log-transformed prior to analysis. Only the 
# Euclidean distance option should be used with log-transformed data. It is meaningless to 
# subject log-transformed data to the {"hellinger", "chord"} transformations 
# available in this function. - One can use either the log(y+1 transformation (log1p() 
# function of {base}), or Anderson et al. (2006) special log transformation available in 
# {vegan}: decostand(mat, "log", logbase=10).
#
# Value --
# A list containing the vector of distances between T1 and T2 for each object and a 
# corresponding vector of p-values. The significant p-values (e.g. p.dist ≤ 0.05) indicate 
# exceptional objects for the difference of their species composition. The p-values should be corrected for multiple testing using function p.adjust() of {stats}. A good general choice is method="holm", which is the default option of the function.
# An output table containing B, C and D.
#
# Author:: Pierre Legendre
# License: GPL-2 
{
  ### Internal functions
  dissim <- function(mat1, mat2, n, method, tr=TRUE, BCD, ref)
  {
    vecD = vector(mode="numeric",length=n)
    if(BCD) { 
      vecB = vector(mode="numeric",length=n)
      vecC = vector(mode="numeric",length=n)
      vecD = vector(mode="numeric",length=n)
    } else { vecB=NA; vecC=NA; vecD=NA }
    #
    # Compute the dissimilarity between T1 and T2 for each object (site)
    # 1. If method = {"hellinger", "chord"}, tr is TRUE
    if(tr) for(i in 1:n) vecD[i] = dist(rbind(mat1[i,], mat2[i,]))
    #
    # 2. Compute the Euclidean distance
    if(method == "euclidean")  
      for(i in 1:n) vecD[i] = dist(rbind(mat1[i,], mat2[i,])) 
      # 3. Compute the Ruzicka or %difference dissimilarity 
      if(method == "ruzicka") dissimil=1       # Quantitative form of Jaccard
      if(method == "%difference") dissimil=2   # Quantitative form of Sørensen
      if(any(method == c("ruzicka", "%difference"))) { 
        for(i in 1:n) {
          tmp = RuzickaD(mat1[i,], mat2[i,], method=method, BCD=BCD, ref=ref) 
          if(BCD) {
            vecB[i] <- tmp$B
            vecC[i] <- tmp$C }
          vecD[i] <- tmp$D
        }
      }
      # Alternative method (not used here) to compute the %difference dissimilarity:
      #	for(i in 1:n) vecD[i] = vegdist(rbind(mat1[i,], mat2[i,]), "bray")         #Slower
      list(vecB=vecB, vecC=vecC, vecD=vecD)
  }
  ###
  transform <- function(mat, method)
  {
    if(method=="hellinger") mat = decostand(mat, "hellinger")
    if(method=="chord")     mat = decostand(mat, "norm")
    mat
  }
  ### End internal functions
  ###
  A <- system.time({
    
    epsilon <- sqrt(.Machine$double.eps)
    method <- match.arg(method, c("euclidean", "hellinger", "chord", "ruzicka", "%difference")) 
    n = nrow(mat1)
    p = ncol(mat1)
    if((nrow(mat2)!=n) | (ncol(mat2)!=p)) stop("The matrices are not of the same size.")
    #
    if(pa.tr) {
      mat1 <- ifelse(mat1>0, 1, 0)
      mat2 <- ifelse(mat2>0, 1, 0) }
    if(any(method == c("hellinger", "chord"))) {
      tr <- TRUE
      require(vegan)
    } else { tr <- FALSE }
    if( (any(method == c("ruzicka", "%difference"))) & BCD) { 
      BCD.mat <- matrix(0,n,3)
      if(method=="ruzicka")    colnames(BCD.mat) <- 
          c("B/(A+B+C)","C/(A+B+C)","D=(B+C)/(A+B+C)")
      if(method=="%difference") colnames(BCD.mat) <- 
          c("B/(2A+B+C)","C/(2A+B+C)","D=(B+C)/(2A+B+C)")
      rownames(BCD.mat) <- paste("Obj",1:n,sep=".")
    } else {
      BCD <- FALSE 
      BCD.mat <- NA }
    ###
    # 1. Compute the reference D for each object from corresponding vectors in the 2 matrices.
    if(tr) { 
      tmp <-dissim(transform(mat1,method), transform(mat2,method),n,method,tr,BCD,ref=FALSE)
    } else { tmp <- dissim(mat1, mat2, n, method, tr, BCD, ref=TRUE) }
    vecD.ref <- tmp$vecD
    if(BCD) { BCD.mat[,1]<-tmp$vecB ; BCD.mat[,2]<-tmp$vecC ; BCD.mat[,3]<-tmp$vecD }
    ###
    if(permute.sp!=3) {   # Permute the data separately in each column.
      # 2. Permutation methods 1 and 2 --
      # Permute *the raw data* by columns. Permute the two matrices in the same way, saving the seed before the two sets of permutations through sample(). 
      # Permutation test for each distance in vector D
      # seed: seed for random number generator, used by the permutation function 
      #       sample(). It is reset to that same value before permuting the values in the  
      #       columns of the second matrix. 
      if(nperm>0) {
        nGE.D = rep(1,n)
        for(iperm in 1:nperm) {
          BCD <- FALSE
          if(permute.sp==1) {    # Permutation methods 1
            seed <- ceiling(runif(1,max=100000))
            # cat("seed =",seed,'\n')
            set.seed(seed)
            mat1.perm <- apply(mat1,2,sample)
            set.seed(seed)
            mat2.perm <- apply(mat2,2,sample)
          } else {  # Permutation methods 2 - Do not force the permutations 
            # to start at the same point in the two matrices.
            mat1.perm <- apply(mat1,2,sample)
            mat2.perm <- apply(mat2,2,sample)
          }
          # 3. Recompute transformations of the matrices and the D values of the paired vectors.
          if(tr) { tmp <- dissim(transform(mat1.perm,method), 
                                 transform(mat2.perm,method), n, method, tr, BCD, ref=FALSE)
          } else { tmp <- dissim(mat1.perm, mat2.perm, n, method, tr, BCD, ref=FALSE) }
          vecD.perm <- tmp$vecD
          ge <- which(vecD.perm+epsilon >= vecD.ref)
          nGE.D[ge] <- nGE.D[ge] + 1
        }
        # 4. Compute the p-value associated with each distance.
        p.dist <- nGE.D/(nperm+1)
      } else { p.dist <- NA }   # if nperm=0
      
    } else if(permute.sp==3) {   
      # 2.bis  Permutation method 3 -- 
      # Permute entire rows in each matrix separately.
      if(nperm>0) {
        seed <- ceiling(runif(1,max=100000))
        set.seed(seed)
        nGE.D = rep(1,n)
        for(iperm in 1:nperm) {
          BCD <- FALSE
          mat1.perm <- mat1[sample(n),]
          mat2.perm <- mat2[sample(n),]
          #
          # 3.bis Recompute the D values of the paired vectors.
          if(tr) { tmp <- dissim(transform(mat1.perm,method), 
                                 transform(mat2.perm,method), n, method, tr, BCD, ref=FALSE)
          } else { tmp <- dissim(mat1.perm, mat2.perm, n, method, tr, BCD, ref=FALSE) }
          vecD.perm <- tmp$vecD
          ge <- which(vecD.perm+epsilon >= vecD.ref)
          nGE.D[ge] <- nGE.D[ge] + 1
        }
        # 4.bis Compute the p-value associated with each distance.
        p.dist <- nGE.D/(nperm+1)
      } else { p.dist <- NA }   # if nperm=0
    }
    p.adj <- p.adjust(p.dist,"holm")
  })
  A[3] <- sprintf("%2f",A[3])
  if(clock) cat("Computation time =",A[3]," sec",'\n')
  #
  list(vecD.ref=vecD.ref, p.dist=p.dist, p.adj=p.adj, BCD.mat=BCD.mat)
}

RuzickaD <- function(vec1, vec2, method="ruzicka", BCD=FALSE, ref=TRUE)
  #
  # Compute the Ruzicka dissimilarity (quantitative form of the Jaccard dissimilarity)
  # or the percentage difference (quantitative form of the Sørensen dissimilarity).
  # A single dissimilarity is computed because there are only two vectors in this function.
  #
  # Arguments --
  # vec1, vec2 : data vectors (species abundance or presence-absence data)
  # method == c("ruzicka", "%difference")
  # BCD=TRUE  : Compute and save the B and C components of the %difference and Ruzicka D.
  #             For the %difference, they are B/(2A+B+C), C/(2A+B+C), D/(2A+B+C).
  #             For the Ruzicka D, they are B/(A+B+C), C/(A+B+C), D/(A+B+C).
# BCD=FALSE : Do not compute the components. BCD=FALSE for D other than %diff and Ruzicka.
# ref=TRUE  : Compute the reference values of D, B and C
#    =FALSE : Under permutation, compute only the value of D. Use separate code (shorter).
#
# License: GPL-2 
# Author:: Pierre Legendre, April 2015
{
  # An algorithm applicable to matrices Y containing two data vectors only
  #
  A <- sum(pmin(vec1, vec2))          # A = sum of minima from comparison of the 2 vectors
  sum.Y <- sum(vec1, vec2)            # Sum of all values in the two vectors, (2A+B+C)
  #
  if(ref) {    # Compute the reference values of statistics D, B and C
    tmp = vec1 - vec2
    B = sum(tmp[tmp>0])                 # Sum of the species losses between T1 and T2
    C = -sum(tmp[tmp<0])                # Sum of the species gains between T1 and T2
    D = B+C                             # Dissimilarity
    
    # Under permutation, compute only the value of D. - Shorter computation time.
  } else { 
    D <- sum.Y-2*A                      # (B+C)
  }
  # Compute the denominator (den) of the Ruzicka or %difference index
  if(method == "ruzicka") { den <-(sum.Y-A)  # den = (A+B+C)
  } else { den <- sum.Y }                # den = (2A+B+C)
  if(!BCD) { B <- NA ; C <- NA }
  list(B=B/den, C=C/den, D=D/den)
}

# Examples -- 
# data(mite)
# res1 = paired.diff(mite[1:10,],mite[61:70,],method="hellinger",nperm=999,permute.sp=1)
# Computation time = 1.971000  sec 
# res2 = paired.diff(mite[1:10,],mite[61:70,],method="hellinger",nperm=999,permute.sp=3)

########End of paired.diff2 function 
##################################################################################################


##################################################################################################
####RAREFYING CLADOCERAN SPECIES RICHNESS                                                        #
##################################################################################################

####Use count data to calculate rarefied species richness for each sample. 
#Note: count data should be verified that it represents most accurate counting/taxonomy scheme.
#Done. 

clad<- read.csv(file.choose()) #schefferville_clad_bd_JULY2015.csv
#File path: C:\Users\Winegardner\Documents\MCGILL\PhD chapters and projects\Schefferville general\Analysis data\Beta diversity data
#Count matrix of Cladocerans from schefferville_lmerdata_Mar2015.xls
#Counts rounded to nearest whole integer.
#Alona OTU2 re-classified as Alona spp.
#Bosmina lumped together into Bosmina spp. 

##Column descriptions
#Sample_ID_Master - Sample ID including core, interval # etc. 
#Sample_ID_T_Order - Numbers the intervals starting at 0 for oldest interval, includes lake code (e.g. DAR_T0, DAR_T1) (for use if comparing 1 interval to the next)
#Counter - Either Katherine Velghe (KV) or Natasha Salter (NS)
#Lake - Dauriat, Knob, Dolly or Denault
#Est_year - Estimated year based on final age models 
#Time_Period - Pre-mining (up to 1953), Mining (1954-1982), Post-mining (1983-present)
#Time_period_grouping - T0_PreM, T1_Min, T2_PosM, for grouping into 3 time periods and using each time period as an individual matrix --> for temporal beta
#Species columns, arranged alphabetically

##Work with DAR and KB only for the time being (Dolly and Denault only have top/bottom samples)
#DAR
clad.dar<- as.data.frame(subset(clad, Lake == "Dauriat ", drop=T))

#KB
clad.kb<- as.data.frame(subset(clad, Lake == "Knob", drop=T))

##Calculate rarefied richness based on the individual lakes (not pooled together)
#DAR
row.names(clad.dar) <- clad.dar[,1]
Srar.dar <- rarefy(clad.dar[,8:42], min(rowSums(clad.dar[,8:42])))
Srar.dar

#Save as dataframe and bind to sample ID column (with lake, year, time period)
Srar.dar<- as.data.frame(Srar.dar)

#KB
row.names(clad.kb) <- clad.kb[,1]
Srar.kb <- rarefy(clad.kb[,8:42], min(rowSums(clad.kb[,8:42])))
Srar.kb

#Save as dataframe and bind to sample ID column (with lake, year, time period)
Srar.kb<- as.data.frame(Srar.kb)

##################################################################################################


##################################################################################################
####ALPHA DIVERSITY (Shannon and Simpson)                                                        #
##################################################################################################

####Calculate alpha diversity- Shannon Weiner and Simpson (evenness) for each intervaal####
##Shannon
#DAR
Shannon.dar<- diversity(clad.dar[,8:42], index = "shannon")
Shannon.dar<- as.data.frame(Shannon.dar)

#KB
Shannon.kb<- diversity(clad.kb[,8:42], index = "shannon")
Shannon.kb<- as.data.frame(Shannon.kb)

##Simpson
#DAR
Simpson.dar<- diversity(clad.dar[,8:42], index = "simpson")
Simpson.dar<- as.data.frame(Simpson.dar)

#KB
Simpson.kb<- diversity(clad.kb[,8:42], index = "simpson")
Simpson.kb<- as.data.frame(Simpson.kb)

##Bind both alpha diversity metrics to dataframe with Srar
#Should now have a matrix with sample/lake information columns as well as Srar, Shannon, Simpson 

#DAR
clad.div.dar<- as.data.frame(cbind(clad.dar, Srar.dar, Shannon.dar, Simpson.dar))

#KB
clad.div.kb<- as.data.frame(cbind(clad.kb, Srar.kb, Shannon.kb, Simpson.kb))

##################################################################################################


##################################################################################################
####TEMPORAL BETA DIVERSITY                                                                      #
##################################################################################################

####Calculate temporal beta diversity for intervals between the different time intervals####

##Dauriat
#So will have to make "matrices" for each individual time sample or better to do in periods? 
#Then compare one time period to next. 
#Calculate total temporal BD for each time period comparison.
#Calculate species gain and loss components of temporal beta for each time period comparison. 
#Calculate which time period comparisons show significant beta (i.e. significant temporal change)

#Subset DAR into the three "Time_period_grouping" 
clad.dar.T0<- as.data.frame(subset(clad.dar, Time_period_grouping == "T0_PreM", drop=T))
clad.dar.T1<- as.data.frame(subset(clad.dar, Time_period_grouping == "T1_Min", drop=T))
clad.dar.T2<- as.data.frame(subset(clad.dar, Time_period_grouping == "T2_PosM", drop=T))

#Melt into long format, so can use ply to get averages 
#T0
clad.dar.T0long<- melt(clad.dar.T0, id.vars = c("Sample_ID_Master", "Sample_ID_T_order", "Counter", "Lake", "Est_year", "Time_period", "Time_period_grouping"))
colnames(clad.dar.T0long) [8]<- 'Taxa'
colnames(clad.dar.T0long) [9]<- 'Count'

clad.darT0.mean<- ddply(clad.dar.T0long, "Taxa", summarize, mean_Count = mean(Count, na.rm=T))
rownames(clad.darT0.mean)<- as.character(clad.darT0.mean[,1])

clad.darT0mean.wide<-as.data.frame(t(clad.darT0.mean)) #but have taxa names as a row- just remove? 
##CONTINUE HERE. 

#T1



#T2






##Rough
#Create an average assemblage for each Time_period_grouping 
clad.darT0.mean<- as.data.frame(colMeans(clad.dar.T0[,8:42]))
colnames(clad.darT0.mean) [1]<- 'Mean_count'

clad.darT1.mean<- as.data.frame(colMeans(clad.dar.T1[,8:42]))
colnames(clad.darT1.mean) [1]<- 'Mean_count'

clad.darT2.mean<- as.data.frame(colMeans(clad.dar.T2[,8:42]))
colnames(clad.darT2.mean) [1]<- 'Mean_count'
  
#Cast these long format back into wide format
#T0
clad.darT0.meanwide<- dcast(clad.darT0.mean, Mean_count)

#T1

#T2
##

#Temporal BD between the different time periods
#T0 to T1
bdtemp.darT0T1<- decompose.D2(clad.dar.T0[,8:42], clad.dar.T1[,8:42], den.type=2) #but need to just have one row for each. 

#T1 to T2

#T0 to T2


#decompose.D2


#Looking at output
temp.D2.quant
temp.D2quant.mat1<- as.data.frame(temp.D2.quant$mat1)
temp.D2quant.mat2<- as.data.frame(temp.D2.quant$mat2)



##Knob 

##################################################################################################


##################################################################################################
####TIME PERIOD GAMMA DIVERSITY                                                                  #
##################################################################################################


##################################################################################################