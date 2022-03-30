## If installing impute or sva for the first time, uncomment the 3 lines below:

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("impute", version = "3.8")
#BiocManager::install("sva")
# BiocManager::install("genefilter")
# n

### Code settigngs:

my_col2<-c("blue",rgb(0,0,1,0.5),"white",rgb(1,0,0,0.5),"red")

# Remove peptides that are more abundant that 10% of the carrier channel from the single cell runs
remove_abundant_peptides<-T

# Use dart-enhanced spectra or spectra
dart_or_spectra<-"dart"


# Figure dimensions in inches
w.t<-5
h.t<-5

my_colors<-c("black","#048ABF","#FF5733")

# Filter to less than X missing data per column or row:
na.col<-0.99
na.row<-0.99

# Imputation
k.t<-3




####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Functions and libraries
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################




#BiocManager::install('flowCore')
# Check for presence of required packages for the import script
req.pkg<-c("plyr","limma","gtools",'Hmisc', "gridExtra", "stringr", "ggpubr", "Rcpp", "reshape2","psych",
           "ggplot2", "gridExtra", "gridGraphics", "grid", "ggridges",
           "matrixStats","patchwork","ggbeeswarm", "scales","stylo", "readr", "tidyverse","rlang","vctrs","tibble","dplyr","dbplyr")

#install.packages("ggbeeswarm")
#devtools::install_github("thomasp85/patchwork")


present.pkg<-installed.packages()[,1]

if(any(!req.pkg%in%present.pkg)){

  install.packages( req.pkg[ !req.pkg%in%present.pkg ] )

}




####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

### load libraries and functions
library(readr)
library(gtools)
library(reshape2)
library(scales)
library(ggplot2)
library(gridGraphics)
library(grid)
library(gridExtra)
library(stringr)
library(ggpubr)
library(ggridges)
library(matrixStats)
library(patchwork)
library(sva)
library(stylo)
library(tidyverse)
library(reticulate)
library(psych)
library(ComplexHeatmap)
library(circlize)
library(limma)
library(Hmisc)
#library("genefilter")
#library(flowCore)

library(ggbeeswarm)

# py_install("pandas")
# py_install('matplotlib')
# py_install('numba')
# py_install('scipy')
source_python('/Users/andrewleduc/Desktop/nPOP_Paper/code/func.py')



####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Column and row normalize

cr_norm<-function(dat){
  
  for(k in 1:ncol(dat)){
    
    dat[,k]<-dat[,k]/median(dat[,k], na.rm = T)
    
    
  }
  
  
  for(k in 1:nrow(dat)){
    
    dat[k,]<-dat[k,]/median(dat[k,], na.rm = T)
    
  }
  
  return(dat)
  
}



####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Column and row normalize on log scale

cr_norm_log<-function(dat){
  
  for(k in 1:ncol(dat)){
    
    dat[,k]<-dat[,k]-median(dat[,k], na.rm = T)
    
    
  }
  
  
  for(k in 1:nrow(dat)){
    
    dat[k,]<-dat[k,]-median(dat[k,], na.rm = T)
    
  }
  
  return(dat)
}


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# K - nearest neighbors imputation

hknn<-function(dat, k){
  
  # Create a copy of the data, NA values to be filled in later
  dat.imp<-dat
  
  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(dat)) )
  #dist.mat<- 1-as.matrix(cor((dat), use="pairwise.complete.obs"))

  #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))
  
  # Column names of the similarity matrix, same as data matrix
  cnames<-colnames(dist.mat)
  
  # For each column in the data... 
  for(X in cnames){
    
    # Find the distances of all other columns to that column 
    distances<-dist.mat[, X]
    
    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]
    
    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered<-dat[ , names(distances.ordered ) ]
    
    # Take the values in the column of interest
    vec<-dat[, X]
    
    # Which entries are missing and need to be imputed...
    na.index<-which( is.na(vec) )
    
    # For each of the missing entries (rows) in column X...
    for(i in na.index){
      
      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      
      #print(length(closest.columns))
      
      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){
        
        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns[1:k] ] )
        
      }
      
      
      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){
        
        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns ] )
        
      }
      
      
    }
    
    # Populate a the matrix with the new, imputed values
    dat.imp[,X]<-vec
    
  }
  
  return(dat.imp)
  
}



####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Calculate coefficient of variation for each and every protein with >= N peptides

prot.cv<-function(mat.t, meta, npeps){
  
  # Get rid of modified peptides
  modified.ind<-grep("(",row.names(mat.t), fixed=T)
  mat.t<-mat.t[-modified.ind,]
  
  # Row normalization
  for(k in 1:nrow(mat.t)){
    
    mat.t[k,]<-mat.t[k,]/mean(mat.t[k,], na.rm = T)
    #mat.t[k,]<-mat.t[k,] - mean(mat.t[k,], na.rm = T)
    
    
  }
  
  prots<-unique(meta$protein)
  
  cv.mat<-matrix(data=NA, nrow=length(prots), ncol=ncol(mat.t))
  
  for(i in 1:nrow(cv.mat)){
    
    peps<-unique( meta$sequence[meta$protein%in%prots[i]] )
    
    if(length(peps)>=npeps){
      
      values.t<-mat.t[rownames(mat.t)%in%peps, ]
      
      # add matrix that is count of peptides going into each cv value
      cvs<-t(sqrt( rowVars(t(values.t), na.rm=T) ) / rowMeans(t(values.t), na.rm=T) )
      
      cv.mat[i,]<-cvs
      
    }
    
    
  }
  
  rownames(cv.mat)<-prots
  colnames(cv.mat)<-colnames(mat.t)
  
  return(cv.mat)
  
}


prot.cvna<-function(mat.t, meta, npeps){
  
  # Get rid of modified peptides
  modified.ind<-grep("(",row.names(mat.t), fixed=T)
  mat.t<-mat.t[-modified.ind,]
  
  # Row normalization
  for(k in 1:nrow(mat.t)){
    
    mat.t[k,]<-mat.t[k,]/mean(mat.t[k,], na.rm = T)
    #mat.t[k,]<-mat.t[k,] - mean(mat.t[k,], na.rm = T)
    
    
  }
  
  prots<-unique(meta$protein)
  
  cv.mat<-matrix(data=NA, nrow=length(prots), ncol=ncol(mat.t))
  
  for(i in 1:nrow(cv.mat)){
    
    peps<-unique( meta$sequence[meta$protein%in%prots[i]] )
    
    if(length(peps)>=npeps){
      
      values.t<-t(mat.t[rownames(mat.t)%in%peps, ])
      
      # add matrix that is count of peptides going into each cv value
      cvs<-ncol(values.t) - na.count(values.t)
      
      cv.mat[i,]<-cvs
      
    }
    
    
  }
  
  rownames(cv.mat)<-prots
  colnames(cv.mat)<-colnames(mat.t)
  
  return(cv.mat)
  
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Calculate coefficient of variation for each and every protein with >= N peptides

prot.cv.log<-function(mat.t, meta, npeps){
  
  # Get rid of modified peptides
  modified.ind<-grep("(",row.names(mat.t), fixed=T)
  mat.t<-mat.t[-modified.ind,]
  
  # Row normalization
  for(k in 1:nrow(mat.t)){
    
    #mat.t[k,]<-mat.t[k,]/mean(mat.t[k,], na.rm = T)
    mat.t[k,]<-mat.t[k,] - mean(mat.t[k,], na.rm = T)
    
    
  }
  
  prots<-unique(ev.melt$protein)
  
  cv.mat<-matrix(data=NA, nrow=length(prots), ncol=ncol(mat.t))
  
  for(i in 1:nrow(cv.mat)){
    
    peps<-unique( ev.melt$sequence[ev.melt$protein%in%prots[i]] )
    
    if(length(peps)>=npeps){
      
      values.t<-mat.t[rownames(mat.t)%in%peps, ]
      
      # add matrix that is count of peptides going into each cv value
      cvs<-t(sqrt( rowVars(t(values.t), na.rm=T) ) / rowMeans(t(values.t), na.rm=T) )
      
      cv.mat[i,]<-cvs
      
    }
    
    
  }
  
  rownames(cv.mat)<-prots
  colnames(cv.mat)<-colnames(mat.t)
  
  return(cv.mat)
  
}



################################################################################################################
# remove.duplicates - remove duplicates across multiple columns
################################################################################################################

# Find unique entries by column1 and column2: 
# Ex remove.duplicates(TMT50,c("Sequence","Charge"))
remove.duplicates<-function(data,Cols){
  
  return(data[!duplicated(data[,Cols]),])
  
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


# Take mean or median, removing NA values
median.na<-function(y){ median(y, na.rm = T) }
mean.na<-function(x){ mean(x, na.rm=T) }


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


# Color palette for cells / control well

my_cell_colors<-c(rgb(0,0,0,0.9),rgb(0,0,1,0.7),rgb(1,0,0,0.7))


# Calculate False Discovery rate from PEP values
calc_fdr <- function(pep) {
  
  return( (cumsum(pep[order(pep)]) / seq(1, length(pep)))[order(order(pep))] )
  
}


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Count the number of NA values in each row of a matrix
na.count<-function(data){
  
  na.v<-c()
  for(i in 1:nrow(data)){
    
    na.v<-c(na.v, length(which(is.na(data[i, ]))) )
    
  }
  
  return(na.v)
  
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
### Matrix missing data filter function
# Filter remove COLUMNS then ROWS that have percent missing data
# (NA) greater than the specified thresholds
# Return the filtered matrix


filt.mat.cr<-function(mat, pct.r,pct.c){
  
  kc<-c()
  for(k in 1:ncol(mat)){
    
    pct.na<-length(which(is.na(mat[,k]))) / length(mat[,k])
    if(pct.na <= pct.c){ kc<-c(kc,k)}
    #print(pct.na)
    
    
  }
  
  mat<-mat[,kc]
  
  kr<-c()
  for(k in 1:nrow(mat)){
    
    pct.na<-length(which(is.na(mat[k,]))) / length(mat[k,])
    if(pct.na <= pct.r){ kr<-c(kr,k)}
    #print(pct.na)
    
    
  }
  
  mat<-mat[kr,]
  
  
  
  return(mat)
  
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
### Matrix missing data filter function
# Filter remove ROWS then COLUMNS that have percent missing data
# (NA) greater than the specified thresholds
# Return the filtered matrix


filt.mat.rc<-function(mat, pct.r,pct.c){
  
  kr<-c()
  for(k in 1:nrow(mat)){
    
    pct.na<-length(which(is.na(mat[k,]))) / length(mat[k,])
    if(pct.na <= pct.r){ kr<-c(kr,k)}
    #print(pct.na)
    
  }
  
  mat<-mat[kr,]
  
  kc<-c()
  for(k in 1:ncol(mat)){
    
    pct.na<-length(which(is.na(mat[,k]))) / length(mat[,k])
    if(pct.na <= pct.c){ kc<-c(kc,k)}
    #print(pct.na)
    
  }
  
  mat<-mat[,kc]
  

  
  
  
  return(mat)
  
}



####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Function to perform collapse, by median

collapse_to_protein<-function(ev.matrix.t, ev.melt, row_norm){
  
  ev.matrix.t.prot<-ev.matrix.t[0,]
  
  if(row_norm){
    
    # Row normalization
    for(k in 1:nrow(ev.matrix.t)){
      
      #ev.matrix.t[k,]<-ev.matrix.t[k,]/median(ev.matrix.t[k,], na.rm = T)
      ev.matrix.t[k,]<-ev.matrix.t[k,] - mean(ev.matrix.t[k,], na.rm = T)
      
    }
  }
  
  prot.v<-c()
  for(X in unique(as.character(ev.melt$protein))){
    
    rows.t<- which(rownames(ev.matrix.t) %in% ev.melt$sequence[ev.melt$protein%in%X])
    
    if(length(rows.t)>0){
      
      #print(length(rows.t))
      mat.t<-matrix(nrow = 1, ncol = ncol(ev.matrix.t))
      if(length(rows.t)>1){   mat.t<-apply(ev.matrix.t[rows.t,],2,median.na) }
      if(length(rows.t)==1){   mat.t<-ev.matrix.t[rows.t,] }
      
      ev.matrix.t.prot<-rbind(ev.matrix.t.prot, mat.t )
      
      prot.v<-c(prot.v, X)
      
    }
    
  }
  
  row.names(ev.matrix.t.prot)<-prot.v
  
  print( dim(ev.matrix.t.prot) )
  
  return( ev.matrix.t.prot )
  
}



##############################################

cv<-function(x){
  
  sd(x, na.rm=T) / mean(x, na.rm=T)
  
}

cvna<-function(x){
  
  sum(!is.na(x))
  
}




####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# Calculate peptide level FDR

pep_fdr_greaterthan<-function(dat, thres){

rmi<-c()

for(X in unique(dat$modseq)){
  
  dp<-dat$dart_PEP[dat$modseq==X]
  
  rawt<-dat$Raw.file[dat$modseq==X]
  
  dpf<-calc_fdr(dp)
  
  rk<-rawt[which(dpf>thres)]
  
  indt<-which(dat$Raw.file%in%rk & dat$modseq==X)
  
  rmi<-c(rmi, indt)
  
}

return(rmi)

}


# Calculate peptide level FDR

pep_fdr_greaterthan_PEP<-function(dat, thres){
  
  rmi<-c()
  
  for(X in unique(dat$modseq)){
    
    dp<-dat$PEP[dat$modseq==X]
    
    rawt<-dat$Raw.file[dat$modseq==X]
    
    dpf<-calc_fdr(dp)
    
    rk<-rawt[which(dpf>thres)]
    
    indt<-which(dat$Raw.file%in%rk & dat$modseq==X)
    
    rmi<-c(rmi, indt)
    
  }
  
  return(rmi)
  
}

Empirical_FDR_dPEP <- function(ev){
  exps <- levels(factor(ev$Raw.file))
  ev$Score_naOmit <- ifelse(is.na(ev$dart_PEP),100,ev$dart_PEP)
  for(i in 1:length(exps)){
    ev_single <- ev %>% dplyr::filter(Raw.file == levels(factor(ev$Raw.file))[i])
    ev_rank <- ev_single %>% dplyr::arrange(Score_naOmit)
    ev_rank$Ind <- seq(1:nrow(ev_rank))
    ev_rank$revCount <- ifelse(ev_rank$Reverse=="+",1,0)
    ev_rank <- ev_rank %>% dplyr::mutate(RevHits = cumsum(revCount),qVal = RevHits/Ind)
    if(i== 1){
      ev_rank_out <- ev_rank
    }
    if(i > 1){
      ev_rank_out <- rbind(ev_rank_out,ev_rank)
    }
  }
  return(ev_rank_out)
}

########################################################################################################################
########################################################################################################################

# Functions for Distance Calculations that give number of pairwise observations

Distance_euclid <- function(prepp, obsThreshold){
  
  # count the number of shared observations between protein (i.e how many single cells have both proteins quantified)
  obser <- pairwiseCount(t(prepp))
  # obtain the upper triangle of shared observation matrix 
  obser[lower.tri(obser, diag = T)] <- 188 
  obser.m <- reshape2::melt(obser)
  obser.m$Var1Var2 <- paste(obser.m$Var1, " ", obser.m$Var2)
  
  # obtain the correlation matrix 
  #ppCor <- cor(t(prepp), method = "spearman", use = "pairwise.complete.obs")
  ppCor <- stats::dist(as.matrix(prepp),method = "euclidean")
  # obtain the upper trianlge of the correlation matrix
  ppCor1 <- (as.matrix(ppCor))
  #ppCor1 <- cor(t(prepp), method = "spearman", use = "pairwise.complete.obs")
  ppCor1[lower.tri(ppCor1, diag = T)] <- 188
  ppCor1.m <- reshape2::melt(ppCor1)
  ppCor1.m$Var1Var2 <- paste(ppCor1.m$Var1, " ",  ppCor1.m$Var2)
  
  # merge the two matrices together 
  ppCor1.m$obser <- obser.m$value[match(ppCor1.m$Var1Var2, obser.m$Var1Var2)]
  # filter for the number of minimum shared observations, and for the upper triangle 
  ppCor1.m.fil <- ppCor1.m %>% dplyr::filter(obser > obsThreshold & value != 188)
  
  #ggplot(data = ppCor1.m.fil, aes(x = value)) + geom_histogram(color = "black", fill = "grey") +
  #  geom_vline(xintercept = 0) + 
  #  theme_classic() + 
  #  labs(title = "Correlations with at least x observations", 
  #       subtitle = "filtered for upper triangle")
  
  return(ppCor1.m.fil)
}
Distance_cor <- function(prepp, obsThreshold){
  
  # count the number of shared observations between protein (i.e how many single cells have both proteins quantified)
  obser <- pairwiseCount(t(prepp))
  # obtain the upper triangle of shared observation matrix 
  obser[lower.tri(obser, diag = T)] <- 188 
  obser.m <- reshape2::melt(obser)
  obser.m$Var1Var2 <- paste(obser.m$Var1, " ", obser.m$Var2)
  
  # obtain the correlation matrix 
  ppCor <- Hmisc::rcorr(t(prepp), type = "spearman")$r
  # obtain the upper trianlge of the correlation matrix
  ppCor1 <- ppCor
  ppCor1[lower.tri(ppCor1, diag = T)] <- 188
  ppCor1.m <- reshape2::melt(ppCor1)
  ppCor1.m$Var1Var2 <- paste(ppCor1.m$Var1, " ",  ppCor1.m$Var2)
  
  # merge the two matrices together 
  ppCor1.m$obser <- obser.m$value[match(ppCor1.m$Var1Var2, obser.m$Var1Var2)]
  # filter for the number of minimum shared observations, and for the upper triangle 
  ppCor1.m.fil <- ppCor1.m %>% dplyr::filter(obser > obsThreshold & value != 188)
  
  ggplot(data = ppCor1.m.fil, aes(x = value)) + geom_histogram(color = "black", fill = "grey") +
    geom_vline(xintercept = 0) + 
    theme_classic() + 
    labs(title = "Correlations with at least x observations", 
         subtitle = "filtered for upper triangle")
  
  return(ppCor1.m.fil)
}
Distance_cosine <- function(prepp, obsThreshold){
  
  # count the number of shared observations between protein (i.e how many single cells have both proteins quantified)
  obser <- pairwiseCount(t(prepp))
  # obtain the upper triangle of shared observation matrix 
  obser[lower.tri(obser, diag = T)] <- 188 
  obser.m <- reshape2::melt(obser)
  obser.m$Var1Var2 <- paste(obser.m$Var1, " ", obser.m$Var2)
  
  # obtain the cosine matrix 
  mmat_p <- !is.na(prepp)
  ppCor <- sparse_cosine_similarity(prepp,mmat_p)
  rownames(ppCor) <- rownames(prepp)
  colnames(ppCor) <- rownames(prepp)
  # obtain the upper trianlge of the correlation matrix
  ppCor1 <- ppCor
  ppCor1[lower.tri(ppCor1, diag = T)] <- 188
  ppCor1.m <- reshape2::melt(ppCor1)
  ppCor1.m$Var1Var2 <- paste(ppCor1.m$Var1, " ",  ppCor1.m$Var2)
  
  # merge the two matrices together 
  ppCor1.m$obser <- obser.m$value[match(ppCor1.m$Var1Var2, obser.m$Var1Var2)]
  # filter for the number of minimum shared observations, and for the upper triangle 
  ppCor1.m.fil <- ppCor1.m %>% dplyr::filter(obser > obsThreshold & value != 188)
  
  ggplot(data = ppCor1.m.fil, aes(x = value)) + geom_histogram(color = "black", fill = "grey") +
    geom_vline(xintercept = 0) + 
    theme_classic() + 
    labs(title = "Correlations with at least x observations", 
         subtitle = "filtered for upper triangle")
  
  return(ppCor1.m.fil)
}
Distance_cosine_weighted <- function(prepp,weights, obsThreshold){
  
  # count the number of shared observations between protein (i.e how many single cells have both proteins quantified)
  obser <- pairwiseCount(t(prepp))
  # obtain the upper triangle of shared observation matrix 
  obser[lower.tri(obser, diag = T)] <- 188 
  obser.m <- reshape2::melt(obser)
  obser.m$Var1Var2 <- paste(obser.m$Var1, " ", obser.m$Var2)
  
  # obtain the cosine matrix 
  mmat_p <- !is.na(prepp)
  weights <- as.matrix(weights)
  ppCor <- sparse_cosine_similarity_weight(prepp,mmat_p,weights)
  rownames(ppCor) <- rownames(prepp)
  colnames(ppCor) <- rownames(prepp)
  # obtain the upper trianlge of the correlation matrix
  ppCor1 <- ppCor
  ppCor1[lower.tri(ppCor1, diag = T)] <- 188
  ppCor1.m <- reshape2::melt(ppCor1)
  ppCor1.m$Var1Var2 <- paste(ppCor1.m$Var1, " ",  ppCor1.m$Var2)
  
  # merge the two matrices together 
  ppCor1.m$obser <- obser.m$value[match(ppCor1.m$Var1Var2, obser.m$Var1Var2)]
  # filter for the number of minimum shared observations, and for the upper triangle 
  ppCor1.m.fil <- ppCor1.m %>% dplyr::filter(obser > obsThreshold & value != 188)
  
  ggplot(data = ppCor1.m.fil, aes(x = value)) + geom_histogram(color = "black", fill = "grey") +
    geom_vline(xintercept = 0) + 
    theme_classic() + 
    labs(title = "Correlations with at least x observations", 
         subtitle = "filtered for upper triangle")
  
  return(ppCor1.m.fil)
}
Distance_dot <- function(prepp, obsThreshold){
  
  # count the number of shared observations between protein (i.e how many single cells have both proteins quantified)
  obser <- pairwiseCount(t(prepp))
  # obtain the upper triangle of shared observation matrix 
  obser[lower.tri(obser, diag = T)] <- 188 
  obser.m <- reshape2::melt(obser)
  obser.m$Var1Var2 <- paste(obser.m$Var1, " ", obser.m$Var2)
  
  # obtain the cosine matrix 
  mmat_p <- !is.na(prepp)
  ppCor <- sparse_dot(prepp,mmat_p)
  rownames(ppCor) <- rownames(prepp)
  colnames(ppCor) <- rownames(prepp)
  # obtain the upper trianlge of the correlation matrix
  ppCor1 <- ppCor
  ppCor1[lower.tri(ppCor1, diag = T)] <- 188
  ppCor1.m <- reshape2::melt(ppCor1)
  ppCor1.m$Var1Var2 <- paste(ppCor1.m$Var1, " ",  ppCor1.m$Var2)
  
  # merge the two matrices together 
  ppCor1.m$obser <- obser.m$value[match(ppCor1.m$Var1Var2, obser.m$Var1Var2)]
  # filter for the number of minimum shared observations, and for the upper triangle 
  ppCor1.m.fil <- ppCor1.m %>% dplyr::filter(obser > obsThreshold & value != 188)
  
  ggplot(data = ppCor1.m.fil, aes(x = value)) + geom_histogram(color = "black", fill = "grey") +
    geom_vline(xintercept = 0) + 
    theme_classic() + 
    labs(title = "Correlations with at least x observations", 
         subtitle = "filtered for upper triangle")
  
  return(ppCor1.m.fil)
}
Distance_euclid_py <- function(prepp, obsThreshold){
  
  # count the number of shared observations between protein (i.e how many single cells have both proteins quantified)
  obser <- pairwiseCount(t(prepp))
  # obtain the upper triangle of shared observation matrix 
  obser[lower.tri(obser, diag = T)] <- 188 
  obser.m <- reshape2::melt(obser)
  obser.m$Var1Var2 <- paste(obser.m$Var1, " ", obser.m$Var2)
  
  # obtain the cosine matrix 
  mmat_p <- !is.na(prepp)
  ppCor <- sparse_euclid(prepp,mmat_p)
  rownames(ppCor) <- rownames(prepp)
  colnames(ppCor) <- rownames(prepp)
  # obtain the upper trianlge of the correlation matrix
  ppCor1 <- ppCor
  ppCor1[lower.tri(ppCor1, diag = T)] <- 188
  ppCor1.m <- reshape2::melt(ppCor1)
  ppCor1.m$Var1Var2 <- paste(ppCor1.m$Var1, " ",  ppCor1.m$Var2)
  
  # merge the two matrices together 
  ppCor1.m$obser <- obser.m$value[match(ppCor1.m$Var1Var2, obser.m$Var1Var2)]
  # filter for the number of minimum shared observations, and for the upper triangle 
  ppCor1.m.fil <- ppCor1.m %>% dplyr::filter(obser > obsThreshold & value != 188)
  
  ggplot(data = ppCor1.m.fil, aes(x = value)) + geom_histogram(color = "black", fill = "grey") +
    geom_vline(xintercept = 0) + 
    theme_classic() + 
    labs(title = "Correlations with at least x observations", 
         subtitle = "filtered for upper triangle")
  
  return(ppCor1.m.fil)
}
Distance_propr_py <- function(prepp, obsThreshold){
  
  # count the number of shared observations between protein (i.e how many single cells have both proteins quantified)
  obser <- pairwiseCount(t(prepp))
  # obtain the upper triangle of shared observation matrix 
  obser[lower.tri(obser, diag = T)] <- 188 
  obser.m <- reshape2::melt(obser)
  obser.m$Var1Var2 <- paste(obser.m$Var1, " ", obser.m$Var2)
  
  # obtain the cosine matrix 
  mmat_p <- !is.na(prepp)
  ppCor <- sparse_propr(prepp,mmat_p)
  rownames(ppCor) <- rownames(prepp)
  colnames(ppCor) <- rownames(prepp)
  # obtain the upper trianlge of the correlation matrix
  ppCor1 <- ppCor
  ppCor1[lower.tri(ppCor1, diag = T)] <- 188
  ppCor1.m <- reshape2::melt(ppCor1)
  ppCor1.m$Var1Var2 <- paste(ppCor1.m$Var1, " ",  ppCor1.m$Var2)
  
  # merge the two matrices together 
  ppCor1.m$obser <- obser.m$value[match(ppCor1.m$Var1Var2, obser.m$Var1Var2)]
  # filter for the number of minimum shared observations, and for the upper triangle 
  ppCor1.m.fil <- ppCor1.m %>% dplyr::filter(obser > obsThreshold & value != 188)
  
  ggplot(data = ppCor1.m.fil, aes(x = value)) + geom_histogram(color = "black", fill = "grey") +
    geom_vline(xintercept = 0) + 
    theme_classic() + 
    labs(title = "Correlations with at least x observations", 
         subtitle = "filtered for upper triangle")
  
  return(ppCor1.m.fil)
}

Distance_cor_sig <- function(prepp, obsThreshold){
  
  # count the number of shared observations between protein (i.e how many single cells have both proteins quantified)
  obser <- pairwiseCount(t(prepp))
  # obtain the upper triangle of shared observation matrix 
  obser[lower.tri(obser, diag = T)] <- 188 
  obser.m <- reshape2::melt(obser)
  obser.m$Var1Var2 <- paste(obser.m$Var1, " ", obser.m$Var2)
  
  # obtain the correlation matrix 
  ppCor <- cor(t(prepp), method = "spearman", use = "pairwise.complete.obs")
  # obtain the upper trianlge of the correlation matrix
  ppCor1 <- ppCor
  ppCor1[lower.tri(ppCor1, diag = T)] <- 188
  ppCor1.m <- reshape2::melt(ppCor1)
  ppCor1.m$Var1Var2 <- paste(ppCor1.m$Var1, " ",  ppCor1.m$Var2)
  
  # merge the two matrices together 
  ppCor1.m$obser <- obser.m$value[match(ppCor1.m$Var1Var2, obser.m$Var1Var2)]
  # filter for the number of minimum shared observations, and for the upper triangle 
  ppCor1.m.fil <- ppCor1.m %>% dplyr::filter(obser > obsThreshold & value != 188)
  
  ggplot(data = ppCor1.m.fil, aes(x = value)) + geom_histogram(color = "black", fill = "grey") +
    geom_vline(xintercept = 0) + 
    theme_classic() + 
    labs(title = "Correlations with at least x observations", 
         subtitle = "filtered for upper triangle")
  
  return(ppCor1.m.fil)
}

########################################################################################################################
########################################################################################################################

#GSEA on different populations

Population_GSEA <- function(data,GO_db,sn){

  if (sn == 2){
    pop1 <- unique(data$Condition)[1]
    pop2 <- unique(data$Condition)[2]
    pops <- c(pop1,pop2)
    lens <- c(length(unique(data$SampID[data$Condition == pop1])),length(unique(data$SampID[data$Condition == pop2])))
    GSEA_output <- data.frame(GO_term = character(), pVal = numeric(),numberOfMatches= numeric(),fractionOfDB_Observed = numeric(),stringsAsFactors = FALSE,
                              Cond1med_int = numeric(),Cond2med_int = numeric())
    
  }
  if (sn==3){
    pop1 <- unique(data$Condition)[1]
    pop2 <- unique(data$Condition)[2]
    pop3 <- unique(data$Condition)[3]
    pops <- c(pop1,pop2,pop3)
    lens <- c(c(length(unique(data$SampID[data$Condition == pop1])),length(unique(data$SampID[data$Condition == pop2])),length(unique(data$SampID[data$Condition == pop3]))))
    GSEA_output <- data.frame(GO_term = character(), pVal = numeric(),numberOfMatches= numeric(),fractionOfDB_Observed = numeric(),stringsAsFactors = FALSE,
                              Cond1med_int = numeric(),Cond2med_int = numeric(),Cond3med_int = numeric())
    
  }
  unique_GO_terms <- unique(GO_db$GO_term_name)

 
  for (i in 1:length(unique_GO_terms)){
    
    GO_term <- unique_GO_terms[i]
    GO_db_lim <- GO_db %>% dplyr::filter(GO_term_name == GO_term)
    data_matches <- data %>% dplyr::filter(Protein %in% GO_db_lim$Uniprot)
    matches_number <- length(unique(data_matches$Protein))
    total_Proteins_in_GOterm <- length(unique(GO_db_lim$Uniprot))
    data_matches <- data_matches %>% filter(is.na(Intensity)==F)
    fractionObserved <- matches_number/total_Proteins_in_GOterm
    dataProt_matches_medInt <- data_matches %>% dplyr::group_by(Condition)  %>% dplyr::summarize(medianInt = median(Intensity,na.rm = T))
    dataProt_matches_medInt$GO_term <- GO_term
    
    
    #make sure majority of samples have observations present
    check <- T
    for (j in 1:length(pops)){
      check <- length(unique(data_matches$SampID[data_matches$Condition == pops[j]])) > (lens[j]*.7)
      if(check == F){
        break
      }
    }
    if(matches_number > 1 && length(unique(data_matches$Condition)) > 1 && check == T ){
      
      AOV_out <- aov(Intensity ~ Condition, data = data_matches)
      AOV_out <- summary(AOV_out)
      AOV_out <- data.frame(AOV_out[[1]])
      pVal <-AOV_out[1,5]
    }
    else{
      pVal = NA
    }
    GSEA_output[i,1] <- GO_term
    GSEA_output[i,2] <- pVal 
    GSEA_output[i,3] <- matches_number
    GSEA_output[i,4] <- fractionObserved
    GSEA_output[i,5]<- dataProt_matches_medInt$medianInt[1]
    GSEA_output[i,6]<- dataProt_matches_medInt$medianInt[2]
    if(sn ==3){
      GSEA_output[i,7]<- dataProt_matches_medInt$medianInt[3]
    }
    
  }
  return(GSEA_output)
}

Population_GSEA_norm <- function(data,GO_db,sn){
  
  unique_GO_terms <- unique(GO_db$GO_term_name)
  
  GSEA_output <- data.frame(GO_term = character(), pVal = numeric(),numberOfMatches= numeric(),fractionOfDB_Observed = numeric(),stringsAsFactors = FALSE,
                            Cond1med_int = numeric(),Cond2med_int = numeric(),Cond3med_int = numeric())
  
  for (i in 1:length(unique_GO_terms)){
    
    GO_term <- unique_GO_terms[i]
    GO_db_lim <- GO_db %>% dplyr::filter(GO_term_name == GO_term)
    data_matches <- data %>% dplyr::filter(Protein %in% GO_db_lim$Uniprot)
    matches_number <- length(unique(data_matches$Protein))
    total_Proteins_in_GOterm <- length(unique(GO_db_lim$Uniprot))
    data_matches <- data_matches %>% filter(is.na(Intensity)==F)
    fractionObserved <- matches_number/total_Proteins_in_GOterm
    dataProt_matches_medInt <- data_matches %>% dplyr::group_by(Condition)  %>% dplyr::summarize(medianInt = median(Intensity,na.rm = T))
    dataProt_matches_medInt$GO_term <- GO_term
    
    # Kruskall Wallis
    #make sure majority of samples have observations present
    
    if(matches_number > 2 && length(unique(data_matches$Condition)) > 2 ){
      
      AOV_out <- aov(Intensity ~ Condition, data = data_matches)
      AOV_out <- summary(AOV_out)
      AOV_out <- data.frame(AOV_out[[1]])
      pVal <-AOV_out[1,5]
      GSEA_output[i,1] <- GO_term
      GSEA_output[i,2] <- pVal 
      GSEA_output[i,3] <- matches_number
      GSEA_output[i,4] <- fractionObserved
      GSEA_output[i,5]<- dataProt_matches_medInt$medianInt[dataProt_matches_medInt$Condition == 'G1']
      GSEA_output[i,6]<- dataProt_matches_medInt$medianInt[dataProt_matches_medInt$Condition == 'S']
      GSEA_output[i,7]<- dataProt_matches_medInt$medianInt[dataProt_matches_medInt$Condition == 'G2']
    }
    else{
      pVal = NA
    }
    GSEA_output[i,1] <- NA
    GSEA_output[i,2] <- NA 
    GSEA_output[i,3] <- NA
    GSEA_output[i,4] <- NA
    GSEA_output[i,5]<- NA
    GSEA_output[i,6]<- NA
    GSEA_output[i,7]<- NA
    
  }
  return(GSEA_output)
}



#Messy other functions (Random)

getDat <- function(evidence){
  evidence$seqCharge <- paste0(evidence$Modified.sequence,evidence$Charge)
  evidence <- evidence %>% dplyr::filter(PEP < .02, Potential.contaminant != '+',Reverse != '+',PIF > .8)
  evidence <- evidence %>% dplyr::select('Raw.file','seqCharge','Intensity','Retention.time','Leading.razor.protein',
                                         'Missed.cleavages',contains('Reporter.intensity.corrected'))
  
  return(evidence)
}
getRI <- function(evidence){
  b <- NULL
  for (i in 1:ncol(evidence)){
    if (median(evidence[,i]==0)){
      a <- colnames(evidence)[i]
      b <- c(b,a)
    }
  }
  evidence <- evidence %>% dplyr::select(-all_of(b))
  return(evidence)
}
row_norm <- function(data,start,norm){
  a <- colnames(data)[1]
  data_cut <- as.matrix(data[,start:ncol(data)])
  
  if(norm == 0 | missing(norm)){
    for (i in 1:nrow(data_cut)){
      data_cut[i,] <- data_cut[i,]/mean(data_cut[i,],na.rm = TRUE)
    }
  } else {
    for (i in 1:nrow(data_cut)){
      data_cut[i,] <- data_cut[i,]-mean(data_cut[i,],na.rm = TRUE)
    }
  }
  data_cut <- as.data.frame(data_cut)
  data<- cbind(data[,1:(start-1)],data_cut)
  if (start == 2){
    colnames(data)[1] <- a
  }
  return(data)
}
col_norm <- function(data,start,norm){
  if(missing(norm) | norm == 0){
    for (i in start:ncol(data)){
      data[,i] <- data[,i]/median(data[,i],na.rm = TRUE)
    }
  } else {
    for (i in start:ncol(data)){
      data[,i] <- data[,i]-median(data[,i],na.rm = TRUE)
    }
  }
  return(data)
}
