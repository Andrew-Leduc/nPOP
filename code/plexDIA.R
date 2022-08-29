source('/Users/andrewleduc/Desktop/GitHub/plexDIA/code/plexDIA_Functions.R')
source('/Users/andrewleduc/Desktop/GitHub/plexDIA/code/functions_parameters.R')


#Failed injections not to include when doing DIA-NN max LFQ
runs <-c('wAL_plexMel65','wAL_plexMel66','wAL_plexMel02','wAL_plexMel04','wAL_plexMel06','wAL_plexMel08','wAL_plexMel10','wAL_plexMel2',
         'wAL_plexMel14','wAL_plexMel16',
         'wAL_plexMel18','wAL_plexMel20','wAL_plexMel22',
         'wAL_plexMel24','wAL_plexMel26','wAL_plexMel28','wAL_plexMel30','wAL_plexMel32')


# load meta data... just Runs X labels
meta<-read.csv('../dat/raw_data/plexDIA/annotation_plexDIA.csv')
meta <- meta %>% filter(!Set %in% runs)
meta <- melt(meta,id.vars = 'Set')
meta$variable <- str_sub(meta$variable,2,2)
meta$Set2 <- meta$Set
meta$Set <- paste0(meta$Set,'.raw.',meta$variable)
meta$variable <- NULL

#Read main search
evnew <- read.delim('../dat/raw_data/plexDIA/report_plexDIA_mel_nPOP.tsv')
evnew <- evnew %>% filter(!Run %in% runs)
#write.table(evnew, file="../dat/raw_data/mymatrix.tsv", row.names=FALSE, col.names=TRUE,sep = '\t',quote=FALSE)

meta_set <- meta %>% filter(Set %in% colnames(Ms1_extract2))


evnew <- evnew %>% filter(Lib.PG.Q.Value < .01)

# Optional... have to test
evnew <- evnew %>% filter(Channel.Q.Value < .01)

# I used the MS1 Extracted for the quant
Ms1_extract<-read.delim('../dat/raw_data/plexDIA/report.pr_matrix_channels_ms1_extracted.tsv')

#For linking Main search that allows you to filter to the MS1 extracted, then filter MS1 extracted
evnew$seqrun <- paste0(evnew$Stripped.Sequence,evnew$Precursor.Charge,evnew$Run,'.raw.',substr(evnew$Precursor.Id[1:nrow(evnew)], 10, 10))
Ms1_extract$seqcharge <- paste0(Ms1_extract$Stripped.Sequence,Ms1_extract$Precursor.Charge)
Ms1_extract <- Ms1_extract[,!grepl('Quality',colnames(Ms1_extract))]
Ms1_extract <- as.data.frame(Ms1_extract) %>% dplyr::select(Protein.Group,seqcharge, contains('D..AL.AL.wAL'))
Ms1_extract <- melt(Ms1_extract, id = c('Protein.Group','seqcharge'))
Ms1_extract$variable <- str_sub(Ms1_extract$variable,10,30)
Ms1_extract$seqrun <- paste0(Ms1_extract$seqcharge,Ms1_extract$variable)
Ms1_extract <- Ms1_extract %>% filter(seqrun %in% evnew$seqrun)
Ms1_extract$seqrun <- NULL

#normalize for computing CVs
Ms1_extract <- dcast(Ms1_extract,Protein.Group+seqcharge~variable,value.var = 'value')
p2p_map <-Ms1_extract %>% dplyr::select(Protein.Group,seqcharge)
rownames(Ms1_extract)<- Ms1_extract$seqcharge

prot_intense <- melt(Ms1_extract)
Ms1_extract$Protein.Group <- NULL
Ms1_extract$seqcharge <- NULL
Ms1_extract <- as.matrix(Ms1_extract)
Ms1_extract <- Ms1_extract[,colnames(Ms1_extract) %in% meta$Set]

#save the un-normalized matrix
Ms1_extract_raw.m <- Ms1_extract
Ms1_extract[Ms1_extract==0] <- NA

#save the un-normalized matrix with 0s NAed
Ms1_extract_raw <- Ms1_extract
for(i in 1:ncol(Ms1_extract)){
  Ms1_extract[,i] <- Ms1_extract[,i]/median(Ms1_extract[,i],na.rm = T)
}
for(i in 1:nrow(Ms1_extract)){
  Ms1_extract[i,] <- Ms1_extract[i,]/mean(Ms1_extract[i,],na.rm = T)
}


CV_mat <- as.data.frame(Ms1_extract)
CV_mat$prot <- p2p_map$Protein.Group
CV_mat <- melt(CV_mat,id.vars = 'prot')


### this is just for plotting CVs vs protein intensity
abs_mat <- as.data.frame(Ms1_extract_raw)
abs_mat$prot <- p2p_map$Protein.Group
abs_mat <- melt(abs_mat,id.vars = 'prot')

abs_mat <- abs_mat %>%
  dplyr::group_by(prot, variable) %>%
  dplyr::summarise(cvq = mean(value,na.rm=T))
###


CV_mat <- CV_mat %>%
  dplyr::group_by(prot, variable) %>%
  dplyr::summarise(cvq = cv(value),ct = sum(is.na(value)==F))

CV_mat2 <- CV_mat
CV_mat2$intense <- abs_mat$cvq
CV_mat_plot <- CV_mat2 %>% filter(ct > 2)


# Plotting intensitys vs CVs
plot(log10(CV_mat_plot$intense),(CV_mat_plot$cvq))
ggplot(CV_mat_plot,aes(x=log10(intense),y=cvq)) + 
  geom_point(alpha = .05,color='red')+ylim(c(0,2))+
  theme_bw() +ylab('Protein CVs')+xlab('log10(Protein intensity)')+
  font("xylab", size=30) +
  font("xy.text", size=20)
  


#CV_mat <- CV_mat %>% filter(is.na(cvq)==F)

CV_mat <- CV_mat %>%
  dplyr::group_by( variable) %>%
  dplyr::summarise(cvq = median(cvq,na.rm = T))


CV_mat <- CV_mat %>% left_join(meta, by = c('variable'='Set'))
CV_mat$SumIntense <- colSums(Ms1_extract_raw, na.rm = T)



CV_plot_plex <- ggplot(data=CV_mat, aes(x=cvq,fill=value)) + geom_density( alpha=0.5,adjust=1.5) + theme_pubr() +
  scale_fill_manual(values=my_col3[c(2,1)]) + 
  xlab("CV") + ylab("Fraction of cells") + rremove("y.ticks") + rremove("y.text") +
  font("xylab", size=30) +
  font("x.text", size=20) +
  coord_cartesian(xlim=c(.1,.8))+
  annotate("text", x=0.2, y= 14, label=paste0("104 cells"), size=10, color=my_col3[c(2)])+
  annotate("text", x=0.64, y= 12, label=paste0("6"," Ctr -"), size=10, color=my_col3[c(1)])+
  annotate("text", x=0.63, y= 14, label=paste0('15'," cells"), size=10, color=my_col3[c(2)])+
  annotate("text", x=0.2, y= 12, label=paste0("0"," Ctr -"), size=10, color=my_col3[c(1)])+
  #annotate("text", x=0.25, y= 3, label="Macrophage-like", size=6) +
  rremove("legend") + geom_vline(xintercept=0.4, lty=2, size=2, color="gray50") + theme(plot.margin = margin(1, 1, 0, 1, "cm"))
CV_plot_plex
#ggsave('/Users/andrewleduc/Desktop/nPOP_Paper/PDFs/CV_plot_plex.pdf',plot = CV_plot_plex,height = unit(6, "in"),width = unit(7, "in"),dpi = 1000)



### CV filter
CV_mat<- CV_mat %>% filter(cvq < .4)


Ms1_extract <- Ms1_extract[,colnames(Ms1_extract) %in% CV_mat$variable]

Ms1_extract.pep_wrt <- Ms1_extract[,colnames(Ms1_extract3_no_imp)]
write.csv(Ms1_extract.pep_wrt, file="../dat/Proccessed_data/plexDIA/peptidesXsc_plexDIA.csv")


# Plotting of negative signal intensity...
pep_and_prot <- as.data.frame(colSums(is.na(Ms1_extract)==F))
pep_and_prot$type <- 'Peptides'
colnames(pep_and_prot)[1]<- 'value'
pep_and_prot <- pep_and_prot %>% filter(value > 250)

frac1 <- (hist(log10(Ms1_extract_raw.m_neg$value),plot=F)$counts) / (sum(hist(log10(Ms1_extract_raw.m_neg$value),plot=F)$counts))
frac2 <- (hist(log10(Ms1_extract_raw.m_pos$value),plot=F)$counts) / (sum(hist(log10(Ms1_extract_raw.m_pos$value),plot=F)$counts))

pep_and_prot1 <- as.data.frame(frac1)
pep_and_prot1$type = ' Negative ctrl'
pep_and_prot1$scaler <- (1:nrow(pep_and_prot1))/2
pep_and_prot1$frac1[6] <- 1
pep_and_prot1$frac1[2:13]= 0
pep_and_prot2<- as.data.frame(frac2)
pep_and_prot2$type <- 'Melanoma'
pep_and_prot2$scaler <- (1:nrow(pep_and_prot2))/2
colnames(pep_and_prot2)[1] <- 'frac1'
pep_and_prot2$frac1[6] <- .148
pep_and_prot3 <- rbind(pep_and_prot1,pep_and_prot2)
pep_and_prot3$frac1 <-round(pep_and_prot3$frac1,2)
neg_sig <- ggplot(pep_and_prot3,aes(y = frac1,x=scaler,fill=type )) + geom_bar(stat="identity",position ="dodge",alpha = 1)+ theme_classic()+
  ylab('Fraction of values') + xlab('log10(Intensity)')+scale_fill_manual(values = c("#000000","#F6BE00"))+
  theme(axis.title.x = element_text(size = 30),axis.title.y = element_text(size = 30),axis.text.x = element_text(size = 25,color = 'black'),
        axis.text.y = element_text(size = 25,color = 'black'),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=0), #change legend title font size
        legend.text = element_text(size=25),
        legend.position=c(.70, .75)) + xlim(c(2.5,7))

#ggsave('/Users/andrewleduc/Desktop/nPOP_Paper/PDFs/neg_sig.pdf',plot = neg_sig,height = unit(6, "in"),width = unit(7, "in"),dpi = 1000)



### Back to regular data proccessing, use raw matrix with only good CV cells for doing DIA-NNmaxLFQ
Ms1_extract <- Ms1_extract[,rownames(pep_and_prot)]
Ms1_extract <- Ms1_extract_raw[rownames(Ms1_extract),colnames(Ms1_extract)]
Ms1_extract <- as.data.frame(Ms1_extract)
Ms1_extract$seqcharge <- rownames(Ms1_extract)
Ms1_extract$seqcharge <- rownames(Ms1_extract)
Ms1_extract <- Ms1_extract %>% left_join(p2p_map, by = c('seqcharge'))

Ms1_extract2 <- reshape2::melt(Ms1_extract,id = c('seqcharge','Protein.Group'))
colnames(Ms1_extract2)[3] <- 'File.Name'
Ms1_extract2 <- diann_maxlfq(Ms1_extract2, id.header = "seqcharge",group.header = "Protein.Group",quantity.header = "value")
Ms1_extract2 <- as.matrix(Ms1_extract2)
Ms1_extract2_noImp <- Ms1_extract2


# Just for ploting total peptides and proteins
hist(colSums(is.na(Ms1_extract2)==F),20)
pep_and_prot2 <- as.data.frame(colSums(is.na(Ms1_extract2)==F))
pep_and_prot2$type <- 'Proteins'
colnames(pep_and_prot2)[1]<- 'value'
pep_and_prot3 <- rbind(pep_and_prot,pep_and_prot2)
median(pep_and_prot2$value)
pep_prot_plex <- ggplot(pep_and_prot3,aes(x = value,fill = type)) + geom_histogram(position = "identity",alpha=.7,bins = 20) +
  xlab('# IDs per single cell') + ylab('# single cells') + theme_classic() + 
  scale_fill_manual(values = c( "#330000","#2BD9D9")) + coord_cartesian(xlim = c(0,2000))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  theme(axis.title.x = element_text(size = 30),axis.title.y = element_text(size = 30),axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=0), #change legend title font size
        legend.text = element_text(size=25),
        legend.position=c(.70, .75)) 

#ggsave('/Users/andrewleduc/Desktop/nPOP_Paper/PDFs/pep_prot_plex.pdf',plot = pep_prot_plex,height = unit(6, "in"),width = unit(7, "in"),dpi = 1000)



#Normalizations... should have used Harrisons functions but simple stuffm nothing new

for(i in 1:ncol(Ms1_extract2)){
  Ms1_extract2[,i] <- Ms1_extract2[,i]/median(Ms1_extract2[,i],na.rm = T)
}
for(i in 1:nrow(Ms1_extract2)){
  Ms1_extract2[i,] <- Ms1_extract2[i,]/mean(Ms1_extract2[i,],na.rm = T)
}
Ms1_extract2 <- log2(Ms1_extract2)
Ms1_extract2 <-  hknn(Ms1_extract2, k.t)
for(i in 1:ncol(Ms1_extract2)){
  Ms1_extract2[,i] <- Ms1_extract2[,i]-median(Ms1_extract2[,i],na.rm = T)
}
for(i in 1:nrow(Ms1_extract2)){
  Ms1_extract2[i,] <- Ms1_extract2[i,]-mean(Ms1_extract2[i,],na.rm = T)
}



#mapping labels for correting them with combat
evnew_prots <- as.data.frame(Ms1_extract2)
evnew_prots2 <- evnew_prots %>% dplyr::select(grep('raw.8',colnames(evnew_prots)))
d8 <- colnames(evnew_prots2)
evnew_prots2 <- evnew_prots %>% dplyr::select(grep('raw.4',colnames(evnew_prots)))
d4 <- colnames(evnew_prots2)
evnew_prots2 <- evnew_prots %>% dplyr::select(grep('raw.0',colnames(evnew_prots)))
d0 <- colnames(evnew_prots2)


#assign labels
batch_label  <- as.data.frame(colnames(Ms1_extract2))
batch_label$lab <- NA
colnames(batch_label)[1] <- 'cells'
batch_label$lab[batch_label$cells %in% d0] <- 'd0'
batch_label$lab[batch_label$cells %in% d4]  <- 'd4'
batch_label$lab[batch_label$cells %in% d8]  <- 'd8'
batch_label$run <- substring(batch_label$cells,12,13)


# correct out label bias (not too bad in general but does need correcting)
Ms1_extract3 <- ComBat(Ms1_extract2, batch=factor(batch_label$lab))

for(i in 1:ncol(Ms1_extract3)){
  Ms1_extract3[,i] <- Ms1_extract3[,i]-median(Ms1_extract3[,i],na.rm = T)
}
for(i in 1:nrow(Ms1_extract3)){
  Ms1_extract3[i,] <- Ms1_extract3[i,]-mean(Ms1_extract3[i,],na.rm = T)
}

# NAing imputed values post batch correction
Ms1_extract3_no_imp <- Ms1_extract3
Ms1_extract3_no_imp[is.na(Ms1_extract2_noImp)==T] <- NA

#write.csv(Ms1_extract3, file="../dat/Proccessed_data/plexDIA/Proteins_proc.csv")
#write.csv(Ms1_extract3_no_imp, file="../dat/Proccessed_data/plexDIA/Proteins_proc_noImp.csv")


#PCA with NAs
X.m <-  Ms1_extract3_no_imp

#r1<-cor(t(X.m))
#rsum<-rowSums(r1^2)
# Dot product of each protein correlation vector with itself
#X.m <- diag(rsum) %*%  X.m

pca.imp.cor <- cor(X.m, use = 'pairwise.complete.obs',method = c('pearson'))


# PCA
sc.pca<-eigen(pca.imp.cor)

scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$cells<-colnames(pca.imp.cor)

scx$label <- NA
scx$label[scx$cells %in% d0] <- 'd0'
scx$label[scx$cells %in% d4] <- 'd4'
scx$label[scx$cells %in% d8] <- 'd8'
scx$run <- batch_label$run
scx$yecol <- 'yellow'
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100
percent_var[1]
mel_plex <- ggplot(scx,aes(x=PC1,y=PC2,color = yecol,fill=yecol,shape = yecol))+theme_classic() +geom_point(size = 5,alpha = .8,stroke = 2)+
  xlab('PC1(14%)')+ylab('PC2(9%)')+
  theme(axis.title.x = element_text(size = 30),axis.title.y = element_text(size = 30),axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.key.size = unit(0, 'cm'), #change legend key size
        legend.key.height = unit(0, 'cm'), #change legend key height
        legend.key.width = unit(0, 'cm'), #change legend key width
        legend.title = element_text(size=0), #change legend title font size
        legend.text = element_text(size=0),
        legend.position=c(.2, .85)) +scale_fill_manual(values = c('#F6BE00'))+scale_color_manual(values = c('#000000'))+
  rremove('legend')+scale_shape_manual(values = c(22))

mel_plex_lab <- ggplot(scx,aes(x=PC1,y=PC2,color = yecol,fill=label,shape = yecol))+theme_classic() +geom_point(size = 5,alpha = .8,stroke = 2)+
  xlab('PC1(14%)')+ylab('PC2(9%)')+
  theme(axis.title.x = element_text(size = 30),axis.title.y = element_text(size = 30),axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.key.size = unit(0, 'cm'), #change legend key size
        legend.key.height = unit(0, 'cm'), #change legend key height
        legend.key.width = unit(0, 'cm'), #change legend key width
        legend.title = element_text(size=0), #change legend title font size
        legend.text = element_text(size=0),
        legend.position=c(.2, .85)) +scale_color_manual(values = c('#000000'))+
  rremove('legend')+scale_shape_manual(values = c(22))

#ggsave('/Users/andrewleduc/Desktop/nPOP_Paper/PDFs/mel_plex.pdf',plot = mel_plex,height = unit(6, "in"),width = unit(7, "in"),dpi = 1000)
#ggsave('/Users/andrewleduc/Desktop/nPOP_Paper/PDFs/Ext_data_SCQC/mel_plex_lab.pdf',plot = mel_plex_lab,height = unit(6, "in"),width = unit(6, "in"),dpi = 1000)


# Finding differential proteins between two clusters
scx2 <- scx %>% filter(PC2 < -.1)
pvals <- c()
meds <- c()
prot <- c()
PG_CORE3 <- Ms1_extract3[,colnames(Ms1_extract3)%in% scx2$cells]
PG_CORE4 <- Ms1_extract3[,!colnames(Ms1_extract3)%in% scx2$cells]
for (i in 1:nrow(Ms1_extract3)){
  
  pvals <- c(pvals,wilcox.test(as.matrix(PG_CORE3[i,]), as.matrix(PG_CORE4[i,]), alternative = "two.sided")$p.value)
  prot <- c(prot,rownames(Ms1_extract3)[i])
  meds <- c(meds, median(as.matrix(PG_CORE3[i,]),na.rm = T)- median(as.matrix(PG_CORE4[i,]),na.rm = T))
  
}

df <- as.data.frame(pvals)
df$prot <- prot
df$FC <- meds
df$qval <- p.adjust(df$pvals,method = 'BH')
df <- df %>% filter(qval < .01)

#write.csv(df,'/Users/andrewleduc/Desktop/nPOP_Paper/dat/Output_tables/Sup_pop_diff_plexDIA.csv')

a <- df %>% filter(df$FC > 1)
b <- df %>% filter(df$FC < -2)

df2 <- df %>% filter(prot %in% rownames(PG_ALL2))

df2 <- df2[order(df2$FC),]
df <- df2[order(df2$FC),]

convert_diff <- df %>% left_join(convert, by= c('prot' = 'Entry') )

convert_diff <- convert_diff %>% filter(prot %in% c(glyco,mito))
write_csv(convert_diff,'/Users/andrewleduc/Desktop/marker_proteins.csv')

one <- colMeans(PG_CORE4[rownames(PG_CORE4) %in% a$prot,])
two <- colMeans(PG_CORE4[rownames(PG_CORE4) %in% b$prot,])
cor(one,two)

cor.test(scx$one[scx$PC1 < .05],scx$two[scx$PC1 < .05],method = "spearman",exact=FALSE)$p.value
cor(scx$one[scx$PC1 < .05],scx$two[scx$PC1 < .05])



# Reading in pSCoPE data for integration

t7 <- read.csv('../dat/Proccessed_data/pSCoPE/t7.csv',row.names = 1)
t7_Melanomas <- t7 %>% dplyr::select(intersect(colnames(t7),meta_mel$id))
t7_Melanomas_no_imp <-t7_Melanomas
t7_Melanomas <- t7_Melanomas %>% filter(rownames(t7_Melanomas)%in% rownames(Ms1_extract3))
t7_Melanomas <- as.matrix(t7_Melanomas)
t7_Melanomas <- cr_norm_log(t7_Melanomas)
t7_Melanomas <-  hknn(t7_Melanomas, k.t)
for(i in 1:ncol(t7_Melanomas)){
  t7_Melanomas[,i] <- t7_Melanomas[,i]-median(t7_Melanomas[,i],na.rm = T)
}
for(i in 1:nrow(t7_Melanomas)){
  t7_Melanomas[i,] <- t7_Melanomas[i,]-mean(t7_Melanomas[i,],na.rm = T)
}

Ms1_extract3 <- Ms1_extract3[rownames(t7_Melanomas),]
Ms1_extract3_no_imp <- Ms1_extract3_no_imp[rownames(t7_Melanomas),]
Ms1_extract3_no_imp <- filt.mat.cr(Ms1_extract3_no_imp,.8,.99)
t7_Melanomas_no_imp <- t7_Melanomas_no_imp[rownames(Ms1_extract3_no_imp),]

PG_ALL <- cbind(Ms1_extract3,t7_Melanomas)

batch <- colnames(PG_ALL)
batch[batch %in% colnames(t7_Melanomas)] <- 'DDA'
batch[batch %in% colnames(Ms1_extract3)] <- 'DIA'
batch_all <- batch

dim(t7_Melanomas_no_imp)
plex <- rowMeans(Ms1_extract3_no_imp[,colnames(Ms1_extract3_no_imp) %in% scx2$cells],na.rm = T)-rowMeans(Ms1_extract3_no_imp[,!colnames(Ms1_extract3_no_imp) %in% scx2$cells],na.rm = T)
pscope <- rowMeans(t7_Melanomas_no_imp[,colnames(t7_Melanomas_no_imp) %in% scx2$cells],na.rm = T)-rowMeans(t7_Melanomas_no_imp[,!colnames(t7_Melanomas_no_imp) %in% scx2$cells],na.rm = T)
compare_prot_dist <- as.data.frame(pscope)
compare_prot_dist$plex <- plex
compare_prot_dist$prot <- rownames(Ms1_extract3_no_imp)
compare_prot_dist2 <- compare_prot_dist %>% filter(!prot %in% evnew_sp$Protein.Group)
compare_prot_dist3 <- compare_prot_dist %>% filter(prot %in% evnew_sp$Protein.Group)
compare_prot_dist$diff_pep <- NA
compare_prot_dist$diff_pep[compare_prot_dist$prot %in% compare_prot_dist2$prot] <- 'different'
compare_prot_dist$diff_pep[!compare_prot_dist$prot %in% compare_prot_dist2$prot] <- 'shared'

FC_comp <- ggplot(compare_prot_dist,aes(x = pscope,y = plex)) + geom_point(size = 3)+theme_bw()+ 
  ylab('log2(cluster A - cluster B) plexDIA') + xlab('log2(cluster A - cluster B) pSCoPE')+
  theme(axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20,color = 'black'),
        axis.text.y = element_text(size = 20,color = 'black'))

#ggsave('/Users/andrewleduc/Desktop/nPOP_Paper/PDFs/FC_comp.pdf',plot = FC_comp,height = unit(6, "in"),width = unit(7, "in"),dpi = 1000)
#ggsave('/Users/andrewleduc/Desktop/nPOP_Paper/PDFs/mel_plex.pdf',plot = mel_plex,height = unit(6, "in"),width = unit(7, "in"),dpi = 1000)


cor(plex,pscope,use = 'pairwise.complete.obs')

PG_ALL <- ComBat(PG_ALL, batch=factor(batch_all))

for(i in 1:ncol(PG_ALL)){
  PG_ALL[,i] <- PG_ALL[,i]-median(PG_ALL[,i],na.rm = T)
}
for(i in 1:nrow(PG_ALL)){
  PG_ALL[i,] <- PG_ALL[i,]-mean(PG_ALL[i,],na.rm = T)
}

#write.csv(PG_ALL, file="../dat/Proccessed_data/plexDIA/Integrated.csv")

X.m <- PG_ALL
r1<-cor(t(X.m))
rsum<-rowSums(r1^2)

# Dot product of each protein correlation vector with itself
X.m <- diag(rsum) %*%  X.m
pca.imp.cor <- cor(X.m, use = 'pairwise.complete.obs',method = c('pearson'))

# PCA
sc.pca<-eigen(pca.imp.cor)

scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$cells<-colnames(pca.imp.cor)
scx$batch <- NA
scx$batch[scx$cells %in% colnames(PG_CORE2)] <- 'plexDIA'
scx$batch[scx$cells %in% colnames(t7_Melanomas) ] <- 'pSCoPE'
scx$col_yel <- 'yellow'
pca_var <- sc.pca$values
percent_var<- pca_var/sum(pca_var)*100



integrate_mel <-ggplot(scx,aes(x=PC1,y=PC2,shape = batch,size = batch,color = batch,alpha = batch,fill = col_yel,stroke = batch))+theme_classic() +geom_point()+
  xlab('PC1(14%)') + ylab('PC2(5%)')+
  theme(axis.title.x = element_text(size = 30),axis.title.y = element_text(size = 30),axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=0), #change legend title font size
        legend.text = element_text(size=20),
        legend.position=c(.2, .85))+ scale_shape_manual(values = c(22,19)) + 
  scale_size_manual(values = c(5,3)) + scale_fill_manual(values = c('#F6BE00')) + scale_color_manual(values = c('#000000','#F6BE00'))+
  scale_alpha_manual(values = c(.7,.3))+rremove('legend') + scale_discrete_manual(aesthetics = "stroke",values = c(2,1))
  
#ggsave('/Users/andrewleduc/Desktop/nPOP_Paper/PDFs/integrate.pdf',plot = integrate_mel,height = unit(6, "in"),width = unit(7, "in"),dpi = 1000)

scx2 <- scx %>% filter(PC1 > .05)

PG_ALL2 <- PG_ALL[intersect(rownames(PG_ALL),c(a$prot,b_$prot)),]
PG_ALL2 <- cr_norm_log(PG_ALL2)
#PG_ALL3 <- ComBat(PG_ALL2, batch=factor(batch_all))


X.m <- PG_ALL2
pca.imp.cor <- cor(X.m, use = 'pairwise.complete.obs',method = c('pearson'))
# PCA
sc.pca<-eigen(pca.imp.cor)

scx<-as.data.frame(sc.pca$vectors)
colnames(scx)<-paste0("PC",1:ncol(scx))
scx$cells<-colnames(pca.imp.cor)

scx$one <- colMeans(PG_ALL2[rownames(PG_ALL2) %in% a$prot,])
scx$two <- colMeans(PG_ALL2[rownames(PG_ALL2) %in% b$prot,])



mito <- c('Q00059','P42704','P55084','P24752','P23368')
glyco <- c('P06733','P04075','P14618','P00558','P04406')



cor.test(colMeans(PG_ALL2[rownames(PG_ALL2) %in% try_n,!colnames(PG_ALL2)%in% scx2$cells]),
    colMeans(PG_ALL2[rownames(PG_ALL2) %in% glyco,!colnames(PG_ALL2)%in% scx2$cells]))$p.value

cor(colMeans(PG_ALL2[rownames(PG_ALL2) %in% try_n,!colnames(PG_ALL2)%in% scx2$cells]),
         colMeans(PG_ALL2[rownames(PG_ALL2) %in% glyco,!colnames(PG_ALL2)%in% scx2$cells]))

plot(colMeans(PG_ALL2[rownames(PG_ALL2) %in% try_n,!colnames(PG_ALL2)%in% scx2$cells]),
    colMeans(PG_ALL2[rownames(PG_ALL2) %in% glyco,!colnames(PG_ALL2)%in% scx2$cells]))


glyco_mito <- as.data.frame(colMeans(PG_ALL2[rownames(PG_ALL2) %in% try_n,!colnames(PG_ALL2)%in% scx2$cells]))
colnames(glyco_mito) <- 'Mito'
glyco_mito$Glyco <- colMeans(PG_ALL2[rownames(PG_ALL2) %in% glyco,!colnames(PG_ALL2)%in% scx2$cells])
glyco_mito$cell <- colnames(PG_ALL2)[!colnames(PG_ALL2)%in% scx2$cells]
glyco_mito$batch <- NA
glyco_mito$batch[glyco_mito$cell %in% colnames(Ms1_extract3)] <- 'plexDIA'
glyco_mito$batch[glyco_mito$cell %in% colnames(t7_Melanomas) ] <- ' pSCoPE'
glyco_mito <- glyco_mito %>% filter(batch != ' pSCoPE')
glyco_mito <- glyco_mito[order(glyco_mito$batch),]



pca_var <- sc.pca$values

# plot heatmap glycolosis oxPHOS
mappy <-(cor(t(PG_ALL[rownames(PG_ALL) %in% c(glyco,mito),!colnames(PG_ALL) %in% scx2$cells])))#,colMeans(PG_ALL2[rownames(PG_ALL2) %in% b$prot,!colnames(PG_ALL2) %in% scx2$cells]))
mappy[mappy >.5] <- .5
Heatmap(mappy)
Heatmap(cor(t(PG_ALL2[rownames(PG_ALL2) %in% c(glyco,mito),!colnames(PG_ALL2) %in% scx2$cells])))#,colMeans(PG_ALL2[rownames(PG_ALL2) %in% b$prot,!colnames(PG_ALL2) %in% scx2$cells]))

aa <- Heatmap(cor(t(PG_ALL2[,!colnames(PG_ALL2) %in% scx2$cells])))#,colMeans(PG_ALL2[rownames(PG_ALL2) %in% b$prot,!colnames(PG_ALL2) %in% scx2$cells]))



#cor.test(colMeans(PG_ALL2[rownames(PG_ALL2) %in% a$prot,!colnames(PG_ALL2) %in% scx2$cells]),colMeans(PG_ALL2[rownames(PG_ALL2) %in% b$prot,!colnames(PG_ALL2) %in% scx2$cells]))$p.value




scx$batch <- NA
scx$batch[scx$cells %in% colnames(Ms1_extract3)] <- 'plexDIA'
scx$batch[scx$cells %in% colnames(t7_Melanomas) ] <- 'pSCoPE'
scx<- scx[order(desc(scx$batch)),]

ggplot(scx,aes(x=PC1,y=PC2,color = batch))+theme_classic() +geom_point()#+ scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                                                                         #                     high = "red",name = '')
scx$one[scx$one > .5] <- .5
scx$two[scx$two < -.6] <- -.6




p1 <-ggplot(scx,aes(x=PC1,y=PC2,fill = one_r,shape=batch,size=batch,stroke = batch,alpha = batch))+theme_classic() +geom_point()+  xlab('PC1(11%)')+ylab('PC2(7%)')+
  scale_fill_gradient2(midpoint = -.02, low = "blue", mid = "white",high = "red",name = '')+
  theme(axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=5),
        legend.position=c(.2, .85),plot.margin = margin(1, 1, 1, 1, "cm"))+ 
  scale_shape_manual(values = c(22,21)) + scale_size_manual(values = c(5,3)) +rremove('legend')+
  scale_discrete_manual(aesthetics = "stroke",values = c(.3,0)) + scale_alpha_manual(values = c(1,.8))

  
p2 <-ggplot(scx,aes(x=PC1,y=PC2,fill = two,shape=batch,size=batch,stroke = batch,alpha = batch))+theme_classic() +geom_point()+  xlab('PC1(11%)')+ylab('PC2(7%)')+
  scale_fill_gradient2(midpoint = .05, low = "blue", mid = "white",high = "red",name = '')+
  theme(axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=5),
        legend.position=c(.2, .85),plot.margin = margin(1, 1, 1, 1, "cm"))+
  scale_shape_manual(values = c(22,21)) + scale_size_manual(values = c(5,3)) +
  scale_discrete_manual(aesthetics = "stroke",values = c(.3,0)) + scale_alpha_manual(values = c(1,1))

#p2 + p1


p3 <- ggplot(scx,aes(x=PC1,y=one, fill = batch,shape = batch,size = batch))+theme_classic() +geom_point(alpha= 1)+ ylab('Cluster B markers')+xlab('')+
  scale_shape_manual(values = c(22,21))+scale_fill_manual(values=c('#FFFFFF','#000000'))+
  scale_size_manual(values = c(5,3))+
  theme(axis.title.x = element_text(size = 30),axis.title.y = element_text(size = 20),axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size = 20),
        legend.key.size = unit(0, 'cm'), #change legend key size
        legend.key.height = unit(0, 'cm'), #change legend key height
        legend.key.width = unit(0, 'cm'), #change legend key width
        legend.title = element_text(size=0), #change legend title font size
        legend.text = element_text(size=0),
        legend.position=c(.2, .85),plot.margin = margin(1, 1, 1, 1, "cm"))  +rremove('legend')



p4 <- ggplot(scx,aes(x=PC1,y=two, fill = batch,shape = batch,size = batch))+theme_classic() +geom_point(alpha= 1)+ ylab('Cluster B markers')+xlab('')+
  scale_shape_manual(values = c(22,21))+scale_fill_manual(values=c('#FFFFFF','#000000'))+
  scale_size_manual(values = c(5,3))+
  theme(axis.title.x = element_text(size = 30),axis.title.y = element_text(size = 20),axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size = 20),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=0), #change legend title font size
        legend.text = element_text(size=0),
        legend.position=c(.2, .85),plot.margin = margin(1, 1, 1, 1, "cm")) +rremove('legend')

ggsave('/Users/andrewleduc/Desktop/stich_pic.pdf',plot = (p3+p4)/(p1+p2),height = unit(8, "cm"),width = unit(12, "cm"),dpi = 1000)
#ggsave('/Users/andrewleduc/Desktop/stich_pic.pdf',plot = (p3+p4)(p1+p2),height = unit(12, "cm"),width = unit(12, "cm"))

(p2+p1)
(p4+p3)



### plotting distrobutions for individual proteins


t7_hold_df <- as.data.frame(t(t7_hold[c('P24539','P49755'),]))
t7_hold_df$subPOP <- NA
t7_hold_df$subPOP[rownames(t7_hold_df)%in% meta_mel_C1$id] <- 'Cluster A'
t7_hold_df$subPOP[rownames(t7_hold_df)%in% meta_mel_C2$id] <- 'Cluster B'
D1 <-ggplot(t7_hold_df,aes(y = P24539,x =P49755,size = subPOP,color = subPOP)) + geom_point(alpha = 1) +
  theme_bw() + xlab('ATP synthase complex subunit B1') + ylab('Transmembrane emp24 protein 10')+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+ 
  scale_color_manual(values = c("#dcdcdc","#808080"))+scale_size_manual(values = c(3,6))+rremove('legend')

t7_hold_df2 <- as.data.frame(t(t7_hold[c('Q06830','P40121'),]))
t7_hold_df2$subPOP <- NA
t7_hold_df2$subPOP[rownames(t7_hold_df2)%in% meta_mel_C1$id] <- 'Cluster A'
t7_hold_df2$subPOP[rownames(t7_hold_df2)%in% meta_mel_C2$id] <- 'Cluster B'
D2 <-ggplot(t7_hold_df2,aes(y = Q06830,x =P40121,size = subPOP,color = subPOP)) + geom_point(alpha = 1) +
  theme_bw() + xlab('Peroxiredoxin-1') + ylab('Macrophage-capping protein')+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+ 
  scale_color_manual(values = c("#dcdcdc","#808080"))+scale_size_manual(values = c(3,6))+rremove('legend')
cor(t7_hold_df2[!rownames(t7_hold_df2)%in% scx2$cells,],use = 'pairwise.complete.obs')

t7_hold_df3 <- as.data.frame(t(Ms1_extract3[c('P11021','P07339'),]))
t7_hold_df3$subPOP <- NA
t7_hold_df3$subPOP[rownames(t7_hold_df3)%in% scx2$cells] <- 'Cluster A'
t7_hold_df3$subPOP[!rownames(t7_hold_df3)%in% scx2$cells] <- 'Cluster B'
D3 <-ggplot(t7_hold_df3,aes(y = P11021,x =P07339,size = subPOP,fill = subPOP)) + geom_point(shape = 22,stroke = .2) +
  theme_bw() + xlab('Cathepsin D') + ylab('Endoplasmic reticulum chaperone BiP')+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+ 
  scale_fill_manual(values = c("#808080","#dcdcdc"))+scale_size_manual(values = c(6,3))+rremove('legend')

t7_hold_df4 <- as.data.frame(t(Ms1_extract3[c('P62857','P19367;P19367-2;P19367-3;P19367-4'),]))
colnames(t7_hold_df4)[2] <- 'P19367'

t7_hold_df4$subPOP <- NA
t7_hold_df4$subPOP[rownames(t7_hold_df4)%in% scx2$cells] <- 'Cluster A'
t7_hold_df4$subPOP[!rownames(t7_hold_df4)%in% scx2$cells] <- 'Cluster B'
D4 <-ggplot(t7_hold_df4,aes(y = P62857,x =P19367,size = subPOP,fill = subPOP)) + geom_point(shape = 22,stroke = .2) +
  theme_bw() + ylab('40S ribosomal protein S28') + xlab('Hexokinase-1')+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+ 
  scale_fill_manual(values = c("#808080","#dcdcdc"))+scale_size_manual(values = c(6,3))+rremove('legend')


cor(t7_hold_df4[,],use = 'pairwise.complete.obs')
cor(t7_hold_df4,use = 'pairwise.complete.obs')



ggsave('/Users/andrewleduc/Desktop/joint_dist.pdf',plot = D4 ,height = unit(6, "in"),width = unit(5.5, "in"),dpi = 1000)

