source("/Users/andrewleduc/Desktop/AnalysisScripts/sourceme.R")
### Read in evidence file
##


ev<-read.delim('/Volumes/GoogleDrive/My Drive/MS/cellenONE/aleduc/Label_efficiency_check/evidence.txt')
unique(ev$Raw.file)
#AL0058
#AL080

ev<-ev[ev$PEP<0.02,]
ev<-ev[ev$PIF > 0.7,]
ev<-ev[ev$Reverse !='+',]
ev<-ev[ev$Potential.contaminant!='+',]


### Calculate number of sites TMT could label (assuming only labeling K, N-terminus)
##
unique(ev$Raw.file)
ev<-num.TMT(ev,'Sequence')
ev<-ev[ev$num.TMT != 0,]
ev$zero <- ev$num.TMT-ev$Acetyl..Protein.N.term.
ev<-ev[ev$zero != 0,]
ev$rat <- (ev$Variable_TMTpro16plex.Lys+ev$Variable_TMTpro16plex.N.term)/(ev$num.TMT-ev$Acetyl..Protein.N.term.)

#ev$Variable_TMTpro16plex.N.term
#ev$TMTPro_Nter_LE


unique(ev$Raw.file)[4]
ev2 <- ev %>% filter(Raw.file == unique(ev$Raw.file)[2])
#ev2 <- ev %>% filter(Raw.file == 'AL078')

### This is youre answer below

mean(ev2$rat, na.rm = T)


