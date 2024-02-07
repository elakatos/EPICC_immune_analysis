source('0.basics.R')

sampleList <- epicc.df$Sample[epicc.df$Type=='Deep']

# Evaluate evolutionary close regions, subclonal escape -------------------

# for each cancer, this file contains burden info, plus a mutation/not column based on the subclonal mutation we explore
measures.files <- list.files('Subclonality/Subclonal_escapes/',pattern='NAmeasures.txt', full.names = T)
measures.total <- data.frame(matrix(vector(),ncol=12))

# define if we only consider deep sequenced samples (or low-pass as well) and if limit to evolutionarily close samples
deep <- T ; close <- T

# go through patient-by-patient and define which biopsies to consider
# then normalise proportional burden and immune dNdS to the considered samples
samp <- 'C518'
measures.df <- read.delim(grep(samp,measures.files,value=T))
measures.df <- subset(measures.df, Region!='C')
if(close) measures.df <- subset(measures.df, Region!='C' & Region!='D')
if(deep) measures.df <- subset(measures.df, Sample %in% sampleList)
measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')] <- apply(measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')],2,function(x) x/mean(x,na.rm=T))
measures.total <- rbind(measures.total, measures.df)

samp <- 'C524'
measures.df <- read.delim(grep(samp,measures.files,value=T))
measures.df$Mut[measures.df$Mut=='fs'] <- 'mutation'; measures.df$Mut[measures.df$Mut=='wt'] <- 'no mut' #define the frameshift as mutation, vs wild-type
if(close) measures.df <- subset(measures.df, Mut!='LOH')
if(deep) measures.df <- subset(measures.df, Sample %in% sampleList)
measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')] <- apply(measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')],2,function(x) x/mean(x,na.rm=T))
measures.total <- rbind(measures.total, measures.df)

samp <- 'C536'
measures.df <- read.delim(grep(samp,measures.files,value=T))
if(deep) measures.df <- subset(measures.df, Sample %in% sampleList)
measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')] <- apply(measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')],2,function(x) x/mean(x,na.rm=T))
measures.total <- rbind(measures.total, measures.df)

samp <- 'C537'
measures.df <- read.delim(grep(samp,measures.files,value=T))
measures.df <- subset(measures.df, Region!='A')
if(close) measures.df <- subset(measures.df, Region=='C')
if(deep) measures.df <- subset(measures.df, Sample %in% sampleList)
measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')] <- apply(measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')],2,function(x) x/mean(x,na.rm=T))
measures.total <- rbind(measures.total, measures.df)

samp <- 'C543'
measures.df <- read.delim(grep(samp,measures.files,value=T))
if(deep) measures.df <- subset(measures.df, Sample %in% sampleList)
measures.df <- subset(measures.df, Sample!='C543_A1_G9')
measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')] <- apply(measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')],2,function(x) x/mean(x,na.rm=T))
measures.total <- rbind(measures.total, measures.df)

samp <- 'C548'
measures.df <- read.delim(grep(samp,measures.files,value=T))
if(deep) measures.df <- subset(measures.df, Sample %in% sampleList)
# since all have mutation, define higher impact or extra mutation as distinguished subclonal immune escape
measures.df$Mut[measures.df$Mut=='mut'] <- 'no mut'; measures.df$Mut[measures.df$Mut=='extra mut' | measures.df$Mut=='higher impact'] <- 'mutation'
if(close)measures.df <- subset(measures.df, Region!='C' & Region!='D')
measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')] <- apply(measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')],2,function(x) x/mean(x,na.rm=T))
measures.total <- rbind(measures.total, measures.df)

samp <- 'C552'
measures.df <- read.delim(grep(samp,measures.files,value=T))
if(deep) measures.df <- subset(measures.df, Sample %in% sampleList)
# since all have mutation, define extra mutation as distinguished subclonal immune escape
measures.df$Mut[measures.df$Mut=='mut'] <- 'no mut'; measures.df$Mut[measures.df$Mut=='extra mut'] <- 'mutation' 
measures.df <- subset(measures.df, Region!='F')
measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')] <- apply(measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')],2,function(x) x/mean(x,na.rm=T))
measures.total <- rbind(measures.total, measures.df)

samp <- 'C554'
measures.df <- read.delim(grep(samp,measures.files,value=T))
if(deep) measures.df <- subset(measures.df, Sample %in% sampleList)
if(close) measures.df <- subset(measures.df, Region=='B')
measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')] <- apply(measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')],2,function(x) x/mean(x,na.rm=T))
measures.total <- rbind(measures.total, measures.df)

samp <- 'C559'
measures.df <- read.delim(grep(samp,measures.files,value=T))
measures.df$Mut[measures.df$Mut=='no LOH'] <- 'no mut'; measures.df$Mut[measures.df$Mut=='full LOH'] <- 'LOH'; #add a category for LOH
if(close) measures.df <- subset(measures.df, Region=='A' | (Region == 'D' & !(Sample %in% c('C559_D1_G5','C559_D1_G9'))))
if(deep) measures.df <- subset(measures.df, Sample %in% sampleList)
measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')] <- apply(measures.df[,c('SNVBurden','PropBurden','SNVBurden_sc','PropBurden_sc','imm_dNdS')],2,function(x) x/mean(x,na.rm=T))
measures.total <- rbind(measures.total, measures.df)

write.table(measures.total, file='Subclonality/Subclonal_escape_normalised.deepWGSclose.txt', quote=F, row.names=F, sep='\t')
