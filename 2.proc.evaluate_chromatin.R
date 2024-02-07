
source('0.basics.R')


# Process chromatin accessibility with neoant -----------------------------

atac.proms <- readRDS('ATAC/perSample_openclose_data.rds')
atac.proms$sample <- gsub('EPICC_','',atac.proms$sample)
samplesATAC <- unique(epicc.df$Patient[epicc.df$Sample %in% atac.proms$sample])

# Generate all-sample matrix for mutations and samples with promoter accessible/not info
atac.muts.total <- data.frame(matrix(vector(),ncol=8))
for (pat in samplesATAC){
  regs <- intersect(atac.proms$sample, epicc.df$Sample[epicc.df$Tissue=='Cancer' & epicc.df$Patient==pat])
  atac.proms.pat <- atac.proms[atac.proms$sample %in% regs, ]
  muts.pat <- read.delim(paste0('Subclonality/',pat,'.total.extended.txt'))
  muts.pat <- subset(muts.pat, !is.na(Neoantigen))
  
  atac.muts.pat <- data.frame(matrix(vector(),ncol=5))
  names(atac.muts.pat) <- c('Patient','numMut','numAccMut','numNotAccMut','mutID')
  for (mut in muts.pat$mutID){
    gene <- muts.pat$GeneID[muts.pat$mutID==mut]
    if (!(gene %in% atac.proms.pat$gene)){next}
    samp.mut <- regs[muts.pat[muts.pat$mutID==mut,regs]==1] #how many samples have it mutated
    acc.atac <- sum(atac.proms.pat$atac_promotor_accessible[atac.proms.pat$gene==gene & atac.proms.pat$potentially_somatic==1 & atac.proms.pat$sample %in% samp.mut]==1, na.rm=T) #how many have it acc. promoter in ATAC
    nonacc.atac <- sum(atac.proms.pat$atac_promotor_accessible[atac.proms.pat$gene==gene & atac.proms.pat$potentially_somatic==1 & atac.proms.pat$sample %in% samp.mut]==0,na.rm=T) #how many do not have acc. promoter in ATAC
    atac.muts.pat[nrow(atac.muts.pat)+1,] <- c(pat, length(samp.mut),acc.atac, nonacc.atac, mut)
  }
  
  atac.muts.pat[,2:4] <- apply(atac.muts.pat[,2:4], 2, as.numeric)
  atac.muts.pat[,c('Category','Neoantigen','MutType')] <- muts.pat[match(atac.muts.pat$mutID, muts.pat$mutID),c('Category','Neoantigen','MutType')]
  atac.muts.total <- rbind(atac.muts.total, atac.muts.pat)
}

write.table(atac.muts.total, file='ATAC/PromoterAccessSomaticForOR_NA_nonNA_allcancers.txt', sep='\t',row.names=F, quote=F)



# Get SCAAlosses ----------------------------------------------------------

# Read in analysis results from SCAA evaluation, carried out clonally for each cancer
scaa <- readRDS('ATAC/analysis_results_raw.rds')

# 1) For each mutation, check if there is a SCAA loss in the promoter of that gene

lossTH <- -1 # threshold for the fold change to call something a loss
scaa.df <- scaa$fc*scaa$sig #modify table to only contain significant fold changes, 0 otherwise
scaa.df <- scaa.df[,!grepl('adenoma',names(scaa.df))]; names(scaa.df) <- sub('.pure','',names(scaa.df))
samplesATAC <- intersect(samples, names(scaa.df))

# Generate all-sample matrix for mutations and samples with promoter SCAA loss
scaa.muts.total <- data.frame(matrix(vector(),ncol=8))
for (pat in samplesATAC){
  scaa.pat <- scaa.df[, pat,drop=F]
  muts.pat <- read.delim(paste0('Subclonality/',pat,'.total.extended.txt'))
  muts.pat <- subset(muts.pat, !is.na(Neoantigen))
  
  scaa.muts.pat <- data.frame(matrix(vector(),ncol=5))
  names(scaa.muts.pat) <- c('Patient','SLossProm','SLossEnh','GeneID','mutID')
  # iterate for each mutation
  for (mut in muts.pat$mutID){
    gene <- muts.pat$GeneID[muts.pat$mutID==mut] # get the ID of the gene the mutation is located in
    geneRowsP <- which(grepl(gene, scaa$flanking_genes) & scaa$set=='proximal') #rows in scaa.df that are for promoters of that gene
    geneRowsE <- which(grepl(gene, scaa$flanking_genes) & scaa$set=='distal') #rows in scaa.df that are for enhancers of that gene
    scaalossP <- NA; scaalossE <- NA
    if (length(geneRowsP)>0){
      scaalossP <- sum(scaa.pat[geneRowsP,]< lossTH) #whether there is a loss SCAA associated with promoter of the gene
    }
    if (length(geneRowsE)>0){
      scaalossE <- sum(scaa.pat[geneRowsE,]< lossTH) #whether there is a loss SCAA associated with enhancer of the gene
    }
    if (is.na(scaalossP) & is.na(scaalossE)){next}
    scaa.muts.pat[nrow(scaa.muts.pat)+1,] <- c(pat, scaalossP, scaalossE, gene, mut)
  }
  
  scaa.muts.pat[,2:3] <- apply(scaa.muts.pat[,2:3], 2, as.numeric)
  scaa.muts.pat[,c('Category','Neoantigen','MutType')] <- muts.pat[match(scaa.muts.pat$mutID, muts.pat$mutID),c('Category','Neoantigen','MutType')]
  scaa.muts.total <- rbind(scaa.muts.total, scaa.muts.pat)
}

write.table(scaa.muts.total, file='ATAC/SCAALoss_NA_nonNA_allcancers.txt', sep='\t',row.names=F, quote=F)


# 2) For each SCAA promoter loss, check if that gene (in that patient) has a mutation/NA attached to it (and if yes, what category?)

scaa <- readRDS('ATAC/analysis_results_raw.rds')

lossTH <- -1 # threshold for the fold change to call something a loss
scaa.df <- scaa$fc*scaa$sig #table only containing significant fold changes, 0 otherwise
scaa.df <- scaa.df[,!grepl('adenoma',names(scaa.df))]; names(scaa.df) <- sub('.pure','',names(scaa.df))
scaa.df <- 1*(scaa.df< lossTH)
samplesATAC <- intersect(samples, colnames(scaa.df))

# Generate all-sample matrix for promoter SCAA loss and if there are genes attached to it
scaa.muts.total <- data.frame(matrix(vector(),ncol=4))
for (pat in samplesATAC){
  scaa.pat <- scaa.df[, pat,drop=F]
  muts.pat <- read.delim(paste0('Subclonality/',pat,'.total.extended.txt'))
  muts.pat <- subset(muts.pat, !is.na(Neoantigen) & MutType=='SNV') #limit to SNVs
  
  scaa.muts.pat <- data.frame(matrix(vector(),ncol=4))
  names(scaa.muts.pat) <- c('Patient','SCAA','Neoantigen','Category')
  for (ind in which(scaa.pat>0)){
    gene <- unlist(scaa$flanking_genes[ind])
    if (scaa$set[ind]=='proximal'){
      lossmuts <- paste0(muts.pat[muts.pat$GeneID %in% gene,'Neoantigen'], collapse=';') #find all mutations and their NA-type associated with this SCAA
      losscats <- paste0(muts.pat[muts.pat$GeneID %in% gene,'Category'], collapse=';') #find all mutations and their subclonality
    }else{next}
    scaa.muts.pat[nrow(scaa.muts.pat)+1,] <- c(pat, row.names(scaa.pat)[ind], lossmuts, losscats)
  }
  scaa.muts.total <- rbind(scaa.muts.total, scaa.muts.pat)
}

# Divide up neo/non-neo categories into non-NA, both NA and non-NA, and SB NA
scaa.muts.total$Neoantigen <- ifelse(scaa.muts.total$Neoantigen=='','',
                                     ifelse(!grepl('SB',scaa.muts.total$Neoantigen),
                                            ifelse(!grepl('Non',scaa.muts.total$Neoantigen),NA,'Non-Neo'),
                                            ifelse(grepl('Non',scaa.muts.total$Neoantigen),'Both','SB Neo' )))
write.table(scaa.muts.total, file='ATAC/SCAALoss_mutsPerSCAA_allcancers.txt', sep='\t',row.names=F, quote=F)

