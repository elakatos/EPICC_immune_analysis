
source('0.basics.R')



# Total burden per biopsy -------------------------------------------------

# use per-patient tables with all mutations to count different categories of (protein-changing) mutations
totalBurden.df <- data.frame(matrix(vector(),ncol=13))
for (samp in samples){
  pat.muts <- read.delim(paste0('Subclonality/',samp,'.total.extended.txt'))
  regs<- epicc.df$Sample[epicc.df$Patient==samp]
  pat.muts$CodingChange[is.na(pat.muts$CodingChange)] <- F
  pat.muts$MutType[is.na(pat.muts$MutType)] <- '-'
  pat.muts$Category[pat.muts$Category=='private_adenoma'] <- 'adenoma_only'
  
  burden.df <- data.frame(Sample = regs,
                          TotalCCSNV = as.numeric(colSums(pat.muts[pat.muts$MutType=='SNV',regs])),
                          TotalCCFS = as.numeric(colSums(pat.muts[pat.muts$MutType=='FS',regs])),
                          TotalWGSTruncal = as.numeric(colSums(pat.muts[pat.muts$Category=='truncal',regs])),
                          TotalWGSTruncalFull = as.numeric(colSums(pat.muts[pat.muts$Category=='truncal_full',regs])),
                          TotalCodingChangeTruncal = as.numeric(colSums(pat.muts[pat.muts$Category=='truncal' & pat.muts$CodingChange,regs])),
                          TotalCodingChangeTruncalFull = as.numeric(colSums(pat.muts[pat.muts$Category=='truncal_full' & pat.muts$CodingChange,regs])),
                          TotalCCSNVTruncal = as.numeric(colSums(pat.muts[pat.muts$Category=='truncal' & pat.muts$MutType=='SNV',regs])),
                          TotalCCSNVTruncalFull = as.numeric(colSums(pat.muts[pat.muts$Category=='truncal_full' & pat.muts$MutType=='SNV',regs])),
                          TotalCCFSTruncal = as.numeric(colSums(pat.muts[pat.muts$Category=='truncal' & pat.muts$MutType=='FS',regs])),
                          TotalCCFSTruncalFull = as.numeric(colSums(pat.muts[pat.muts$Category=='truncal_full' & pat.muts$MutType=='FS',regs])),
                          TotalCCSNVAdonly = as.numeric(colSums(pat.muts[pat.muts$Category=='adenoma_only' & pat.muts$MutType=='SNV',regs])),
                          TotalCCFSAdonly = as.numeric(colSums(pat.muts[pat.muts$Category=='adenoma_only' & pat.muts$MutType=='FS',regs])),
                          TotalCCSNVPrivate = as.numeric(colSums(pat.muts[pat.muts$Category=='private' & pat.muts$MutType=='SNV',regs])),
                          TotalCCFSPrivate = as.numeric(colSums(pat.muts[pat.muts$Category=='private' & pat.muts$MutType=='FS',regs]))
  )
  totalBurden.df <- rbind(totalBurden.df, burden.df)
}

write.table(totalBurden.df, file='Burden/Total_burden_clonality_allsample.txt', sep='\t', row.names=F, quote=F)

# Compute per-biopsy neoantigen burden ------------------------------------

# Compute burden for each biopsy (according to various filtering/selection criteria)
totalBurden.df <- read.delim('Burden/Total_burden_clonality_allsample.txt')
hlaAlt.df <- read.delim('Immune_escape/EPICC_escape_master_file.txt')

burdenPerReg.df <- data.frame(matrix(vector(),ncol=4))
for (samp in samples){
  #read in both SNV and FS neoantigen tables and filter for novel and SB neoantigens and compute (NA mutation) burden
  #(note that these files follow NeoPredPipe output format, but have been processed to remove HLA-related info, and add clonality information)
  epTable <- read.table(paste0('Burden/Neopred_raw/',samp,'.neoantigens.extended.txt'),
                        sep='\t', stringsAsFactors = F, header=T)
  epTable.indel <- read.table(paste0('Burden/Neopred_raw/',samp,'.neoantigens.Indels.extended.txt'),
                              sep='\t', stringsAsFactors = F, header=T)
  epTable <- subset(epTable, Novelty==1); epTable.indel <- subset(epTable.indel, Novelty==1)
  epTable <- subset(epTable, BindLevel=='SB'); epTable.indel <- subset(epTable.indel, BindLevel=='SB')
  epTable <- epTable[!duplicated(epTable$LineID),] # deduplicate for unique mutations
  epTable.indel <- epTable.indel[!duplicated(epTable.indel$LineID),]
  regs <- names(epTable[startsWith(names(epTable), samp)])
  # optionally filter for subclonality category when comparing clonal/subclonal mutations
  #epTable <- subset(epTable, (Category %in% c('truncal','truncal_full'))); epTable.indel <- subset(epTable.indel, (Category %in% c('truncal','truncal_full')))
  
  burdenPerReg.df <- rbind(burdenPerReg.df, data.frame(Patient = samp,
                                                       Sample=regs,
                                                       SNVBurden=as.numeric(colSums(epTable[,regs,drop=F], na.rm=T)),
                                                       IndelBurden=as.numeric(colSums(epTable.indel[,regs,drop=F], na.rm=T))))
  
}

burdenPerReg.df <- subset(burdenPerReg.df, Sample %in% epicc.df$Sample) #subset to only include biopsies passing QC
# add sample characteristics
burdenPerReg.df[,c('Purity','Tissue','SampleType')] <- epicc.df[match(burdenPerReg.df$Sample, epicc.df$Sample),c('Purity','Tissue','SampleType')]
burdenPerReg.df$MSI <- ifelse((burdenPerReg.df$Patient %in% msiList) & (burdenPerReg.df$Patient=='C516' | burdenPerReg.df$Tissue=='Cancer'), 'MSI','MSS')
# add total number of protein-changing SNV/FS mutations - or the same but with clonal/subclonal
burdenPerReg.df[,c('TotalSNV','TotalIndel')] <- totalBurden.df[match(burdenPerReg.df$Sample, totalBurden.df$Sample),c('TotalCCSNV','TotalCCFS')]
#burdenPerReg.df[,c('TotalSNV','TotalIndel')] <- totalBurden.df[match(burdenPerReg.df$Sample, totalBurden.df$Sample),c('TotalCCSNVTruncal','TotalCCFSTruncal')]


# Process immune dNdS -----------------------------------------------------
# read in raw output from SOPRANO and extract important information

dnds.df <-read.delim("Burden/Immune_dNdS_raw/ALLEPICC_SUMMARY_PASS.txt",header=F,
                     col.names = c("Barcode","coverage","ON_dnds","ON_lowci",
                                   "ON_highci","ON_muts","OFF_dnds",
                                   "OFF_lowci","OFF_highci","OFF_muts",
                                   "pval","ON_na","ON_Na","ON_ns","ON_Ns",
                                   "OFF_na","OFF_Na","OFF_ns","OFF_Ns",'V20','V21','V22')) %>% dplyr::select(where(~sum(!is.na(.x)) > 0))
dnds.df$Sample <- epicc.df$Sample[match(dnds.df$Barcode,epicc.df$Barcode)]
dnds.df$imm_dnds <- dnds.df$ON_dnds/dnds.df$OFF_dnds
dnds.df$imm_lowci <- dnds.df$ON_lowci/dnds.df$OFF_lowci
dnds.df$imm_highci <- dnds.df$ON_highci/dnds.df$OFF_highci
dnds.df <- dnds.df[!duplicated(dnds.df$Barcode),]

# Compute VAF -------------------------------------------------------------

getVAFmatrix <- function(tumGT, mutcalls){
  tumNR <- apply(tumGT,2, function(z) sapply(z, function(x) as.numeric(unlist(strsplit(x,':'))[5])))
  tumNV <- apply(tumGT,2, function(z) sapply(z, function(x) as.numeric(unlist(strsplit(x,':'))[6])))
  tumVAF <- tumNV/tumNR
  tumVAF <- tumVAF * mutcalls # update to 0 for mutations that were indicated non-present after phylo analysis
  return(tumVAF)
}

# Get VAF tables for each patient separately - use mutation tables to obtain annotation of mutations
for (samp in samples){
  samp.df <- epicc.df[epicc.df$Patient==samp,]
  
  pat.muts <- read.delim(paste0('Subclonality/All_mutations/',samp,'.total.extended.txt'))
  regs <- samp.df$Sample
  
  # VCF files are available from https://data.mendeley.com/datasets/7wx3chtsxx/2
  v <- paste0('VCF/',samp,'_GRCh38_platypus_filtered_annotated_notdbsnp.vcf.bgz.vcf')
  pat.vcf <- read.vcfR(v, verbose=F)
  vcf.regs <- colnames(pat.vcf@gt)[-1]
  vcf.regs <- gsub('EPICC_','',vcf.regs) # format name of VCF GT columns to match sample names
  vcf.regs <- sub('_D1$','',vcf.regs)
  colnames(pat.vcf@gt)[-1] <- vcf.regs
  vcf.mutID <- paste0(pat.vcf@fix[,'CHROM'],':',pat.vcf@fix[,'POS'],'_',pat.vcf@fix[,'REF'],'/',pat.vcf@fix[,'ALT'])
  
  vafinfo.matrix <- pat.vcf@gt[match(pat.muts$mutID, vcf.mutID) ,regs]
  muts.matrix <- as.matrix(pat.muts[,regs])
  # get VAF matrix
  vaf.df <- as.data.frame(getVAFmatrix(vafinfo.matrix, muts.matrix))
  # correct for sample purity
  vaf.df <- t(t(vaf.df)/samp.df$Purity)
  vaf.df <- cbind(vaf.df, pat.muts[,!grepl(samp, names(pat.muts))])
  
  write.table(vaf.df, paste0('Subclonality/VAF_raw/',samp,'.total.VAFcorrected.txt'),
              sep='\t', row.names=F, quote=F)
}

# Combine all samples and patients together for Neo/Non-Neo (only protein changing mutations) 1/VAF values
vaf.df.total <- data.frame(matrix(vector(),ncol=8))
for (samp in samples){
  samp.df <- sampleList.df[sampleList.df$Patient==samp,]
  regs <- samp.df$Sample
  
  pat.muts <- read.delim(paste0('Subclonality/',samp,'.total.extended.txt'))
  pat.vaf <- read.delim(paste0('Subclonality/VAF_raw/',samp,'.total.VAFcorrected.txt'))
  # Carry out filtering for mutation type, etc.
  pat.vaf <- subset(pat.vaf, !is.na(Neoantigen))
  vaf.df.m <- reshape2::melt(pat.vaf[,c(regs,'Neoantigen', 'Category', 'MutType','mutID')], id=c('Neoantigen', 'Category','MutType','mutID') )
  vaf.df.m <- subset(vaf.df.m, value>0)
  vaf.df.m$invVAF <- 1/vaf.df.m$value
  vaf.df.m$Patient <- samp 
  
  vaf.df.total <- rbind(vaf.df.total, vaf.df.m)
}


# Evaluate (sub)clonality -------------------------------------------------

categoriseClonality <- function(x, samp.df){
  # input information: number of biopsies (samples) and regions (larger cancer areas) a mutation is present in
  category <- 'other'
  regs <- samp.df$Sample
  # mutations present in a single area of the tumour are regional. if only one area is sequenced, we cannot distinguish
  if (x$numRegions==1 & length(unique(samp.df$Region))>1){
    category <- 'regional'
  }
  # mutations present in a single biopsy are definitely regional, and private if we can prove it is a single sample from that area
  if (x$numSamples==1){
    withEpi <- regs[x[,regs]==1]
    areaCount <- sum(samp.df$Region==samp.df$Region[samp.df$Sample==withEpi]) # number of biopsies from that region in total
    # if we have >1 sample from the area but only 1 with antigen, it's categorised private (otherwise stays regional)
    if (areaCount>1){
      category <- 'private'
    }
  }
  # mutations present in at least two areas but not in all areas (with at least two samples with absence) are shared subclonal
  if (x$numRegions>1 & x$numRegions<length(unique(samp.df$Region)) & x$numSamples<(length(regs)-1)){
    category <- 'shared_subclonal'
  }
  # mutations present in only the adenoma are assigned a specific category of this
  if (x$numSamplesCancer==0){
    category <- 'adenoma_only'
    if (x$numSamples==1){
      category <- 'private_adenoma'
    }
  }
  # mutations present in all cancer samples (even if absent in adenoma) are noted to be truncal in cancer
  if (x$numSamplesCancer == sum(samp.df$Tissue=='Cancer')){
    category <- 'truncal'
  }
  # mutations present in all samples (even on top of cancer samples) are fully truncal
  if (x$numSamples == length(regs) & x$numSamples > x$numSamplesCancer){
    category <- 'truncal_full'
  }
  return(category)
}


# Annotate tables of all mutations with clonality and neoantigen information
for (samp in samples){
  samp.df <- sampleList.df[sampleList.df$Patient==samp,]
  
  regs <- epicc.df$Sample[epicc.df$Patient==samp]
  pat.muts <- read.delim(paste0('Subclonality/',samp,'.total.extended.txt'))
  pat.muts$numSamples <- rowSums(pat.muts[,intersect(regs, epicc.df$Sample)])
  pat.muts <- subset(pat.muts, numSamples>0)
  
  # compute number of regions mutation is present in
  for (area in unique(samp.df$Region)){
    regs.area <- samp.df$Sample[samp.df$Region==area]
    pat.muts[,paste0('Region_',area)] <- as.numeric(rowSums(pat.muts[,regs.area,drop=F])>0)
  }
  pat.muts$numRegions <- rowSums(pat.muts[,startsWith(names(pat.muts),'Region_'),drop=F])
  regs.canc <- samp.df$Sample[samp.df$Tissue=='Cancer']
  pat.muts$numSamplesCancer <- rowSums(pat.muts[,regs.canc])
  
  pat.muts$Category <- NA
  for (i in 1:nrow(pat.muts)){
    pat.muts$Category[i] <- categoriseClonality(pat.muts[i,], samp.df)
  }
  
  # check if mutation is neoantigen
  epTable <- read.table(paste0('Burden/Neopred_raw/',samp,'.neoantigens.extended.txt'),
                        sep='\t', stringsAsFactors = F, header=T)
  epTable.indel <- read.table(paste0('Burden/Neopred_raw/',samp,'.neoantigens.Indels.extended.txt'),
                              sep='\t', stringsAsFactors = F, header=T)
  epTable <- rbind(epTable, epTable.indel)
  # different levels: non-antigen, weak antigen, strong antigen, strong with specific A & R value, strong with RecoPo>=0.1
  # when strong neoantigens are evaluated, the last three categories are considered
  pat.muts$Neoantigen <- ifelse(pat.muts$CodingChange, 'Non-Neo', NA)
  pat.muts$Neoantigen[pat.muts$mutID %in% epTable$mutID] <- 'Neo'
  pat.muts$Neoantigen[pat.muts$mutID %in% epTable$mutID[epTable$BindLevel=='SB']] <- 'SB Neo'
  pat.muts$Neoantigen[pat.muts$mutID %in% epTable$mutID[epTable$BindLevel=='SB' & !is.na(epTable$RecoPo) & (epTable$RP_A > 7.5 | epTable$RP_R > 1e-6)]] <- 'SBRP Neo'
  pat.muts$Neoantigen[pat.muts$mutID %in% epTable$mutID[epTable$BindLevel=='SB' & epTable$RecoPo>=0.1 & !is.na(epTable$RecoPo)]] <- 'SBRecopo Neo'
  
  write.table(pat.muts, paste0('Subclonality/',samp,'.total.extended.txt'),
              sep='\t', row.names=F, quote=F)
}
