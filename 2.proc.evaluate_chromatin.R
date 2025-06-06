
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



# Transcription factor and APG-SCAAloss -----------------------------------

# Function to read in bed file of transcription factor binding sites into GRanges
read_bed_file = function(f) {
  require(GenomicRanges)
  
  stopifnot(file.exists(f))
  
  fl = readr::read_tsv(f, col_names=FALSE, col_types = readr::cols(), n_max=2)
  if (!suppressWarnings(is.na(as.numeric(fl[1,2])))) { # asssume no header
    if (NROW(fl) == 0) return(NULL)
    bed_data = readr::read_tsv(f, col_names = FALSE, col_types = "cii")
  } else { # assume header
    if (NROW(fl) == 1) return(NULL)
    stopifnot(!suppressWarnings(is.na(as.numeric(fl[2,2])))) 
    bed_data = readr::read_tsv(f, col_names = TRUE, col_types = "cii")
  }
  colnames(bed_data) = c("chr","start","end")
  bed_data = as(bed_data, "GRanges")
  return(bed_data)
}

# get Granges object for peaks associated with at least X SCAloss in APG
apg.genes <- read.delim('Immune_escape/APG_list.txt')
scaa <- readRDS('ATAC/analysis_results_raw.rds')
scaa.ind.list <- c()
scaa.df <- scaa$fc*scaa$sig #table only containing significant fold changes, 0 otherwise
names(scaa.df) <- gsub('.pure','',names(scaa.df))
for (gene in apg.genes$GeneName){
  scaa.ind <- which(sapply(scaa$flanking_symbol, function(z) gene %in% unlist(z)))
  #scaa.loss <- rowSums(scaa.df[scaa.ind,]< -1) > 1 # require loss to be present in >1 sample
  scaa.loss <- rowSums(scaa.df[scaa.ind,]< -1) > 0 # require loss to be present in any sample
  if (sum(scaa.loss)>0){ # only add if we detected at least one SCAAloss associated with that peak
    scaa.ind.list <- c(scaa.ind.list, scaa.ind[scaa.loss])
  }
}
peak_pos_scaloss <- scaa$peak_annotation@anno[scaa.ind.list]

# get TF overlap with peaks
bed_folder <- 'ATAC/TFs/motif_matches_v3.1_beds/'
bed_files <- list.files(bed_folder, pattern='*.bed.gz')
tf.scaloss <- data.frame(matrix(vector(),ncol=6))

for (i in 1:length(bed_files)){
  print(i)
  tf_bind <- read_bed_file(paste0(bed_folder,bed_files[i]))
  tf_dist <- as.data.frame(distanceToNearest(tf_bind-1000, peak_pos_scaloss-250))
  tf_overlap <- tf_dist[tf_dist$distance<1000,] # if distance is less than 1000, TF binding site overlaps peak
  tf_overlap$name <- peak_pos_scaloss$SYMBOL[tf_overlap$subjectHits]
  tf_overlap$annotation <- peak_pos_scaloss$annotation[tf_overlap$subjectHits]
  if (nrow(tf_overlap)>0){
    tf_overlap$tf_name <- gsub('.bed.gz','',bed_files[i])
    tf.scaloss <- rbind(tf.scaloss, tf_overlap)
  }
}
# shorten name of TF
tf.scaloss$tf_name <- sapply(tf.scaloss$tf_name, function(z) unlist(strsplit(z,'_'))[[3]])
# saved as TF_peakoverlap_APG_recurscaloss.txt or TF_peakoverlap_APG_scaloss.txt


# TF binding in general vs APG-SCAAloss loci ------------------------------

# All predicted TF binding sites vs APG-SCAAloss binding, higlighting TFs we identified
bed_folder <- 'ATAC/TFs/motif_matches_v3.1_beds/'
bed_files <- list.files(bed_folder, pattern='*.bed.gz')

binding_site_count <- data.frame(matrix(vector(), ncol=2))
names(binding_site_count) <- c('count','tf_name')
for (i in 1:length(bed_files)){
  print(i)
  tf_bind <- read_bed_file(paste0(bed_folder,bed_files[i]))
  binding_site_count[nrow(binding_site_count)+1,] <- c(length(tf_bind), gsub('.bed.gz','',bed_files[i]))
}
binding_site_count$tf_name <- sapply(binding_site_count$tf_name, function(z) unlist(strsplit(z,'_'))[[3]])
binding_site_count$count <- as.numeric(binding_site_count$count)

# Annotate with count of binding sites overlapping (recurr) APG-SCAAloss peaks
tf.scaloss <- read.delim('ATAC/TFs/TF_peakoverlap_APG_recurscaloss.txt')
tf.scaloss.sum <- as.data.frame(table(tf.scaloss$name, tf.scaloss$tf_name))
tf.scaloss.pertf <- aggregate(tf.scaloss.sum$Freq, by=list(tf.scaloss.sum$Var2), sum); names(tf.scaloss.pertf) <- c('tf_name','APG_count')
tf.scaloss.pertf <- tf.scaloss.pertf[order(tf.scaloss.pertf$APG_count, decreasing = T),]

binding_site_count$apgscaa_bind <- tf.scaloss.pertf$APG_count[match(binding_site_count$tf_name, tf.scaloss.pertf$tf_name)]
binding_site_count$apgscaa_bind[is.na(binding_site_count$apgscaa_bind)] <- 0

# Annotate with number of unique APGs affected by these sites
tf.scaloss.sum$Freq <- 1*(tf.scaloss.sum$Freq>0) # summarise for each TF
tf.scaloss.pertf <- aggregate(tf.scaloss.sum$Freq, by=list(tf.scaloss.sum$Var2), sum); names(tf.scaloss.pertf) <- c('tf_name','APG_count')
binding_site_count$apgs_affected <- tf.scaloss.pertf$APG_count[match(binding_site_count$tf_name, tf.scaloss.pertf$tf_name)]
binding_site_count$apgs_affected[is.na(binding_site_count$apgs_affected)] <- 0
tf.selected_recur <- tf.scaloss.pertf[tf.scaloss.pertf$APG_count>2,] #take ones binding >2 different APGs' promoters
binding_site_count$in_recurr <- binding_site_count$tf_name %in% tf.selected_recur$tf_name
ggplot(binding_site_count, aes(x=count, y=apgscaa_bind, colour=in_recurr, alpha=in_recurr)) +
  geom_point(size=2) + scale_x_log10() +
  scale_colour_manual(values=setNames(c('grey30','darkgreen'),c('FALSE','TRUE'))) +
  scale_alpha_manual(values=c(0.25, 1)) +
  theme_mypub() +
  guides(colour='none', alpha='none') + labs(x='# of binding sites', y='# of sites in APG SCAAloss prom.') +
  geom_text(data = subset(binding_site_count, in_recurr ),
            aes(label = tf_name), color = "darkgreen", nudge_x = -0.25, nudge_y = 1, size=3)

# Test the two groups (in_recurr and not) and their APG-SCAAloss binding sites vs other sites
apg_bind.sum <- aggregate(binding_site_count$apgscaa_bind, by=list(binding_site_count$in_recurr), sum)
names(apg_bind.sum) <- c('APG_SCAA_TF', 'apgscaa_bind')
apg_bind.sum$bind <- aggregate(binding_site_count$count, by=list(binding_site_count$in_recurr), sum)$x
apg_bind.sum$bind <- apg_bind.sum$bind - apg_bind.sum$apgscaa_bind
apg_bind.sum
fisher.test(apg_bind.sum[,2:3], alternative = 'less')

# Annotate with binding site count in SCAAloss promoters (any gene)
scaa <- readRDS('ATAC/analysis_results_raw.rds')
scaa.df <- scaa$fc*scaa$sig #table only containing significant fold changes, 0 otherwise
names(scaa.df) <- gsub('.pure','',names(scaa.df))
scaa.loss.ind <- which(rowSums(scaa.df< -1) > 0) # all SCAA losses regardless of location
peak_pos_lossprom <- subset(scaa$peak_annotation@anno[scaa.loss.ind], grepl('Promoter',annotation))

lossprom_binding_site_count <- data.frame(matrix(vector(), ncol=2))
names(lossprom_binding_site_count) <- c('count','tf_name')
for (i in 1:length(bed_files)){
  print(i)
  tf_bind <- read_bed_file(paste0(bed_folder,bed_files[i]))
  tf_dist <- as.data.frame(distanceToNearest(tf_bind-1000, peak_pos_lossprom-250))
  tf_overlap <- tf_dist[tf_dist$distance<1000,] # if distance is less than 1000, TF binding site overlaps peak
  lossprom_binding_site_count[nrow(lossprom_binding_site_count)+1,] <- c(nrow(tf_overlap), gsub('.bed.gz','',bed_files[i]))
}
lossprom_binding_site_count$tf_name <- sapply(lossprom_binding_site_count$tf_name, function(z) unlist(strsplit(z,'_'))[[3]])
lossprom_binding_site_count$count <- as.numeric(lossprom_binding_site_count$count)

binding_site_count$promscaa_bind <- lossprom_binding_site_count$count[match(binding_site_count$tf_name, lossprom_binding_site_count$tf_name)]
ggplot(binding_site_count, aes(x=count, y=promscaa_bind, colour=in_recurr, alpha=in_recurr)) +
  geom_point(size=2) + scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values=setNames(c('grey30','darkgreen'),c('FALSE','TRUE'))) +
  scale_alpha_manual(values=c(0.25, 1)) +
  theme_mypub() +
  guides(colour='none', alpha='none') + labs(x='# of binding sites', y='# of sites in SCAAloss prom.') +
  geom_text(data = subset(binding_site_count, in_recurr ),
            aes(label = tf_name), color = "darkgreen", nudge_x = -0.2, nudge_y = 0.2, size=3)

# Test the two groups (in_recurr and not) and their APG-SCAAloss binding sites vs other promoter loss sites sites
apg_bind.sum <- aggregate(binding_site_count$apgscaa_bind, by=list(binding_site_count$in_recurr), sum)
names(apg_bind.sum) <- c('APG_SCAA_TF', 'apgscaa_bind')
apg_bind.sum$bind <- aggregate(binding_site_count$promscaa_bind, by=list(binding_site_count$in_recurr), sum)$x
apg_bind.sum$bind <- apg_bind.sum$bind - apg_bind.sum$apgscaa_bind
apg_bind.sum
fisher.test(apg_bind.sum[,2:3], alternative = 'less')

# Highlight those with enough binding sites overlapping with APG-SCAAloss peaks
bind_tfs <- subset(binding_site_count, apgscaa_bind>=3)
bind_tfs$apgs_affected_names <- sapply(bind_tfs$tf_name, function(tf)  paste0(tf.scaloss.sum[tf.scaloss.sum$Var2==tf & tf.scaloss.sum$Freq>0,'Var1'], collapse=';'))
