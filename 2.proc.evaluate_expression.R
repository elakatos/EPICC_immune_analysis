
source('0.basics.R')



# Compute phylogenetic signal ---------------------------------------------

# Compute the phylogenetic signal of APGs/ICPs by comparing expression and evolutionary trees using Pagel's lambda
# we use a set of 100 pruned trees where only biopsies with shared WGS&RNA are kept and branch lengths for lpWGS are bootstrapped
# for more details, see Househam & Heide, Nature, 2022.

treelist <- readRDS('RNA/pgls.trees.rds')
row.names(expTable.master) <- expTable.master$GeneID

# we compute T-cell score based on the combined expression of these 12 genes
tcag <- c('ENSG00000108691', 'ENSG00000277632', 'ENSG00000275302', 'ENSG00000138755', 'ENSG00000169245',
          'ENSG00000153563', 'ENSG00000241106', 'ENSG00000242574', 'ENSG00000204252', 'ENSG00000113088',
          'ENSG00000163600', 'ENSG00000125347')
exprtcag <- sapply(expTable.master[(expTable.master$GeneID %in% tcag),2:ncol(expTable.master)], function(x) log10(x+1))
expTable.master['Tcell_score',] <- c('Tcell_score',apply(exprtcag, 2, mean))
expTable.master[,2:ncol(expTable.master)] <- apply(expTable.master[,2:ncol(expTable.master)],2,as.numeric)

# pd-l1, ctla-4, tcellscore, apgs, hla-a, -b, -c
apg.genes <- read.delim('Immune_escape/APG_list.txt')$GeneName
genes.tocheck.df <- gene.table[(gene.table$HGNC.symbol %in% c(apg.genes,'PD-L1','CTLA-4')) & (gene.table$Chromosome.scaffold.name %in% c(1:22,'X','Y')), ]
genes.tocheck.df <- genes.tocheck.df[!duplicated(genes.tocheck.df$Gene.stable.ID),]
genes.tocheck <- c('Tcell_score',genes.tocheck.df$Gene.stable.ID)

gene_lambda_list <- list()
for(pat in names(treelist)) {
  print(pat)
  trimtree <- treelist[[pat]][[1]]
  treesam <- paste0(pat,'_',trimtree$tip.label)
  rna.measure <- as.data.frame(t(expTable.master[genes.tocheck,treesam])) # take expression from full expression table
  row.names(rna.measure) <- gsub('C\\d+_(\\S+)','\\1',row.names(rna.measure))
  
  siglist <- list()
  for(i in c(1:length(treelist[[pat]]))) {
    trimtree <- treelist[[pat]][[as.character(i)]]
    
    sigphylo <- data.frame(Gene=colnames(rna.measure),Lambda=0,Lpval=0,Kstat=0,Kpval=0)
    for(j in c(1:ncol(rna.measure))) {
      curgen <- rna.measure[,j];names(curgen) <- row.names(rna.measure)
      if(var(curgen)>0) {
        res <- phylosig(trimtree,curgen,method='lambda',test=T)
        sigphylo[j,'Lambda'] <- res$lambda;sigphylo[j,'Lpval'] <- res$P
      } else {
        res <- phylosig(trimtree,curgen,method='lambda',test=F)
        sigphylo[j,'Lambda'] <- res$lambda;sigphylo[j,'Lpval'] <- 1
      }
    }
    siglist[[as.character(i)]] <- sigphylo
    print(paste0(pat,': ',signif((i/length(treelist[[pat]])*100),2),'%'))
  }
  
  # get all lambda and p- values and compute medians
  lam <- do.call(cbind, lapply(siglist,function(x) { x$Lambda }))
  lam_pval <- do.call(cbind, lapply(siglist,function(x) { x$Lpval }))
  resdf <- data.frame(MedLambda=rowMedians(lam),
                      LamPval=rowMedians(lam_pval))
  row.names(resdf) <- colnames(rna.measure)
  
  gene_lambda_list[[pat]] <- resdf
}
# transform into data frame
exp.lambda.df <- data.frame(matrix(vector(),ncol=4))
for (i in 1:length(gene_lambda_list)){
  lambda.tmp <- gene_lambda_list[[i]]
  lambda.tmp$Patient <- names(gene_lambda_list)[i]
  lambda.tmp$Measure <- row.names(lambda.tmp)
  exp.lambda.df <- rbind(exp.lambda.df, lambda.tmp)
}

saveRDS(exp.lambda.df, file='RNA/PhyloLambda_rna_measures.allRNA.rds')


# Gather all prot-changing mutations with gene info -----------------------

mutTable.total <- matrix(vector(),ncol=5)
for (samp in samples){
  mutTable <- read.delim(paste0('Subclonality/',samp,'.total.extended.txt'))
  mutTable <- subset(mutTable, !is.na(mutTable$Neoantigen))
  mutTable$Sample <- samp
  mutTable.total <- rbind(mutTable.total, mutTable[,c('Sample','GeneID','Category','MutType','Neoantigen')])
}

# Process BRC files -------------------------------------------------------

# list all bam-readcount output files
rnaList <- sub('.BRC.txt','',list.files('RNA/BamReadCount_raw/', pattern='.BRC.txt'))

for (region in rnaList){
  maxCol <- max(count.fields(paste0('RNA/BamReadCount_raw/',region,'.BRC.txt'), sep = "\t"))
  brc.df <- data.frame(matrix(vector(),ncol=8))
  names(brc.df)=c('contig','position','refAllele','totalCount','A','C','G','T')
  # the file might be empty in which case we have -Inf columns
  if(maxCol!=-Inf){
    altCol <- maxCol-10
    # read in A, C, T, G info from bam-readcount file (ignore FS)
    brc.tmp <- read.table(paste0('RNA/BamReadCount_raw/',region,'.BRC.txt'),
                          col.names=c('contig','position','refAllele','totalCount','X','A','C','G','T','N',
                                      paste0('altAllele',1:(altCol))),
                          colClasses=c('character','numeric',rep('character',maxCol-2)),
                          fill=T, sep='\t', stringsAsFactors = F)
    brc.tmp <- subset(brc.tmp, totalCount>1)
    if (nrow(brc.tmp)>0){
      for (base in c('A','C','G','T')){
        # number of reads with a given base is the second : separated field
        brc.tmp[,base] <- sapply(brc.tmp[,base], function(z) as.numeric(unlist(strsplit(z,':'))[2]))
      }
      brc.df <- rbind(brc.df, brc.tmp[,c(1:4,6:9)])
    }
  }
  write.table(brc.df, paste0('RNA/BamReadCount_raw/',region,'.BRC.proc.txt'),
              sep='\t', row.names=F, quote=F)
}




# Process RNA muts into expr/not matrix -----------------------------------

samplesRNA <- unique(epicc.df$Patient[epicc.df$MatchRNA=='Yes' & epicc.df$Type=='Deep'])
rna.bam.df <- read.delim('RNA/EPICC_RNA_bamList.filtered.txt', header=F)

for (pat in samplesRNA){
  total.muts.pat <- read.delim(paste0('Subclonality/',pat,'.total.extended.txt'))
  total.muts.pat <- subset(total.muts.pat, CodingChange & MutType=='SNV')
  row.names(total.muts.pat) <- total.muts.pat$mutID
  rna.muts.pat <- total.muts.pat[,colnames(total.muts.pat) %in% rna.bam.df$V4, drop=F]
  
  # for each sample with both WGS and RNAseq, we read in bam-readcount evidence
  # note that for a handful of samples there were two batches of RNA processed which should be merged together
  for (samp in rna.bam.df$V4[rna.bam.df$V2==pat]){
    samp.bc <- rna.bam.df$V3[rna.bam.df$V4==samp]
    j <- samp.bc[1]
    brc.df <- read.delim(paste0('RNA/BamReadCount_raw/',j,'.BRC.proc.txt'))
    if (nrow(brc.df)==0){
      rna.muts.pat[,samp] <- NA
      next
    }
    brc.df$Loc <- paste0(brc.df$contig,':',brc.df$position)
    # get mutation coordinate and the alternate allele we want to count the evidence for
    muts.loc <- sapply(row.names(total.muts.pat), function(i) unlist(strsplit(i,'_'))[1])
    muts.alt <- sapply(row.names(total.muts.pat), function(i) unlist(strsplit(i,'/'))[2])
    # find total reads and alt reads mapping to that location; if there's no info in BRC file, leave as NA
    rna.totalcount <- sapply(1:nrow(total.muts.pat), function(i) ifelse(muts.loc[i] %in% brc.df$Loc,
                                                                        brc.df[brc.df$Loc==muts.loc[i],'totalCount'],
                                                                        NA))
    rna.altcount <- sapply(1:nrow(total.muts.pat), function(i) ifelse((muts.loc[i] %in% brc.df$Loc) & (nchar(muts.alt[i])==1), #ignore FS/multi-allele
                                                                      brc.df[brc.df$Loc==muts.loc[i],muts.alt[i]],
                                                                      NA))
    #in cases with 2 RNAseq files, we add total and alt count from the two files
    if (length(samp.bc)>1){
      j <- samp.bc[2]
      brc.df <- read.delim(paste0('RNA/BamReadCount_raw/',j,'.BRC.proc.txt'))
      brc.df$Loc <- paste0(brc.df$contig,':',brc.df$position)
      muts.loc <- sapply(row.names(total.muts.pat), function(i) unlist(strsplit(i,'_'))[1])
      muts.alt <- sapply(row.names(total.muts.pat), function(i) unlist(strsplit(i,'/'))[2])
      
      rna.totalcount.tmp <- sapply(1:nrow(total.muts.pat), function(i) ifelse(muts.loc[i] %in% brc.df$Loc,
                                                                              brc.df[brc.df$Loc==muts.loc[i],'totalCount'],
                                                                              NA))
      rna.altcount.tmp <- sapply(1:nrow(total.muts.pat), function(i) ifelse(muts.loc[i] %in% brc.df$Loc,
                                                                            brc.df[brc.df$Loc==muts.loc[i],muts.alt[i]],
                                                                            NA))
      # note the manual handling of NA: if both are NA, we keep NA, otherwise the sum with NA=0
      rna.totalcount <- sapply(1:length(rna.totalcount), function(i) ifelse(is.na(rna.totalcount[i]),
                                                                            ifelse(is.na(rna.totalcount.tmp[i]), NA, rna.totalcount.tmp[i]),
                                                                            ifelse(is.na(rna.totalcount.tmp[i]), rna.totalcount[i], rna.totalcount.tmp[i] + rna.totalcount[i])))
      rna.altcount <- sapply(1:length(rna.altcount), function(i) ifelse(is.na(rna.altcount[i]),
                                                                        ifelse(is.na(rna.altcount.tmp[i]), NA, rna.altcount.tmp[i]),
                                                                        ifelse(is.na(rna.altcount.tmp[i]), rna.altcount[i], rna.altcount.tmp[i] + rna.altcount[i])))
      
    }
    # evaluate the mutation based on total and alt counts
    # we require at least 3 mutated reads to detect the mutation
    # and we require over 10 total reads with 0 alt reads at the position to confidently "not detect" the mutation
    rna.muts.pat[,samp] <- ifelse(rna.altcount>=3 , 1, ifelse(rna.altcount==0 & rna.totalcount>10, 0, NA))
  }
  write.table(rna.muts.pat, file=paste0('RNA/BamReadCount_raw/BRC_mutation_matrix_',pat,'.txt'), sep='\t', quote=F)
  
}

# Evaluate cancer-wise (number of biopsies with/without mut) information for each mutation
exp.muts.total <- data.frame(matrix(vector(),ncol=7))
for (pat in samplesRNA){
  rna.muts.pat <- read.delim(paste0('RNA/BamReadCount_raw/BRC_mutation_matrix_',pat,'.txt'))
  regs <- intersect(names(rna.muts.pat), epicc.df$Sample[epicc.df$Tissue=='Cancer'])
  muts.pat.extended <- read.delim(paste0('Subclonality/',pat,'.total.extended.txt'))
  row.names(muts.pat.extended) <- muts.pat.extended$mutID
  muts.pat <- muts.pat.extended[row.names(rna.muts.pat),regs,drop=F]
  
  exp.muts.pat <- data.frame(matrix(vector(),ncol=4))
  # columns: cancer/patient, mutated in WGS, mutation expressed in RNA, mutation not expressed in RNA
  # note that we only consider biopsies with WGS&RNA, and #expressed + #notexpressed might not equal mut (due to NAs)
  names(exp.muts.pat) <- c('Patient','numMut','numExprMut','numNotExprMut')
  for (mut in row.names(rna.muts.pat)){
    samp.mut <- sum(muts.pat[mut,]==1) #how many samples have it mutated in WGS
    mut.rna <- sum(rna.muts.pat[mut, samp.mut, drop=F]==1, na.rm=T) #how many have it mutated in RNA
    wt.rna <- sum(rna.muts.pat[mut, samp.mut, drop=F]==0, na.rm=T) #how many have it non-mutated in the RNA
    exp.muts.pat[nrow(exp.muts.pat)+1,] <- c(pat, samp.mut,mut.rna, wt.rna)
  }
  exp.muts.pat$mutID <- row.names(rna.muts.pat)
  
  exp.muts.pat[,2:4] <- apply(exp.muts.pat[,2:4], 2, as.numeric)
  # add annotation regarding (sub)clonality category and Neo
  exp.muts.pat[,c('Category','Neoantigen')] <- muts.pat.extended[match(exp.muts.pat$mutID, muts.pat.extended$mutID),c('Category','Neoantigen')]
  exp.muts.total <- rbind(exp.muts.total, exp.muts.pat)
}

write.table(exp.muts.total, file='RNA/Expression_NA_nonNA_allcancers.txt', sep='\t',row.names=F, quote=F)
