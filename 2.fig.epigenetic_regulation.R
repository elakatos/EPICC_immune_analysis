source('0.basics.R')

# Expression information
expTable.master <- read.delim('RNA/All_EPICC_tpm.txt')
vst.deseq <- readRDS('RNA/allgenes.vsd.ensembl_withC516adenoma.rds')
geneexp <- as.data.frame(assay(vst.deseq))

# SCAA of APGs ------------------------------------------------------------

# Plot heatmap of SCAAs in promoters of APGs
apg.scaa.df <- read.delim('Immune_escape/EPICC_APG_SCAA_master_file.txt')

# get a per-cancer mutation matrix - if any biopsy in a cancer/adenoma has a mutation, we denote it mutated
apgmut.master <- read.delim('Immune_escape/EPICC_APG_MUT_master_file.txt')
apgmut.simple <- data.frame(matrix(vector(),ncol=34))
names(apgmut.simple) <- c('Patient',names(apgmut.master)[2:34])
for (pat in unique(apgmut.master$Patient_AdCar)){
  apgmut.pat <- apgmut.master[apgmut.master$Patient_AdCar==pat,]
  apgmut.simple[nrow(apgmut.simple)+1,] <- c(pat,apply(apgmut.pat[,2:34],2,max))
}

# make a long-format dataframe and plot heatmap
apg.scaa.m <- reshape2::melt(apg.scaa.df, id='Sample')
apg.scaa.m$mut <- sapply(1:nrow(apg.scaa.m), function(i) ifelse(as.character(apg.scaa.m$variable[i]) %in% names(apgmut.simple),
                                                                apgmut.simple[apgmut.simple$Patient==apg.scaa.m$Sample[i],as.character(apg.scaa.m$variable[i])],
                                                                NA))
ggplot(apg.scaa.m, aes(x=Sample, y=variable)) +
  geom_tile(aes(fill=as.factor(sign(value)))) + theme_mypub() +
  scale_fill_manual(values=setNames(c('skyblue3','white','coral3'),c(-1,0,1)), labels=c('loss','','gain')) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), axis.title=element_blank()) +
  geom_point(aes(shape=as.factor(mut))) +
  scale_shape_manual(values=setNames(c(18,NA,NA),c(1,0,NA))) +
  labs(fill='SCAA') + guides(shape='none')


# Compare proportion of losses to randomised genes
apg.scaa.count <- table(sign(as.matrix(apg.scaa.df[,-1])))[c(1,3)]
rand.scaa.counts <- read.delim('ATAC/SCAA_randomised_genes.txt')
prop <- rand.scaa.counts$losses/(rand.scaa.counts$losses+rand.scaa.counts$gains)
table(prop>(apg.scaa.count[1]/(apg.scaa.count[1]+apg.scaa.count[2])))


# Plot normalised expression for gene-sample combinations with SCAA
apg <- read.delim('Immune_escape/APG_list.txt')
apg.cna.df <- read.delim('Immune_escape/EPICC_APG_CN_master_file.txt')
apg.scaa.df <- read.delim('Immune_escape/EPICC_APG_SCAA_master_file.txt')
#expTable.apg <- expTable.master[expTable.master$GeneID %in% apg$EnsID[apg$Use],]
expTable.apg <- geneexp[apg$EnsID[apg$Use],]

# get mean of all normal samples
normTable.m <- as.data.frame(matrix(vector(),ncol=4))
norm <- rna.df$Sample[rna.df$Tissue=='Normal']
apg.norm <- rowMeans(expTable.apg[,norm,drop=F])
#names(apg.norm) <- expTable.apg$GeneID
# for each patient get genes with SCAA and then compute normalised expression compared to all normals
for (pat in intersect(apg.scaa.df$Sample, rna.df$Patient)){
  scaa.genes <- names(apg.scaa.df)[apg.scaa.df[apg.scaa.df$Sample==pat,]!=0]
  scaa.geneid <- apg$EnsID[apg$GeneName %in% scaa.genes]
  if (length(scaa.geneid)==0){next}
  samps <- rna.df$Sample[rna.df$Patient==pat & rna.df$Tissue!='Normal']; samps <- intersect(samps, names(geneexp))
  #apg.tmp <- expTable.apg[expTable.apg$GeneID %in% scaa.geneid,samps]
  #apg.tmp <- apg.tmp/apg.norm[expTable.apg$GeneID[expTable.apg$GeneID %in% scaa.geneid]]
  #apg.tmp$GeneID <- expTable.apg$GeneID[expTable.apg$GeneID %in% scaa.geneid]
  apg.tmp <- expTable.apg[scaa.geneid, samps]
  apg.tmp <- apg.tmp/apg.norm[scaa.geneid]
  apg.tmp$GeneID <- scaa.geneid
  apg.tmp$GeneName <- apg$GeneName[match(apg.tmp$GeneID, apg$EnsID)]
  apg.tmp.m <- reshape2::melt(apg.tmp, id=c('GeneID','GeneName'))
  normTable.m <- rbind(normTable.m, apg.tmp.m)
}

# add CN information 
normTable.m$value[normTable.m$value==Inf] <- NA
normTable.m$CN <- sapply(1:nrow(normTable.m), function(i) ifelse(as.character(normTable.m$variable[i]) %in% apg.cna.df$Sample,
                                                                 apg.cna.df[apg.cna.df$Sample==as.character(normTable.m$variable[i]),normTable.m$GeneName[i]],
                                                                 NA))
normTable.m$Patient <- rna.df$Patient[match(normTable.m$variable,rna.df$Sample)]
normTable.m$variable <- sapply(as.character(normTable.m$variable), function(z) substr(z, 6,nchar(z)))
normTable.sub <- subset(normTable.m, !is.na(CN)) #exclude samples/genes with NA CN
# ggplot(normTable.sub, aes(x=variable, y=GeneName, fill=log2(value+1))) + 
#   geom_tile() + theme_mypub() +
#   scale_fill_gradientn(colours=c('#2166ac','#67a9cf','white','firebrick'), values=c(0,0.5/5, 1/5, 1),
#                        limits=c(0, 5), na.value='grey70', breaks=c(0,1,2,3,4), labels=c(0, 1, 2, 4,8)) +
#   theme(axis.text.x = element_text(angle=90, hjust=0.5, size=10), axis.text.y=element_text(size=10)) +
#   theme(strip.background = element_blank(), axis.title=element_blank()) +
#   # geom_point(aes(shape=as.factor(CN))) + #uncomment to include an additional marker to indicate CN loss/neutral/gain
#   scale_shape_manual(values=setNames(c(15,18,NA,NA,NA,NA,NA),c(1,2,3,4,5,6,NA))) +
#   facet_wrap(.~Patient, scales='free') + labs(fill='Normalised\nexpression')

ggplot(na.omit(normTable.sub), aes(x=variable, y=GeneName, fill=log2(value))) + 
  geom_tile() + theme_mypub() +
  scale_fill_gradientn(colours=c('#2166ac','#67a9cf','white','firebrick'), values=c(0,0.4,0.5, 1),
                       limits=c(log2(0.5), log2(2)), na.value='grey70') +
  theme(axis.text.x = element_text(angle=90, hjust=0.5, size=10), axis.text.y=element_text(size=10)) +
  theme(strip.background = element_blank(), axis.title=element_blank()) +
  # geom_point(aes(shape=as.factor(CN))) + #uncomment to include an additional marker to indicate CN loss/neutral/gain
  scale_shape_manual(values=setNames(c(15,18,NA,NA,NA,NA,NA),c(1,2,3,4,5,6,NA))) +
  facet_wrap(.~Patient, scales='free') + labs(fill='Normalised\nexpression\n(log2)')


# Plot expresion for SCAAloss vs wild-type tumour regions
gene <- 'ERAP2'
ensid <- apg$EnsID[apg$GeneName==gene]
gene.scaa.regions <- apg.scaa.df$Sample[apg.scaa.df[,gene]<0]
gene.cn.regions <- apg.cna.df$Sample[apg.cna.df[,gene]==2] #add info about neutral/not CN

gene.df <- epicc.df[epicc.df$MatchRNA=='Yes',c('Sample','Patient','Region','Type')]
gene.df$SCAA <- (gene.df$Patient %in% gene.scaa.regions)
gene.df$CN2 <- gene.df$Sample %in% gene.cn.regions
gene.df$Purity <- epicc.df$Purity[match(gene.df$Sample, epicc.df$Sample)]
#gene.df$Exp <- as.numeric(expTable.master[expTable.master$GeneID==ensid,gene.df$Sample])
gene.df <- subset(gene.df, Sample %in% names(geneexp))
gene.df$Exp <- as.numeric(geneexp[ensid,gene.df$Sample])

ggplot(gene.df, aes(x=paste0(SCAA), y=Exp)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_mypub() +
  geom_jitter(height=0, width=0.2, aes(colour=Type, shape=Type)) +
  stat_compare_means(comparisons=list(c('FALSE','TRUE'))) +
  scale_x_discrete(labels=c('wild-type','SCAA loss')) +
  labs(x=gene,y=paste0(gene,' expression (VST)')) +
  scale_colour_manual(values=c('darkblue','skyblue'),labels=c('deep WGS','low-pass WGS')) +
  scale_shape_discrete(labels=c('deep WGS','low-pass WGS')) +
  scale_y_log10(breaks=c(6, 10, 16)) #+ facet_wrap(.~CN2) #uncomment to test neutral/non-neutral CN separately



# Phylogenetic signal and ATAC peaks --------------------------------------

# Plot ATAC tracks for promoter of phylogenetic APG example
counts <- readRDS("ATAC/atacseq_counts_per_peak.rds")
getRegionATACCoverage <- function(samples_to_count, counts){
  # read in all samples and merge into one, as if a mega-ATAC file was observed and peaks accummulated
  gr <- unlist(GRangesList(lapply(1:length(samples_to_count), function(s) {
    gr <-import(paste0('ATAC/Examples/EPICC_',samples_to_count[s],'_C1_nucleosome_free-GRCh38-proper_pairs-shifted.bed.gz'),
                format="BED",
                genome="hg38",
                which=region)
    gr <- resize(gr, width = 100, fix = "center")
    return(gr)
  })))
  measured_samples<-paste0('EPICC_',samples_to_count,'_C1')
  # normalise to peak counts of sample to eliminate inter-sample variability
  samples_counts <- colSums(counts[,measured_samples,drop=F]) / 1000000 #divides all the reads of each sample by 1M 
  samples_cov <- coverage(gr)/sum(samples_counts)
  bins <- unlist(tile(region,width=1))
  samples_cov.overlap <- samples_cov[names(samples_cov) %in% levels(seqnames(bins))]
  samples_cov_1b <- GenomicRanges::binnedAverage(bins, samples_cov.overlap, "score")
  return(samples_cov_1b)
}

scaa.loc <- 'chr6:44246124-44247274' #HSP90AB1 promoter regions
region <- GRanges(seqnames = unlist(strsplit(scaa.loc,':'))[1],
                  ranges = unlist(strsplit(scaa.loc,':'))[2])
pat <- 'C559'
regsClade1 <- c('D','C')
regsClade2 <- c('B')
samples_to_count <- atac.df$Sample[atac.df$Patient==pat & atac.df$Region %in% regsClade1]
samples_cov_1b <- getRegionATACCoverage(samples_to_count, counts)
samples_cov.tmp <- data.frame(Position=samples_cov_1b@ranges@start, Score=samples_cov_1b$score, Sample='C&D')
samples_cov.df <- samples_cov.tmp
samples_to_count <- atac.df$Sample[atac.df$Patient==pat & atac.df$Region %in% regsClade2]
samples_cov_1b <- getRegionATACCoverage(samples_to_count, counts)
samples_cov.tmp <- data.frame(Position=samples_cov_1b@ranges@start, Score=samples_cov_1b$score, Sample='B')
samples_cov.df <- rbind(samples_cov.df, samples_cov.tmp)

ggplot(samples_cov.df, aes(x=Position, y=Score, colour=Sample)) + geom_line() +
  theme_mypub() +
  labs(title=paste0(scaa.loc,'\nHSP90AB1 promoter'), y='Peak coverage in C559', colour='Cancer\nregion') +
  scale_colour_manual(values=c('orange','brown'))


# Closed promoter and Neo co-occurrence -----------------------------------

# Compute and plot odds ratio for neoantigen being (somaticly) non-accessible
muttype <- 'SNV'
neoIdentifier <- c('SB Neo','SBRP Neo','SBRecopo Neo')
nonneoIdentifier <- c('Non-Neo')
atac.muts.total <- read.delim('ATAC/PromoterAccessSomaticForOR_NA_nonNA_allcancers.txt')
atac.muts.total <- subset(atac.muts.total, (numAccMut + numNotAccMut)>1 & MutType==muttype)

# go through all, clonal and subclonal mutations and compute a Fisher test
# the column numNotAccMut is the number of biopsies within a cancer with that mutation being somatically non-accessible
atac.muts.sub <- subset(atac.muts.total, (Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & !(Category %in% c('adenoma_only','private_adenoma')))
f.all <- fisher.test(atac.muts.sub$Neoantigen %in% neoIdentifier, atac.muts.sub$numNotAccMut>0)
atac.muts.sub <- subset(atac.muts.total, (Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & Category %in% c('truncal','truncal_full'))
f.cl <- fisher.test(atac.muts.sub$Neoantigen %in% neoIdentifier, atac.muts.sub$numNotAccMut>0)
atac.muts.sub <- subset(atac.muts.total, (Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & Category %in% c('private','regional','shared_subclonal','other'))
f.sc <- fisher.test(atac.muts.sub$Neoantigen %in% neoIdentifier,atac.muts.sub$numNotAccMut>0)

# summarise results and plot as rectangular markers with errorbars
f.df <- data.frame(Mutation=c('All','Clonal','Subclonal'), OR=c(f.all$estimate,f.cl$estimate, f.sc$estimate),
                   Conf_low=c(f.all$conf.int[1],f.cl$conf.int[1],f.sc$conf.int[1]), Conf_high=c(f.all$conf.int[2],f.cl$conf.int[2],f.sc$conf.int[2]),
                   P=c(f.all$p.value,f.cl$p.value, f.sc$p.value))
ggplot(f.df) +
  geom_rect(aes(ymin=OR-0.025, ymax = OR+0.025, xmin=as.numeric(as.factor(Mutation))-0.25, xmax=as.numeric(as.factor(Mutation))+0.25, fill=Mutation)) +
  scale_fill_manual(values=c('grey50',colorRampPalette(brewer.pal(8, "PuOr"))(8)[c(7,2)]),
                    labels=c('All','Clonal','Subclonal')) +
  scale_colour_manual(values=c('grey50',colorRampPalette(brewer.pal(8, "PuOr"))(8)[c(7,2)]),
                      labels=c('All','Clonal','Subclonal')) +
  geom_errorbar(aes(x=as.numeric(as.factor(Mutation)), ymin=Conf_low, ymax=Conf_high, colour=Mutation), width=0.1) +
  theme_mypub() + 
  geom_hline(yintercept = 1, linetype='dashed') +
  labs(y='Odds ratio') + guides(fill='none', colour='none') +
  geom_point(data=f.df[f.df$P<0.05,], aes(x=as.numeric(as.factor(Mutation)), y=3), shape=8) +
  scale_x_continuous(breaks=c(1,2,3), labels=f.df$Mutation) + theme(axis.title.x = element_blank())


# repeat the same for MMRp/d cancers only
atac.muts.sub <- subset(atac.muts.total, (Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & !(Category %in% c('adenoma_only','private_adenoma')))
chisq.test(paste0(atac.muts.sub$Neoantigen %in% neoIdentifier,atac.muts.sub$numNotAccMut>0), atac.muts.sub$Patient %in% msiList) # test MMRd/p proportions
atac.muts.sub <- subset(atac.muts.sub, (Patient %in% msiList)) #in msiList for MMRd, ! in msiList for MMRd
f.all.sub <- fisher.test(atac.muts.sub$Neoantigen %in% neoIdentifier, atac.muts.sub$numNotAccMut>0)
atac.muts.sub <- subset(atac.muts.total, (Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & Category %in% c('truncal','truncal_full'))
chisq.test(paste0(atac.muts.sub$Neoantigen %in% neoIdentifier,atac.muts.sub$numNotAccMut>0), atac.muts.sub$Patient %in% msiList) # test MMRd/p proportions
atac.muts.sub <- subset(atac.muts.sub, (Patient %in% msiList))
f.cl.sub <- fisher.test(atac.muts.sub$Neoantigen %in% neoIdentifier, atac.muts.sub$numNotAccMut>0)
atac.muts.sub <- subset(atac.muts.total, (Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & Category %in% c('private','regional','shared_subclonal','other'))
chisq.test(paste0(atac.muts.sub$Neoantigen %in% neoIdentifier,atac.muts.sub$numNotAccMut>0), atac.muts.sub$Patient %in% msiList) # test MMRd/p proportions
atac.muts.sub <- subset(atac.muts.sub, (Patient %in% msiList))
f.sc.sub <- fisher.test(atac.muts.sub$Neoantigen %in% neoIdentifier,atac.muts.sub$numNotAccMut>0)

f.df.sub <- data.frame(Mutation=c('All','Clonal','Subclonal'), OR=c(f.all.sub$estimate,f.cl.sub$estimate, f.sc.sub$estimate),
                   Conf_low=c(f.all.sub$conf.int[1],f.cl.sub$conf.int[1],f.sc.sub$conf.int[1]), Conf_high=c(f.all.sub$conf.int[2],f.cl.sub$conf.int[2],f.sc.sub$conf.int[2]),
                   P=c(f.all.sub$p.value,f.cl.sub$p.value, f.sc.sub$p.value))
f.df.sub



# Proportion of SCAAloss/Neo co-occurrence --------------------------------
muttype <- 'SNV'
neoIdentifier <- c('SB Neo','SBRP Neo','SBRecopo Neo')
nonneoIdentifier <- c('Non-Neo')

# Compute and plot the proportion of mutations that are affected by SCAAloss (for each cancer)
scaa.muts.total <- read.delim('ATAC/SCAALoss_NA_nonNA_allcancers.txt')
scaa.muts.total <- subset(scaa.muts.total, !is.na(SLossProm) & MutType==muttype & !(Category %in% c('adenoma_only','private_adenoma')) )

# denote SCAAloss with binary 1/0, filter for only SB Neo/Non-Neo and get proportion by average of 1/0
scaa.muts.total$SCAALoss <- 1*(scaa.muts.total$SLossProm>0)
atac.muts.sub <- subset(scaa.muts.total, (Neoantigen %in% c(neoIdentifier, nonneoIdentifier)))
notexp.pat <- aggregate(atac.muts.sub$SCAALoss, by=list(atac.muts.sub$Patient, atac.muts.sub$Neoantigen %in% neoIdentifier), mean)
notexp.pat$Escape <- patientEsc.df$Escape[match(notexp.pat$Group.1, patientEsc.df$Patient)]
notexp.pat$MSI <- ifelse(notexp.pat$Group.1 %in% msiList, 'MMRd','MMRp')
# exclude patients (if any) where there is only a proportion value for Neo or Non-Neo; plot with log-scale
pat.count <- table(notexp.pat$Group.1)
ggplot(notexp.pat[!(notexp.pat$Group.1 %in% (names(pat.count)[pat.count<2])),], aes(x=Group.2, y=log10(x+0.01)+2, colour=Escape)) +
  geom_point(size=2) + geom_line(aes(group=Group.1), alpha=0.7, size=1.2) +
  theme_mypub() + theme(axis.title.x=element_blank()) +
  scale_colour_manual(values=c('grey40','goldenrod','brown')) +
  scale_x_discrete(expand=c(0.1,0.1), labels=c('Non-Neo','SB Neo')) +
  scale_y_continuous(breaks=c(0, 0.544,1.0414), labels=c('0%','2.5%','10%')) +
  stat_compare_means(paired=T, comparisons=list(c('TRUE','FALSE'))) +
  labs(y='Proportion of mutations\nwith SCAAloss') #+ facet_wrap(.~MSI)

# Compute and plot the number of SCAAlosses that are preceding a neoantigen/non-neoantigen (for each cancer)
scaa.muts.total <- read.delim('ATAC/SCAALoss_mutsPerSCAA_allcancers.txt')
scaa.muts.total <- subset(scaa.muts.total, !is.na(Neoantigen) & !(Category %in% c('adenoma_only','private_adenoma')) )

# denote neoantigen SCAAloss with binary 1/0, filter out other mutations
scaa.muts.total$InNA <- 1*(scaa.muts.total$Neoantigen!='Non-Neo' )
atac.muts.sub <- subset(scaa.muts.total, (Neoantigen!=''))
loss.pat <- as.data.frame(table(atac.muts.sub$Patient, atac.muts.sub$InNA))
loss.pat$Escape <-patientEsc.df$Escape[match(loss.pat$Var1, patientEsc.df$Patient)]
loss.pat$MSI <- ifelse(loss.pat$Var1 %in% msiList, 'MMRd','MMRp')
ggplot(loss.pat, aes(x=Var2, y=log10(Freq+0.1)+1, colour=Escape)) +
  geom_point(size=2) + geom_line(aes(group=Var1), alpha=0.5, size=1.2) +
  theme_mypub() + theme(axis.title.x=element_blank()) +
  scale_colour_manual(values=c('grey40','goldenrod','brown')) +
  scale_x_discrete(expand=c(0.1,0.1), labels=c('Before Non-Neo','Before Neo')) +
  scale_y_continuous(breaks=c(0, 1.0414, 1.70757, 2.6031) , labels=c(0, 1, 5, 40))+
  stat_compare_means(paired=T, comparisons=list(c('0','1'))) +
  labs(y='Number of promoter SCAAloss') #+ facet_wrap(.~MSI)


