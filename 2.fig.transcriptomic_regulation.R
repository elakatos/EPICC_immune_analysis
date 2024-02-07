
source('0.basics.R')

# RNA information
apg <- read.delim('Immune_escape/APG_list.txt')
geneGroup.df <- readRDS('RNA/gene_clustering_and_id_conversion.rds')
expTable.master <- read.delim('RNA/All_EPICC_tpm.txt')

# Phylogenetic signal of APGs ---------------------------------------------

# Plot phylogenetic signal for APGs
exp.lambda.df <- readRDS('RNA/PhyloLambda_rna_measures.allRNA.rds')
exp.lambda.df$Name <- apg$GeneName[match(exp.lambda.df$Measure, apg$EnsID)]
# omit T-cell score as it is not an APG but a composite of T-cell associated genes
exp.lambda.df <- exp.lambda.df[exp.lambda.df$Measure!='Tcell_score',]

ggplot(exp.lambda.df, aes(x=Patient, y=Name, fill=-log10(LamPval))) +
  geom_tile() +
  theme_mypub() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), axis.title=element_blank()) +
  scale_fill_gradientn(colours=c('#ffffe6','#fed976', '#e31a1c'),
                       na.value='grey70') +
  geom_point(data = subset(exp.lambda.df, LamPval<0.05), aes(x = Patient, y = Name), 
             shape = 8) +
  scale_x_discrete(expand=c(0,0)) +
  labs(fill='-log10(p)')


# Gene group and Neo co-occurrence ----------------------------------------

muttype <- 'SNV' # choose between SNV/FS
neoIdentifier <- c('SBRecopo Neo', 'SBRP Neo','SB Neo')
nonneoIdentifier <- c('Non-Neo')
mutTable.total <- read.delim('Subclonality/NA_nonNA_allcancers.txt')
mutTable.total$GeneGroup <- geneGroup.df$Group[match(mutTable.total$GeneID, geneGroup.df$ensembl_gene_id)]
mutTable.total <- mutTable.total[!is.na(mutTable.total$GeneGroup),]

# Plot neoantigen/not for clonal mutation of each gene group 
mutTable.total.sub <- mutTable.total[mutTable.total$MutType==muttype & (mutTable.total$Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & mutTable.total$Category%in% c('truncal','truncal_full'),]
mutTable.total.sub$Neoantigen[mutTable.total.sub$Neoantigen %in% neoIdentifier] <- 'SB Neo'; mutTable.total.sub$Neoantigen[mutTable.total.sub$Neoantigen %in% nonneoIdentifier] <- 'Non-Neo'
ct <- chisq.test(mutTable.total.sub$GeneGroup, mutTable.total.sub$Neoantigen)
ggplot(mutTable.total.sub, aes(x=GeneGroup, fill=Neoantigen)) +
  geom_bar(position='fill', width=0.75) +
  #scale_fill_manual(values=setNames(c("#b16710",'grey80'),c('SB Neo','Non-Neo')), labels=c('Not neoant','Neoant')) +
  scale_fill_manual(values=setNames(c("#5a4990",'grey80'),c('SB Neo','Non-Neo')), labels=c('Not neoant','Neoant')) +
  scale_y_continuous(labels=percent_format(), expand=c(0,0), limits=c(0, 1.1)) +
  theme_mypub() + labs(x='Gene group', y=paste0('Clonal ',muttype,'s'), fill='') +
  geom_signif(aes(y=1), y_position=c(1.03), xmin=c(1), xmax=c(4),
              annotation=scientific(ct$p.value,digits=3), tip_length=0, size=0.5)

# Plot neoantigen/not for subclonal mutation of each gene group 
mutTable.total.sub <- mutTable.total[mutTable.total$MutType==muttype & (mutTable.total$Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & mutTable.total$Category%in% c( 'private','regional','shared_subclonal','other'),]
mutTable.total.sub$Neoantigen[mutTable.total.sub$Neoantigen %in% neoIdentifier] <- 'SB Neo'; mutTable.total.sub$Neoantigen[mutTable.total.sub$Neoantigen %in% nonneoIdentifier] <- 'Non-Neo'
ct <- chisq.test(mutTable.total.sub$GeneGroup, mutTable.total.sub$Neoantigen)
ggplot(mutTable.total.sub, aes(x=GeneGroup, fill=Neoantigen)) +
  geom_bar(position='fill', width=0.75) +
  scale_fill_manual(values=setNames(c("#b16710",'grey80'),c('SB Neo','Non-Neo')), labels=c('Not neoant','Neoant')) +
  scale_y_continuous(labels=percent_format(), expand=c(0,0), limits=c(0, 1.1)) +
  theme_mypub() + labs(x='Gene group', y=paste0('Subclonal ',muttype,'s'), fill='') +
  geom_signif(aes(y=1), y_position=c(1.03), xmin=c(1), xmax=c(4),
              annotation=scientific(ct$p.value,digits=3), tip_length=0, size=0.5)


# Compute and plot odds ratio for neoantigen being located in each gene group
f.df.total <- data.frame(matrix(vector(),ncol=6))
for (g in 1:4){
  group.genes <- geneGroup.df$ensembl_gene_id[geneGroup.df$Group==g]
  mutTable.total.sub <- mutTable.total[mutTable.total$MutType==muttype & (mutTable.total$Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & !(mutTable.total$Category %in% c('adenoma_only','private_adenoma','truncal_full')),]
  f.all <- fisher.test(mutTable.total.sub$GeneID %in% group.genes, mutTable.total.sub$Neoantigen %in% neoIdentifier)
  mutTable.total.sub <- mutTable.total[mutTable.total$MutType==muttype & (mutTable.total$Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & mutTable.total$Category%in% c('truncal'),]
  f.cl <- fisher.test(mutTable.total.sub$GeneID %in% group.genes, mutTable.total.sub$Neoantigen %in% neoIdentifier)
  mutTable.total.sub <- mutTable.total[mutTable.total$MutType==muttype & (mutTable.total$Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & mutTable.total$Category %in% c('private','regional','shared_subclonal','other'),]
  f.sc <- fisher.test(mutTable.total.sub$GeneID %in% group.genes, mutTable.total.sub$Neoantigen %in% neoIdentifier)
  f.df <- data.frame(Mutation=c('All','Clonal','Subclonal'), OR=c(f.all$estimate,f.cl$estimate, f.sc$estimate),
                     Conf_low=c(f.all$conf.int[1],f.cl$conf.int[1],f.sc$conf.int[1]), Conf_high=c(f.all$conf.int[2],f.cl$conf.int[2],f.sc$conf.int[2]),
                     P=c(f.all$p.value,f.cl$p.value, f.sc$p.value),
                     Group=g)
  f.df.total <- rbind(f.df.total, f.df)
}
f.df.total$Xloc <- as.numeric(as.factor(f.df.total$Group))
f.df.total$Xloc <- ifelse(f.df.total$Mutation=='All', f.df.total$Xloc-0.27,
                          ifelse(f.df.total$Mutation=='Subclonal', f.df.total$Xloc+0.27,f.df.total$Xloc  ))
# plot result as rectangular indicators with error bars
ggplot(f.df.total) +
  geom_rect(aes(ymin=OR-0.015, ymax = OR+0.015, xmin=Xloc-0.1, xmax=Xloc+0.1, fill=Mutation)) +
  scale_fill_manual(values=c('grey50',colorRampPalette(brewer.pal(8, "PuOr"))(8)[c(7,2)]),
                    labels=c('All','Clonal','Subclonal')) +
  scale_colour_manual(values=c('grey50',colorRampPalette(brewer.pal(8, "PuOr"))(8)[c(7,2)]),
                      labels=c('All','Clonal','Subclonal')) +
  geom_errorbar(aes(x=Xloc, ymin=Conf_low, ymax=Conf_high, colour=Mutation), width=0.1) +
  theme_mypub() + 
  geom_hline(yintercept = 1, linetype='dashed') +
  labs(y='Odds ratio', x='Gene group', fill=paste0(muttype,'s'), colour=paste0(muttype,'s')) +
  geom_point(data=f.df.total[f.df.total$P<0.05,], aes(x=Xloc, y=1.6), shape=8)


# Consistently expressed gene and Neo co-occurrence -----------------------

# Compute and plot odds ratio for neoantigen being located in a consistently expressed gene

# establish the list of consistently expressed genes, use the criteria for lung: >1TPM for 95% of cases
expTable.master.sub <- expTable.master[,names(expTable.master) %in% c(rna.df[rna.df$Region!='E',]$Sample, 'GeneID')]
const.exp <- expTable.master.sub$GeneID[rowSums(expTable.master.sub[,-1]>1)>276]

muttype <- 'SNV' # choose from SNV/FS
neoIdentifier <- c('SBRecopo Neo','SBRP Neo','SB Neo')
nonneoIdentifier <- c('Non-Neo')
mutTable.total <- read.delim('Subclonality/NA_nonNA_allcancers.txt')

mutTable.total.sub <- mutTable.total[mutTable.total$MutType==muttype & (mutTable.total$Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & !(mutTable.total$Category) %in% c('adenoma_only','private_adenoma','truncal_full'),]
f.all <- fisher.test(mutTable.total.sub$GeneID %in% const.exp, mutTable.total.sub$Neoantigen %in% neoIdentifier)
mutTable.total.sub <- mutTable.total[mutTable.total$MutType==muttype & (mutTable.total$Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & mutTable.total$Category%in% c('truncal'),]
f.cl <- fisher.test(mutTable.total.sub$GeneID %in% const.exp, mutTable.total.sub$Neoantigen %in% neoIdentifier)
mutTable.total.sub <- mutTable.total[mutTable.total$MutType==muttype & (mutTable.total$Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & mutTable.total$Category %in% c('private','regional', 'other', 'shared_subclonal'),]
f.sc <- fisher.test(mutTable.total.sub$GeneID %in% const.exp, mutTable.total.sub$Neoantigen %in% neoIdentifier)
f.df <- data.frame(Mutation=paste0(c('All ','Clonal ','Subclonal '),muttype,'s'), OR=c(f.all$estimate,f.cl$estimate, f.sc$estimate),
                   Conf_low=c(f.all$conf.int[1],f.cl$conf.int[1],f.sc$conf.int[1]), Conf_high=c(f.all$conf.int[2],f.cl$conf.int[2],f.sc$conf.int[2]),
                   P=c(f.all$p.value,f.cl$p.value, f.sc$p.value))
# plot combined result as rectangular markers with error bars
ggplot(f.df) +
  geom_rect(aes(ymin=OR-0.015, ymax = OR+0.015, xmin=as.numeric(as.factor(Mutation))-0.25, xmax=as.numeric(as.factor(Mutation))+0.25, fill=Mutation)) +
  scale_fill_manual(values=c('grey50',colorRampPalette(brewer.pal(8, "PuOr"))(8)[c(7,2)]),
                    labels=c('All','Clonal','Subclonal')) +
  scale_colour_manual(values=c('grey50',colorRampPalette(brewer.pal(8, "PuOr"))(8)[c(7,2)]),
                      labels=c('All','Clonal','Subclonal')) +
  geom_errorbar(aes(x=as.numeric(as.factor(Mutation)), ymin=Conf_low, ymax=Conf_high, colour=Mutation), width=0.1) +
  theme_mypub() + 
  geom_hline(yintercept = 1, linetype='dashed') +
  labs(y='Odds ratio') + guides(fill='none', colour='none') +
  geom_point(data=f.df[f.df$P<0.05,], aes(x=as.numeric(as.factor(Mutation)), y=1.2), shape=8) +
  scale_x_continuous(breaks=c(1,2,3), labels=f.df$Mutation) + theme(axis.title.x = element_blank())


# Mutation RNA depletion and Neo co-occurrence -------------------------------

# Compute and plot odds ratio for neoantigen (present in DNA) being depleted at RNA level
neoIdentifier <- c('SB Neo','SBRP Neo','SBRecopo Neo')
nonneoIdentifier <- c('Non-Neo')
exp.muts.total <- read.delim('RNA/Expression_NA_nonNA_allcancers.txt')
exp.muts.total <- subset(exp.muts.total, (numExprMut + numNotExprMut)>1)

# go through all, clonal and subclonal and compute a Fisher test
# column numNotExprMut is the number of biopsies within a cancer with that mutation being transcriptionally depleted
exp.muts.sub <- subset(exp.muts.total, (Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & !(Category %in% c('adenoma_only','private_adenoma')))
f.all <- fisher.test(exp.muts.sub$Neoantigen %in% neoIdentifier,exp.muts.sub$numNotExprMut>0)
exp.muts.sub <- subset(exp.muts.total, (Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & Category %in% c('truncal','truncal_full'))
f.cl <- fisher.test(exp.muts.sub$Neoantigen %in% neoIdentifier,exp.muts.sub$numNotExprMut>0)
exp.muts.sub <- subset(exp.muts.total, (Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & Category %in% c('private','regional','shared_subclonal','other'))
f.sc <- fisher.test(exp.muts.sub$Neoantigen %in% neoIdentifier,exp.muts.sub$numNotExprMut>0)
# summarise results and plot as rectangular markers with errorbars
f.df <- data.frame(Mutation=c('All','Clonal','Subclonal'), OR=c(f.all$estimate,f.cl$estimate, f.sc$estimate),
                   Conf_low=c(f.all$conf.int[1],f.cl$conf.int[1],f.sc$conf.int[1]), Conf_high=c(f.all$conf.int[2],f.cl$conf.int[2],f.sc$conf.int[2]),
                   P=c(f.all$p.value,f.cl$p.value, f.sc$p.value))
ggplot(f.df) +
  geom_rect(aes(ymin=OR-0.04, ymax = OR+0.04, xmin=as.numeric(as.factor(Mutation))-0.25, xmax=as.numeric(as.factor(Mutation))+0.25, fill=Mutation)) +
  scale_fill_manual(values=c('grey50',colorRampPalette(brewer.pal(8, "PuOr"))(8)[c(7,2)]),
                    labels=c('All','Clonal','Subclonal')) +
  scale_colour_manual(values=c('grey50',colorRampPalette(brewer.pal(8, "PuOr"))(8)[c(7,2)]),
                      labels=c('All','Clonal','Subclonal')) +
  geom_errorbar(aes(x=as.numeric(as.factor(Mutation)), ymin=Conf_low, ymax=Conf_high, colour=Mutation), width=0.1) +
  theme_mypub() + 
  geom_hline(yintercept = 1, linetype='dashed') +
  labs(y='Odds ratio') + guides(fill='none', colour='none') +
  geom_point(data=f.df[f.df$P<0.05,], aes(x=as.numeric(as.factor(Mutation)), y=4), shape=8) +
  scale_x_continuous(breaks=c(1,2,3), labels=f.df$Mutation) + theme(axis.title.x = element_blank())


# Compute and plot the proportion of mutations that are depleted (for each cancer)
exp.muts.total <- read.delim('RNA/Expression_NA_nonNA_allcancers.txt')
exp.muts.total <- subset(exp.muts.total, (numExprMut + numNotExprMut)>1)

# denote depletion in at least one biopsy (evidenced by numNotExprMut>0) as 1/0, compute proportion as average
exp.muts.total$NotExp <- 1*(exp.muts.total$numNotExprMut>0)
exp.muts.sub <- subset(exp.muts.total, (Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & Category %in% c('truncal','truncal_full'))
notexp.pat <- aggregate(exp.muts.sub$NotExp, by=list(exp.muts.sub$Patient, exp.muts.sub$Neoantigen %in% neoIdentifier), mean)
notexp.pat$Escape <- patientEsc.df$Escape[match(notexp.pat$Group.1, patientEsc.df$Patient)]
# exclude patients with only a value for Neo or Non-Neo as no comparison can be made
pat.count <- table(notexp.pat$Group.1)
ggplot(notexp.pat[!(notexp.pat$Group.1 %in% (names(pat.count)[pat.count<2])),], aes(x=Group.2, y=x, colour=Escape)) +
  geom_point(size=2) + geom_line(aes(group=Group.1), alpha=0.7, size=1.2) +
  theme_mypub() + theme(axis.title.x=element_blank()) +
  scale_colour_manual(values=c('grey40','goldenrod','brown')) +
  scale_x_discrete(expand=c(0.1,0.1)) +
  scale_y_continuous(labels=percent_format(), limit=c(0, 1.1)) +
  stat_compare_means(paired=T, comparisons=list(c('TRUE','FALSE'))) +
  labs(y='Proportion of mutations\nnot expressed')


# Compute and plot the average proportion of biopsies within a cancer that have a clonal mutations depleted
exp.muts.total <- read.delim('RNA/Expression_NA_nonNA_allcancers.txt')
exp.muts.total <- subset(exp.muts.total, (numExprMut + numNotExprMut)>1)

# compute percent of biopsies with depletion as number with depletion divided by total number with sufficient reads in RNA
exp.muts.total$NotExprPc <- exp.muts.total$numNotExprMut/(exp.muts.total$numExprMut + exp.muts.total$numNotExprMut)
exp.muts.sub <- subset(exp.muts.total, (Neoantigen %in% c(neoIdentifier, nonneoIdentifier)) & Category %in% c('truncal','truncal_full'))
notexppc.pat <- aggregate(exp.muts.sub$NotExprPc, by=list(exp.muts.sub$Patient, exp.muts.sub$Neoantigen %in% neoIdentifier), mean)
notexppc.pat$Escape <- patientEsc.df$Escape[match(notexppc.pat$Group.1, patientEsc.df$Patient)]
# exclude patients with only a value for Neo or Non-Neo as no comparison can be made
pat.count <- table(notexppc.pat$Group.1)
ggplot(notexppc.pat[!(notexppc.pat$Group.1 %in% (names(pat.count)[pat.count<2])),], aes(x=Group.2, y=x, colour=Escape)) +
  geom_point(size=2) + geom_line(aes(group=Group.1),alpha=0.7, size=1.2) +
  theme_mypub() + theme(axis.title.x=element_blank()) +
  scale_x_discrete(expand=c(0.1,0.1), labels=c('Non-Neo','Neoantigen')) +
  scale_y_continuous(labels=percent_format(), limit=c(0, 0.75)) +
  scale_colour_manual(values=c('grey40','goldenrod','brown')) +
  stat_compare_means(paired=T, comparisons=list(c('TRUE','FALSE'))) +
  labs(y='Proportion of mutated samples\nnot expressing mutation')



