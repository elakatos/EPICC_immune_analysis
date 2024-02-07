
source('0.basics.R')


# APG mutation and expression ---------------------------------------------

expTable.master <- read.delim('RNA/All_EPICC_tpm.txt')
apg <- read.delim('Immune_escape/APG_list.txt')
apgmut.master <- read.delim('Immmune_escape/EPICC_APGMUT_master_file.txt')

# Plot expression of a given gene in wild-type and mutated biopsies
gene <- 'RFXAP' #or B2M or NLRC5
ensid <- apg$EnsID[apg$GeneName==gene]
gene.muts.regions <- apgmut.master$Sample[apgmut.master[,gene]>0]
gene.df <- epicc.df[epicc.df$MatchRNA=='Yes',c('Sample','Patient','Region','Type')]
gene.df$Mut <- gene.df$Sample %in% gene.muts.regions
gene.df$Exp <- as.numeric(expTable.master[expTable.master$GeneID==ensid,gene.df$Sample])
ggplot(gene.df, aes(x=Mut, y=Exp)) +
  geom_boxplot(outlier.shape=NA) + theme_mypub() +
  geom_jitter(height=0, width=0.1, aes(colour=Type, shape=Type)) +
  stat_compare_means(comparisons=list(c('FALSE','TRUE'))) +
  scale_x_discrete(labels=c('wild-type','mutated')) + labs(x=gene,y=paste0(gene,' expression (TPM)')) +
  scale_colour_manual(values=c('darkblue','skyblue'), labels=c('deep WGS','low-pass WGS'))

# Plot expression of the sum of all HLA type-I in biopsies mutated/wild-type for an APG
gene.df$HLAA <- as.numeric(expTable.master[expTable.master$GeneID=="ENSG00000206503",gene.df$Sample])
gene.df$HLAB <- as.numeric(expTable.master[expTable.master$GeneID=="ENSG00000234745",gene.df$Sample])
gene.df$HLAC <- as.numeric(expTable.master[expTable.master$GeneID=="ENSG00000204525",gene.df$Sample])
ggplot(gene.df, aes(x=Mut, y=HLAA+HLAB+HLAC)) +
  geom_boxplot(outlier.shape=NA) + theme_mypub() +
  geom_jitter(height=0, width=0.2, aes(colour=Type, shape=Type)) +
  stat_compare_means(comparisons=list(c('FALSE','TRUE'))) +
  scale_x_discrete(labels=c('wild-type','mutated')) + labs(x=gene,y=('HLA class-I expression (TPM)')) +
  scale_colour_manual(values=c('darkblue','skyblue'),labels=c('deep WGS','low-pass WGS')) + scale_y_log10() +
  scale_shape_discrete(labels=c('deep WGS','low-pass WGS'))

# Adenoma vs carcinoma ----------------------------------------------------

burden.df <- read.delim('Burden/Burden_master_table.allsample.txt')

# subset to only take into account patients who have an adenoma
burden.sub <- subset(burden.df, Patient %in% adenomaList)
ggplot(burden.sub, aes(x=Tissue, y=PropBurden, fill=Tissue)) +
  geom_boxplot(outlier.shape=NA) + geom_jitter(height=0, width=0.1) +
  theme_mypub() +
  stat_compare_means(label.y.npc=0.95, comparisons=list(c('Adenoma','Cancer'))) + 
  scale_fill_manual(values=setNames(c('grey80','grey40'),c('Adenoma','Cancer'))) +
  guides(fill='none') + labs(y='Proportional neoantigen burden')



# Burden association with cancer-level escape -----------------------------

burden.df <- read.delim('Burden/Burden_master_table.allsample.txt')

# Proportional neoantigen burden
ggplot(burden.df[burden.df$Tissue=='Cancer' & burden.df$MSI=='MSS',], aes(x=patEscape, y=SNVBurden/TotalSNV, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.2, height=0) +
  theme_mypub() + guides(fill='none') +
  scale_fill_manual(values=c('grey40','goldenrod','brown')) +
  stat_compare_means(comparisons=list(c('No','Yes'),c('Partial','Yes'),c('No','Partial')), label.y = c(0.635, 0.615, 0.6)) +
  labs(y='Proportional neoantigen burden', x='Cancer escape')

# Immune dNdS with confidence intervals
burden.df <- subset(burden.df, !is.na(imm_dNdS) & imm_dNdS<Inf)
burden.sub <- burden.df[burden.df$Tissue=='Cancer' & burden.df$MSI=='MSS',]
median.imm <- aggregate(burden.sub$imm_dNdS, by=list(burden.sub$patEscape), function(x) median(x, na.rm=T))
# note that coord_cartesian omits the top of some large CIs but does not affect test results
ggplot(burden.sub, aes(x=patEscape, y=imm_dNdS, colour=patEscape)) +
  geom_hline(yintercept = 1, linetype='dashed') +
  geom_point(aes(group=Sample),position=position_dodge(width=0.5), size=2) +
  geom_linerange(aes(ymin=imm_lowci,ymax=imm_highci, group=Sample), position=position_dodge(width=0.5), alpha=0.5) +
  theme_mypub() + guides(colour='none') +
  scale_colour_manual(values=c('grey40','goldenrod','brown')) +
  labs(y='Immune dNdS', x='Cancer escape') +
  coord_cartesian(ylim=c(0,5)) +
  geom_crossbar(data=median.imm, aes(x = Group.1, colour=Group.1, y=x, ymin = x, ymax = x), width = 0.5) +
  stat_compare_means(comparisons=list(c('No','Yes'),c('Partial','Yes'),c('No','Partial')), label.y=c(4.25, 3.75, 3.25))

# Immune dNdS comparison of FF-WGS and FFPE-PS
total.dnds.df <- read.delim('Burden/Immdnds_comparison_table.allsharedsample.txt')

median.imm <- aggregate(total.dnds.df$imm_dnds, by=list(total.dnds.df$patEscape, total.dnds.df$Type), function(x) median(x, na.rm=T))
names(median.imm) <- c('patEscape','Type','imm_dnds')
ggplot(total.dnds.df[total.dnds.df$patEscape!='Yes',], aes(y=imm_dnds, x=Type, colour=patEscape)) +
  geom_hline(yintercept = 1, linetype='dashed') +
  geom_point(aes(group=Sample),position=position_dodge(width=0.5), size=2) +
  geom_linerange(aes(ymin=imm_lowci,ymax=imm_highci, group=Sample), position=position_dodge(width=0.5), alpha=0.5) +
  theme_mypub() + guides(colour='none') +
  scale_colour_manual(values=c('grey40','goldenrod','brown')) +
  labs(y='Immune dNdS', x='Cancer escape') + facet_wrap(.~paste0(patEscape,' cancer escape'))+ 
  stat_compare_means(comparison=list(c('FF_WGS','FFPE_PS')), label.y = c(3.2, 3.2)) +
  coord_cartesian(ylim =c(0, 3.5)) +
  geom_crossbar(data=median.imm[median.imm$patEscape!='Yes',], aes(ymin = imm_dnds, ymax = imm_dnds), width = 0.5) +
  theme(strip.background = element_blank(), axis.title.x=element_blank())


# VAF depletion association with escape -----------------------------------

burden.df <- read.delim('Burden/Burden_master_table.allsample.txt')
vaf.df.total <- read.delim('Subclonality/NA_nonNA_VAF.perSample.allcancers.txt')

# VAF association in escaped cancers
# set filters for: SNVs, fully escaped MMRp cancers, non-clonal mutations
typeFilt <- 'SNV'
patientFilt <- setdiff(patientEsc.df$Patient[patientEsc.df$Escape=='Yes'], msiList)
categoryFilt <- c('private','regional','shared_subclonal','other')
sampleFilt <- burden.df$Sample[burden.df$Escape=='Yes']
neoIdentifier <- c('SB Neo','SBRecopo Neo','SBRP Neo') # combine all categories with strong neoantigens
nonneoIdentifier <- 'Non-Neo'

vaf.df.total <- subset(vaf.df.total, variable %in% epicc.df$Sample[epicc.df$Tissue=='Cancer']) #exclude adenomas
vaf.df.m <- subset(vaf.df.total, (Category %in% categoryFilt) & MutType==typeFilt & (Patient %in% patientFilt) & (variable %in% sampleFilt) & Neoantigen %in% c(neoIdentifier, nonneoIdentifier))

# Comparison of the proportion of mutations with VAF<0.3 (biopsy-by-biopsy) between SB Neo & non-Neo
vaf.df.tab <- as.data.frame(table(vaf.df.m$Neoantigen %in% neoIdentifier, vaf.df.m$value<0.3, vaf.df.m$variable))
vaf.df.denom <- as.data.frame(table(vaf.df.m$Neoantigen %in% neoIdentifier, vaf.df.m$variable))
vaf.df.tab$Total <- vaf.df.denom$Freq[match(paste0(vaf.df.tab$Var1, vaf.df.tab$Var3), paste0(vaf.df.denom$Var1, vaf.df.denom$Var2))]
vaf.df.tab$Proportion <- vaf.df.tab$Freq/vaf.df.tab$Total

ggplot(vaf.df.tab[vaf.df.tab$Var2!=F,], aes(x=Var1, y=Proportion*100)) +
  geom_violin() + geom_point(size=2, colour='grey40') + geom_line(aes(group=Var3), alpha=0.5) +
  theme_mypub() + theme(axis.title.x = element_blank()) +
  stat_compare_means(paired=T, comparisons=list(c('FALSE','TRUE'))) +
  scale_x_discrete(labels=c('Non-Neo','SB Neo')) + labs(y='Mutations VAF<0.3 [%]')

# Cumulative VAF curve of all biopsies
ks <- ks.test(vaf.df.m$invVAF[vaf.df.m$Neoantigen %in% neoIdentifier], vaf.df.m$invVAF[vaf.df.m$Neoantigen==nonneoIdentifier])
vafseq <- seq(0, max(vaf.df.m$invVAF)+1,by=1)
vaf.cum.df <- data.frame(invVAF=vafseq,
                         cumNeo=sapply(vafseq, function(z) sum(vaf.df.m$invVAF<=z & (vaf.df.m$Neoantigen %in% neoIdentifier), na.rm=T)/sum(vaf.df.m$Neoantigen %in% neoIdentifier,na.rm=T)),
                         cumNonNeo=sapply(vafseq, function(z) sum(vaf.df.m$invVAF<=z & vaf.df.m$Neoantigen==nonneoIdentifier, na.rm=T)/sum(vaf.df.m$Neoantigen==nonneoIdentifier,na.rm=T)),
                         cumAllNeo=sapply(vafseq, function(z) sum(vaf.df.m$invVAF<=z & (vaf.df.m$Neoantigen!=nonneoIdentifier), na.rm=T)/sum(vaf.df.m$Neoantigen!=nonneoIdentifier,na.rm=T))
)

ggplot(vaf.cum.df, aes(x=invVAF,y=cumNeo)) +
  geom_point(colour='firebrick') + geom_line(colour='firebrick') +
  geom_line(aes(y=cumNonNeo)) +
  theme_mypub() + annotate('text', x=5, y=0.7, label=paste0('KS-test p=',round(ks$p.value,3))) +
  labs(x='1/VAF',y='Cumulative count')

# VAF association in partially escaped cancers
# set filters for: SNVs, Escape biopsies within partially escaped MMRp cancers, non-clonal mutations
typeFilt <- 'SNV'
patientFilt <- setdiff(patientEsc.df$Patient[patientEsc.df$Escape=='Partial'], msiList)
categoryFilt <- c('private','regional','shared_subclonal','other')
sampleFilt <- burden.df$Sample[burden.df$Escape=='Yes']
neoIdentifier <- c('SB Neo','SBRecopo Neo','SBRP Neo') # combine all categories with strong neoantigens
nonneoIdentifier <- 'Non-Neo'

vaf.df.total <- subset(vaf.df.total, variable %in% epicc.df$Sample[epicc.df$Tissue=='Cancer']) #exclude adenomas
vaf.df.m <- subset(vaf.df.total, (Category %in% categoryFilt) & MutType==typeFilt & (Patient %in% patientFilt) & (variable %in% sampleFilt) & Neoantigen %in% c(neoIdentifier, nonneoIdentifier))

# Comparison of the proportion of mutations with VAF<0.3 (biopsy-by-biopsy) between SB Neo & non-Neo
vaf.df.tab <- as.data.frame(table(vaf.df.m$Neoantigen %in% neoIdentifier, vaf.df.m$value<0.3, vaf.df.m$variable))
vaf.df.denom <- as.data.frame(table(vaf.df.m$Neoantigen %in% neoIdentifier, vaf.df.m$variable))
vaf.df.tab$Total <- vaf.df.denom$Freq[match(paste0(vaf.df.tab$Var1, vaf.df.tab$Var3), paste0(vaf.df.denom$Var1, vaf.df.denom$Var2))]
vaf.df.tab$Proportion <- vaf.df.tab$Freq/vaf.df.tab$Total

ggplot(vaf.df.tab[vaf.df.tab$Var2!=F,], aes(x=Var1, y=Proportion*100)) +
  geom_violin() + geom_point(size=2, colour='grey40') + geom_line(aes(group=Var3), alpha=0.5) +
  theme_mypub() + theme(axis.title.x = element_blank()) +
  stat_compare_means(paired=T, comparisons=list(c('FALSE','TRUE'))) +
  scale_x_discrete(labels=c('Non-Neo','SB Neo')) + labs(y='Mutations VAF<0.3 [%]')

# Cumulative VAF curve of all biopsies
ks <- ks.test(vaf.df.m$invVAF[vaf.df.m$Neoantigen %in% neoIdentifier], vaf.df.m$invVAF[vaf.df.m$Neoantigen==nonneoIdentifier])
vafseq <- seq(0, max(vaf.df.m$invVAF)+1,by=1)
vaf.cum.df <- data.frame(invVAF=vafseq,
                         cumNeo=sapply(vafseq, function(z) sum(vaf.df.m$invVAF<=z & (vaf.df.m$Neoantigen %in% neoIdentifier), na.rm=T)/sum(vaf.df.m$Neoantigen %in% neoIdentifier,na.rm=T)),
                         cumNonNeo=sapply(vafseq, function(z) sum(vaf.df.m$invVAF<=z & vaf.df.m$Neoantigen==nonneoIdentifier, na.rm=T)/sum(vaf.df.m$Neoantigen==nonneoIdentifier,na.rm=T)),
                         cumAllNeo=sapply(vafseq, function(z) sum(vaf.df.m$invVAF<=z & (vaf.df.m$Neoantigen!=nonneoIdentifier), na.rm=T)/sum(vaf.df.m$Neoantigen!=nonneoIdentifier,na.rm=T))
)

ggplot(vaf.cum.df, aes(x=invVAF,y=cumNeo)) +
  geom_point(colour='firebrick') + geom_line(colour='firebrick') +
  geom_line(aes(y=cumNonNeo)) +
  theme_mypub() + annotate('text', x=15, y=0.7, label=paste0('KS-test p=',round(ks$p.value,3))) +
  labs(x='1/VAF',y='Cumulative count')


# Immune infiltrate association with escape -------------------------------

ip.inf.df <- readRDS('FFPE_samples/cycif_summary_with_normal.rds')
ip.inf.df$sample_type <- factor(ip.inf.df$sample_type, levels=c('Normal','Superficial','Invasive','Node'))
ip.inf.sub <- subset(ip.inf.df, !is.na(patEscape) & sample_type!='Normal')

# PD-1+ cells
ggplot(ip.inf.sub, aes(y=per_epi_cell.PD1plus, x=patEscape, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1, height=0) +
  theme_mypub() +
  scale_fill_manual(values=c('grey40','goldenrod','brown')) +
  stat_compare_means(comparisons=list(c('No','Partial'),c('Partial','Yes'),c('No','Yes'))) +
  guides(fill='none') +
  labs(x='Cancer escape', y='PD-1+ cells [per epi cell]') 

# VISTA+ cells
# note that coord_cartesian omits outliers but does not affect test results
ggplot(ip.inf.sub, aes(y=per_epi_cell.VISTAplus, x=patEscape, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1, height=0) +
  coord_cartesian(ylim=c(0, 0.35)) + # comment out to see full results
  theme_mypub() +
  scale_fill_manual(values=c('grey40','goldenrod','brown')) +
  stat_compare_means(comparisons=list(c('No','Partial'),c('No','Yes'),c('Partial','Yes')),
                     label.y = c(0.1, 0.15, 0.2)) +
  guides(fill='none') +
  labs(x='Cancer escape', y='VISTA+ cells [per epi cell]')

# Stromal cells
ggplot(ip.inf.sub, aes(y=per_epi_cell.Stromal_cells, x=patEscape, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1, height=0) +
  theme_mypub() +
  scale_fill_manual(values=c('grey40','goldenrod','brown')) +
  stat_compare_means(comparisons=list(c('No','Partial'),c('Partial','Yes'),c('No','Yes'))) +
  guides(fill='none') +
  labs(x='Cancer escape', y='Stromal cells [per epi cell]') 

# CTLs
ggplot(ip.inf.sub, aes(y=per_epi_cell.Cytotoxic_T_cell, x=patEscape, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1, height=0) +
  theme_mypub() +
  scale_fill_manual(values=c('grey40','goldenrod','brown')) +
  stat_compare_means(comparisons=list(c('No','Partial'),c('No','Yes'),c('Partial','Yes')),
                     label.y=c(0.4, 0.5, 0.6)) +
  guides(fill='none') +
  labs(x='Cancer escape', y='CTLs [per epi cell]') 

# CTLA-4 Tregs
ggplot(ip.inf.sub, aes(y=per_epi_cell.CD3plusCD4plusCTLA4plusCD45ROplusFOXP3plus, x=patEscape, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1, height=0) +
  theme_mypub() +
  scale_fill_manual(values=c('grey40','goldenrod','brown')) +
  stat_compare_means(comparisons=list(c('No','Partial'),c('No','Yes'),c('Partial','Yes')),
                     label.y=c(0.3, 0.35, 0.4)) +
  guides(fill='none') +
  labs(x='Cancer escape', y='CTLA-4+ Tregs [per epi cell]') 




