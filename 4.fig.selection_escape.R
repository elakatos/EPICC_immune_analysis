
source('0.basics.R')


# APG mutation and expression ---------------------------------------------

expTable.master <- read.delim('RNA/All_EPICC_tpm.txt')
vst.deseq <- readRDS('RNA/allgenes.vsd.ensembl_withC516adenoma.rds')
geneexp <- as.data.frame(assay(vst.deseq))
apg <- read.delim('Immune_escape/APG_list.txt')
apgmut.master <- read.delim('Immune_escape/EPICC_APG_MUT_master_file.txt')

# Plot expression of a given gene in wild-type and mutated biopsies
gene <- 'B2M' # B2M or NLRC5
ensid <- apg$EnsID[apg$GeneName==gene]
gene.muts.regions <- apgmut.master$Sample[apgmut.master[,gene]>0]
gene.df <- epicc.df[epicc.df$MatchRNA=='Yes',c('Sample','Patient','Region','Type')]
gene.df$Mut <- gene.df$Sample %in% gene.muts.regions
#gene.df$Exp <- as.numeric(expTable.master[expTable.master$GeneID==ensid,gene.df$Sample])
gene.df <- subset(gene.df, Sample %in% names(geneexp))
gene.df$Exp <- as.numeric(geneexp[ensid,gene.df$Sample]) # updated computation based on VST values
ggplot(gene.df, aes(x=Mut, y=Exp)) +
  geom_boxplot(outlier.shape=NA) + theme_mypub() +
  geom_jitter(height=0, width=0.1, aes(colour=Type, shape=Type)) +
  stat_compare_means(comparisons=list(c('FALSE','TRUE'))) +
  scale_x_discrete(labels=c('wild-type','mutated')) + labs(x=gene,y=paste0(gene,' expression (VST)')) +
  scale_shape_discrete(labels=c('deep WGS','low-pass WGS')) +
  scale_colour_manual(values=c('darkblue','skyblue'), labels=c('deep WGS','low-pass WGS'))

# Plot expression of the sum of all HLA type-I in biopsies mutated/wild-type for an APG
gene.df$HLAA <- as.numeric(geneexp["ENSG00000206503",gene.df$Sample])
gene.df$HLAB <- as.numeric(geneexp["ENSG00000234745",gene.df$Sample])
gene.df$HLAC <- as.numeric(geneexp["ENSG00000204525",gene.df$Sample])
ggplot(gene.df, aes(x=Mut, y=HLAA+HLAB+HLAC)) +
  geom_boxplot(outlier.shape=NA) + theme_mypub() +
  geom_jitter(height=0, width=0.2, aes(colour=Type, shape=Type)) +
  stat_compare_means(comparisons=list(c('FALSE','TRUE'))) +
  scale_x_discrete(labels=c('wild-type','mutated')) + labs(x=gene,y=('HLA class-I expression (VST)')) +
  scale_colour_manual(values=c('darkblue','skyblue'),labels=c('deep WGS','low-pass WGS')) + #scale_y_log10() +
  scale_shape_discrete(labels=c('deep WGS','low-pass WGS'))

# Adenoma vs carcinoma ----------------------------------------------------

burden.df <- read.delim('Burden/Burden_master_table.allsample.txt')

# subset to only take into account patients who have an adenoma
burden.sub <- subset(burden.df, Patient %in% adenomaList)
# optional: filter out C530, as its adenomas are of low quality
# burden.sub <- subset(burden.sub, !(Patient %in% c('C530')))
burden.other <- lmer(PropBurden ~ 1 + (1 | Patient), data=burden.sub)
burden.lmer <- lmer(PropBurden ~ Tissue + (1 | Patient), data=burden.sub)
summary(burden.lmer)
p.lmer <- anova(burden.lmer, burden.other)$`Pr(>Chisq)`[2]
ggplot(burden.sub, aes(x=Tissue, y=PropBurden, fill=Tissue)) +
  geom_boxplot(outlier.shape=NA) + geom_jitter(height=0, width=0.1) +
  theme_mypub() +
  scale_fill_manual(values=setNames(c('grey80','grey40'),c('Adenoma','Cancer'))) +
  stat_compare_means(label.y.npc=0.95, comparisons=list(c('Adenoma','Cancer')), label.y=0.65) + 
  guides(fill='none') + labs(y='Proportional neoantigen burden') +
  annotate('text',x=1, y=0.75, label=paste0('p(fixed effect)=',scientific(p.lmer,2)))

# subset to MMRp
burden.sub <- burden.sub[burden.sub$MSI=='MSS',]
burden.other <- lmer(PropBurden ~ 1 + (1 | Patient), data=burden.sub)
burden.lmer <- lmer(PropBurden ~ Tissue + (1 | Patient), data=burden.sub)
summary(burden.lmer)
p.lmer <- anova(burden.lmer, burden.other)$`Pr(>Chisq)`[2]
ggplot(burden.sub, aes(x=Tissue, y=PropBurden, fill=Tissue)) +
  geom_boxplot(outlier.shape=NA) + geom_jitter(height=0, width=0.1) +
  theme_mypub() +
  stat_compare_means(label.y.npc=0.95, comparisons=list(c('Adenoma','Cancer')), label.y=0.65) + 
  scale_fill_manual(values=setNames(c('grey80','grey40'),c('Adenoma','Cancer')) ) +
  guides(fill='none') + labs(y='Proportional neoantigen burden') +
  annotate('text',x=1, y=0.75, label=paste0('p(fixed effect)=',scientific(p.lmer,2)))


# Plot association between total mutation burden and neoantigen burden
ggplot(burden.df, aes(x=TotalWGS, y=SNVBurden, colour=MSI)) + geom_point() +
  theme_mypub() + scale_x_log10() + scale_y_log10() + facet_wrap(.~Tissue) +
  scale_colour_manual(values=setNames(c('#bc432d','#568b7e'),c('MSI','MSS'))) +
  labs(x='Total mut burden (WGS)', y='Total neoantigen burden', colour='') +
  theme(strip.background = element_blank())
ggplot(burden.df, aes(x=TotalWGS, y=PropBurden, colour=MSI)) + geom_point() +
  theme_mypub() + scale_x_log10() + facet_wrap(.~Tissue) +
  scale_colour_manual(values=setNames(c('#bc432d','#568b7e'),c('MSI','MSS'))) +
  labs(x='Total mut burden (WGS)', y='Prop neoantigen burden', colour='') +
  theme(strip.background = element_blank())

burden.sub <- subset(burden.df, Patient %in% adenomaList & MSI=='MSS')
ggplot(burden.sub, aes(x=Tissue, y=TotalWGS, fill=Tissue)) +
  geom_boxplot(outlier.shape=NA) + geom_jitter(height=0, width=0.1) +
  theme_mypub() +
  stat_compare_means(label.y.npc=0.95, comparisons=list(c('Adenoma','Cancer'))) + 
  scale_fill_manual(values=setNames(c('grey80','grey40'),c('Adenoma','Cancer'))) +
  guides(fill='none') + labs(y='Total mut burden')

# Burden association with cancer-level escape -----------------------------

burden.df <- read.delim('Burden/Burden_master_table.allsample.txt')
burden.df$patEscape <- patientEsc.df$EscapeWithEpi[match(burden.df$Patient, patientEsc.df$Patient)] #add epigenetic escape in too
burden.df$patEscape <- factor(burden.df$patEscape, levels=c('No','Epigenetic','Partial','Yes'))

# Proportional neoantigen burden
burden.sub <- subset(burden.df, Tissue=='Cancer' & MSI=='MSS')
burden.other <- lmer(PropBurden ~ 1 + (1 | Patient), data=burden.sub)
burden.lmer <- lmer(PropBurden ~ patEscape + (1 | Patient), data=burden.sub)
summary(burden.lmer)
p.lmer <- anova(burden.lmer, burden.other)$`Pr(>Chisq)`[2]
ggplot(burden.df[burden.df$Tissue=='Cancer' & burden.df$MSI=='MSS',], aes(x=patEscape, y=SNVBurden/TotalSNV, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.2, height=0) +
  theme_mypub() + guides(fill='none') +
  scale_fill_manual(values=c('grey40','yellowgreen','goldenrod','brown')) +
  stat_compare_means(comparisons=list(c('No','Yes'),c('Epigenetic','Partial'),c('Epigenetic','Yes')), label.y = c(0.65, 0.6, 0.625)) +
  labs(y='Proportional neoantigen burden', x='Cancer escape') +
  annotate('text',x=1.5, y=0.7, label=paste0('p(fixed effect)=',scientific(p.lmer,2)))


# Immune dNdS with confidence intervals
burden.df <- subset(burden.df, !is.na(imm_dNdS) & imm_dNdS<Inf)
burden.sub <- burden.df[burden.df$Tissue=='Cancer' & burden.df$MSI=='MSS' & burden.df$patEscape!='Epigenetic',]
median.imm <- aggregate(burden.sub$imm_dNdS, by=list(burden.sub$patEscape), function(x) median(x, na.rm=T))
# note that coord_cartesian omits the top of some large CIs but does not affect test results
ggplot(burden.sub, aes(x=patEscape, y=imm_dNdS, colour=patEscape)) +
  geom_hline(yintercept = 1, linetype='dashed') +
  geom_point(aes(group=Sample),position=position_dodge(width=0.5), size=2) +
  geom_linerange(aes(ymin=imm_lowci,ymax=imm_highci, group=Sample), position=position_dodge(width=0.5), alpha=0.5) +
  theme_mypub() + guides(colour='none') +
  scale_colour_manual(values=c('grey40','goldenrod','brown')) +
  labs(y='Immune dNdS', x='Cancer escape') +
 # coord_cartesian(ylim=c(0,5)) +
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
# add in epigenetic escape information that is originally not present in dataset
# since all epigenetically escaped cases are without genetic escape, we add this escape as a new category, Epi
ip.inf.sub$epiEscape <- patientEsc.df$EpigeneticEscape[match(ip.inf.sub$PatientEPICC, patientEsc.df$Patient)]
ip.inf.sub$patEscape[ip.inf.sub$patEscape=='No' & ip.inf.sub$epiEscape!='No' & !is.na(ip.inf.sub$epiEscape)] <- ip.inf.sub$epiEscape[ip.inf.sub$patEscape=='No' & ip.inf.sub$epiEscape!='No' & !is.na(ip.inf.sub$epiEscape)]
ip.inf.sub$patEscape[ip.inf.sub$patEscape=='Weak'] <- 'Epi'
ip.inf.sub$patEscape <- factor(ip.inf.sub$patEscape, levels=c('No','Epi','Partial','Yes'))

# PD-1+ cells
inf.lmer <- lmer(per_epi_cell.PD1plus ~ patEscape + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(per_epi_cell.PD1plus ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
p.lmer <- anova(inf.lmer, inf.other)$`Pr(>Chisq)`[2] # comparison of models with and without fixed effect
ggplot(ip.inf.sub, aes(x=patEscape, y=per_epi_cell.PD1plus, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(height=0, width=0.1) +
  theme_mypub() +
  labs(y='PD-1+ cells [per epi cell]', fill='Immune\nescape', x='Cancer escape') + guides(fill='none') +
  scale_fill_manual(values=c('grey40','yellowgreen','goldenrod','brown'))  +
  stat_compare_means(comparisons=list(c('No','Partial'),c('Partial','Yes'),c('No','Yes')), label.y = c(0.06,0.075, 0.09)) +
  #annotate('text', x=1.5, y=0.11, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  scale_x_discrete(labels=c('No','Epigenetic','Partial gen.','Genetic'))

# VISTA+ cells
# note that coord_cartesian omits outliers but does not affect test results
inf.lmer <- lmer(per_epi_cell.VISTAplus ~ patEscape + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(per_epi_cell.VISTAplus ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
p.lmer <- anova(inf.lmer, inf.other)$`Pr(>Chisq)`[2]
ggplot(ip.inf.sub, aes(x=patEscape, y=per_epi_cell.VISTAplus, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(height=0, width=0.1) +
  theme_mypub() +
  labs(y='VISTA+ cells [per epi cell]', fill='Immune\nescape', x='Cancer escape') + guides(fill='none') +
  scale_fill_manual(values=c('grey40','yellowgreen','goldenrod','brown'))  +
  stat_compare_means(comparisons=list(c('No','Partial'),c('Partial','Yes'),c('Epi','Partial')), label.y = c(0.04,0.06, 0.08)) +
  #coord_cartesian(ylim=c(0, 0.25)) + # comment out to see full results
  #annotate('text', x=1.5, y=0.25, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  scale_x_discrete(labels=c('No','Epigenetic','Partial gen.','Genetic'))

# Stromal cells
inf.lmer <- lmer(per_epi_cell.Stromal_cells ~ patEscape + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(per_epi_cell.Stromal_cells ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
p.lmer <- anova(inf.lmer, inf.other)$`Pr(>Chisq)`[2]
ggplot(ip.inf.sub, aes(x=patEscape, y=per_epi_cell.Stromal_cells, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(height=0, width=0.1) +
  theme_mypub() +
  labs(y='Stromal cells [per epi cell]', fill='Immune\nescape', x='Cancer escape') + guides(fill='none') +
  scale_fill_manual(values=c('grey40', 'yellowgreen','goldenrod','brown')) +
  stat_compare_means(comparisons=list(c('No','Yes'),c('Epi','Yes'),c('Partial','Yes')), label.y = c(0.52, 0.575, 0.62)) +
  annotate('text', x=1.5, y=0.65, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  scale_x_discrete(labels=c('No','Epigenetic','Partial gen.','Genetic'))

# CD163+ cells
# note that coord_cartesian omits outliers but does not affect test results
inf.lmer <- lmer(per_epi_cell.CD163plus ~ patEscape + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(per_epi_cell.CD163plus ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
p.lmer <- anova(inf.lmer, inf.other)$`Pr(>Chisq)`[2]
ggplot(ip.inf.sub, aes(x=patEscape, y=per_epi_cell.CD163plus, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(height=0, width=0.1) +
  theme_mypub() +
  labs(y='CD163+ cells [per epi cell]', fill='Immune\nescape', x='Cancer escape') + guides(fill='none') +
  scale_fill_manual(values=c('grey40','yellowgreen','goldenrod','brown'))  +
  stat_compare_means(comparisons=list(c('No','Partial')), label.y = c(0.16)) +
  coord_cartesian(ylim=c(0, 0.22)) +  # comment out to see full results
  annotate('text', x=1.5, y=0.21, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  scale_x_discrete(labels=c('No','Epigenetic','Partial gen.','Genetic'))
  

# CTLs
inf.lmer <- lmer(per_epi_cell.Cytotoxic_T_cell ~ patEscape + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(per_epi_cell.Cytotoxic_T_cell ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
p.lmer <- anova(inf.lmer, inf.other)$`Pr(>Chisq)`[2]
ggplot(ip.inf.sub, aes(x=patEscape, y=per_epi_cell.Cytotoxic_T_cell, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(height=0, width=0.1) +
  theme_mypub() +
  labs(y='CTLs [per epi cell]', fill='Immune\nescape', x='Cancer escape') + guides(fill='none') +
  scale_fill_manual(values=c('grey40','yellowgreen','goldenrod','brown')) +
  stat_compare_means(comparisons=list(c('Partial','No'),c('Epi','Yes'),c('No','Yes')), label.y = c(0.68, 0.75, 0.8)) +
  annotate('text', x=1.5, y=0.95, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  scale_x_discrete(labels=c('No','Epigenetic','Partial gen.','Genetic'))

# CTLA-4 Tregs
inf.lmer <- lmer(per_epi_cell.CD3plusCD4plusCTLA4plusCD45ROplusFOXP3plus ~ patEscape + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(per_epi_cell.CD3plusCD4plusCTLA4plusCD45ROplusFOXP3plus ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
p.lmer <- anova(inf.lmer, inf.other)$`Pr(>Chisq)`[2]
ggplot(ip.inf.sub, aes(x=patEscape, y=per_epi_cell.CD3plusCD4plusCTLA4plusCD45ROplusFOXP3plus, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(height=0, width=0.1) +
  theme_mypub() + 
  labs(y='CTLA-4+ Tregs [per epi cell]', fill='Immune\nescape', x='Cancer escape') + guides(fill='none') +
  scale_fill_manual(values=c('grey40','yellowgreen','goldenrod','brown')) +
  stat_compare_means(comparisons=list(c('No','Partial'),c('No','Yes')), label.y = c(0.325, 0.375)) +
  annotate('text', x=1.5, y=0.45, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  scale_x_discrete(labels=c('No','Epigenetic','Partial gen.','Genetic'))

# Comparison of non-escaped and epigenetically escaped cancers only
ip.inf.sub <- subset(ip.inf.sub, (!is.na(epiEscape) & epiEscape!='No') | patEscape=='No')

# CD68+ cells
inf.lmer <- lmer(per_epi_cell.CD68plus ~ epiEscape + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(per_epi_cell.CD68plus ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(per_epi_cell.CD3plusCD4plusCTLA4plusCD45ROplusFOXP3plus ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
p.lmer <- anova(inf.lmer, inf.other)$`Pr(>Chisq)`[2]
ggplot(ip.inf.sub[!is.na(ip.inf.sub$epiEscape),], aes(x=epiEscape, y=per_epi_cell.CD68plus, fill=epiEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(height=0, width=0.1) +
  theme_mypub() +
  stat_compare_means(comparisons=list(c('No','Weak')), label.y=c(0.2)) +
  labs(y='CD68+ cells [per epi cell]', x='Epi. escape') + guides(fill='none') +
  scale_fill_manual(values=c('grey40','yellowgreen')) +
  annotate('text', x=1, y=0.25, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  scale_x_discrete(labels=c('No','Yes'))

# CD45RO+ cells
inf.lmer <- lmer(per_epi_cell.CD45ROplus ~ epiEscape + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(per_epi_cell.CD45ROplus ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
p.lmer <- anova(inf.lmer, inf.other)$`Pr(>Chisq)`[2]
ggplot(ip.inf.sub[!is.na(ip.inf.sub$epiEscape),], aes(x=epiEscape, y=per_epi_cell.CD45ROplus, fill=epiEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(height=0, width=0.1) +
  theme_mypub() +
  stat_compare_means(comparisons=list(c('No','Weak')), label.y=c(0.14)) +
  labs(y='CD45RO+ cells [per epi cell]', x='Epi. escape') + guides(fill='none') +
  scale_fill_manual(values=c('grey40','yellowgreen')) +
  annotate('text', x=1, y=0.24, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  scale_x_discrete(labels=c('No','Yes'))


# TGFbeta signalling and escape status ------------------------------------

vst.deseq <- readRDS('RNA/allgenes.vsd.ensembl_withC516adenoma.rds')
geneexp <- as.data.frame(assay(vst.deseq))
burden.df <- read.delim('Burden/Burden_master_table.allsample.txt')
burden.df$patEscape <- patientEsc.df$EscapeWithEpi[match(burden.df$Patient,patientEsc.df$Patient)] # update with epigenetic escape
burden.df$patEscape <- factor(burden.df$patEscape, levels = c('No','Epigenetic','Partial','Yes'))

# F-TBRS signature, following Mariathasan et al., 2018
genes.ftbrs <- c('ACTA2', 'ACTG2', 'ADAM12', 'ADAM19', 'CNN1', 'COL4A1', 'CCN2', 'CTPS1', 'RFLNB', 'FSTL3', 'HSPB1', 'IGFBP3', 'PXDC1', 'SEMA7A', 'SH3PXD2A', 'TAGLN', 'TGFBI', 'TNS1', 'TPM1')
ensid.ftbrs <- c("ENSG00000120708","ENSG00000183688","ENSG00000168994","ENSG00000187498","ENSG00000070404","ENSG00000149591","ENSG00000138623","ENSG00000107957",
                 "ENSG00000118523","ENSG00000135074","ENSG00000163017","ENSG00000171793","ENSG00000107796","ENSG00000130176","ENSG00000140416","ENSG00000146674",
                 "ENSG00000106211","ENSG00000079308","ENSG00000148848")
exp.ftbrs <- geneexp[ensid.ftbrs,]
z.ftbrs <- t(scale(t(exp.ftbrs))) # transform to z score
z.pca <- prcomp(t(z.ftbrs), scale. = F, center = F); sig.ftbrs <- z.pca$x[,1, drop=F] # take first principal component as signature

burden.df$FTBRS <- sig.ftbrs[match(burden.df$Sample, row.names(sig.ftbrs))]
burden.sub <- subset(burden.df, Tissue=='Cancer' & MSI=='MSS') #limit to MMRp CRCs
burden.other <- lmer(FTBRS ~ 1 + (1 | Patient), data=burden.sub)
burden.lmer <- lmer(FTBRS ~ patEscape + (1 | Patient), data=burden.sub)
p.lmer <- anova(burden.lmer, burden.other)$`Pr(>Chisq)`[2]
ggplot(burden.df[burden.df$Tissue=='Cancer' & burden.df$MSI=='MSS',], aes(x=patEscape, y=FTBRS, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.2, height=0) +
  theme_mypub() + guides(fill='none') +
  scale_fill_manual(values=c('grey40','yellowgreen', 'goldenrod','brown')) +
  stat_compare_means(comparisons=list(c('No','Epigenetic'),c('Epigenetic','Yes'), c('No','Partial')),
                     label.y = c(9.5, 8, 6.5)) +
  annotate('text', x=1.5, y=12, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  labs(y='F-TBRS', x='Cancer escape')

# TGBFBR2 expression
gene <- 'TGFBR2'
ensid <- 'ENSG00000163513'
exp.tmp <- t(geneexp[ensid, ])

burden.df$Gene.Exp <- exp.tmp[match(burden.df$Sample,row.names(exp.tmp)),]
burden.sub <- subset(burden.df, Tissue=='Cancer' & MSI=='MSS')
burden.other <- lmer(Gene.Exp ~ 1 + (1 | Patient), data=burden.sub)
burden.lmer <- lmer(Gene.Exp ~ patEscape + (1 | Patient), data=burden.sub)
p.lmer <- anova(burden.lmer, burden.other)$`Pr(>Chisq)`[2]
ggplot(burden.df[burden.df$Tissue=='Cancer' & burden.df$MSI=='MSS',], aes(x=patEscape, y=Gene.Exp, fill=patEscape)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.2, height=0) +
  theme_mypub() + guides(fill='none') +
  scale_fill_manual(values=c('grey40','yellowgreen', 'goldenrod','brown')) +
  stat_compare_means(comparisons=list(c('No','Yes'),c('Partial','Yes'),c('Epigenetic','Yes')),
                     label.y=c(18.5,15.5, 17)) +
  annotate('text', x=1.5, y=20.5, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  labs(y=paste0(gene, ' expression (VST)'), x='Cancer escape')

