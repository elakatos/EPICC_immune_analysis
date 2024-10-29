
source('0.basics.R')


# Multivariable regression (FF-WGS) ---------------------------------------

burden.df <- read.delim('Burden/Burden_master_table.allsample.txt')
burden.df$MSI <- factor(burden.df$MSI, levels=c('MSS','MSI'))
burden.df$Tissue <- relevel(as.factor(burden.df$Tissue), ref='Cancer')
burden.df$Patient <- paste0('=',burden.df$Patient); burden.df$Patient <- relevel(as.factor(burden.df$Patient), ref='=C531')

# Beta regression for proportional burden of FF-WGS
# note that Patient is not plotted in main figure, but included below
all.test.mv <- betareg(PropBurden ~ SampleType + Tissue + Patient + Escape + Purity, data=burden.df)
summary(all.test.mv)
plot_summs(all.test.mv, coefs=c('Type=Gland'='SampleTypeGland','Tissue=Adenoma'='TissueAdenoma',
                                'Potential escape'='EscapeWeak','Escape'='EscapeYes','Purity'='Purity'),
           colors=c('purple4')) +
  theme_mypub() + theme(axis.title.y = element_blank())+ labs(x='Change in prop neoantigen burden')
# include patient (smallest and largest estimated coeff) as well to show in comparison to other variables
which.min(all.test.mv$coefficients$mean[startsWith(names(all.test.mv$coefficients$mean), 'Patient')])
which.max(all.test.mv$coefficients$mean[startsWith(names(all.test.mv$coefficients$mean), 'Patient')])
plot_summs(all.test.mv, coefs=c('Type=Gland'='SampleTypeGland','Tissue=Adenoma'='TissueAdenoma',
                                'Potential escape'='EscapeWeak','Escape'='EscapeYes','Purity'='Purity',
                                'Patient (min/max)'='Patient=C518', 'Patient (min/max)'='Patient=C519'),
           colors=c('purple4')) +
  theme_mypub() + theme(axis.title.y = element_blank())+ labs(x='Change in prop neoantigen burden')

# full plot with all inter-tumour effects as well
plot_summs(all.test.mv, colors=c('purple4'), omit.coefs = c('(phi)','(Intercept)')) +
  theme_mypub() +
  theme(axis.title.y=element_blank()) + labs(x='Change in prop neoantigen burden')


# Linear regression on immune dNdS of FF-WGS
burden.sub <- subset(burden.df, !is.na(imm_dNdS) & imm_dNdS<Inf)
dnds.test.mv <- lm(imm_dNdS~ SampleType + Tissue + Patient + Escape + Purity, data=burden.sub)
summary(dnds.test.mv)
plot_summs(dnds.test.mv, coefs=c('Type=Gland'='SampleTypeGland','Tissue=Adenoma'='TissueAdenoma',
                                 'Potential escape'='EscapeWeak','Escape'='EscapeYes','Purity'='Purity'),
           colors=c('darkblue')) +
  theme_mypub() + theme(axis.title.y = element_blank())+ labs(x='Change in immune dNdS')
# full plot with all inter-tumour effects as well
dnds.test.mv <- lm(imm_dNdS~ SampleType + Tissue + Patient + Escape + Purity, data=burden.sub)
plot_summs(dnds.test.mv, omit.coefs = c('(phi)','(Intercept)'),
           colors=c('darkblue')) +
  theme_mypub() + theme(axis.title.y = element_blank())+ labs(x='Change in immune dNdS')



# Multivariable regression (FFPE-PS) --------------------------------------

ip.burden.df <- read.delim('Burden/Burden_master_table.allsample.FFPE.txt')
ip.burden.df$SampleType <- relevel(as.factor(ip.burden.df$SampleType), ref='superficial')
ip.burden.df$PatientEPICC <- paste0('=',ip.burden.df$PatientEPICC); ip.burden.df$Patient <- relevel(as.factor(ip.burden.df$PatientEPICC),ref='=C531')

# Beta regression on proportional burden of FFPE-PS
all.test.mv <- betareg(PropBurden ~ SampleType + Patient + Purity, data=ip.burden.df)
summary(all.test.mv)
plot_summs(all.test.mv, coefs=c('Type=Node'='SampleTypenode','Type=Invasive'='SampleTypeinvasive',
                                'Purity'='Purity'),
           colors=c('purple4')) +
  theme_mypub() + theme(axis.title.y = element_blank())+ labs(x='Change in prop neoantigen burden')

# include patient (smallest and largest estimated coeff) as well to show in comparison to other variables
which.min(all.test.mv$coefficients$mean[startsWith(names(all.test.mv$coefficients$mean), 'Patient')])
which.max(all.test.mv$coefficients$mean[startsWith(names(all.test.mv$coefficients$mean), 'Patient')])
plot_summs(all.test.mv, coefs=c('Type=Node'='SampleTypenode','Type=Invasive'='SampleTypeinvasive',
                                'Purity'='Purity', 'Patient (min/max)'='Patient=C537','Patient (min/max)'='Patient=C550'),
           colors=c('purple4')) +
  theme_mypub() + theme(axis.title.y = element_blank())+ labs(x='Change in prop neoantigen burden')
# full plot with all inter-tumour effects as well
plot_summs(all.test.mv, omit.coefs = c('(phi)','(Intercept)'),
           colors=c('purple4')) +
  theme_mypub() + theme(axis.title.y = element_blank())+ labs(x='Change in prop neoantigen burden')


# Linear regression on immune dNdS of FFPE-PS
ip.burden.sub <- subset(ip.burden.df, imm_dnds<Inf & !is.na(imm_dnds) & !is.nan(imm_dnds))
all.test.mv.dnds <- lm(imm_dnds ~ SampleType + Patient + Purity, data=ip.burden.sub)
summary(all.test.mv.dnds)
plot_summs(all.test.mv.dnds, coefs=c('Type=Node'='SampleTypenode','Type=Invasive'='SampleTypeinvasive',
                                     'Purity'='Purity'),
           colors=c('darkblue')) +
  theme_mypub() + theme(axis.title.y = element_blank())+ labs(x='Change in immune dNdS')
# full plot
plot_summs(all.test.mv.dnds,exp=T, omit.coefs = c('(phi)','(Intercept)'),
           colors=c('darkblue')) +
  theme_mypub() + theme(axis.title.y = element_blank())+ labs(x='Change in immune dNdS')



# Clonal vs subclonal burden ----------------------------------------------

# (Mean) clonal vs subclonal proportional neoantigen burden in FF-WGS
burdenPerReg.df.comp <- read.delim('Burden/Neoantigen_burden_deepWGS.SBmut.comparetruncal.txt')
# exclude adenomas and biopsies with too few mutations
burdenPerReg.df.comp.sub <- burdenPerReg.df.comp[burdenPerReg.df.comp$Tissue=='Cancer' & burdenPerReg.df.comp$TotalSNV>10,]
burdenPerReg.comp.m <- reshape2::melt(burdenPerReg.df.comp.sub[,c('Sample','Patient','patEscape','PropBurden','PropBurdenTruncal')], id=c('Patient','patEscape','Sample'))
burdenPerReg.comp.perpatient <- aggregate(burdenPerReg.comp.m$value, by=list(burdenPerReg.comp.m$Patient, burdenPerReg.comp.m$patEscape,
                                                                             burdenPerReg.comp.m$variable), mean)
names(burdenPerReg.comp.perpatient) <- c('Patient','patEscape','variable','value')
burdenPerReg.comp.perpatient$variable <- factor(burdenPerReg.comp.perpatient$variable, levels=c('PropBurdenTruncal','PropBurden'))
ggplot(burdenPerReg.comp.perpatient, aes(x=variable, y=value, colour=patEscape)) +
  geom_point(size=2) + geom_line(aes(group=Patient), alpha=0.65, linewidth=1.5) +
  scale_colour_manual(values=c('grey40','goldenrod','brown')) +
  theme_mypub() + theme(axis.title.x=element_blank()) +
  scale_x_discrete(labels=c('Clonal SNVs','Subclonal SNVs'), expand = c(0.1,0.1)) +
  stat_compare_means(comparisons=list(c('PropBurden','PropBurdenTruncal')), paired=T) +
  labs(y='Prop. neoantigen burden', colour='Escape')

# (Mean) clonal vs subclonal proportional neoantigen burden in FFPE-PS
burdenPerReg.df.comp <- read.delim('Burden/Immunopanel_burden_all.SBmut.comparetruncal.txt')
burdenPerReg.df.comp.sub <- burdenPerReg.df.comp[burdenPerReg.df.comp$TotalSNV>10,]
burdenPerReg.comp.m <- reshape2::melt(burdenPerReg.df.comp.sub[,c('Sample','Patient','patEscape','PropBurden','PropBurdenTruncal')], id=c('Patient','patEscape','Sample'))
burdenPerReg.comp.perpatient <- aggregate(burdenPerReg.comp.m$value, by=list(burdenPerReg.comp.m$Patient, burdenPerReg.comp.m$patEscape,
                                                                             burdenPerReg.comp.m$variable), mean)
names(burdenPerReg.comp.perpatient) <- c('Patient','patEscape','variable','value')
burdenPerReg.comp.perpatient$variable <- factor(burdenPerReg.comp.perpatient$variable, levels=c('PropBurdenTruncal','PropBurden'))
ggplot(burdenPerReg.comp.perpatient, aes(x=variable, y=value, colour=patEscape)) +
  geom_point(size=2) + geom_line(aes(group=Patient), alpha=0.65, linewidth=1.5) +
  scale_colour_manual(values=c('grey40','goldenrod','brown')) +
  theme_mypub() + theme(axis.title.x=element_blank()) +
  scale_x_discrete(labels=c('Clonal SNVs','Subclonal SNVs'), expand = c(0.1,0.1)) +
  stat_compare_means(comparisons=list(c('PropBurden','PropBurdenTruncal')), paired=T) +
  labs(y='Prop. neoantigen burden', colour='Escape')


# Burden distance between biopsies ----------------------------------------

burden.df <- read.delim('Burden/Burden_master_table.allsample.txt')
# limit to cancer samples
burden.df <- subset(burden.df, Tissue=='Cancer')
row.names(burden.df) <- burden.df$Sample

# Distance of proportional burden computed for all comparison and then divided up into three categories
var <- 'PropBurden'
dist.matrix <- dist(burden.df[,var,drop=F])
dist.df <- reshape2::melt(as.matrix(dist.matrix))
dist.df <- subset(dist.df, Var1!=Var2)
dist.df$CompID <- apply(dist.df[,1:2], 1, function(z) paste0(sort(z), collapse = ','))
dist.df <- dist.df[!duplicated(dist.df$CompID),]
dist.df$Patient1 <- burden.df$Patient[match(dist.df$Var1, burden.df$Sample)]
dist.df$Patient2 <- burden.df$Patient[match(dist.df$Var2, burden.df$Sample)]
dist.df$Region1 <- epicc.df$Region[match(dist.df$Var1, epicc.df$Sample)]
dist.df$Region2 <- epicc.df$Region[match(dist.df$Var2, epicc.df$Sample)]

dist.df$Category <- ifelse(dist.df$Patient1!=dist.df$Patient2, 'between-patient',
                           ifelse(dist.df$Region1==dist.df$Region2, 'within-region', 'between-region'))

ggplot(dist.df, aes(x=Category, y=value, fill=Category)) + geom_violin() +
  geom_boxplot(width=0.075, fill='grey80') +
  labs(x='Comparison type', y=paste0('Prop neoantiben burden','\ndistance')) +
  scale_fill_manual(values=c('#7570b3','#1b9e77','#d95f02')) +
  theme_mypub() + guides(fill='none') +
  stat_compare_means(comparisons=list(c('between-region','within-region'),c('between-patient','between-region')),
                     label.y = c(0.15,0.2))


# Distance of immune dNdS computed for all comparison and then divided up into three categories
var <- 'imm_dNdS'
dist.matrix <- dist(burden.df[,var,drop=F])
dist.df <- reshape2::melt(as.matrix(dist.matrix))
dist.df <- subset(dist.df, Var1!=Var2)
dist.df$CompID <- apply(dist.df[,1:2], 1, function(z) paste0(sort(z), collapse = ','))
dist.df <- dist.df[!duplicated(dist.df$CompID),]
dist.df$Patient1 <- burden.df$Patient[match(dist.df$Var1, burden.df$Sample)]
dist.df$Patient2 <- burden.df$Patient[match(dist.df$Var2, burden.df$Sample)]
dist.df$Region1 <- epicc.df$Region[match(dist.df$Var1, epicc.df$Sample)]
dist.df$Region2 <- epicc.df$Region[match(dist.df$Var2, epicc.df$Sample)]

dist.df$Category <- ifelse(dist.df$Patient1!=dist.df$Patient2, 'between-patient',
                           ifelse(dist.df$Region1==dist.df$Region2, 'within-region', 'between-region'))

ggplot(dist.df, aes(x=Category, y=value, fill=Category)) + geom_violin() +
  geom_boxplot(width=0.075, fill='grey80') +
  labs(x='Comparison type', y=paste0('Immune dNdS','\ndistance')) +
  scale_fill_manual(values=c('#7570b3','#1b9e77','#d95f02')) +
  theme_mypub() + guides(fill='none') +
  stat_compare_means(comparisons=list(c('between-region','within-region'),c('between-patient','between-region')),
                     label.y = c(2,2.5))



# VAF distribution for subclonal ------------------------------------------

burden.df <- read.delim('Burden/Burden_master_table.allsample.txt')
vaf.df.total <- read.delim('Subclonality/NA_nonNA_VAF.perSample.allcancers.txt')

# VAF association for subclonal NAs in all non-escaped cancers/biopsies
# set filters for: SNVs, Non-escape biopsies within non-escaped cancers, subclonal mutations
typeFilt <- 'SNV'
patientFilt <- patientEsc.df$Patient[patientEsc.df$Escape=='No']
categoryFilt <- c('private','regional')
sampleFilt <- burden.df$Sample[burden.df$Escape=='No']
neoIdentifier <- c('SB Neo','SBRecopo Neo','SBRP Neo') # combine all categories with strong neoantigens
nonneoIdentifier <- 'Non-Neo'

vaf.df.total <- subset(vaf.df.total, variable %in% epicc.df$Sample[epicc.df$Tissue=='Cancer']) #exclude adenomas
vaf.df.m <- subset(vaf.df.total, (Category %in% categoryFilt) & MutType==typeFilt & (Patient %in% patientFilt) & (variable %in% sampleFilt) & Neoantigen %in% c(neoIdentifier, nonneoIdentifier))

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

# Distribution of mutation sharedness -------------------------------------

# define what is neoantigen and what is not
neoIdentifier <- c('SB Neo','SBRecopo Neo', 'SBRP Neo')
nonneoIdentifier <- c('Non-Neo')
snvInd <- 'SNV'

for (samp in samples){
  # this mutation table contains mutations in each row and has a column (numSamplesCancer) with the number of biopsies that mutation is present in
  # set which mutation table to use to plot for FF-WGS or FFPE-PS
  mutTable <- read.delim(paste0('Subclonality/',samp,'.total.extended.txt')) #FF-WGS
 # mutTable <- read.delim(paste0('Subclonality/',samp,'.total.extended.FFPE.txt')) #FFPE-PS
  
  # this is all mutations, which I annotated for coding change, so I first filter it to be coding and SNV
  mutTable <- subset(mutTable, !is.na(CodingChange) & !is.na(MutType) & MutType==snvInd & numSamplesCancer>0)
  # additionally filtered to only contain annotated non-neoantigen
  mutTable.mut <- subset(mutTable, Neoantigen %in% nonneoIdentifier)
  if(nrow(mutTable.mut)==0){next}
  count.tab.mut <- as.data.frame(table(mutTable.mut$numSamplesCancer)); count.tab.mut$Type <- 'Non-Neo'
  count.tab.mut$CumFreq <- cumsum(count.tab.mut$Freq)/sum(count.tab.mut$Freq)
  
  # similarly, for neoantigens filter to only contain those that fit the definition
  mutTable.ep <- subset(mutTable, Neoantigen %in% neoIdentifier)
  count.tab.ep <- as.data.frame(table(mutTable.ep$numSamplesCancer))
  count.tab.ep$CumFreq <- cumsum(count.tab.ep$Freq)/sum(count.tab.ep$Freq)
  if (nrow(count.tab.ep)==0){next}
  count.tab.ep$Type <- 'SB NA'
  
  count.tab <- rbind(count.tab.mut, count.tab.ep)
  count.tab$Var1 <- as.numeric(as.character(count.tab$Var1))
  count.tab <- rbind(count.tab, data.frame(Var1=c(0,0),Freq=c(0,0),Type=c('Non-Neo','SB NA'),CumFreq=c(0,0)))
  
  # Kolmogorov-Smirnov test between the set of biopsy-count values
  ptest <- ks.test(mutTable.ep$numSamplesCancer,
                   mutTable.mut$numSamplesCancer)
  # plot cumulative distribution with K-S test p-value
  p <- ggplot(count.tab, aes(x=Var1, y=CumFreq, colour=Type)) + 
    geom_point(size=2, aes(shape=Type)) +
    geom_step(aes(group=Type)) +
    theme_mypub() + scale_colour_manual(values=c('grey50','firebrick')) +
    labs(x='Biopsies with mutation', y='Cumulative count',
         title=paste0(samp)) +
    annotate('text', label=paste0('p=',round(ptest$p.value,2)), x=2, y=0.75) +
    scale_x_continuous(limits=c(0, max(count.tab$Var1)))
  print(p)
}


# Normalised burden with subclonal escape ---------------------------------

measures.total <- read.delim('Subclonality/Subclonal_escape_normalised.deepWGSclose.txt')
measures.total$Patient <- as.factor(measures.total$Patient)

# Normalised proportional burden in subclonally escaped regions, with t-test comparing to mu=1
mutt.prop <- t.test(measures.total$PropBurden[measures.total$Mut=='mutation'], mu=1)
ggplot(measures.total[measures.total$Mut=='mutation',], aes(y=PropBurden,x=Mut)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.25, height=0, aes(colour=Patient)) +
  theme_mypub() + theme(axis.title.x=element_blank()) +
  geom_hline(yintercept = 1, linetype='dashed') +
  labs(y='Normalised proportional \nneoantigen burden') +
  annotate('text', y=1.15, x=1, label=paste0('t-test p=',scientific(mutt.prop$p.value,digits=3))) +
  scale_colour_brewer(palette='Set1')

# Normalised immune dNdS in subclonally escaped regions, with t-test comparing to mu=1
mutt.dnds <- t.test(measures.total$imm_dNdS[measures.total$Mut=='mutation'], mu=1)
ggplot(measures.total[measures.total$Mut=='mutation',], aes(y=imm_dNdS,x=Mut)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.25, height=0, aes(colour=Patient)) +
  theme_mypub() + theme(axis.title.x=element_blank()) +
  geom_hline(yintercept = 1, linetype='dashed') +
  labs(y='Normalised \nimmune dNdS') +
  annotate('text', y=1.4, x=1, label=paste0('t-test p=',scientific(mutt.dnds$p.value,digits=3))) +
  scale_colour_brewer(palette='Set1')



# VAF distribution of sample types in FFPE --------------------------------

ip.vaf.df <- read.delim('FFPE_samples/IgNA_nonIgNA_VAF.perSample.allcancers.FFPE.txt')
neoIdentifier <- 'Na'
nonneoIdentifier <- 'nonIgNa'

# VAF association (depletion signal) within Invasive margin
# Compare immunogenic and non-immunogenic mutation VAF
typeFilter <- 'invasive margin'
ip.vaf.df.m <- subset(ip.vaf.df, sample_type==typeFilter) #only include samples from specific biopsy types

# Plot cumulative 1/VAF curve of all biopsies and compute KS-test to detect depletion
ks <- ks.test(ip.vaf.df.m$inverse_VAF[ip.vaf.df.m$con %in% neoIdentifier], ip.vaf.df.m$inverse_VAF[ip.vaf.df.m$con==nonneoIdentifier])
vafseq <- seq(0, max(ip.vaf.df.m$inverse_VAF)+1,by=1)
vaf.cum.df <- data.frame(invVAF=vafseq,
                         cumNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & (ip.vaf.df.m$con %in% neoIdentifier), na.rm=T)/sum(ip.vaf.df.m$con %in% neoIdentifier,na.rm=T)),
                         cumNonNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & ip.vaf.df.m$con==nonneoIdentifier, na.rm=T)/sum(ip.vaf.df.m$con==nonneoIdentifier,na.rm=T)),
                         cumAllNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & (ip.vaf.df.m$con!=nonneoIdentifier), na.rm=T)/sum(ip.vaf.df.m$con!=nonneoIdentifier,na.rm=T))
)

ggplot(vaf.cum.df, aes(x=invVAF,y=cumNeo)) +
  geom_point(colour='firebrick', size=1) + geom_line(colour='firebrick') +
  geom_line(aes(y=cumNonNeo)) +
  theme_mypub() + annotate('text', x=150, y=0.7, label=paste0('KS-test p=',round(ks$p.value,4))) +
  labs(x='1/VAF',y='Cumulative count')

# Comparison of the proportion of mutations with VAF in specified range (biopsy-by-biopsy) between imm. Neo & non-imm.
# Using range for small subclones, 0.05<VAF<0.1
vaf_limit1 <-0.05
vaf_limit2 <-0.1

vaf.df.tab <- as.data.frame(table(ip.vaf.df.m$con %in% neoIdentifier, ip.vaf.df.m$VAF<vaf_limit2 & ip.vaf.df.m$VAF>=vaf_limit1, ip.vaf.df.m$tumour_region))
vaf.df.denom <- as.data.frame(table(ip.vaf.df.m$con %in% neoIdentifier, ip.vaf.df.m$tumour_region))
vaf.df.tab$Total <- vaf.df.denom$Freq[match(paste0(vaf.df.tab$Var1, vaf.df.tab$Var3), paste0(vaf.df.denom$Var1, vaf.df.denom$Var2))]
vaf.df.tab$Proportion <- vaf.df.tab$Freq/vaf.df.tab$Total

ggplot(vaf.df.tab[vaf.df.tab$Var2!=F,], aes(x=Var1, y=Proportion*100)) +
  geom_violin() + geom_point(size=2, colour='grey40') + geom_line(aes(group=Var3), alpha=0.5) +
  theme_mypub() + theme(axis.title.x = element_blank()) +
  stat_compare_means(paired=T, comparisons=list(c('FALSE','TRUE'))) +
  scale_x_discrete(labels=c('Non. imm.','Imm. Neo')) + labs(y='Mutations 0.05<VAF<0.1 [%]')


# VAF association (depletion signal) within Superficial samples
# Compare immunogenic and non-immunogenic mutation VAF
typeFilter <- 'superficial'
ip.vaf.df.m <- subset(ip.vaf.df, sample_type==typeFilter) #only include samples from specific biopsy types

# Plot cumulative 1/VAF curve of all biopsies and compute KS-test to detect depletion
ks <- ks.test(ip.vaf.df.m$inverse_VAF[ip.vaf.df.m$con %in% neoIdentifier], ip.vaf.df.m$inverse_VAF[ip.vaf.df.m$con==nonneoIdentifier])
vafseq <- seq(0, max(ip.vaf.df.m$inverse_VAF)+1,by=1)
vaf.cum.df <- data.frame(invVAF=vafseq,
                         cumNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & (ip.vaf.df.m$con %in% neoIdentifier), na.rm=T)/sum(ip.vaf.df.m$con %in% neoIdentifier,na.rm=T)),
                         cumNonNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & ip.vaf.df.m$con==nonneoIdentifier, na.rm=T)/sum(ip.vaf.df.m$con==nonneoIdentifier,na.rm=T)),
                         cumAllNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & (ip.vaf.df.m$con!=nonneoIdentifier), na.rm=T)/sum(ip.vaf.df.m$con!=nonneoIdentifier,na.rm=T))
)
ggplot(vaf.cum.df, aes(x=invVAF,y=cumNeo)) +
  geom_point(colour='firebrick', size=1) + geom_line(colour='firebrick') +
  geom_line(aes(y=cumNonNeo)) +
  theme_mypub() + annotate('text', x=150, y=0.7, label=paste0('KS-test p=',round(ks$p.value,4))) +
  labs(x='1/VAF',y='Cumulative count')


# VAF association (depletion signal) within Node samples
# Compare immunogenic and non-immunogenic mutation VAF
typeFilter <- 'node'
ip.vaf.df.m <- subset(ip.vaf.df, sample_type==typeFilter) #only include samples from specific biopsy types

# Plot cumulative 1/VAF curve of all biopsies and compute KS-test to detect depletion
ks <- ks.test(ip.vaf.df.m$inverse_VAF[ip.vaf.df.m$con %in% neoIdentifier], ip.vaf.df.m$inverse_VAF[ip.vaf.df.m$con==nonneoIdentifier])
vafseq <- seq(0, max(ip.vaf.df.m$inverse_VAF)+1,by=1)
vaf.cum.df <- data.frame(invVAF=vafseq,
                         cumNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & (ip.vaf.df.m$con %in% neoIdentifier), na.rm=T)/sum(ip.vaf.df.m$con %in% neoIdentifier,na.rm=T)),
                         cumNonNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & ip.vaf.df.m$con==nonneoIdentifier, na.rm=T)/sum(ip.vaf.df.m$con==nonneoIdentifier,na.rm=T)),
                         cumAllNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & (ip.vaf.df.m$con!=nonneoIdentifier), na.rm=T)/sum(ip.vaf.df.m$con!=nonneoIdentifier,na.rm=T))
)
ggplot(vaf.cum.df, aes(x=invVAF,y=cumNeo)) +
  geom_point(colour='firebrick', size=1) + geom_line(colour='firebrick') +
  geom_line(aes(y=cumNonNeo)) +
  theme_mypub() + annotate('text', x=150, y=0.7, label=paste0('KS-test p=',round(ks$p.value,4))) +
  labs(x='1/VAF',y='Cumulative count')


# VAF association (depletion signal) within Invasive samples, using weaker immun. definition
# Compare SB Neo and non-immunogenic mutation VAF
ip.vaf.df <- read.delim('FFPE_samples/SBNA_nonIgNA_VAF.perSample.allcancers.FFPE.txt')
neoIdentifier <- 'SBNa'
nonneoIdentifier <- 'nonIgNa'
typeFilter <- 'invasive margin'
ip.vaf.df.m <- subset(ip.vaf.df, sample_type==typeFilter) #only include samples from specific biopsy types

# Plot cumulative 1/VAF curve of all biopsies and compute KS-test to detect depletion
ks <- ks.test(ip.vaf.df.m$inverse_VAF[ip.vaf.df.m$con %in% neoIdentifier], ip.vaf.df.m$inverse_VAF[ip.vaf.df.m$con==nonneoIdentifier])
vafseq <- seq(0, max(ip.vaf.df.m$inverse_VAF)+1,by=1)
vaf.cum.df <- data.frame(invVAF=vafseq,
                         cumNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & (ip.vaf.df.m$con %in% neoIdentifier), na.rm=T)/sum(ip.vaf.df.m$con %in% neoIdentifier,na.rm=T)),
                         cumNonNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & ip.vaf.df.m$con==nonneoIdentifier, na.rm=T)/sum(ip.vaf.df.m$con==nonneoIdentifier,na.rm=T)),
                         cumAllNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & (ip.vaf.df.m$con!=nonneoIdentifier), na.rm=T)/sum(ip.vaf.df.m$con!=nonneoIdentifier,na.rm=T))
)
ggplot(vaf.cum.df, aes(x=invVAF,y=cumNeo)) +
  geom_point(colour='firebrick', size=1) + geom_line(colour='firebrick') +
  geom_line(aes(y=cumNonNeo)) +
  theme_mypub() + annotate('text', x=150, y=0.7, label=paste0('KS-test p=',round(ks$p.value,4))) +
  labs(x='1/VAF',y='Cumulative count')

# VAF association (depletion signal) within Superficial samples, using weaker immun. definition
# Compare SB Neo and non-immunogenic mutation VAF
typeFilter <- 'superficial'
ip.vaf.df.m <- subset(ip.vaf.df, sample_type==typeFilter) #only include samples from specific biopsy types

# Plot cumulative 1/VAF curve of all biopsies and compute KS-test to detect depletion
ks <- ks.test(ip.vaf.df.m$inverse_VAF[ip.vaf.df.m$con %in% neoIdentifier], ip.vaf.df.m$inverse_VAF[ip.vaf.df.m$con==nonneoIdentifier])
vafseq <- seq(0, max(ip.vaf.df.m$inverse_VAF)+1,by=1)
vaf.cum.df <- data.frame(invVAF=vafseq,
                         cumNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & (ip.vaf.df.m$con %in% neoIdentifier), na.rm=T)/sum(ip.vaf.df.m$con %in% neoIdentifier,na.rm=T)),
                         cumNonNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & ip.vaf.df.m$con==nonneoIdentifier, na.rm=T)/sum(ip.vaf.df.m$con==nonneoIdentifier,na.rm=T)),
                         cumAllNeo=sapply(vafseq, function(z) sum(ip.vaf.df.m$inverse_VAF<=z & (ip.vaf.df.m$con!=nonneoIdentifier), na.rm=T)/sum(ip.vaf.df.m$con!=nonneoIdentifier,na.rm=T))
)
ggplot(vaf.cum.df, aes(x=invVAF,y=cumNeo)) +
  geom_point(colour='firebrick', size=1) + geom_line(colour='firebrick') +
  geom_line(aes(y=cumNonNeo)) +
  theme_mypub() + annotate('text', x=150, y=0.7, label=paste0('KS-test p=',round(ks$p.value,4))) +
  labs(x='1/VAF',y='Cumulative count')
