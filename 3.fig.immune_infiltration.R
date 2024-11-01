
source('0.basics.R')


# Infiltrating cells (CyCIF) ----------------------------------------------

ip.inf.df <- readRDS('FFPE_samples/cycif_summary_with_normal.rds')
names(ip.inf.df) <- gsub('\\+','plus',gsub(' ','_',names(ip.inf.df)))
ip.inf.df$sample_type <- factor(ip.inf.df$sample_type, levels=c('Normal','Superficial','Invasive','Node'))

ip.inf.sub <- ip.inf.df

inf.lmer <- lmer(per_epi_cell.Cytotoxic_T_cell ~ sample_type + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(per_epi_cell.Cytotoxic_T_cell ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
p.lmer <- anova(inf.lmer, inf.other)$`Pr(>Chisq)`[2] # comparison of models with and without fixed effect
# CTLs
# note: scale limit omits two outliers from plotting - it does not affect the results
ggplot(ip.inf.sub, aes(y=per_epi_cell.Cytotoxic_T_cell, x=sample_type, fill=sample_type)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1, height=0) + theme_mypub() +
  stat_compare_means(comparisons=list(
    c('Normal','Superficial'),c('Normal','Invasive'),c('Normal','Node')),
    label.y = c(0.2, 0.24, 0.275)) +
  scale_fill_manual(values=setNames(c('skyblue','#6d6d4e','plum2','darkorange4'),c('Normal','Superficial','Invasive','Node'))) +
  scale_y_continuous(limits=c(0, 0.36)) + # comment this out to see full results
  annotate('text', x=1.5, y=0.36, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  guides(fill='none') + labs(y='Cytotoxic T-cells\n[per epi. cell]', x='Sample type')

# CTLA-4 Tregs
inf.lmer <- lmer(per_epi_cell.CD3plusCD4plusCTLA4plusCD45ROplusFOXP3plus ~ sample_type + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(per_epi_cell.CD3plusCD4plusCTLA4plusCD45ROplusFOXP3plus ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
p.lmer <- anova(inf.lmer, inf.other)$`Pr(>Chisq)`[2] # comparison of models with and without fixed effect
ggplot(ip.inf.sub, aes(y=per_epi_cell.CD3plusCD4plusCTLA4plusCD45ROplusFOXP3plus, x=sample_type, fill=sample_type)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1, height=0) + theme_mypub() +
  stat_compare_means(comparisons=list(
    c('Normal','Superficial'),c('Normal','Invasive'),c('Normal','Node')),
    label.y=c(0.19, 0.23, 0.27)) +
  scale_fill_manual(values=setNames(c('skyblue','#6d6d4e','plum2','darkorange4'),c('Normal','Superficial','Invasive','Node'))) +
  annotate('text', x=1.5, y=0.35, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  guides(fill='none') + labs(y='CTLA4+ Tregs\n[per epi. cell]', x='Sample type')

# PD-L1 epithelial cell fraction compared between tumour regions
ip.inf.sub <- subset(ip.inf.df,sample_type!='Normal')

inf.lmer <- lmer(fraction.PDL1plus_epithelial_cells ~ sample_type + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(fraction.PDL1plus_epithelial_cells ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
p.lmer <- anova(inf.lmer, inf.other)$`Pr(>Chisq)`[2] # comparison of models with and without fixed effect
ggplot(ip.inf.sub, aes(y=fraction.PDL1plus_epithelial_cells, x=sample_type, fill=sample_type)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1, height=0) + theme_mypub() +
  stat_compare_means(comparisons=list(c('Superficial','Invasive'),c('Superficial','Node'),c('Node','Invasive')),
                     label.y=c(0.05, 0.07, 0.09)) +
  scale_fill_manual(values=setNames(c('skyblue','#6d6d4e','plum2','darkorange4'),c('Normal','Superficial','Invasive','Node'))) +
  annotate('text', x=1.5, y=0.11, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  guides(fill='none') + labs( y='Fraction of PD-L1+ epi. cells',x='Sample type')

# Other PD-L1+ cells across regions
ip.inf.sub <- subset(ip.inf.df,sample_type!='Normal')
inf.lmer <- lmer(fraction.PDL1plus ~ sample_type + (1 | PatientEPICC), data=ip.inf.sub)
inf.other <- lmer(fraction.PDL1plus ~ 1 + (1 | PatientEPICC), data=ip.inf.sub)
p.lmer <- anova(inf.lmer, inf.other)$`Pr(>Chisq)`[2] # comparison of models with and without fixed effect
ggplot(ip.inf.sub, aes(y=per_epi_cell.PDL1plus, x=sample_type, fill=sample_type)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1, height=0) + theme_mypub() +
  stat_compare_means(comparisons=list(c('Superficial','Invasive'),c('Superficial','Node'),c('Node','Invasive')),
                     label.y=c(0.1, 0.13, 0.16)) +
  scale_fill_manual(values=setNames(c('skyblue','#6d6d4e','plum2','darkorange4'),c('Normal','Superficial','Invasive','Node'))) +
  annotate('text', x=1.5, y=0.22, label=paste0('p(fixed effect)=',scientific(p.lmer,2))) +
  guides(fill='none') + labs( y='PD-L1+ non-epi. cells [per epi cell]',x='Sample type')


# Regression between burden and infiltrates -------------------------------

ip.inf.df <- readRDS('FFPE_samples/cycif_summary_with_normal.rds')
ip.inf.df$sample_type <- factor(ip.inf.df$sample_type, levels=c('Normal','Superficial','Invasive','Node'))
ip.inf.sub <- subset(ip.inf.df, !is.na(patEscape) & sample_type!='Normal')
burden.df <- read.delim('Burden/Burden_master_table.allsample.FFPE.txt')
burden.df$SampleType <- factor(burden.df$SampleType, levels=c('superficial','invasive','node'))
ip.inf.sub$PropBurden <- burden.df$PropBurden[match(ip.inf.sub$roi, burden.df$Sample)]

mv.test <- lmer(PropBurden ~ sample_type + per_epi_cell.VISTAplus + per_epi_cell.Stromal_cells + fraction.PDL1plus_epithelial_cells + per_epi_cell.Cytotoxic_T_cell + per_epi_cell.CD3plusCD4plusCTLA4plusCD45ROplusFOXP3plus +per_epi_cell.CD68plus + per_epi_cell.CD45ROplus + per_epi_cell.CD163plus + per_epi_cell.PD1plus + (1 | PatientEPICC),
                data = ip.inf.sub)

plot_summs(mv.test, coefs=c('VISTA+ cells'='per_epi_cell.VISTAplus','Stromal cells'='per_epi_cell.Stromal_cells',
                            'PD-L1+ epi cells'='fraction.PDL1plus_epithelial_cells','CTLs'='per_epi_cell.Cytotoxic_T_cell',
                            'CTLA-4+ Tregs'='per_epi_cell.CD3plusCD4plusCTLA4plusCD45ROplusFOXP3plus',
                            'PD1+ cells'='per_epi_cell.PD1plus','CD163+ cells'='per_epi_cell.CD163plus',
                            'CD68+ cells'='per_epi_cell.CD68plus', 'CD45RO+ cells' = 'per_epi_cell.CD45ROplus'
),
colors=c('purple4')) +
  theme_mypub() + theme(axis.title.y = element_blank())+ labs(x='Change in prop neoantigen burden')

# Infiltrating lymphocytes (H&E) ------------------------------------------

# Total Lymphocyte count
ip.inf.df <- readRDS('FFPE_samples/classifier_all.rds')
ip.inf.df$sample_type <- factor(ip.inf.df$sample_type, levels=c('Normal','Adj. normal','Superficial','Invasive','Node'))

# Lymphocytes per epi cell
# note: coord_cartesian 1 outlier from plotting to help visualisation, but does not affect test results
ggplot(ip.inf.df, aes(y = lym_to_all_epi, x = sample_type, fill=sample_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2, height=0) +
  coord_cartesian(ylim=c(0,1.65)) +
  scale_fill_manual(values=setNames(c('skyblue','#105e83','#6d6d4e','plum2','darkorange4'),
                                    c('Normal','Adj. normal','Superficial','Invasive','Node'))) +
  stat_compare_means(comparisons=list(c('Normal','Superficial'),c('Adj. normal','Superficial'),
                                      c('Normal','Invasive'),c('Adj. normal','Invasive'),
                                      c('Normal','Node'),c('Adj. normal','Node')),
                     label.y = c(0.7,0.5,1.1,0.9,1.5,1.3))+
  theme_mypub() + labs(y='Lymphocytes [per epi. cell]', x='Sample type') +
  guides(fill='none')

# Lymphocyte fraction
# note: coord_cartesian 1 outlier from plotting to help visualisation, but does not affect test results
ggplot(ip.inf.df, aes(y = lym_fraction, x = sample_type, fill=sample_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2, height=0) +
  #coord_cartesian(ylim=c(0,1.65)) +
  scale_fill_manual(values=setNames(c('skyblue','#105e83','#6d6d4e','plum2','darkorange4'),
                                    c('Normal','Adj. normal','Superficial','Invasive','Node'))) +
  stat_compare_means(comparisons=list(c('Normal','Superficial'),c('Adj. normal','Superficial'),
                                      c('Normal','Invasive'),c('Adj. normal','Invasive'),
                                      c('Normal','Node'),c('Adj. normal','Node')),
                     label.y = c(0.25,0.2,0.35,0.3,0.45,0.4))+
  theme_mypub() + labs(y='Lymphocyte fraction', x='Sample type') +
  guides(fill='none')


# Lymphocytes classified as Infiltrating Lymphocytes
ip.inf.df <- read.delim('FFPE_samples/classifier_itlr.txt')
ip.inf.sub <- ip.inf.df
ip.inf.sub$sample_type <- factor(ip.inf.sub$sample_type, levels=c('Superficial','Invasive','Node'))

ggplot(ip.inf.sub, aes(y=ITLR, x=sample_type, fill=sample_type)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1, height=0) +
  theme_mypub() +
  scale_fill_manual(values=setNames(c('#6d6d4e','plum2','darkorange4'),c('Superficial','Invasive','Node'))) +
  stat_compare_means(comparisons=list(c('Superficial','Invasive'),c('Superficial','Node'),c('Node','Invasive')), label.y = c(0.2, 0.35, 0.4)) +
  guides(fill='none') +
  labs(x='Sample type', y='Infiltrating tumour lymphocyte ratio')


# Distance to CTLs --------------------------------------------------------

# Normal compared to tumour regions
dist.df <- readRDS('FFPE_samples/epicell_CTL_distances.rds')

# note: coord_cartesian omits outliers from plotting to help visualisation, but does not affect test results
ggplot(dist.df, aes(y=min_micron, x=sample_type, fill=sample_type)) +
  geom_boxplot(outlier.shape=NA) + theme_mypub() +
  stat_compare_means(comparisons=list(
    c('Normal','Superficial'),c('Normal','Invasive'),c('Normal','Node')),
    label.y=c(300,350,400)) +
  coord_cartesian(ylim=c(0,450)) +
  scale_fill_manual(values=setNames(c('skyblue','#6d6d4e','plum2','darkorange4'),c('Normal','Superficial','Invasive','Node'))) +
  guides(fill='none') + labs(y=expression('Min. distance to CTL ['~mu~ 'm]'), x='Sample type')


# PD-L1+ compared to other epithelial cells, in different regions
# note: coord_cartesian omits outliers from plotting to help visualisation, but does not affect test results
dist.df <- readRDS('FFPE_samples/epicell_PDL1plus_CTL_distances.rds')

# Distance in Invasive margin
ggplot(dist.df[dist.df$sample_type=='Invasive',], aes(y=min_micron, x=factor, fill=factor)) +
  geom_boxplot(outlier.shape=NA) + theme_mypub() +
  stat_compare_means(comparisons=list(
    c('PDL1+ epithelial cells','Epithelial cells')),
    label.y=c(1150)) +
  coord_cartesian(ylim=c(0,1400)) +
  scale_fill_manual(values=setNames(c('#821b82','plum2'),c('PDL1+ epithelial cells','Epithelial cells'))) +
  guides(fill='none') + labs(y=expression('Min. distance to CTL ['~mu~ 'm]'), x='Cell type') +
  scale_x_discrete(labels=c('Epithelial cells\n(CN: 0, 1, 5, 14)','PDL1+ epithelial cells\n(CN: 3, 10)'))

# Distance in All cancer samples
ggplot(dist.df, aes(y=min_micron, x=factor, fill=factor)) +
  geom_boxplot(outlier.shape=NA) + theme_mypub() +
  stat_compare_means(comparisons=list(
    c('PDL1+ epithelial cells','Epithelial cells')),
    label.y=c(1250)) +
  coord_cartesian(ylim=c(0,1400)) +
  scale_fill_manual(values=setNames(c('#821b82','plum2'),c('PDL1+ epithelial cells','Epithelial cells'))) +
  guides(fill='none') + labs(y=expression('Min. distance to CTL ['~mu~ 'm]'), x='Cell type') +
  scale_x_discrete(labels=c('Epithelial cells\n(CN: 0, 1, 5, 14)','PDL1+ epithelial cells\n(CN: 3, 10)'))

# Distance of only PD-L1+ in different regions
dist.df$sample_type <- factor(dist.df$sample_type, levels=c('Superficial','Invasive','Node'))
ggplot(dist.df[dist.df$factor=='PDL1+ epithelial cells',], aes(y=min_micron, x=sample_type, fill=sample_type)) +
  geom_boxplot(outlier.shape=NA) + theme_mypub() +
  stat_compare_means(comparisons=list(c('Superficial','Invasive'),c('Superficial','Node'),c('Node','Invasive')),
                     label.y = c(600,700,800)) +
  coord_cartesian(ylim=c(0,1000)) +
  scale_fill_manual(values=setNames(c('skyblue','#6d6d4e','plum2','darkorange4'),c('Normal','Superficial','Invasive','Node'))) +
  guides(fill='none') + labs(y=expression('Min. distance to CTL ['~mu~ 'm]'), x='Sample type') +
  scale_x_discrete()


# Cellular neighbourhood distribution -------------------------------------

cn.df <- readRDS('FFPE_samples/cycif_cn_summary.rds')

# Proportion of cells in CN3 in each roi
cn.sub <- cn.df
# compute the proportion of cells in CN3 as the mean of 1 (in CN3) and 0 (not in CN3) indicators
cn.sub$CN3 <- 1*(cn.sub$neighborhood10=='3')
cn3.roi <- aggregate(cn.sub$CN3, by=list(cn.sub$roi, cn.sub$sample_type), mean)
names(cn3.roi) <- c('roi','sample_type','x')
cn3.roi$sample_type <- factor(cn3.roi$sample_type, levels=c('Normal','Superficial','Invasive','Node'))
ggplot(cn3.roi, aes(x=sample_type, y=x, fill=sample_type)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1, height=0) +
  theme_mypub() +
  scale_fill_manual(values=setNames(c('skyblue','#6d6d4e','plum2','darkorange4'),c('Normal','Superficial','Invasive','Node'))) +
  stat_compare_means(comparisons=list(c('Normal','Superficial'),c('Superficial','Invasive'),c('Superficial','Node'),c('Invasive','Node')),
                     label.y = c(0.05, 0.075, 0.1, 0.15)) +
  guides(fill='none') + labs(x='Sample type', y='Fraction of cells in CN3')

# Proportion of cells in CN10 in each roi
cn.sub <- cn.df
cn.sub$CN10 <- 1*(cn.sub$neighborhood10=='10')
cn10.roi <- aggregate(cn.sub$CN10, by=list(cn.sub$roi, cn.sub$sample_type), mean)
names(cn10.roi) <- c('roi','sample_type','x')
cn10.roi$sample_type <- factor(cn10.roi$sample_type, levels=c('Normal','Superficial','Invasive','Node'))
ggplot(cn10.roi, aes(x=sample_type, y=x, fill=sample_type)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.1, height=0) +
  theme_mypub() +
  scale_fill_manual(values=setNames(c('skyblue','#6d6d4e','plum2','darkorange4'),c('Normal','Superficial','Invasive','Node'))) +
  stat_compare_means(comparisons=list(c('Normal','Superficial'),c('Normal','Invasive'),c('Superficial','Invasive'),c('Superficial','Node'))) +
  guides(fill='none') + labs(x='Sample type', y='Fraction of cells in CN10')


# Ripley's H function -----------------------------------------------------

# Spatial distribution of PDL1+ and PDL1 epithelial cells
# note that for this analysis, ki67+ have been clubbed together with general pheno
ripley.df <- readRDS("FFPE_samples/spatialstats_PDL1plus_PD1plus.rds")

normal <- "inc_normal"
stat_to_use <- "RHFunc_PDL1_PD1"

if (normal=="inc_normal"){
  ripley.df<- ripley.df %>%
    mutate(sample_type = fct_relevel(sample_type, 
                                     "Normal", "Superficial", "Invasive", "Node"))
}else{
  ripley.df <- ripley.df[!ripley.df$sample_type=="Normal",]
  ripley.df<- ripley.df %>%
    mutate(sample_type = fct_relevel(sample_type, 
                                     "Superficial", "Invasive", "Node"))
}

ggplot(ripley.df, aes(y = ripley.df[,stat_to_use], x = sample_type)) +
  geom_boxplot(outlier.shape = NA, aes(fill=sample_type)) +
  geom_jitter(width=0.2, height=0)+
  scale_fill_manual(values=setNames(c('skyblue','#6d6d4e','plum2','darkorange4'),c('Normal','Superficial','Invasive','Node'))) +
  stat_compare_means(comparisons=list(
                                      c('Superficial','Invasive'),c('Superficial','Node'))) +
  labs(x='Sample type', y="Ripley's H function") +
  theme_mypub() + guides(fill='none') +
  geom_hline(yintercept = 0, linetype='dashed', alpha=0.5)


