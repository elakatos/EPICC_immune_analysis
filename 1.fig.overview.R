
source('0.basics.R')

# Overview figure of FF-WGS samples ---------------------------------------

burden.df <- read.delim('Burden/Burden_master_table.allsample.txt')

apg.df <- read.delim('Immune_escape/EPICC_APG_MUT_master_file.txt')
apg.df <- apg.df[,1:34]
apg.scaa.df <- read.delim('Immune_escape/EPICC_APG_SCAA_master_file.perBiopsy.txt')
apg.scaa.df[,-1] <- 2*as.numeric(apg.scaa.df[,-1]<0)
apg.df$TAP2 <- 0 # add tap2 to fill out to full list
apg.scaa.df <- apg.scaa.df[match(apg.df$Sample, apg.scaa.df$Sample), names(apg.df)[-1]]
apg.scaa.df[is.na(apg.scaa.df)] <- 0
apg.df[,-1] <- apg.df[,-1] + apg.scaa.df
apg.df <- apg.df[,c(T,colSums(apg.df[,-1])>0)]
apg.df <- apg.df[,c(1,1+order(colSums(apg.df[,-1]),decreasing=T))]

# reorder patients as MSI-MSS and high-low proportional burden
burden.cancer <- subset(burden.df, Tissue=='Cancer')
xPat.df <- aggregate(burden.cancer$PropBurden, by=list(burden.cancer$Patient, burden.cancer$MSI), mean); names(xPat.df)[1:2] <- c('Patient','MSI')
xPat.df <- xPat.df[order(xPat.df$MSI, (-xPat.df$x)),]
xPat.df$PosShift <- 3*(0:28)

#xPat.df <- data.frame(Patient=unique(burden.df$Patient), PosShift=3*(0:28))
burden.df <- burden.df[order( xPat.df$PosShift[match(burden.df$Patient, xPat.df$Patient)] ),]
burden.df$Pos <- 1:nrow(burden.df); burden.df$Pos <- burden.df$Pos + xPat.df$PosShift[match(burden.df$Patient, xPat.df$Patient)]

apg.df <- apg.df[match(burden.df$Sample,apg.df$Sample),]
apg.df[,-1] <- apply(apg.df[,-1], 2, as.character)
burden.df <- cbind(burden.df, apg.df[,-1])
burden.df$imm_dNdS[burden.df$imm_dNdS==Inf] <- NA
burden.df$imm_highci[burden.df$imm_highci==Inf] <- NA
burden.df$imm_highci[burden.df$imm_highci>3.5] <- 3.5

# Colour code for immune escape and sample type categories
totalCols <- setNames(c('grey60','goldenrod','goldenrod','brown', '#baba9e','#6d6d4e', '#bc432d','#568b7e','darkblue','skyblue', 'tan4','skyblue3',
                        'deepskyblue4','deepskyblue4','grey80', '#d582a0','#b63d4d','#894487','grey80','grey80'),
                      c('No','Weak','Partial','Yes','Adenoma','Cancer', 'MSI','MSS','Deep','LP','1','2',
                        'weak_evidence_LOH','LOH','AI','weak_evidence_mut','HLA_snv','HLA_fs','FALSE','0'))

# Panel for burdens and categorical variables
sh <- 7 # offset to shift different burdens along the y axis
ggplot(burden.df) + theme_mypub() +
  geom_bar(aes(x=Pos, y=log10(TotalSyn)*2), stat='identity') +
  geom_point(aes(x=Pos, y=6+PropBurden*10), size=1.5) +
  geom_point(aes(x=Pos, y=9+sh+(imm_dNdS-1)*2), shape=15, colour='darkblue') +
  geom_errorbar(aes(x=Pos, ymin=9+sh+(imm_lowci-1)*2, ymax=9+sh+(imm_highci-1)*2), width=0.1, alpha=0.8, colour='#494977')+
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=16.5+sh, ymax=17.5+sh, fill=MSI)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=15.5+sh, ymax=16.5+sh, fill=Tissue)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=14.5+sh, ymax=15.5+sh, fill=Type)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-2.5, ymax=-1.5, fill=HLAsnv)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-3.5, ymax=-2.5, fill=HLAfs)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-4.5, ymax=-3.5, fill=HLAloh)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-5.8, ymax=-5.5, fill=B2M)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-6.1, ymax=-5.8, fill=CALR)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-6.4, ymax=-6.1, fill=CIITA)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-6.7, ymax=-6.4, fill=CREB1)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-7, ymax=-6.7, fill=ERAP1)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-7.3, ymax=-7, fill=ERAP2)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-7.6, ymax=-7.3, fill=TBK1)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-7.9, ymax=-7.6, fill=NLRC5)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-8.2, ymax=-7.9, fill=HSP90AA1)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-8.5, ymax=-8.2, fill=HSPA2)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-8.8, ymax=-8.5, fill=HSPA4)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-9.1, ymax=-8.8, fill=HSPA5)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-9.4, ymax=-9.1, fill=HSPA6)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-9.7, ymax=-9.4, fill=PSMA7)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-10, ymax=-9.7, fill=PSME2)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-10.3, ymax=-10, fill=NFYA)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-10.6, ymax=-10.3, fill=NFYB)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-10.9, ymax=-10.6, fill=NFYC)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-11.2, ymax=-10.9, fill=RFX5)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-11.5, ymax=-11.2, fill=RFXAP)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-11.8, ymax=-11.5, fill=RFXANK)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-12.1, ymax=-11.8, fill=TAP2)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-12.4, ymax=-12.1, fill=TAPBP)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-12.7, ymax=-12.4, fill=TBK1)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-14.5, ymax=-13.5, fill=Escape)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-15.5, ymax=-14.5, fill=PatEscape)) +
  #geom_hline(yintercept = c(0,5, 9,8)) +
  geom_hline(yintercept = c(0,4,6,11,(9+sh),(8+sh),(11+sh))) +
  scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0)) +
  scale_fill_manual(values=totalCols) + theme(axis.title = element_blank())

# Separate checkpoint expression panel (continuous fill scale)
ggplot(burden.df) + theme_mypub() +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=0, ymax=1, fill=PDL1)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-1, ymax=0, fill=CTLA4)) +
  scale_fill_gradient(low='papayawhip',high='darkred',na.value='grey80' ,limits=c(0,2.431)) +
  scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0)) +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = 'top') +
  labs(fill='IC expr.')



# Overview of FFPE-PS samples ---------------------------------------------

ip.burden.df <- read.delim('Burden/Burden_master_table.allsample.FFPE.txt')

ip.inf.df <- readRDS('FFPE_samples/cycif_summary_with_normal.rds')
names(ip.inf.df) <- gsub('\\+','plus',gsub(' ','_',names(ip.inf.df)))
ip.inf.plot <- ip.inf.df[match(ip.burden.df$Sample, ip.inf.df$roi),c(paste0('per_epi_cell.',c('Cytotoxic_T_cell','PD1plus','CTLA4plus','CD3plusCD4plusCTLA4plusCD45ROplusFOXP3plus',
                                                                                              'VISTAplus','CD163plus','CD68plus', 'Stromal_cells')),
                                                                     'fraction.PDL1plus_epithelial_cells')]
ip.burden.df <- cbind(ip.burden.df, ip.inf.plot)
ip.burden.df$SampleType <- factor(ip.burden.df$SampleType, levels=c('normal','superficial','invasive','node'))
ip.burden.df <- ip.burden.df[order(ip.burden.df$PatientEPICC, ip.burden.df$SampleType),]

# reorder patients as high-low proportional burden
xPat.df <- aggregate(ip.burden.df$PropBurden, by=list(ip.burden.df$PatientEPICC), function(x) mean(x, na.rm=T)); names(xPat.df) <- c('PatientEPICC','x')
xPat.df <- xPat.df[order(xPat.df$x, decreasing = T),]
xPat.df$PosShift <- 3*(0:10)

ip.burden.df <- ip.burden.df[order( xPat.df$PosShift[match(ip.burden.df$PatientEPICC, xPat.df$PatientEPICC) ]),]
ip.burden.df$imm_highci[ip.burden.df$imm_highci>3.5] <- 3.5

#xPat.df <- data.frame(Patient=unique(ip.burden.df$PatientEPICC), PosShift=3*(0:10))
ip.burden.df$Pos <- 1:nrow(ip.burden.df); ip.burden.df$Pos <- ip.burden.df$Pos + xPat.df$PosShift[match(ip.burden.df$PatientEPICC, xPat.df$Patient)]

# Plot burdens and infiltration (continouous fill scale)
sh <- 9 #offset to shift different burdens along the y axis
ggplot(ip.burden.df) + theme_mypub() +
  geom_point(aes(x=Pos, y=0+PropBurden*10), size=1.5, colour='purple4') +
  geom_point(aes(x=Pos, y=3+sh+(imm_dnds-1)*2), shape=15, colour='darkblue') +
  geom_errorbar(aes(x=Pos, ymin=3+sh+(imm_lowci-1)*2, ymax=3+sh+(imm_highci-1)*2), width=0.1, alpha=0.8, colour='#494977') +
  #geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=16.5+sh, ymax=17.5+sh, fill=SampleType)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-2.5, ymax=-1.5, fill=per_epi_cell.Cytotoxic_T_cell)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-3.5, ymax=-2.5, fill=per_epi_cell.CD3plusCD4plusCTLA4plusCD45ROplusFOXP3plus)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-4.5, ymax=-3.5, fill=per_epi_cell.VISTAplus)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-5.5, ymax=-4.5, fill=per_epi_cell.CD163plus)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-6.5, ymax=-5.5, fill=per_epi_cell.CD68plus)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-7.5, ymax=-6.5, fill=per_epi_cell.PD1plus)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-8.5, ymax=-7.5, fill=per_epi_cell.CTLA4plus)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-9.5, ymax=-8.5, fill=per_epi_cell.Stromal_cells)) +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=-10.5, ymax=-9.5, fill=fraction.PDL1plus_epithelial_cells)) +
  scale_fill_gradient2(low='papayawhip',mid='darkred',high='#290000',na.value='grey80' ,trans='log') +
  scale_y_continuous(expand=c(0,0), breaks=c(0)) + scale_x_continuous(expand=c(0,0)) +
  geom_hline(yintercept = c(0,5,(3+sh),(2+sh),(5+sh))) +
  labs(fill='inf') +
  theme(axis.title = element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_blank())

# Plot burdens and categorical variable
ggplot(ip.burden.df) + theme_mypub() +
  geom_point(aes(x=Pos, y=0+PropBurden*10), size=1.5, colour='purple4') +
  geom_point(aes(x=Pos, y=3+sh+(imm_dnds-1)*2), shape=15, colour='darkblue') +
  geom_rect(aes(xmin=Pos-0.5, xmax=Pos+0.5,ymin=7+sh, ymax=8+sh, fill=SampleType)) +
  scale_fill_manual(values=setNames(c('skyblue','#6d6d4e','plum2','darkorange4'),c('normal','superficial','invasive','node'))) +
  geom_hline(yintercept = c(0,5,(3+sh),(2+sh),(5+sh))) +
  scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0)) +
  theme(axis.title = element_blank(), axis.ticks.x=element_blank(), axis.text.x = element_blank())


# Proportional burden vs immune dNdS --------------------------------------

# In FF-WGS samples, faceted by MMRp/d status
burden.df <- read.delim('Burden/Burden_master_table.allsample.txt')
burden.df$patEscape <- patientEsc.df$EscapeWithEpi[match(burden.df$Patient, patientEsc.df$Patient)]
burden.df$patEscape <- factor(burden.df$patEscape, levels=c('No','Epigenetic','Partial','Yes'))

burden.sub <- subset(burden.df, !is.na(imm_dNdS) & imm_dNdS<Inf)
ggplot(burden.sub, aes(x=PropBurden, y=imm_dNdS, colour=patEscape)) +
  geom_hline(yintercept = 1, colour='grey50', linetype='dashed') +
  geom_point(size=2, aes(shape=Tissue)) +
  scale_shape_manual(values=c(17,16)) +
  scale_colour_manual(values=c('grey40','yellowgreen','goldenrod','brown')) +
  geom_linerange(aes(ymin=imm_lowci,ymax=imm_highci, group=Sample), alpha=0.5) +
  theme_mypub() + stat_cor(label.x.npc = 0.525, label.y.npc=0.35) +
  coord_cartesian(ylim=c(0,5)) + # comment out to see full CIs
  facet_wrap(.~MSI) + labs(x='Proportional neoantigen burden',y='Immune dNdS', colour='Cancer\nescape') +
  theme(strip.background = element_blank())


# In FFPE-PS samples
burden.df <- read.delim('Burden/Burden_master_table.allsample.FFPE.txt')
burden.df$patEscape <- patientEsc.df$EscapeWithEpi[match(burden.df$PatientEPICC, patientEsc.df$Patient)]
burden.df$patEscape <- factor(burden.df$patEscape, levels=c('No','Epigenetic','Partial','Yes'))

burden.sub <- subset(burden.df, !is.na(imm_dnds) & imm_dnds<Inf)
ggplot(burden.sub, aes(x=PropBurden, y=imm_dnds, colour=patEscape)) +
  geom_hline(yintercept = 1, colour='grey50', linetype='dashed') +
  geom_point(size=2, aes(shape=SampleType)) +
  scale_shape_manual(values=c(18,8,15)) +
  scale_colour_manual(values=c('grey40','yellowgreen', 'goldenrod','brown')) +
  geom_linerange(aes(ymin=imm_lowci,ymax=imm_highci, group=Sample), alpha=0.5) +
  theme_mypub() + stat_cor() +
  coord_cartesian(ylim=c(0,5)) +
  labs(x='Proportional neoantigen burden',y='Immune dNdS', colour='Cancer\nescape') +
  theme(strip.background = element_blank())
