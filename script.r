#requirements-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(grid)
library(ggplot2)
library(lemon) # for fig2
library(venn) # for fig1
library(plyr)
library(reshape2)
library(RColorBrewer)
library(ggsankey) # for fig3A
library(tidyverse) # for fig3A
library(dplyr) # for fig3A

#workdir-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#git clone and cd ALS/  
setwd("data/")

source('../EWCE/bootstrap_enrichment_test.r')
source('../EWCE/cell_list_dist.r')
source('../EWCE/generate_controlled_bootstrap_geneset.r')
source('../EWCE/get_summed_proportions.r')

ALS_P<-c('ALS2','ANG','ANXA11','ATXN2','C9orf72','CCNF','CHCHD10','CHMP2B','CYLD','DCTN1','ERBB4','FIG4','FUS','HNRNPA1','KIF5A','MATR3','NEFH','NEK1','OPTN','PFN1','PRPH','SETX','SIGMAR1','SOD1','SPG11','SQSTM1','TARDBP','TBK1','TUBA4A','UBQLN2','VAPB','VCP')
ALS_S<-c('ACSL5','ALCAM','ATXN1','ATXN3','C21orf2','C9orf72','CABIN1','CAMK1G','CNTN4','DPP6','FGD4','FGGY','FNBP1','GGNBP2','GPR133','GPX3','INPP5B','ITGA9','ITPR2','KIF5A','MOBKL2B','MOBP','MYO18B','NEK1','OPCML','PFKP','RAPGEF5','SARM1','SCFD1','SOD1','SUSD2','TBK1','TNIP1','TYW3','UNC13A','IFNK','CRYZ','ERGIC1','PTPRN2','COG3','SPATA2','SCL9A8','CFAP410','G2E3','CLCN3','ZDHHC6','B4GALNT1','MOB3B')
HMN<-c('BSCL2','DCTN1','DNAJB2','FBXO38','GARS','HSPB1','HSPB3','HSPB8','IGHMBP2','PLEKHG5','REEP1','SIGMAR1','SLC5A7','TRPV4')
SA<-c('AFG3L2','ATXN1','ATXN10','ATXN2','ATXN3','ATXN7','ATXN8','ATXN8OS','BEAN','CACNA1A','CACNA1G','CCDC88C','DAB1','EEF2','ELOVL4','ELOVL5','FAT2','FGF14','GRM1','ITPR1','KCNC3','KCND3','MME','NOP56','PDYN','PLD3','PPP2R2B','PRKCG','PUM1','SPTBN2','STUB1','TBP','TGM6','TMEM240','TRPC3','TTBK2')
SMA<-c('AR','ASAH1','ASCC1','ATP7A','BICD2','CHCHD10','DNAJB2','DYNC1H1','GARS1','IGHMBP2','PLEKHG5','SIGMAR1','SMN1','SMN2','TRIP4','TRPV4','UBA1','VAPBC')
SPG<-c('ACP33','ALDH18A1','AMPD2','AP4B1','AP4E1','AP4M1','AP4S1','AP5Z1','ARL6IP1','ATL1','ATP13A2','B4GALNT1','BSCL2','C19orf12','CAPN1','CPT1C','CYP2U1','CYP7B1','DDHD1','DDHD2','DSTYK','ENTPD1','ERLIN1','ERLIN2','FA2H','FARS2','GBA2','GJC2','HPDL','HSPD1','IBA57','KIF1A','KIF5A','L1CAM','MAG','MTRFR','NIPA1','NT5C2','PCYT2','PGN','PLP1','PNPLA6','REEP1','RTN2','SELENOI','SLC33A1','SPAST','SPG11','SPG20','TECPR2','TFG','UBAP1','UCHL1','VPS37A','WASHC5','ZFYVE26','ZFYVE27')

GeneList<-list('ALS-pathogenicity'=ALS_P,'ALS-susceptibility'=ALS_S,'HMN'=HMN,'SA'=SA,'SMA'=SMA,'SPG'=SPG)
ALSGeneList<-list(ALS_P,ALS_S)
NDsGeneList<-list(unique(ALS_P,ALS_S),HMN,SA,SMA,SPG)

#Fig1-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf('Fig1A.pdf',width=8,height=8)
venn(ALSGeneList,snames ='ALS-Pathogenicity,ALS-Susceptibility',zcolor = brewer.pal(5,"Set2")[c(1,1)],opacity =0.5,ilcs=2,sncs=2,ggplot=F,box=F)
dev.off()
pdf('Fig1B.pdf',width=8,height=8)
venn(NDsGeneList,snames ='ALS,HMN,SA,SMA,SPG',zcolor = brewer.pal(5,"Set2"),opacity =0.5,ilcs=2,sncs=2,ggplot=F,box=F)
dev.off()

#Fig2-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('M2H.rda')
m2h<-unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])
bg.mouse<-unique(m2h$MGI.symbol)

ctdFile<-data.frame(species=c('mouse',
                              'mouse',
                              'human',
                              'mouse',
                              'mouse',
                              'human',
                              'mouse'),
                    region=c('mouse brain',
                             'mouse spinal cord',
                             'human muscle',
                             'mouse spinal cord (Milich dataset)',
                             'mouse spinal cord (Sathyamurthy dataset)',
                             'human spinal cord',
                             'mouse brain'),
                    ctdName=c('discovery_mouse_brain_celltype.rda',
                              'discovery_mouse_spinal_cord_celltype.rda',
                              'discovery_human_muscle_celltype.rda',
                              'replication_Milich.rda',
                              'replication_Sathyamurthy.rda',
                              'validation_human_spinal_cord_celltype.rda',
                              'replication_mouse_brain_celltype.rda'))

odf<-data.frame(celltype=NA,FDR=NA,gene=NA,source=NA,p=NA)
for(i in 1:6){
  x<-GeneList[[i]]
  for(j in 1:3){
    load(ctdFile[j,3])

    if(ctdFile[j,1]=='mouse'){
      bg<-bg.mouse
      homology<-unique(m2h[m2h$HGNC.symbol %in% x,"MGI.symbol"])
      hits<-homology[homology %in% attr(ctd[[1]]$specificity,'dimnames')[[1]]]
    }else{
      bg<-attr(ctd[[1]]$specificity,'dimnames')[[1]]
      hits<-x[x %in% bg]
    }

    set.seed(2022)
    rdf<-bootstrap_enrichment_test(sct_data=ctd,
                                   hits=hits,
                                   bg=bg,
                                   reps=10000,
                                   genelistSpecies=ctdFile[j,1],
                                   sctSpecies=ctdFile[j,1],
                                   annotLevel=1)
    rdf$results$celltype<-row.names(rdf$results)
    rdf$results$FDR<-p.adjust(rdf$results$p,method='fdr')
    rdf$results$gene<-names(GeneList)[i]
    rdf$results$source<-ctdFile[j,2]
    rdf<-rdf$results[,c('celltype','FDR','gene','source','p')]
    odf<-rbind(odf,rdf)
  }
}
odf<-odf[-1,]
row.names(odf)<-1:nrow(odf)
odf$gene<-factor(odf$gene,levels=c("ALS-susceptibility","ALS-pathogenicity","HMN","SA","SMA","SPG"),order=T)

fig2A<-ggplot(odf[odf$source=='mouse brain',],aes(celltype,-log10(FDR)))+
  geom_histogram(stat = 'identity',position = "dodge",aes(fill = gene))+
  geom_hline(yintercept = -log10(0.05),colour='red')+
  xlab('')+
  ylab(expression(paste("-",log[10]," FDR P-value")))+
  ylim(c(0,2))+
  ggtitle('A (mouse brain)')+
  theme_classic()+
  scale_fill_brewer(palette = 'Set2')+
  coord_flip()+
  theme(legend.title = element_blank(),
        axis.title = element_text(colour='black'),
        axis.text = element_text(colour='black'))

fig2B<-ggplot(odf[odf$source=='mouse spinal cord',],aes(celltype,-log10(FDR)))+
  geom_histogram(stat = 'identity',position = "dodge",aes(fill = gene))+
  geom_hline(yintercept = -log10(0.05),colour='red')+
  xlab('')+
  ylab(expression(paste("-",log[10]," FDR P-value")))+
  ylim(c(0,2.5))+
  ggtitle('B (mouse spinal cord)')+
  theme_classic()+
  scale_fill_brewer(palette = 'Set2')+
  coord_flip()+
  theme(legend.title = element_blank(),
        axis.title = element_text(colour='black'),
        axis.text = element_text(colour='black'))

fig2C<-ggplot(odf[odf$source=='human muscle',],aes(celltype,-log10(FDR)))+
  geom_histogram(stat = 'identity',position = "dodge",aes(fill = gene))+
  geom_hline(yintercept = -log10(0.05),colour='red')+
  xlab('')+
  ylab(expression(paste("-",log[10]," FDR P-value")))+
  ylim(c(0,1.75))+
  ggtitle('C (human muscle)')+
  theme_classic()+
  scale_fill_brewer(palette = 'Set2')+
  coord_flip()+
  theme(legend.title = element_blank(),
        axis.title = element_text(colour='black'),
        axis.text = element_text(colour='black'))

pdf('Fig2.pdf',width=15,height=7)
grid_arrange_shared_legend(fig2A,fig2B,fig2C,ncol=3,position='top')
dev.off()

#Fig3-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fig3A
sdf<-odf[odf$FDR<0.07,]
ldf<-sdf %>% make_long(gene,celltype,source)
fig3A<-ggplot(ldf, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
  geom_sankey(flow.alpha = .6,node.color = "gray30") +
  geom_sankey_label(size = 3, color = "black") +
  theme_sankey(base_size = 18) +
  labs(x = NULL)+ggtitle("A")+
  theme(legend.position = "none",
        plot.title = element_text(size=20,hjust=.15),
        axis.text = element_blank(),
        axis.title = element_text(size=18,colour="black"),
        strip.text = element_text(size=15,colour="black"))


#supplement table 6-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
odfrep<-data.frame(celltype=NA,p=NA,gene=NA,source=NA)
for(i in 1:6){
  x<-GeneList[[i]]
  for(j in 4:7){
    load(ctdFile[j,3])
    
    if(ctdFile[j,1]=='mouse'){
      bg<-bg.mouse
      homology<-unique(m2h[m2h$HGNC.symbol %in% x,"MGI.symbol"])
      hits<-homology[homology %in% attr(ctd[[1]]$specificity,'dimnames')[[1]]]
    }else{
      bg<-attr(ctd[[1]]$specificity,'dimnames')[[1]]
      hits<-x[x %in% bg]
    }
    
    set.seed(2022)
    rdf<-bootstrap_enrichment_test(sct_data=ctd,
                                   hits=hits,
                                   bg=bg,
                                   reps=10000,
                                   genelistSpecies=ctdFile[j,1],
                                   sctSpecies=ctdFile[j,1],
                                   annotLevel=1)
    rdf$results$celltype<-row.names(rdf$results)
    rdf$results$gene<-names(GeneList)[i]
    rdf$results$source<-ctdFile[j,2]
    
    rdf<-rdf$results[,c('celltype','gene','source','p')]
    odfrep<-rbind(odfrep,rdf)
  }
}
odfrep<-odfrep[-1,]
row.names(odfrep)<-1:nrow(odfrep)
mndf<-odfrep[odfrep$celltype %in% c('Alpha.motor.neurons','Gamma.motor.neurons','alpha_mn','gamma_mn'),]

mndf$celltype<-as.character(factor(mndf$celltype,
                                    levels=c('Alpha.motor.neurons','alpha_mn','Gamma.motor.neurons','gamma_mn'),
                                    labels=c('Alpha motor neurons','Alpha motor neurons','Gamma motor neurons','Gamma motor neurons')))

pdf<-odf[(odf$gene=='ALS-pathogenicity' | odf$gene=='ALS-susceptibility') & odf$source=='mouse spinal cord' & (odf$celltype=='Alpha motor neurons' | odf$celltype=='Gamma motor neurons'),-2]

mndf<-rbind(pdf,mndf)
write.table(mndf,file='TableS6_1.xls',sep='\t',quote=F,row.names=F,col.names=T)
meta<-function(p){z<-qnorm(p/2);w<-rep(1,length(p));zscore<-sum(z*w)/(sum(w^2))^0.5;metapvalue<-2*(1-pnorm(abs(-zscore)));return(metapvalue)}
metadf<-ddply(mndf,.(celltype,gene),summarise,metapvalue=meta(p))
write.table(metadf,file='TableS6_2.xls',sep='\t',quote=F,row.names=F,col.names=T)

#Fig4-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fig4A-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df<-read.table("brain_and_SP.sd.txt",sep="\t",row.names = 1,head=T)
fig4A<-ggplot(df,aes(brain))+
  geom_histogram(bins=100,fill= brewer.pal(5,"Set2")[5])+
  xlab('Strictness')+ylab('Counts')+ylim(c(-150,1300))+
  theme_classic()+ggtitle('A')+
  theme(axis.title = element_text(colour='black'),
        axis.text = element_text(colour='black'))+
  annotate(geom = 'point',x = c(1/47,1/1.177),y = c(-25,-25),shape=24,fill="red",colour='red',size=3)+
  annotate(geom = 'text',x = c(1/47,1/1.177),y = c(0,0),label = c('Fam187b','Ank1'), fontface = 'italic',vjust=2,hjust=0.2)

#funtion-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
strictness_test<-function(sddf=NA,geneList=NA){
  results<-apply(sddf,2,function(x){
    bg.genes<-row.names(sddf[!is.na(x),])
    x<-x[!is.na(x)]

    hits<-geneList[geneList %in% bg.genes]
    sizes<-length(hits)

    target.mean<-mean(x[bg.genes %in% hits],na.rm=T)
    bg.distribution<-sapply(rep(1,10000),function(y){
      mean(sample(x=x,size=sizes),na.rm=T)
    })
    sd_value<-sd(bg.distribution,na.rm=T)
    mean_value<-mean(bg.distribution,na.rm=T)
    p_value<-1-pnorm(target.mean,mean=mean_value,sd=sd_value)
    p_value
  })
  return(results)
}

#Fig4B-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
HI<-c('ANK2','APC','BCL11B','CD2AP','COMT','CSF2RA','DMPK','DYRK1A','EGR1','EHMT1','FGF10','FOXP1','FOXP2','GCH1','GHRL','GTF2I','HOXD13','IGF1','KCNQ2','LMX1B','MAPT','MC4R','NF1','NKX2-5','NLGN4X','NLRP3','NPAS3','NSD1','PARK2','PAX6','PIK3R1','PRODH','PTEN','RELN','SATB2','SCN1A','SEMA5A','SHFM1','SHMT1','SNCA','SPR','ST7','TBX1','TCF4','TGFB1','TPM1','TSC1','TSC2','WWOX')
LOFT<-c('A2M','ABCA10','ABCA8','ABCC11','ABCC12','ABHD12B','ABHD14B','ACSBG2','ACSM2A','ACSM2B','ACSM3','ADAM2','ADPRHL1','ADSSL1','AHNAK2','AKAP3','ALDH1B1','ANKRD30A','ANKRD35','ANO5','AP1G2','APIP','APOBEC3A','APOBEC3B','ASB15','ASPSCR1','ATP10B','ATP11A','ATP12A','ATP2C2','BPHL','BPIFA3','BTN3A3','BTNL2','BTNL8','BTNL9','C10orf53','C11orf40','C1orf127','C2orf40','C3orf14','C4orf46','C4orf50','C6','CABYR','CAPN9','CARD6','CCDC121','CCDC13','CCDC60','CCDC66','CD180','CD36','CD96','CDH19','CDK11A','CDKL2','CELA1','CEP72','CES1','CES5A','CFHR1','CFHR2','CFHR3','CHD1L','CHIT1','CHPF2','CLCN1','CLYBL','CNKSR1','COL16A1','COL6A2','COL6A5','CPXM2','CROT','CRYGN','CRYZ','CSH1','CTSE','CYP2A6','CYP2C8','CYP2D6','CYP2F1','CYP3A5','CYP4B1','DCDC2B','DCHS2','DDX60','DEFB126','DHDH','DMBT1','DNAH7','DQX1','DUOX2','ECT2L','EFCAB13','EFCAB3','EFCAB5','EFCAB6','ENOSF1','ENPEP','EPPK1','EPX','ERAP1','ERV3-1','EXO5','FAM129C','FAM151A','FAM187B','FAM45A','FAM81B','FCGBP','FCGR2A','FCN3','FLG','FLG2','FMO2','FRG2B','FUT2','FUT6','GADL1','GBGT1','GBP3','GBP4','GCFC2','GH2','GJB4','GLB1L2','GMPR','GOLGA8S','GP6','GPATCH2L','GRIN3B','GRK7','GYPB','HELB','HK3','HLA-B','HLA-DPA1','HPSE','HRG','HRNR','IDI2','IFIH1','IFNK','IL17RC','IL3RA','IQCH','ITIH1','KIAA0753','KIAA1257','KIAA1586','KIR3DL1','KLK14','KLK3','KRT4','KRT77','KRT83','LMF2','LMO7','LPA','LRRC39','LRTM1','MANEA','MAP3K4','MAZ','MCF2L','MCOLN3','MFSD9','MGAM','MLANA','MMP10','MOGAT1','MOK','MOXD1','MS4A6A','MST1','MUC17','MUC6','MUTYH','MYBBP1A','MYH1','MYH13','MYH8','MYO1A','MYOC','MYOF','NAALAD2','NBPF14','NBPF15','NEIL1','NLRP13','NLRP9','NOP16','NUDT8','OARD1','OBSCN','OCEL1','OR8S1','PAPLN','PDE11A','PDIA2','PGPEP1L','PHRF1','PKD1L2','PKHD1L1','PLA2G2C','PLA2G4D','PLA2R1','PLEKHG7','PLIN4','PNLIPRP3','POLM','POTEH','PPEF2','PPL','PPP1R3A','PRAMEF2','PRB1','PRB2','PRB4','PSG1','PSG11','PSG4','PSG9','PTCHD3','PTGDR','PXDNL','PZP','RAI1','RERGL','RETSAT','RFPL1','RGPD4','RGS11','RHD','RNF32','ROPN1B','RP1L1','RPTN','RTKN2','RTP1','SAMD11','SEMG2','SERHL2','SERPINA10','SERPINA9','SERPINB3','SFI1','SIGLEC1','SIGLEC5','SLC17A9','SLC22A10','SLC22A14','SLC22A25','SLC26A10','SLC5A4','SLCO1B1','SLFN13','SPATA31A6','SPATA4','SPATC1','SPNS3','SULT1A2','SULT1C4','SYNM','SYTL2','TAF6','TCF3','TCHHL1','TEKT3','TGM4','THBS4','THEM5','TIGD6','TLR10','TLR5','TMC2','TMEM82','TMIE','TMPRSS7','TNN','TRIM22','TRIM45','TRIM48','TRIM59','TRMT10B','TRMT2A','TTC38','TTN','UGT2B10','UGT2B17','UGT2B28','UMODL1','UNC93A','UPB1','UPK3A','UPP2','USP45','USP6','VILL','VWA3B','VWA7','WDR27','WDR90','XIRP1','XRRA1','ZAN','ZNF223','ZNF229','ZNF257','ZNF30','ZNF343','ZNF396','ZNF417','ZNF486','ZNF528','ZNF544','ZNF587','ZNF599','ZNF611','ZNF790','ZNF83','ZNF831','ZNF844','ZNF846','ZNF860','ZNF878','ZNF92','ZRANB3')
NNN<-c('KCNJ10','KCNJ11','KCNJ2','KCNJ5','KCNQ1','KCNQ1OT1','KCNQ4','KCTD7','KIAA0196','KIF1B','KIF5A','KIF7','KIT','KLF1','KLF11','KLF6','KLHDC8B','KLK4','KRAS','KRT1','KRT10','KRT13','KRT16','KRT17','KRT4','KRT6A','KRT6B','KRT81','KRT83','KRT86','L2HGDH','LAMA2','LAMP2','LARGE','LBR','LCA5','LCAT','LDB3','LEMD3','LEPRE1','LFNG','LHFPL5','LIG4','LIPI','LITAF','LMBR1','LMBRD1','LOXHD1','LOXL1','LPL','LPP','LRAT','LRP5','LRRC8A','LRTOMT','LTBP2','LTBP3','LYST','LYZ','LZTS1','MAMLD1','MAN2B1','MANBA','MAP3K1','MAP3K8','MAPK10','MAPK8IP1','MARVELD2','MC1R','MCCC2','MCOLN1','MED25','MEFV','MEN1','MESP2','MFN2','MFSD8','MGAT2','MINPP1','MITF','MKKS','MKS1','MLH1','MLH3','MLLT10','MMAA','MMAB','MMACHC','MMADHC','MMP13','MMP20','MN1','MOGS','MPDU1','MPI','MPL','MPLKIP','MPZ','MS4A1','MSH2','MSH6','MSX1','MSX2','MTMR2','MTTP','MUC5B','MUTYH','MVK','MXI1','MYBPC3','MYC','MYH11','MYH14','MYH3','MYH6','MYH7','MYH9','MYL2','MYL3','MYLK2','MYO15A','MYO1A','MYO3A','MYO6','MYO7A','MYOC','MYOT','NBN','NCF1','NCF2','NCOA4','NCSTN','NDP','NDRG1','NDUFA10','NDUFA12','NDUFA13','NDUFA2','NDUFA9','NDUFS3','NDUFS4','NDUFS7','NDUFS8','NEFH','NEFL','NEUROD1','NEXN','NF2','NHLRC1','NIPA1','NIPBL','NKX2-6','NLRP7','NME1','NOBOX','NODAL','NOG','NOTCH1','NOTCH2','NPC2','NPHP1','NPHP3','NPHP4','NPM1','NPPA','NPR2','NR0B1','NR0B2','NR2E3','NR4A3','NR5A1','NRL','NT5C3A','NUP214','NUP62','NYX','OAS1','OAT','OCA2','OCLN','OCRL','ODC1','OFD1','OGG1','OPCML','OPN1LW','OPN1MW','OPTN','OSMR','OSTM1','OTC','OTOA','OTOF','OTX2','PABPN1','PADI4','PALB2','PANK2','PARK7','PAX3','PAX4','PAX7','PAX8','PAX9','PC','PCCA',
       'PCCB','PCDH15','PCM1','PDE6B','PDE6C','PDE6G','PDE8B','PDGFRA','PDGFRL','PDHA1','PDX1','PDYN','PDZD7','PEX1','PEX10','PEX13','PEX14','PEX19','PEX2','PEX26','PEX3','PEX5','PFKM','PHEX','PHF6','PHGDH','PHKA2','PHKB','PHKG2','PHYH','PINK1','PITPNM3','PITX2','PKD1','PKD2','PKHD1','PLA2G2A','PLA2G6','PLAG1','PLCB1','PLOD1','PLP1','PMM2','PMP22','PNKP','PNP','PNPLA6','POLG','POLH','POLR1C','POLR1D','POMT1','POMT2','PORCN','POU3F4','POU4F3','PPARG','PPARGC1B','PPIB','PPOX','PPP1R3A','PPP2R1B','PPT1','PRCC','PRCD','PRF1','PRICKLE2','PRKAG2','PRKAR1A','PRKCG','PRKCH','PRNP','PROK2','PROKR2','PROM1','PRPF3','PRPF31','PRPF8','PRPH','PRPH2','PRPS1','PRX','PSAT1','PSENEN','PSPH','PTCH1','PTCH2','PTGIS','PTHLH','PTPN1','PTPN11','PTPN22','PTPRJ','PTPRQ','PTRF','PYGM','RAB7A','RAD51','RAD54L','RAX','RAX2','RB1','RB1CC1','RBM20','RBM8A','RD3','RDH12','RDH5','RDX','RECQL4','REEP1','REN','RET','RFT1','RFX5','RFXANK','RFXAP','RGR','RHAG','RHO','RLBP1','RMRP','RNASEH2A','RNASEH2B','RNASEH2C','RNF139','RNF6','ROBO2','ROR2','RP1','RP1L1','RP2','RP9','RPE65','RPGR','RPGRIP1','RPL11','RPL35A','RPL5','RPS10','RPS17','RPS19','RPS24','RPS7','RSPH4A','RSPH9','RSPO4','RUNX1','RUNX2','RXFP2','RYR1','RYR2','SAG','SAMHD1','SBF2','SCN1B','SCN3B','SCN4A','SCN4B','SCN5A','SCN9A','SCNN1B','SCNN1G','SCO2','SDC3','SDHA','SDHB','SDHC','SEC23B','SEMA4A','SEPN1','SERPINA1','SERPINB6','SERPING1','SETX','SFTPA1','SFTPA2','SGCA','SGCB','SGCD','SGCG','SGSH','SH2B3','SH3BP2','SH3TC2','SHH','SIM1','SIX1','SIX5','SLC11A2','SLC12A1','SLC12A3','SLC16A1','SLC16A2','SLC17A5','SLC17A8','SLC19A3','SLC22A18','SLC22A4','SLC22A5','SLC25A15','SLC25A22','SLC25A38','SLC25A4','SLC26A2','SLC26A4','SLC2A10','SLC2A4','SLC33A1','SLC34A1','SLC34A2','SLC35A1','SLC35C1','SLC37A4','SLC39A4','SLC3A1','SLC40A1','SLC45A2','SLC4A11','SLC52A3','SLC5A2','SLC5A5','SLC6A19','SLC7A9','SMAD9','SMARCB1','SMN1','SMN2','SMPD1','SNAI2','SNCB','SNRNP200','SNTA1','SOD1','SOS1','SOX10','SOX2','SOX9','SPATA7','SPG11','SPG20','SPG7','SPINK5','SPTA1','SPTAN1','SPTBN2','SPTLC1','SPTLC2','SQSTM1','SRD5A3','SRY','STAT3','STIL','STK11','STOX1','STRA6','STRC','STX11','STX16','STXBP2','SUFU','SUMF1','SUMO1','SURF1','TAP1','TAP2','TAPBP','TARDBP','TAT','TAZ','TBP','TBX20','TCAP','TCIRG1','TCOF1','TDP1','TECTA','TERC',
       'TERT','TF','TFAP2A','TFR2','TGFB3','TGFBR1','TGFBR2','TGIF1','TGM1','THPO','THRB','TLL1','TLR2','TLR4','TMC1','TMC6','TMC8','TMEM43','TMEM67','TMIE','TMPRSS3','TMPRSS6','TNFRSF11A','TNFRSF13B','TNFRSF13C','TNFSF11','TNNC1','TNNI2','TNNI3','TNNT2','TNNT3','TOPORS','TP63','TPM2','TPO','TPP1','TPRN','TREX1','TRIM24','TRIM37','TRIOBP','TRIP11','TRPC6','TRPM1','TRPS1','TRPV4','TSG101','TTBK2','TTC8','TTN','TTR','TUBA1A','TULP1','TWIST1','TYRP1','UBR1','UCHL1','UCP1','UCP3','UGT1A1','UMOD','UNC13D','UPK3A','UROD','UROS','USH1C','USH1G','USH2A','USP9Y','VAPB','VCL','VHL','VSX2','VWF','WAS','WDPCP','WDR35','WDR36','WDR72','WHSC1L1','WNK1','WRN','WT1','XPA','XPC','XRCC3','XYLT1','XYLT2','ZEB1','ZFHX3','ZFYVE27','ZIC2','ZIC3','ZNF469','ZNF513','ZNF592')
HI.mouse<-unique(m2h[m2h$HGNC.symbol %in% HI,"MGI.symbol"])
LOFT.mouse<-unique(m2h[m2h$HGNC.symbol %in% LOFT,"MGI.symbol"])
NNN.mouse<-unique(m2h[m2h$HGNC.symbol %in% NNN,"MGI.symbol"])

result1<-strictness_test(df,HI.mouse)
result2<-strictness_test(df,LOFT.mouse)
result3<-strictness_test(df,NNN.mouse)
sdf<-rbind(as.data.frame(t(result1)),
           as.data.frame(t(result2)),
           as.data.frame(t(result3)))

sdf$gene<-c('HI','LOFT','NNN')
mdf<-melt(sdf,id='gene')
mdf$value<-p.adjust(mdf$value,method='fdr')
fig4B<-ggplot(mdf,aes(gene,-log10(value)))+
  geom_histogram(width=0.75,stat = 'identity',position = "dodge",aes(fill = gene))+
  geom_hline(yintercept = -log10(0.05),colour='red')+
  xlab('')+ylab(expression(paste('-',log[10],"FDR P-value")))+
  facet_wrap(~variable,nrow=1)+scale_fill_brewer(palette = 'Set2')+
  theme_bw()+ggtitle('B')+coord_flip()+
  theme(panel.grid  = element_blank(),
        axis.title = element_text(colour='black'),
        axis.text = element_text(colour='black'),
        legend.title = element_blank(),
        legend.position ='none',
        legend.background  = element_rect(colour = 'black'))

#Fig4C-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ALS_S.mouse<-unique(m2h[m2h$HGNC.symbol %in% ALS_S,"MGI.symbol"])
ALS_P.mouse<-unique(m2h[m2h$HGNC.symbol %in% ALS_P,"MGI.symbol"])
HMN.mouse<-unique(m2h[m2h$HGNC.symbol %in% HMN,"MGI.symbol"])
SA.mouse<-unique(m2h[m2h$HGNC.symbol %in% SA,"MGI.symbol"])
SMA.mouse<-unique(m2h[m2h$HGNC.symbol %in% SMA,"MGI.symbol"])
SPG.mouse<-unique(m2h[m2h$HGNC.symbol %in% SPG,"MGI.symbol"])

df<-read.table("SP.sd.txt",sep="\t",row.names = 1,head=T)

result1<-strictness_test(df,ALS_S.mouse)
result2<-strictness_test(df,ALS_P.mouse)
result3<-strictness_test(df,HMN.mouse)
result4<-strictness_test(df,SA.mouse)
result5<-strictness_test(df,SMA.mouse)
result6<-strictness_test(df,SPG.mouse)

result1<-p.adjust(result1,method='fdr')
result2<-p.adjust(result2,method='fdr')
result3<-p.adjust(result3,method='fdr')
result4<-p.adjust(result4,method='fdr')
result5<-p.adjust(result5,method='fdr')
result6<-p.adjust(result6,method='fdr')

sdf<-rbind(as.data.frame(t(result1)),
           as.data.frame(t(result2)),
           as.data.frame(t(result3)),
           as.data.frame(t(result4)),
           as.data.frame(t(result5)),
           as.data.frame(t(result6)))

sdf$gene<-c("ALS-susceptibility","ALS-pathogenicity","HMN","SA","SMA","SPG")
mdf<-melt(sdf,id='gene')
mdf$gene<-factor(mdf$gene,levels=unique(mdf$gene),order=T)

fig4C<-ggplot(mdf,aes(variable,-log10(value)))+
  geom_histogram(stat = 'identity',position = "dodge",aes(fill = gene))+
  geom_hline(yintercept = -log10(0.05),colour='red')+
  xlab('')+scale_y_continuous(expression(paste('-',log[10],"FDR P-value")))+
  facet_grid(~gene,scale='free')+ggtitle('C')+
  theme_bw()+coord_flip()+scale_fill_brewer(palette = 'Set2')+
  theme(panel.grid  = element_blank(),
        axis.title = element_text(colour='black'),
        axis.text = element_text(colour='black'),
        legend.position = 'none')

#Fig4D-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ALS_P_lof<-c('ALS2','ATXN2','C9orf72','CYLD','ERBB4','FUS','NEK1','OPTN','PFN1','PRPH','SIGMAR1','TARDBP','TBK1')
ALS_P_gof<-c('FIG4','NEFH','SOD1','VAPB')
ALS_P_unknown<-c('ANG','ANXA11','CCNF','CHCHD10','CHMP2B','DCTN1','HNRNPA1','KIF5A','MATR3','SETX','SPG11','SQSTM1','TUBA4A','UBQLN2','VCP')
ALS_P_lof.mouse<-unique(m2h[m2h$HGNC.symbol %in% ALS_P_lof,"MGI.symbol"])
ALS_P_gof.mouse<-unique(m2h[m2h$HGNC.symbol %in% ALS_P_gof,"MGI.symbol"])
ALS_P_unknown.mouse<-unique(m2h[m2h$HGNC.symbol %in% ALS_P_unknown,"MGI.symbol"])

result1<-strictness_test(df,ALS_P_lof.mouse)
result2<-strictness_test(df,ALS_P_gof.mouse)
result3<-strictness_test(df,ALS_P_unknown.mouse)

result1<-p.adjust(result1,method='fdr')
result2<-p.adjust(result2,method='fdr')
result3<-p.adjust(result3,method='fdr')

sdf<-rbind(as.data.frame(t(result1)),
           as.data.frame(t(result2)),
           as.data.frame(t(result3)))

sdf$gene<-c("LoF","GoF","unknown")
mdf<-melt(sdf,id='gene')
mdf$gene<-factor(mdf$gene,levels=unique(mdf$gene),order=T)
mdf$P<-'ALS-pathogenicity'

fig4D<-ggplot(mdf,aes(variable,-log10(value)))+
  geom_histogram(stat = 'identity',position = "dodge",aes(fill = gene))+
  geom_hline(yintercept = -log10(0.05),colour='red')+
  xlab('')+scale_y_continuous(expression(paste('-',log[10],"FDR P-value")))+
  facet_grid(~P+gene,scale='free')+ggtitle('D')+
  theme_bw()+coord_flip()+scale_fill_manual(values = c("#FC8D62","#FC8D6299","#FC8D6280"))+
  #scale_fill_brewer(palette = 'Set2')+
  theme(panel.grid  = element_blank(),
        axis.title = element_text(colour='black'),
        axis.text = element_text(colour='black'),
        legend.position = 'none')#,

#Fig4E-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tmp1<-df[row.names(df) %in% ALS_P_lof.mouse,'Gamma.motor.neurons',drop=F]
tmp2<-df[row.names(df) %in% ALS_P_gof.mouse,'Gamma.motor.neurons',drop=F]
tmp3<-df[row.names(df) %in% ALS_P_unknown.mouse,'Gamma.motor.neurons',drop=F]
tmp1$gene<-'LoF'
tmp2$gene<-'GoF'
tmp3$gene<-'unknown'
gdf<-rbind(tmp1,tmp2,tmp3)

fig4E<-ggplot(gdf,aes(gene,Gamma.motor.neurons))+
  geom_boxplot(aes(fill=gene))+
  annotate(geom = 'point',x = 'GoF',y =  gdf[row.names(gdf)=='Sod1','Gamma.motor.neurons'],colour='red',size=3)+
  ylab('Strictness in gamma MNs')+
  xlab('ALS-pathogenicity genes')+
  ggtitle('E')+
  theme_classic()+
  coord_flip()+
  scale_fill_manual(values = c("#FC8D62","#FC8D6299","#FC8D6250"))+
  #scale_fill_brewer(palette = 'Set2')+
  theme(legend.position = 'none',axis.title = element_text(colour='black'),
        axis.text = element_text(colour='black'))


pdf('Fig4.pdf',width=10,height=10)
grid.newpage()
pushViewport(viewport(layout=grid.layout(12,100)))
print(fig4A,vp=viewport(layout.pos.row=1:3,layout.pos.col=22:60))
print(fig4B,vp=viewport(layout.pos.row=1:3,layout.pos.col=61:100))
print(fig4C,vp=viewport(layout.pos.row=4:8,layout.pos.col=1:100))
print(fig4D,vp=viewport(layout.pos.row=9:12,layout.pos.col=1:65))
print(fig4E,vp=viewport(layout.pos.row=9:12,layout.pos.col=66:100))
dev.off()

#Fig3-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fig3B-C-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#function-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
circos_function<-function(r=1){
  x<-seq(-r,r,length=1000)
  y<-(r^2-x^2)^0.5
  return(data.frame(x=c(x,rev(x)),y=c(y,-y)))
}

circos_rect<-function(xita1,xita2,r,width){
  xita<-seq(xita1,xita2,length=100)
  arc_out_x<-cos(xita)*r
  arc_out_y<-sin(xita)*r
  arc_in_x<-cos(rev(xita))*(r-width)
  arc_in_y<-sin(rev(xita))*(r-width)
  
  r_out_in<-seq(r,r-width,length=10)
  
  width_out_in_x<-cos(xita2)*r_out_in
  width_out_in_y<-sin(xita2)*r_out_in
  
  width_in_out_x<-cos(xita1)*rev(r_out_in)
  width_in_out_y<-sin(xita1)*rev(r_out_in)
  
  return(data.frame(x=c(arc_out_x,width_out_in_x,arc_in_x,width_in_out_x),
                    y=c(arc_out_y,width_out_in_y,arc_in_y,width_in_out_y)))
}

point_function<-function(xita,r=1){
  x<-cos(xita)*r
  y<-sin(xita)*r
  return(data.frame(x=x,y=y))
}

cos_connect<-function(x1,x2,dis,r){
  if(dis<0.001*r){
    dis<-0.001*r
  }
  
  x<-seq(x1,x2,length=100)
  y<- 0.5*dis/r*cos(x*pi/abs(x1-x2))-dis
  return(data.frame(x=x,y=y))
}

change_coord<-function(xita,x,y){
  xita<-xita-1.5*pi
  xb=x*cos(xita)-y*sin(xita)
  yb=x*sin(xita)+y*cos(xita)
  return(data.frame(x=xb,y=yb))
}

p2p<-function(xita1,xita2,r){
  if(xita1<0){
    xita1<-xita1+2*pi
  }
  if(xita2<0){
    xita2<-xita2+2*pi
  }
  if(abs(xita1-xita2)>pi){
    if(xita1>xita2){
      xita1<-xita1-2*pi
    }else{
      xita2<-xita2-2*pi
    }
  }
  
  midpoint_r<-abs(cos((xita1-xita2)/2))*r
  
  chrod <- sin(abs((xita1-xita2)/2))*r
  xita<-(xita1+xita2)/2
  
  cdf<-cos_connect(-chrod,chrod,midpoint_r,r)
  bdf<-change_coord(xita,cdf[,1],cdf[,2])
  return(bdf)
}

arc_connect<-function(xita1,xita2,r){
  
  if(xita1<0){
    xita1<-xita1+2*pi
  }
  if(xita2<0){
    xita2<-xita2+2*pi
  }
  if(abs(xita1-xita2)>pi){
    if(xita1>xita2){
      xita1<-xita1-2*pi
    }else{
      xita2<-xita2-2*pi
    }
  }
  
  xita<-seq(xita2,xita1,length=100)
  x<-cos(xita)*r
  y<-sin(xita)*r
  return(data.frame(x=x,y=y))
}


xita2xy<-function(xita,r){
  x<-cos(xita)*r
  y<-sin(xita)*r
  return(data.frame(x=x,y=y))
}

connect_ploy<-function(df1,df2){
  p1_s<-df1[1,]
  p1_e<-df1[nrow(df1),]
  p2_s<-df2[1,]
  p2_e<-df2[nrow(df2),]
  d1<-sum((p1_s-p2_s)^2)^0.5
  d2<-sum((p1_s-p2_e)^2)^0.5
  d3<-sum((p1_e-p2_s)^2)^0.5
  d4<-sum((p1_e-p2_e)^2)^0.5
  
  if(sum(d1<=c(d2,d3,d4))==3){
    x<-c(rev(df1[,1]),df2[,1])
    y<-c(rev(df1[,2]),df2[,2])
  }else if(sum(d2<=c(d1,d3,d4))==3){
    x<-c(df2[,1],df1[,1])
    y<-c(df2[,2],df1[,2])
  }else if(sum(d3<=c(d2,d1,d4))==3){
    x<-c(df1[,1],df2[,1])
    y<-c(df1[,2],df2[,2])
  }else if(sum(d4<=c(d2,d3,d1))==3){
    x<-c(df2[,1],rev(df1[,1]))
    y<-c(df2[,2],rev(df1[,2]))
  }
  
  return(data.frame(x=x,y=y))
}

#data-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set0<-c(brewer.pal(8, 'Set2'),brewer.pal(5, 'Set1'))
load('discovery_mouse_spinal_cord_celltype.rda')
bg.genes<-attr(ctd[[1]]$specificity,'dimnames')[[1]]
ALS_P_exp<-data.frame(ctd[[1]]$specificity[bg.genes %in% ALS_P.mouse,])
ALS_S_exp<-data.frame(ctd[[1]]$specificity[bg.genes %in% ALS_S.mouse,])


#Fig3B-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ALS_P_exp$Gene<-row.names(ALS_P_exp)
tmpdf<-melt(ALS_P_exp,id='Gene')
g2g<-tmpdf[tmpdf[,3]>0.2,c(2,1,3)]
g2gRaw<-g2g
g2g<-g2g[,1:2]
names(g2g)<-c("go","gene")
g2g[,2]<-factor(g2g[,2])
g2g<-g2g[order(g2g[,2]),]
goNum<-ddply(g2g,.(go),c("nrow"))
goNum<-goNum[order(-goNum[,2]),]
tagL<-sum(goNum[,2])
tag<-unique(goNum[,1])
tagN<-length(tag)-1
breakUnit<-(1.49-0.51)*0.1*pi/tagN
arcUnit<-(1.49-0.51)*0.9*pi/tagL
start0<-0.51*pi
rdf<-data.frame()
tdf<-data.frame()
for(i in tag){
  startP<-start0
  endP<-startP+goNum[goNum[,1]==i,2]*arcUnit
  start0<-endP+breakUnit
  tmp<-circos_rect(startP,endP,1.1,0.09)
  tmp$group<-i
  rdf<-rbind(rdf,tmp)
  midP<-(startP+endP)/2
  tmp2<-data.frame(tag=i,xita=midP)
  tdf<-rbind(tdf,tmp2)
  marker<-seq(startP,endP,by=arcUnit)
  g2g$gene2goStart[g2g[,1]==i]<-rev(marker[-length(marker)])
  g2g$gene2goEnd[g2g[,1]==i]<-rev(marker[-1])
}

geneNum<-ddply(g2g,.(gene),c("nrow"))
#tagL<-sum(geneNum[,2])
tag<-unique(g2g[,2])
tagN<-length(tag)
breakUnit<-(1.49-0.51)*0.1*pi/(tagN-1)
arcUnit<-(1.49-0.51)*0.9*pi/tagN

start0<- -0.49*pi
rdf2<-data.frame()
g2g[,1]<-factor(g2g[,1],levels=goNum[,1],order=T)
g2g<-g2g[order(g2g[,1]),]

for(i in tag){
  startP<-start0
  endP<-startP+arcUnit
  start0<-endP+breakUnit
  tmp<-circos_rect(startP,endP,1.1,0.045)
  tmp$group<-i
  rdf2<-rbind(rdf2,tmp)
  midP<-(startP+endP)/2
  tmp2<-data.frame(tag=i,xita=midP)
  tdf<-rbind(tdf,tmp2)
  weightGene<-cumsum(g2gRaw[g2gRaw[,2]==i,3])
  #marker<-seq(startP,endP,length=geneNum[geneNum[,1]==i,2]+1)#by=arcUnit)
  marker<-startP+(endP-startP)*weightGene/weightGene[length(weightGene)]
  marker<-c(startP,marker)
  g2g$go2geneStart[g2g[,2]==i]<-rev(marker[-length(marker)])
  g2g$go2geneEnd[g2g[,2]==i]<-rev(marker[-1])
}

ndf_use<-data.frame()
for(i in seq(nrow(g2g))){
  #for(i in row.names(g2g[g2g[,2]=="Gene20",])){
  xita<-c(g2g[i,3],g2g[i,4],g2g[i,5],g2g[i,6])
  xita<-xita[order(xita)]
  ndf<-p2p(xita[1],xita[4],1)
  ndf1 <-arc_connect(xita[1],xita[2],1)
  ndf2<-p2p(xita[2],xita[3],1)
  ndf3 <-arc_connect(xita[3],xita[4],1)
  
  ndf_all<-connect_ploy(ndf,ndf1)
  ndf_all<-connect_ploy(ndf_all,ndf2)
  ndf_all<-connect_ploy(ndf_all,ndf3)
  ndf_all$connect<-paste0(g2g[i,1],"_",g2g[i,2])
  ndf_all$group<-g2g[i,1]
  
  ndf_use<-rbind(ndf_use,ndf_all)
}

cdf<-circos_function(1)

tdf[,c(3,4)]<-xita2xy(tdf[,2],1.105)

fig3B<-ggplot()+
  geom_polygon(data=rdf,aes(x,y,group=group,fill=group))+
  geom_polygon(data=rdf2,aes(x,y,group=group),fill="black")+
  geom_polygon(data=ndf_use,aes(x,y,group=connect,fill=group),alpha=0.5)+
  geom_text(data=tdf,aes(x,y,label=tag,angle=xita*180/pi),size=2.5,hjust=0)+
  xlim(c(-2,2))+ylim(c(-2,2))+ggtitle('B')+scale_fill_manual(values=set0)+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #panel.background = element_blank(),
        panel.background=element_rect(fill="transparent",colour=NA),
        plot.background=element_rect(fill="transparent",colour=NA))


#Fig3C-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ALS_S_exp$Gene<-row.names(ALS_S_exp)
tmpdf<-melt(ALS_S_exp,id='Gene')
g2g<-tmpdf[tmpdf[,3]>0.2,c(2,1,3)]
g2gRaw<-g2g
g2g<-g2g[,1:2]
names(g2g)<-c("go","gene")
g2g[,2]<-factor(g2g[,2])
g2g<-g2g[order(g2g[,2]),]
goNum<-ddply(g2g,.(go),c("nrow"))
goNum<-goNum[order(-goNum[,2]),]
tagL<-sum(goNum[,2])
tag<-unique(goNum[,1])
tagN<-length(tag)-1
breakUnit<-(1.49-0.51)*0.1*pi/tagN
arcUnit<-(1.49-0.51)*0.9*pi/tagL
start0<-0.51*pi
rdf<-data.frame()
tdf<-data.frame()
for(i in tag){
  startP<-start0
  endP<-startP+goNum[goNum[,1]==i,2]*arcUnit
  start0<-endP+breakUnit
  tmp<-circos_rect(startP,endP,1.1,0.09)
  tmp$group<-i
  rdf<-rbind(rdf,tmp)
  midP<-(startP+endP)/2
  tmp2<-data.frame(tag=i,xita=midP)
  tdf<-rbind(tdf,tmp2)
  marker<-seq(startP,endP,by=arcUnit)
  g2g$gene2goStart[g2g[,1]==i]<-rev(marker[-length(marker)])
  g2g$gene2goEnd[g2g[,1]==i]<-rev(marker[-1])
}

geneNum<-ddply(g2g,.(gene),c("nrow"))
#tagL<-sum(geneNum[,2])
tag<-unique(g2g[,2])
tagN<-length(tag)
breakUnit<-(1.49-0.51)*0.1*pi/(tagN-1)
arcUnit<-(1.49-0.51)*0.9*pi/tagN

start0<- -0.49*pi
rdf2<-data.frame()
g2g[,1]<-factor(g2g[,1],levels=goNum[,1],order=T)
g2g<-g2g[order(g2g[,1]),]

for(i in tag){
  startP<-start0
  endP<-startP+arcUnit
  start0<-endP+breakUnit
  tmp<-circos_rect(startP,endP,1.1,0.045)
  tmp$group<-i
  rdf2<-rbind(rdf2,tmp)
  midP<-(startP+endP)/2
  tmp2<-data.frame(tag=i,xita=midP)
  tdf<-rbind(tdf,tmp2)
  weightGene<-cumsum(g2gRaw[g2gRaw[,2]==i,3])
  #marker<-seq(startP,endP,length=geneNum[geneNum[,1]==i,2]+1)#by=arcUnit)
  marker<-startP+(endP-startP)*weightGene/weightGene[length(weightGene)]
  marker<-c(startP,marker)
  g2g$go2geneStart[g2g[,2]==i]<-rev(marker[-length(marker)])
  g2g$go2geneEnd[g2g[,2]==i]<-rev(marker[-1])
}

ndf_use<-data.frame()
for(i in seq(nrow(g2g))){
  #for(i in row.names(g2g[g2g[,2]=="Gene20",])){
  xita<-c(g2g[i,3],g2g[i,4],g2g[i,5],g2g[i,6])
  xita<-xita[order(xita)]
  ndf<-p2p(xita[1],xita[4],1)
  ndf1 <-arc_connect(xita[1],xita[2],1)
  ndf2<-p2p(xita[2],xita[3],1)
  ndf3 <-arc_connect(xita[3],xita[4],1)
  
  ndf_all<-connect_ploy(ndf,ndf1)
  ndf_all<-connect_ploy(ndf_all,ndf2)
  ndf_all<-connect_ploy(ndf_all,ndf3)
  ndf_all$connect<-paste0(g2g[i,1],"_",g2g[i,2])
  ndf_all$group<-g2g[i,1]
  
  ndf_use<-rbind(ndf_use,ndf_all)
}

cdf<-circos_function(1)

tdf[,c(3,4)]<-xita2xy(tdf[,2],1.105)

fig3C<-ggplot()+
  geom_polygon(data=rdf,aes(x,y,group=group,fill=group))+
  geom_polygon(data=rdf2,aes(x,y,group=group),fill="black")+
  geom_polygon(data=ndf_use,aes(x,y,group=connect,fill=group),alpha=0.5)+
  geom_text(data=tdf,aes(x,y,label=tag,angle=xita*180/pi),size=2.5,hjust=0)+
  xlim(c(-2,2))+ylim(c(-2,2))+ggtitle('C')+scale_fill_manual(values=set0)+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #panel.background = element_blank(),
        panel.background=element_rect(fill="transparent",colour=NA),
        plot.background=element_rect(fill="transparent",colour=NA))


pdf('Fig3.pdf',width=10,height=10)
grid.newpage()
pushViewport(viewport(layout=grid.layout(100,100)))
print(fig3A,vp=viewport(layout.pos.row=1:50,layout.pos.col=1:100))
print(fig3B,vp=viewport(layout.pos.row=51:100,layout.pos.col=1:50))
print(fig3C,vp=viewport(layout.pos.row=51:100,layout.pos.col=51:100))
dev.off()
