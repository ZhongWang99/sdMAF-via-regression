
source("../1000_genome/chr7/Semilog Manhattan Plot-copy.R")
source("../1000_genome/chr7/sdMAF comparison.R")
source("../sdMAF test.R")

ispar38 <- function(pos)
{
  ifelse(pos<=2781479,"PAR1",ifelse(pos>=155701383,"PAR2",ifelse(pos>=89145000 & pos<=92745001,"PAR3","NPR")))
}



########################################## replication part 4 ##########################################
require(dplyr)
# after chrXgnomad.(pop).SNP.FreqTable are imported and matched to sdMAFsumORIGIN_chrX_1kgenome$POS

sdMAFsumORIGIN_chrX_1kgenome = read.csv("../1000_genome/chrX/sdMAFsumORIGIN2_chrX.csv",stringsAsFactors = F)
pos_1kgenome = sdMAFsumORIGIN_chrX_1kgenome$POS

chrXgnomad.nfe.SNP.FreqTable = read.csv("chrXgnomad/chrXgnomad.nfe.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.afr.SNP.FreqTable = read.csv("chrXgnomad/chrXgnomad.afr.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.amr.SNP.FreqTable = read.csv("chrXgnomad/chrXgnomad.amr.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.eas.SNP.FreqTable = read.csv("chrXgnomad/chrXgnomad.eas.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.sas.SNP.FreqTable = read.csv("chrXgnomad/chrXgnomad.sas.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.nfe.AFtest = read.csv("chrXgnomad.nfe.AFtest.csv",stringsAsFactors = F)
chrXgnomad.afr.AFtest = read.csv("chrXgnomad.afr.AFtest.csv",stringsAsFactors = F)
chrXgnomad.amr.AFtest = read.csv("chrXgnomad.amr.AFtest.csv",stringsAsFactors = F)
chrXgnomad.eas.AFtest = read.csv("chrXgnomad.eas.AFtest.csv",stringsAsFactors = F)
chrXgnomad.sas.AFtest = read.csv("chrXgnomad.sas.AFtest.csv",stringsAsFactors = F)

pos_nfe = chrXgnomad.nfe.SNP.FreqTable$POS
pos_afr = chrXgnomad.afr.SNP.FreqTable$POS
pos_amr = chrXgnomad.amr.SNP.FreqTable$POS
pos_eas = chrXgnomad.eas.SNP.FreqTable$POS
pos_sas = chrXgnomad.sas.SNP.FreqTable$POS

ind_nfe = match(pos_1kgenome,pos_nfe)
ind_afr = match(pos_1kgenome,pos_afr)
ind_amr = match(pos_1kgenome,pos_amr)
ind_eas = match(pos_1kgenome,pos_eas)
ind_sas = match(pos_1kgenome,pos_sas)
na_ind = is.na(ind_nfe)|is.na(ind_afr)|is.na(ind_amr)|is.na(ind_eas)|is.na(ind_sas)
ind_nfe = ind_nfe[!na_ind]
ind_afr = ind_afr[!na_ind]
ind_amr = ind_amr[!na_ind]
ind_eas = ind_eas[!na_ind]
ind_sas = ind_sas[!na_ind]

ind_target = which(!na_ind)

chrXgnomad.nfe.SNP.FreqTable = chrXgnomad.nfe.SNP.FreqTable[ind_nfe,]
chrXgnomad.afr.SNP.FreqTable = chrXgnomad.afr.SNP.FreqTable[ind_afr,]
chrXgnomad.amr.SNP.FreqTable = chrXgnomad.amr.SNP.FreqTable[ind_amr,]
chrXgnomad.eas.SNP.FreqTable = chrXgnomad.eas.SNP.FreqTable[ind_eas,]
chrXgnomad.sas.SNP.FreqTable = chrXgnomad.sas.SNP.FreqTable[ind_sas,]
chrXgnomad.nfe.AFtest = chrXgnomad.nfe.AFtest[ind_nfe,]
chrXgnomad.afr.AFtest = chrXgnomad.afr.AFtest[ind_afr,]
chrXgnomad.amr.AFtest = chrXgnomad.amr.AFtest[ind_amr,]
chrXgnomad.eas.AFtest = chrXgnomad.eas.AFtest[ind_eas,]
chrXgnomad.sas.AFtest = chrXgnomad.sas.AFtest[ind_sas,]


# replicate the results (signif in only 1df test results, 1000 genomes) in gnomad data

chrX_p3hc.AFtest_temp = read.csv("../1000_genome/chrX/chrX_p3hc.AFtest.csv",stringsAsFactors = F)
chrX_p3hc.SNP.FreqTable_temp = read.csv('../1000_genome/chrX/chrX_p3hc.SNP.FreqTable.csv',stringsAsFactors = F)
sdMAFsumORIGIN_chrX_1kgenome = read.csv("../1000_genome/chrX/sdMAFsumORIGIN2_chrX.csv",stringsAsFactors = F)
MetaORIGIN_chrX_1kgenome = read.csv("../1000_genome/chrX/MetaORIGIN2_chrX.csv",stringsAsFactors = F)

signif_pos1 = chrX_p3hc.AFtest_temp$BP[which(chrX_p3hc.SNP.FreqTable_temp$A_RAF>=0.05 & 
                                               chrX_p3hc.SNP.FreqTable_temp$A_RAF<=0.95 &
                                               chrX_p3hc.AFtest_temp$WALD1df.HWD>=(-log10(5e-8)))]
signif_pos2 = MetaORIGIN_chrX_1kgenome$POS[MetaORIGIN_chrX_1kgenome$pvalue>=(-log10(5e-8))]
signif_pos = unique(c(signif_pos1,signif_pos2))

crp_pval = c()
for(p in signif_pos)
{
  crp_pval = c(crp_pval, sdMAFsumORIGIN_chrX_1kgenome$pvalue[which(sdMAFsumORIGIN_chrX_1kgenome$POS==p)[1]])
} 
# target snps
snps = signif_pos[which(crp_pval<6.4)]
snps = snps[-which(snps==155787669)]
snps = snps[-which(snps==7843247)]
snps = sort(snps)


sdMAFcomparisonORIGIN_chrXgnomad = read.csv("sdMAFcomparisonORIGIN_chrXgnomad_match1kgenomeHC.csv",stringsAsFactors = F)
sdMAFsumORIGIN_chrXgnomad = read.csv("sdMAFsumORIGIN_chrXgnomad_match1kgenomeHC.csv",stringsAsFactors = F)

ind_gnomad = match(snps,sdMAFsumORIGIN_chrXgnomad$POS)
length(na.omit(ind_gnomad)) # results in 8 snps matched
snps = snps[which(!is.na(ind_gnomad))]


basic = cbind.data.frame(pop = c("nfe", "eas", "amr", "afr", "sas"),
                         samplesize = NA)
for(p in basic$pop) 
{
  chrX_p3hc.SNP.FreqTable = get(paste0("chrXgnomad.",p,".SNP.FreqTable"))
  chrX_p3hc.AFtest = get(paste0("chrXgnomad.",p,".AFtest"))
  assign(paste0(p,".FreqTable"), chrX_p3hc.SNP.FreqTable %>% filter(POS%in%snps) %>% as.data.frame)
  assign(paste0(p,".AFtest"), chrX_p3hc.AFtest %>% filter(BP%in%snps) %>% as.data.frame)
}
targ.dt = data.frame(POP=basic$pop, sdMAF=NA, sdMAF.p=NA)

save(nfe.FreqTable, eas.FreqTable, amr.FreqTable, afr.FreqTable, sas.FreqTable,
     nfe.AFtest, eas.AFtest, amr.AFtest, afr.AFtest, sas.AFtest,
     basic, targ.dt, snps, file="gnomad_replication_part4.RData")
load("gnomad_replication_part4.RData")

# plots
library(ggplot2)
library(ggpubr)
require(dplyr)
library(ggExtra)
library(gridExtra)
for(k in 1:length(snps))
{
  snp = snps[k]
  par = ispar38(snp)
  #dt.all.FreqTable = chrX_p3hc.SNP.FreqTable_temp %>% filter(POS==snp) %>% as.data.frame
  #dt.all.AFtest = chrX_p3hc.AFtest_temp %>% filter(BP==snp) %>% as.data.frame
  
  # choose nfe as reference dataset
  A1 = nfe.FreqTable$A1[which(nfe.FreqTable$POS==snp)[1]]
  dt.1df = numeric(6)
  for(i in 1:nrow(targ.dt))
  {
    p = targ.dt$POP[i]
    dt.FreqTable = get(paste0(p,".FreqTable"))
    dt.AFtest = get(paste0(p,".AFtest"))
    index1 = grep("F_A1A1", colnames(dt.FreqTable))
    
    dt.FreqTable = dt.FreqTable %>% filter(POS==snp) %>% as.data.frame
    dt.AFtest = dt.AFtest %>% filter(BP==snp) %>% as.data.frame
    
    basic$samplesize[i] = sum(dt.FreqTable[1,index1:(index1+5)])
    targ.dt$sdMAF[i] = ifelse(A1==dt.FreqTable$A1[1],dt.FreqTable$F.M_RAF[1],-dt.FreqTable$F.M_RAF[1])
    targ.dt$sdMAF.p[i] = dt.AFtest$WALD1df.HWD[1]
    if(A1==dt.FreqTable$A1[1])
    {
      dt.1df = dt.1df + dt.FreqTable[1,index1:(index1+5)]
    }
    else
    {
      dt.1df = dt.1df + dt.FreqTable[1,c(index1+2,index1+1,index1,index1+5,index1+4,index1+3)]
    }
  }
  pval.set = targ.dt$sdMAF.p
  sdMAF.set = targ.dt$sdMAF
  meta.p = meta(c(pval.set,sdMAF.set),basic$samplesize)
  if(par=="PAR1"|par=="PAR2") pval.1df = wald.1df.hwd.auto(as.numeric(dt.1df))
  else pval.1df = wald.1df.hwd.xchr(as.numeric(dt.1df))
  
  ## change the scale of y values for the purpose of visualization:
  ##        1. truncate -log10(p)>50 to 50, <1 to 1
  ##        2. then change to log(-log10(p))
  targ.dt = targ.dt %>% mutate(sdMAF.p = case_when(sdMAF.p<1 ~ 1,
                                                   sdMAF.p>50 ~ 50,
                                                   TRUE ~ sdMAF.p)) %>% 
    mutate(sdMAF.p = 0.5/log10(50)*log10(sdMAF.p) + 0.25) %>% as.data.frame
  
  
  ## Stacked plot
  gt = paste0(paste0("POS:",snp," (",par,")"),
              " \n1df -log10(p):",round(pval.1df,2),
              "; Meta -log10(p):",round(meta.p,2),
              ";\n5df -log10(p):",round(sdMAFsumORIGIN_chrXgnomad$pvalue[which(sdMAFsumORIGIN_chrXgnomad$POS==snp)],2))
  print(gt)
  assign(paste0("g",k), ggplot(targ.dt) +
           geom_hline(yintercept = 0.25, linetype="solid", size=1, color = "grey") +
           geom_hline(yintercept = 0.5/log10(50)*log10(-log10(5e-8)) + 0.25, linetype="dashed", size=0.5, color = "red") +
           geom_hline(yintercept = 0, size=0.5, color = "black") +
           geom_point(aes(x=POP,y=sdMAF), size=2, col = "blue", shape=16) +
           geom_point(aes(x=POP,y=sdMAF.p), size=2, col="black", shape=17) +
           ggtitle(gt) +
           xlab('Population') + 
           theme(plot.title = element_text(size=9, face="bold"),
                 axis.text.x = element_text(size=9),
                 axis.text.y.left = element_text(size=9),
                 axis.text.y.right = element_text(size=9),
                 axis.title.x = element_text(size=15),
                 axis.title.y.left = element_text(size=11,hjust=0.2),
                 axis.title.y.right = element_text(size=11,hjust=0.2),
                 legend.text = element_text(size=10),
                 plot.margin = unit(c(1,1,1,1), "cm")) +
           scale_y_continuous(name = 'Female - Male MAF', 
                              breaks = seq(-0.2,0.2,0.1), 
                              labels = seq(-0.2,0.2,0.1),
                              limits = c(-0.25,0.8),
                              sec.axis = sec_axis(trans = ~(.-0.25)*log10(50)/0.5, name = expression(-log[10](italic(p))), breaks = c(0,log10(5),log10(10),log10(50)),
                                                  labels = c(1,5,10,50))) +
           geom_vline(xintercept = c(7.5,11.5,16.5,21.5), linetype="dashed", color = "grey")
  )
}


# select a representative subset
tiff("chrX gnomAD 1df signif vs 5df non-signif (match1kgenome).tiff",width=3500,height=2400,res=300)
ggarrange(g1,g2,g3,g4,g5,
          labels = c("A", "B", "C", "D", "E"),
          nrow=2,ncol=3)
#plot_grid(gx1,gx2,labels=c('A','B'),ncol=1,nrow=2)
dev.off()




########################################## replication part 5 ##########################################

# For comparison of extreme SNPs between gnomAD and 1k genome HC (to generate a 3 by 2 plots)

load("gnomad_replication_part4.RData")
snps = c(120067941,141578688,149768848)  # obtained via manual comparison

chrX_p3hc.SNP.FreqTable_temp = read.csv('/Users/MarkWang/StatGen/1000_genome/chrX_p3hc_EUR.SNP.FreqTable.csv',stringsAsFactors = F)
chrX_p3hc.AFtest_temp = read.csv("/Users/MarkWang/StatGen/1000_genome/chrX_p3hc.AFtest.csv",stringsAsFactors = F)
basic_1kg = cbind.data.frame(pop = c("EUR", "EAS", "AMR", "AFR", "SAS"),
                             samplesize = c(503,504,347,661,489))
for(p in basic_1kg$pop) 
{
  chrX_p3hc.SNP.FreqTable = read.csv(file=paste0('/Users/MarkWang/StatGen/1000_genome/chrX_p3hc_',p,'.SNP.FreqTable.csv'),stringsAsFactors = F)
  chrX_p3hc.AFtest = read.csv(file=paste0('/Users/MarkWang/StatGen/1000_genome/chrX_p3hc_',p,'.AFtest.csv'),stringsAsFactors = F)
  assign(paste0(p,".FreqTable_1kg"), chrX_p3hc.SNP.FreqTable %>% filter(POS%in%snps) %>% as.data.frame)
  assign(paste0(p,".AFtest_1kg"), chrX_p3hc.AFtest %>% filter(BP%in%snps) %>% as.data.frame)
}
targ.dt_1kg = data.frame(POP=basic_1kg$pop, sdMAF=NA, sdMAF.p=NA)

# plots
library(ggplot2)
library(ggpubr)
require(dplyr)
library(ggExtra)
library(gridExtra)
for(k in 1:length(snps))
{
  snp = snps[k]
  par = ispar38(snp)
  
  ########################## 1k genome HC ##################################
  dt.all.FreqTable = chrX_p3hc.SNP.FreqTable_temp %>% filter(POS==snp) %>% as.data.frame
  dt.all.AFtest = chrX_p3hc.AFtest_temp %>% filter(BP==snp) %>% as.data.frame
  A1 = dt.all.FreqTable$A1[1]
  A2 = dt.all.FreqTable$A2[1]
  print(paste("1k genome HC",A1))
  print(paste("1k genome HC",A2))
  targ.rs = as.character(dt.all.FreqTable$ID)[1]
  for(i in 1:nrow(targ.dt_1kg))
  {
    p = targ.dt_1kg$POP[i]
    dt.FreqTable = get(paste0(p,".FreqTable_1kg"))
    dt.AFtest = get(paste0(p,".AFtest_1kg"))
    dt.FreqTable = dt.FreqTable %>% filter(POS==snp) %>% as.data.frame
    dt.AFtest = dt.AFtest %>% filter(BP==snp) %>% as.data.frame
    targ.dt_1kg$sdMAF[i] = ifelse(A1==dt.FreqTable$A1[1],dt.FreqTable$F.M_RAF[1],-dt.FreqTable$F.M_RAF[1])
    targ.dt_1kg$sdMAF.p[i] = dt.AFtest$WALD1df.HWD[1]
  }
  pval.set = targ.dt_1kg$sdMAF.p
  sdMAF.set = targ.dt_1kg$sdMAF
  meta.p = meta_iv(c(pval.set,sdMAF.set),basic_1kg$samplesize)
  
  ## change the scale of y values for the purpose of visualization:
  ##        1. truncate -log10(p)>500 to 500, <1 to 1
  ##        2. then change to log(-log10(p))
  targ.dt_1kg = targ.dt_1kg %>% mutate(sdMAF.p = case_when(sdMAF.p<1 ~ 1,
                                                           sdMAF.p>50 ~ 50,
                                                           TRUE ~ sdMAF.p)) %>% 
    mutate(sdMAF.p = 0.5/log10(50)*log10(sdMAF.p) + 0.25) %>% as.data.frame
  
  
  ## Stacked plot
  gt = paste0("1k Genome High Coverage\n",
              ifelse(targ.rs==".",paste0("POS:",snp," (",ispar38(snp),")"),
                     paste0("POS:",snp," (",ispar38(snp),") ",targ.rs)), 
              " \n1df -log10(p):",round(dt.all.AFtest$WALD1df.HWD[1],2),
              " (Mega), ",round(meta.p,2), " (Meta)",
              ";\n5df -log10(p):",round(sdMAFsumORIGIN_chrX$pvalue[which(sdMAFsumORIGIN_chrX$POS==snp)],2))
  print(gt)
  assign(paste0("g1",k), ggplot(targ.dt_1kg) +
           geom_hline(yintercept = 0.25, linetype="solid", size=1, color = "grey") +
           geom_hline(yintercept = 0.5/log10(50)*log10(-log10(5e-8)) + 0.25, linetype="dashed", size=0.5, color = "red") +
           geom_hline(yintercept = 0, size=0.5, color = "black") + 
           geom_point(aes(x=POP,y=sdMAF), size=2, col = "blue", shape=16) +
           geom_point(aes(x=POP,y=sdMAF.p), size=2, col="black", shape=17) +
           ggtitle(gt) +
           xlab('Population') + 
           theme(plot.title = element_text(size=9, face="bold"),
                 axis.text.x = element_text(size=9),
                 axis.text.y.left = element_text(size=9),
                 axis.text.y.right = element_text(size=9),
                 axis.title.x = element_text(size=15),
                 axis.title.y.left = element_text(size=11,hjust=0.2),
                 axis.title.y.right = element_text(size=11,hjust=0.2),
                 legend.text = element_text(size=10),
                 plot.margin = unit(c(1,1,1,1), "cm")) +
           scale_y_continuous(name = 'Female - Male MAF', 
                              breaks = seq(-0.2,0.2,0.1), 
                              labels = seq(-0.2,0.2,0.1),
                              limits = c(-0.25,0.8),
                              sec.axis = sec_axis(trans = ~(.-0.25)*log10(50)/0.5, name = expression(-log[10](italic(p))), breaks = c(0,log10(5),log10(10),log10(50)),
                                                  labels = c(1,5,10,50))) +
           geom_vline(xintercept = c(7.5,11.5,16.5,21.5), linetype="dashed", color = "grey")
  )
  
  ########################## gnomAD ##################################
  
  # choose nfe as reference dataset
  dt.1df = numeric(6)
  for(i in 1:nrow(targ.dt))
  {
    p = targ.dt$POP[i]
    dt.FreqTable = get(paste0(p,".FreqTable"))
    dt.AFtest = get(paste0(p,".AFtest"))
    index1 = grep("F_A1A1", colnames(dt.FreqTable))
    
    dt.FreqTable = dt.FreqTable %>% filter(POS==snp) %>% as.data.frame
    dt.AFtest = dt.AFtest %>% filter(BP==snp) %>% as.data.frame
    
    basic$samplesize[i] = sum(dt.FreqTable[1,index1:(index1+5)])
    targ.dt$sdMAF[i] = ifelse(A1!=dt.FreqTable$A1[1],dt.FreqTable$F.M_RAF[1],-dt.FreqTable$F.M_RAF[1])
    targ.dt$sdMAF.p[i] = dt.AFtest$WALD1df.HWD[1]
    if(A1==dt.FreqTable$A1[1])
    {
      dt.1df = dt.1df + dt.FreqTable[1,index1:(index1+5)]
    }
    else
    {
      dt.1df = dt.1df + dt.FreqTable[1,c(index1+2,index1+1,index1,index1+5,index1+4,index1+3)]
    }
  }
  pval.set = targ.dt$sdMAF.p
  sdMAF.set = targ.dt$sdMAF
  meta.p = meta_iv(c(pval.set,sdMAF.set),basic$samplesize)
  if(par=="PAR1"|par=="PAR2") pval.1df = wald.1df.hwd.auto(as.numeric(dt.1df))
  else pval.1df = wald.1df.hwd.xchr(as.numeric(dt.1df))
  
  ## change the scale of y values for the purpose of visualization:
  ##        1. truncate -log10(p)>50 to 50, <1 to 1
  ##        2. then change to log(-log10(p))
  targ.dt = targ.dt %>% mutate(sdMAF.p = case_when(sdMAF.p<1 ~ 1,
                                                   sdMAF.p>50 ~ 50,
                                                   TRUE ~ sdMAF.p)) %>% 
    mutate(sdMAF.p = 0.5/log10(50)*log10(sdMAF.p) + 0.25) %>% as.data.frame
  
  
  ## Stacked plot
  gt = paste0("gnomAD v3.1.2\n",
              paste0("POS:",snp," (",par,")"),
              " \n1df -log10(p):",round(pval.1df,2),
              " (Mega), ",round(meta.p,2), " (Meta)",
              ";\n5df -log10(p):",round(sdMAFsumORIGIN_chrXgnomad$pvalue[which(sdMAFsumORIGIN_chrXgnomad$POS==snp)],2))
  print(gt)
  assign(paste0("g2",k), ggplot(targ.dt) +
           geom_hline(yintercept = 0.25, linetype="solid", size=1, color = "grey") +
           geom_hline(yintercept = 0.5/log10(50)*log10(-log10(5e-8)) + 0.25, linetype="dashed", size=0.5, color = "red") +
           geom_hline(yintercept = 0, size=0.5, color = "black") +
           geom_point(aes(x=toupper(POP),y=sdMAF), size=2, col = "blue", shape=16) +
           geom_point(aes(x=toupper(POP),y=sdMAF.p), size=2, col="black", shape=17) +
           ggtitle(gt) +
           xlab('Population') + 
           theme(plot.title = element_text(size=9, face="bold"),
                 axis.text.x = element_text(size=9),
                 axis.text.y.left = element_text(size=9),
                 axis.text.y.right = element_text(size=9),
                 axis.title.x = element_text(size=15),
                 axis.title.y.left = element_text(size=11,hjust=0.2),
                 axis.title.y.right = element_text(size=11,hjust=0.2),
                 legend.text = element_text(size=10),
                 plot.margin = unit(c(1,1,1,1), "cm")) +
           scale_y_continuous(name = 'Female - Male MAF', 
                              breaks = seq(-0.2,0.2,0.1), 
                              labels = seq(-0.2,0.2,0.1),
                              limits = c(-0.25,0.8),
                              sec.axis = sec_axis(trans = ~(.-0.25)*log10(50)/0.5, name = expression(-log[10](italic(p))), breaks = c(0,log10(5),log10(10),log10(50)),
                                                  labels = c(1,5,10,50))) +
           geom_vline(xintercept = c(7.5,11.5,16.5,21.5), linetype="dashed", color = "grey")
  )
  
}


# select a representative subset
tiff("chrX gnomadVS1kgenomeHC_comparison 1df signif vs 5df non-signif.tiff",width=2400,height=3500,res=300)
ggarrange(g11,g21,g12,g22,g13,g23,
          labels = c("A1", "A2", "B1", "B2", "C1", "C2"),
          nrow=3,ncol=2)
#plot_grid(gx1,gx2,labels=c('A','B'),ncol=1,nrow=2)
dev.off()
