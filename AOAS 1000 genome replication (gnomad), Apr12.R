
source("../1000_genome/chr7/Semilog Manhattan Plot-copy.R")
source("../1000_genome/chr7/sdMAF comparison.R")
source("../sdMAF test.R")

ispar38 <- function(pos)
{
  ifelse(pos<=2781479,"PAR1",ifelse(pos>=155701383,"PAR2",ifelse(pos>=89145000 & pos<=92745001,"PAR3","NPR")))
}


########################################## replication part 1 ##########################################
# folder named chrX2
sdMAFsumORIGIN_chrX_1kgenome = read.csv("../1000_genome/chrX/sdMAFsumORIGIN2_chrX.csv",stringsAsFactors = F)
pos_1kgenome = sdMAFsumORIGIN_chrX_1kgenome$POS

chrXgnomad.nfe.SNP.FreqTable = read.csv("chrXgnomad.nfe.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.afr.SNP.FreqTable = read.csv("chrXgnomad.afr.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.amr.SNP.FreqTable = read.csv("chrXgnomad.amr.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.eas.SNP.FreqTable = read.csv("chrXgnomad.eas.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.sas.SNP.FreqTable = read.csv("chrXgnomad.sas.SNP.FreqTable.csv",stringsAsFactors = F)

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

par = ispar38(chrXgnomad.nfe.SNP.FreqTable$POS)
target_ind1 = par=="NPR"|par=="PAR3"
target_ind2 = par=="PAR1"|par=="PAR2"



ind1length = length(which(target_ind1))
ind2length = length(which(target_ind2))
target_pos = c(chrXgnomad.nfe.SNP.FreqTable$POS[target_ind1],chrXgnomad.nfe.SNP.FreqTable$POS[target_ind2])

pop = c("nfe", "eas", "amr", "afr", "sas")
for(p in pop) 
{
  chrXgnomad.SNP.FreqTable = get(paste0("chrXgnomad.",p,".SNP.FreqTable"))
  ind1 = grep("F_A1A1", colnames(chrXgnomad.SNP.FreqTable))
  
  females = chrXgnomad.SNP.FreqTable[c(which(target_ind1),which(target_ind2)),ind1:(ind1+2)]
  males = chrXgnomad.SNP.FreqTable[c(which(target_ind1),which(target_ind2)),(ind1+3):(ind1+5)]
  
  assign(paste0('males_',p),males)
  assign(paste0('females_',p),females)
  print(paste0("Done with ",p," phase1..."))
}

pop = c("eas", "amr", "afr", "sas")
for(p in pop) 
{
  males = get(paste0('males_',p))
  females = get(paste0('females_',p))
  doc_p = numeric(ind1length+ind2length)
  for(i in 1:ind1length)
  {
    doc_p[i] = sdMAFcomparison_X(list(c(as.numeric(females_nfe[i,]),as.numeric(males_nfe[i,])),
                                      c(as.numeric(females[i,]),as.numeric(males[i,]))))
  }
  for(i in (ind1length+1):(ind1length+ind2length))
  {
    doc_p[i] = sdMAFcomparison_A(list(c(as.numeric(females_nfe[i,]),as.numeric(males_nfe[i,])),
                                      c(as.numeric(females[i,]),as.numeric(males[i,]))))
  }
  assign(paste0('pval_',p),doc_p)
  print(paste0("Done with ",p," phase2..."))
}

new_pvalue = c()
for(p in pop) new_pvalue = c(new_pvalue, get(paste0("pval_",p)))

write.csv(data.frame(POS=rep(target_pos,length(pop)),pvalue=new_pvalue,pop=rep(pop,each=ind1length+ind2length)),
          file="sdMAFcomparisonORIGIN_chrXgnomad_match1kgenomeHC.csv",row.names = F)



# Plots for application session (chrX sdMAF sum)
target_pos = c(chrXgnomad.nfe.SNP.FreqTable$POS[target_ind1],chrXgnomad.nfe.SNP.FreqTable$POS[target_ind2])
pop = c("nfe", "eas", "amr", "afr", "sas")
for(p in pop) 
{
  chrXgnomad.SNP.FreqTable1 = get(paste0("chrXgnomad.",p,".SNP.FreqTable"))[which(target_ind1),]
  chrXgnomad.SNP.FreqTable2 = get(paste0("chrXgnomad.",p,".SNP.FreqTable"))[which(target_ind2),]
  ind1 = grep("F_A1A1", colnames(chrXgnomad.SNP.FreqTable1))
  
  dt1 = chrXgnomad.SNP.FreqTable1[,ind1:(ind1+5)]
  dt2 = chrXgnomad.SNP.FreqTable2[,ind1:(ind1+5)]
  
  pval1 = apply(dt1,MARGIN = 1,FUN = wald.1df.hwd.xchr)
  pval2 = apply(dt2,MARGIN = 1,FUN = wald.1df.hwd.auto)
  
  assign(paste0('stat1_',p),qchisq(10^(-pval1),df=1,lower.tail = F))
  assign(paste0('stat2_',p),qchisq(10^(-pval2),df=1,lower.tail = F))
  print(paste0("Done with ",p,"..."))
}

new_stat1 = numeric(ind1length)
for(p in pop) new_stat1 = new_stat1 + get(paste0("stat1_",p))
new_stat2 = numeric(ind2length)
for(p in pop) new_stat2 = new_stat2 + get(paste0("stat2_",p))
new_pvalue1 = -pchisq(new_stat1,df=length(pop),lower.tail = F,log.p=T)/log(10)
new_pvalue2 = -pchisq(new_stat2,df=length(pop),lower.tail = F,log.p=T)/log(10)

write.csv(data.frame(POS=target_pos,pvalue=c(new_pvalue1,new_pvalue2)),
          file="sdMAFsumORIGIN_chrXgnomad_match1kgenomeHC.csv",row.names = F)



# generate PLOT for analysis
sdMAFcomparisonORIGIN_chrXgnomad = read.csv("sdMAFcomparisonORIGIN_chrXgnomad_match1kgenomeHC.csv",stringsAsFactors = F)
sdMAFsumORIGIN_chrXgnomad = read.csv("sdMAFsumORIGIN_chrXgnomad_match1kgenomeHC.csv",stringsAsFactors = F)


tiff("chrX_gnomad.sdMAFcomparison_NFEvsothers_Manhattan(match1kgenomeHC).tiff",width=3000,height=4000,res=300)
outer = FALSE
line = -0.2
cex = 2
adj  = 0.025
nsnp = ind1length + ind2length
char = c("A","B","C","D")
ind = 0
pop = c("eas", "amr", "afr", "sas")
par(mfrow=c(length(pop)+1,1))
ylim = 1000
#ylim = max(10,ceiling(max(na.omit(sdMAFcomparisonORIGIN_chrXgnomad$pvalue)))*2)
#ylim = 400
for(p in pop)
{
  ind = ind + 1
  pos = sdMAFcomparisonORIGIN_chrXgnomad$POS[sdMAFcomparisonORIGIN_chrXgnomad$pop==p]
  pval = sdMAFcomparisonORIGIN_chrXgnomad$pvalue[sdMAFcomparisonORIGIN_chrXgnomad$pop==p]
  
  # sdMAF comparison test (Combined RAF>=0.05)
  m.data.temp = cbind.data.frame(CHR=23,BP=pos,
                                 SNP=NA,
                                 LOGP=pval)
  
  lgt = nrow(m.data.temp)
  snp.sp = as.character(sample(0:100000000,size=lgt,replace = F))
  m.data.temp$SNP = snp.sp
  snphighlight = m.data.temp$SNP[which(m.data.temp$BP<=2781479|m.data.temp$BP>=155701383|(m.data.temp$BP>=89145000 & m.data.temp$BP<=92745001))]
  m.data.temp$BP = m.data.temp$BP/1000000
  m.data.temp$LOGP[which(m.data.temp$LOGP<1)] = 1
  
  manhattan.semilog(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='X chromosome position (Mb)',ylab=expression(-log[10](italic(p))),
                    highlight = snphighlight,xaxt="n",
                    ylim = c(1, ylim),
                    cex.axis=1,cex.lab=1.2,cex=0.5)
  axis(1, at = seq(0, 160, by = 20), cex.axis=1,cex.lab=1.2,cex=0.5)
  abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
  title(outer=outer,adj=adj,main=char[ind],cex.main=cex,col="black",font=2,line=line)
  title(main=paste(toupper(p),"vs","NFE sdMAF comparison"),outer=outer,cex.main=cex,col="black",font=2,line=line)
  
  print(paste0("Done with ",p,"..."))
}

m.data.temp2 = cbind.data.frame(CHR=23,BP=pos,
                                SNP=NA,
                                LOGP=sdMAFsumORIGIN_chrXgnomad$pvalue)

lgt = nrow(m.data.temp2)
snp.sp = as.character(sample(0:100000000,size=lgt,replace = F))
m.data.temp2$SNP = snp.sp
snphighlight2 = m.data.temp2$SNP[which(m.data.temp2$BP<=2781479|m.data.temp2$BP>=155701383|(m.data.temp2$BP>=89145000 & m.data.temp2$BP<=92745001))]
m.data.temp2$BP = m.data.temp2$BP/1000000
m.data.temp2$LOGP[m.data.temp2$LOGP<1] = 1
m.data.temp2$LOGP[m.data.temp2$LOGP>1000] = 1000

manhattan.semilog(m.data.temp2, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='X chromosome position (Mb)',ylab=expression(-log[10](italic(p))),
                  highlight = snphighlight2,xaxt="n",
                  ylim = c(1, 1000),
                  cex.axis=1,cex.lab=1.2,cex=0.5)
axis(1, at = seq(0, 160, by = 20), cex.axis=1,cex.lab=1.2,cex=0.5)
abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
title(outer=outer,adj=adj,main="E",cex.main=cex,col="black",font=2,line=line)
title(main="Multi-population sdMAF",outer=outer,cex.main=cex,col="black",font=2,line=line)
dev.off()


# Pval-Pval comparison plot (run on Mac)
require(scales)
sdMAFsumORIGIN_chrX_1kgenome = read.csv("sdMAFsumORIGIN2_chrX.csv",stringsAsFactors = F)
sdMAFsumORIGIN_chrXgnomad = read.csv("sdMAFsumORIGIN_chrXgnomad_match1kgenomeHC.csv",stringsAsFactors = F)

save(ind_target,file="ind_target.RData")
load("ind_target.RData")

tiff("PPplot(sdMAFsum_1kgenomeHCvsgnomAD).tiff",width=2400,height=2000,res=300)
data1 = cbind.data.frame(p_1kgenome=sdMAFsumORIGIN_chrX_1kgenome$pvalue[ind_target], 
                         p_gnomAD=sdMAFsumORIGIN_chrXgnomad$pvalue)
data1$p_1kgenome[which(data1$p_1kgenome==Inf)] = 1000
data1$p_gnomAD[which(data1$p_gnomAD==Inf)] = 1000
ggplot(data1, aes(p_1kgenome, p_gnomAD)) + 
  xlab(expression(-log[10](italic(p['1k genome HC'])))) +
  ylab(expression(-log[10](italic(p['gnomAD'])))) +
  scale_x_continuous(trans='log10',breaks = c(1e-4,1e-2,1e0,1e2),labels = c('0.0001','0.01','1','100'), limits=c(0.00001,1100)) +
  scale_y_continuous(trans='log10',breaks = c(1e-4,1e-2,1e0,1e2),labels = c('0.0001','0.01','1','100'), limits=c(0.00001,1100)) +
  geom_hex(bins = 100) + 
  geom_abline(slope=1, intercept = 0, linetype="dashed", col='grey') +
  geom_hline(yintercept = -log10(5e-8), linetype="dotted", col='red') + 
  geom_vline(xintercept = -log10(5e-8), linetype="dotted", col='red')
dev.off()



########################################## replication part 2 ##########################################
require(dplyr)
# after chrXgnomad.(pop).SNP.FreqTable are imported and matched to sdMAFsumORIGIN_chrX_1kgenome$POS

sdMAFsumORIGIN_chrX_1kgenome = read.csv("../1000_genome/chrX/sdMAFsumORIGIN2_chrX.csv",stringsAsFactors = F)
pos_1kgenome = sdMAFsumORIGIN_chrX_1kgenome$POS

chrXgnomad.nfe.SNP.FreqTable = read.csv("chrXgnomad.nfe.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.afr.SNP.FreqTable = read.csv("chrXgnomad.afr.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.amr.SNP.FreqTable = read.csv("chrXgnomad.amr.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.eas.SNP.FreqTable = read.csv("chrXgnomad.eas.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.sas.SNP.FreqTable = read.csv("chrXgnomad.sas.SNP.FreqTable.csv",stringsAsFactors = F)
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


# replicate the results (signif in only 5df test results, 1000 genomes) in gnomad data

chrX_p3hc.AFtest_temp = read.csv("../1000_genome/chrX/chrX_p3hc.AFtest.csv",stringsAsFactors = F)
sdMAFsumORIGIN_chrX_1kgenome = read.csv("../1000_genome/chrX/sdMAFsumORIGIN2_chrX.csv",stringsAsFactors = F)

signif_pos = sdMAFsumORIGIN_chrX_1kgenome$POS[which(sdMAFsumORIGIN_chrX_1kgenome$pvalue>=-log10(5e-8))]
crp_pval = c()
for(p in signif_pos)
{
  crp_pval = c(crp_pval, chrX_p3hc.AFtest_temp$WALD1df.HWD[which(chrX_p3hc.AFtest_temp$BP==p)[1]])
}  

#snps = signif_pos[which(crp_pval<(-log10(5e-8)))]
snps = signif_pos[which(crp_pval<3)]
#snps = signif_pos[which(crp_pval<2)]
snps = snps[-which(snps==74874311)]
snps = snps[-which(snps==155787669)]
snps = sort(snps)

snps = c(2749353, 2773342, 2774420)
sdMAFcomparisonORIGIN_chrXgnomad = read.csv("sdMAFcomparisonORIGIN_chrXgnomad_match1kgenomeHC.csv",stringsAsFactors = F)
sdMAFsumORIGIN_chrXgnomad = read.csv("sdMAFsumORIGIN_chrXgnomad_match1kgenomeHC.csv",stringsAsFactors = F)

ind_gnomad = match(snps,sdMAFsumORIGIN_chrXgnomad$POS)
length(na.omit(ind_gnomad)) # results in 2 snps matched
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
     basic, targ.dt, snps, file="gnomad_replication_part2.RData")

load("gnomad_replication_part2.RData")

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
tiff("chrX gnomAD 1df non-signif vs 5df signif (match1kgenome).tiff",width=3500,height=3500,res=300)
ggarrange(g1,g2,g3,g4,g5,g6,g7,g8,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          nrow=3,ncol=3)
#plot_grid(gx1,gx2,labels=c('A','B'),ncol=1,nrow=2)
dev.off()



########################################## replication part 3 ##########################################

# For comparison of extreme SNPs between gnomAD and 1k genome HC (to generate a 3 by 2 plots)

load("gnomad_replication_part2.RData")
snps = c(2749353, 2773342, 2774420)

chrX_p3hc.SNP.FreqTable_temp = read.csv('/Users/MarkWang/StatGen/1000_genome/chrX_p3hc.SNP.FreqTable.csv',stringsAsFactors = F)
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
              " (Mega), ",round(meta.p,2)," (Meta)",
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
tiff("chrX gnomadVS1kgenomeHC_comparison 1df non-signif vs 5df signif.tiff",width=2400,height=3500,res=300)
ggarrange(g11,g21,g12,g22,g13,g23,
          labels = c("A1", "A2", "B1", "B2", "C1", "C2"),
          nrow=3,ncol=2)
#plot_grid(gx1,gx2,labels=c('A','B'),ncol=1,nrow=2)
dev.off()




##############################################################################
# Chromosome-wide plot
##############################################################################

source("Semilog Manhattan Plot-copy.R")

chr = c(1:22)

tiff("AllAuto_gnomad.5Manhattan(match1kgenomeHC).tiff",width=3000,height=5000,res=300)
outer = FALSE
line = 0
cex = 2
adj  = 0.025
par(mfrow=c(5,1))
pop = c("eas", "amr", "afr", "sas")
char = c("A","B","C","D")
ind = 0
for(p in pop) 
{
  ind = ind + 1
  CHR = NULL
  BP = NULL
  LOGP = NULL
  for(ch in chr)
  {
    sdMAFcomparison_auto = read.csv(paste0('chr',ch,'_2/','sdMAFcomparisonORIGIN_chr',ch,'gnomad_match1kgenomeHC.csv'))
    nsnp = nrow(sdMAFcomparison_auto)/4
    pos = sdMAFcomparison_auto$POS[((ind-1)*nsnp+1):(ind*nsnp)]
    pval = sdMAFcomparison_auto$pvalue[((ind-1)*nsnp+1):(ind*nsnp)]
    
    CHR = append(CHR,rep(ch,nsnp))
    BP = append(BP,pos)
    LOGP = append(LOGP,pval)
    print(ch)
  }
  
  m.data.temp = cbind.data.frame(CHR=CHR,BP=BP,
                                 SNP=NA,
                                 LOGP=LOGP)
  m.data.temp$BP = m.data.temp$BP/1000000
  m.data.temp$LOGP[which(m.data.temp$LOGP<1)] = 1
  m.data.temp$LOGP[which(m.data.temp$LOGP>400)] = 400
  
  manhattan.semilog(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='Chromosome',ylab=expression(-log[10](italic(p))),
                    cex.axis=1,cex.lab=1.2,cex=0.5,ylim=c(1,401))
  abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
  title(outer=outer,adj=adj,main=char[ind],cex.main=cex,col="black",font=2,line=line)
  title(main=paste(toupper(p),"vs","NFE sdMAF comparison"),outer=outer,cex.main=cex,col="black",font=2,line=line)
  
  print(paste0("Done with ",p,"..."))
  
}

CHR = NULL
BP = NULL
LOGP = NULL
for(c in chr)
{
  sdMAFsum_chrA = read.csv(paste0("chr",c,"_2/sdMAFsumORIGIN_chr",c,"gnomad_match1kgenomeHC.csv"))
  num_snp = nrow(sdMAFsum_chrA)
  CHR = append(CHR,rep(c,num_snp))
  BP = append(BP,sdMAFsum_chrA$POS)
  LOGP = append(LOGP,sdMAFsum_chrA$pvalue)
  print(c)
}

m.data.temp = cbind.data.frame(CHR=CHR,BP=BP,SNP=NA,LOGP=LOGP)
m.data.temp$BP = m.data.temp$BP/1000000
m.data.temp$LOGP[which(m.data.temp$LOGP<1)] = 1
m.data.temp$LOGP[which(m.data.temp$LOGP>400)] = 400

manhattan.semilog(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='Chromosome',ylab=expression(-log[10](italic(p))),
                  cex.axis=0.6,cex.lab=1.2,cex=0.5,ylim=c(1,401))
abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
title(outer=outer,adj=adj,main="E",cex.main=cex,col="black",font=2,line=line)
title(main="Multi-population sdMAF",outer=outer,cex.main=cex,col="black",font=2,line=line)

dev.off()








