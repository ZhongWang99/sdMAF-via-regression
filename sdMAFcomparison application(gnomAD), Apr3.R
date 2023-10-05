
source("../1000_genome/chr7/Semilog Manhattan Plot-copy.R")
source("../1000_genome/chr7/sdMAF comparison.R")

# Plots for application session (chr7 sdMAF comparison)

chr7gnomad.nfe.SNP.FreqTable = read.csv("chr7gnomad.nfe.SNP.FreqTable.csv",stringsAsFactors = F)
chr7gnomad.afr.SNP.FreqTable = read.csv("chr7gnomad.afr.SNP.FreqTable.csv",stringsAsFactors = F)
chr7gnomad.amr.SNP.FreqTable = read.csv("chr7gnomad.amr.SNP.FreqTable.csv",stringsAsFactors = F)
chr7gnomad.eas.SNP.FreqTable = read.csv("chr7gnomad.eas.SNP.FreqTable.csv",stringsAsFactors = F)
chr7gnomad.sas.SNP.FreqTable = read.csv("chr7gnomad.sas.SNP.FreqTable.csv",stringsAsFactors = F)

# match the POS from two data

pos_nfe = chr7gnomad.nfe.SNP.FreqTable$POS
pos_afr = chr7gnomad.afr.SNP.FreqTable$POS
pos_amr = chr7gnomad.amr.SNP.FreqTable$POS
pos_eas = chr7gnomad.eas.SNP.FreqTable$POS
pos_sas = chr7gnomad.sas.SNP.FreqTable$POS

ind_afr = match(pos_nfe,pos_afr)
ind_amr = match(pos_nfe,pos_amr)
ind_eas = match(pos_nfe,pos_eas)
ind_sas = match(pos_nfe,pos_sas)
na_ind = is.na(ind_afr)|is.na(ind_amr)|is.na(ind_eas)|is.na(ind_sas)
ind_nfe = which(!na_ind)
ind_afr = ind_afr[!na_ind]
ind_amr = ind_amr[!na_ind]
ind_eas = ind_eas[!na_ind]
ind_sas = ind_sas[!na_ind]

chr7gnomad.nfe.SNP.FreqTable = chr7gnomad.nfe.SNP.FreqTable[ind_nfe,]
chr7gnomad.afr.SNP.FreqTable = chr7gnomad.afr.SNP.FreqTable[ind_afr,]
chr7gnomad.amr.SNP.FreqTable = chr7gnomad.amr.SNP.FreqTable[ind_amr,]
chr7gnomad.eas.SNP.FreqTable = chr7gnomad.eas.SNP.FreqTable[ind_eas,]
chr7gnomad.sas.SNP.FreqTable = chr7gnomad.sas.SNP.FreqTable[ind_sas,]


target_ind = chr7gnomad.nfe.SNP.FreqTable$A_RAF>=0.05 & chr7gnomad.nfe.SNP.FreqTable$A_RAF<=0.95 &
                chr7gnomad.afr.SNP.FreqTable$A_RAF>=0.05 & chr7gnomad.afr.SNP.FreqTable$A_RAF<=0.95 &
                chr7gnomad.amr.SNP.FreqTable$A_RAF>=0.05 & chr7gnomad.amr.SNP.FreqTable$A_RAF<=0.95 &
                chr7gnomad.eas.SNP.FreqTable$A_RAF>=0.05 & chr7gnomad.eas.SNP.FreqTable$A_RAF<=0.95 &
                chr7gnomad.sas.SNP.FreqTable$A_RAF>=0.05 & chr7gnomad.sas.SNP.FreqTable$A_RAF<=0.95

pop = c("nfe", "afr", "amr", "eas", "sas")
for(p in pop) 
{
  chr7gnomad.SNP.FreqTable = get(paste0("chr7gnomad.",p,".SNP.FreqTable"))[target_ind,]
  ind1 = grep("F_A1A1", colnames(chr7gnomad.SNP.FreqTable))
  
  males = chr7gnomad.SNP.FreqTable[,(ind1+3):(ind1+5)]
  females = chr7gnomad.SNP.FreqTable[,ind1:(ind1+2)]
  
  assign(paste0('males_',p),males)
  assign(paste0('females_',p),females)
  print(paste0("Done with ",p," phase1..."))
}


pop = c("nfe", "afr", "amr", "eas", "sas")

for(p in pop) 
{
  chr7gnomad.AFtest = read.csv(file=paste0('chr7gnomad.',p,'.AFtest.csv'))[get(paste0("ind_",p)),][target_ind,]
  assign(paste0('stat_',p),qchisq(10^(-chr7gnomad.AFtest$WALD1df.HWD),df=1,lower.tail = F))
  print(p)
}

new_stat = stat_nfe + stat_afr + stat_amr + stat_eas + stat_sas
new_pvalue = -pchisq(new_stat,df=5,lower.tail = F,log.p=T)/log(10)


write.csv(data.frame(POS=chr7gnomad.AFtest$BP,pvalue=new_pvalue),
          file="sdMAFsum_chr7gnomad.csv",row.names = F)



pop = c("afr", "amr", "eas", "sas")
for(p in pop) 
{
  males = get(paste0('males_',p))
  females = get(paste0('females_',p))
  doc_p = numeric(nrow(males))
  for(i in 1:nrow(males))
  {
    doc_p[i] = sdMAFcomparison_A(list(c(as.numeric(females_nfe[i,]),as.numeric(males_nfe[i,])),
                                      c(as.numeric(females[i,]),as.numeric(males[i,]))))
  }
  assign(paste0('pval_',p),doc_p)
  print(paste0("Done with ",p," phase2..."))
}

new_pvalue = c(pval_afr, pval_amr, pval_eas, pval_sas)

write.csv(data.frame(POS=rep(chr7gnomad.nfe.SNP.FreqTable[target_ind,]$POS,length(pop)),pvalue=new_pvalue),
          file="sdMAFcomparison_chr7gnomad.csv",row.names = F)

sdMAFcomparison_chr7gnomad = read.csv("sdMAFcomparison_chr7gnomad.csv")
sdMAFsum_chr7gnomad = read.csv("sdMAFsum_chr7gnomad.csv")





tiff("chr7_gnomad.sdMAFcomparison_NFEvsothers_Manhattan(RAF_0.05).tiff",width=3000,height=5000,res=300)
outer = FALSE
line = -2
cex = 2
adj  = 0.025
par(mfrow=c(5,1))
char = c("A","B","C","D")
ind = 0
nsnp = nrow(sdMAFcomparison_chr7gnomad)/4
pop = c("afr", "amr", "eas", "sas")
ylim1 = max(10,ceiling(max(na.omit(sdMAFcomparison_chr7gnomad$pvalue)))*2)
ylim2 = max(10,ceiling(max(na.omit(sdMAFsum_chr7gnomad$pvalue)))*2)
for(p in pop) 
{
  ind = ind + 1
  pos = sdMAFcomparison_chr7gnomad$POS[((ind-1)*nsnp+1):(ind*nsnp)]
  pval = sdMAFcomparison_chr7gnomad$pvalue[((ind-1)*nsnp+1):(ind*nsnp)]
  
  m.data.temp = cbind.data.frame(CHR=7,BP=pos,
                                 SNP=NA,
                                 LOGP=pval)
  m.data.temp$BP = m.data.temp$BP/1000000
  m.data.temp$LOGP[which(m.data.temp$LOGP<1)] = 1
  
  manhattan.semilog(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='Chromosome 7 position (Mb)',ylab=expression(-log[10](italic(p))),
                    cex.axis=1,cex.lab=1.2,cex=0.5,ylim = c(1, ylim1))
  abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
  title(outer=outer,adj=adj,main=char[ind],cex.main=cex,col="black",font=2,line=line)
  title(main=paste(toupper(p),"vs","NFE sdMAF comparison"),outer=outer,cex.main=cex,col="black",font=2,line=line)
  
  print(paste0("Done with ",p,"..."))
}
m.data.temp2 = cbind.data.frame(CHR=7,BP=pos,
                                SNP=NA,
                                LOGP=sdMAFsum_chr7gnomad$pvalue)
m.data.temp2$BP = m.data.temp2$BP/1000000
m.data.temp2$LOGP[which(m.data.temp2$LOGP<1)] = 1

manhattan.semilog(m.data.temp2, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='Chromosome 7 position (Mb)',ylab=expression(-log[10](italic(p))),
                  cex.axis=1,cex.lab=1.2,cex=0.5,ylim = c(1, ylim2))
abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
title(outer=outer,adj=adj,main="E",cex.main=cex,col="black",font=2,line=line)
title(main="Multi-population sdMAF",outer=outer,cex.main=cex,col="black",font=2,line=line)
dev.off()




################################################### chrX ###################################################

# Exclude SAS

source("../1000_genome/chr7/Semilog Manhattan Plot-copy.R")
source("../1000_genome/chr7/sdMAF comparison.R")
source("../1000_genome/chr7/Score&Wald, Dec 21-copy.R")


ispar38 <- function(pos)
{
  ifelse(pos<=2781479,"PAR1",ifelse(pos>=155701383,"PAR2",ifelse(pos>=89145000 & pos<=92745001,"PAR3","NPR")))
}

chrXgnomad.nfe.SNP.FreqTable = read.csv("chrXgnomad.nfe.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.afr.SNP.FreqTable = read.csv("chrXgnomad.afr.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.amr.SNP.FreqTable = read.csv("chrXgnomad.amr.SNP.FreqTable.csv",stringsAsFactors = F)
chrXgnomad.eas.SNP.FreqTable = read.csv("chrXgnomad.eas.SNP.FreqTable.csv",stringsAsFactors = F)

# match the POS from two data

pos_nfe = chrXgnomad.nfe.SNP.FreqTable$POS
pos_afr = chrXgnomad.afr.SNP.FreqTable$POS
pos_amr = chrXgnomad.amr.SNP.FreqTable$POS
pos_eas = chrXgnomad.eas.SNP.FreqTable$POS

ind_afr = match(pos_nfe,pos_afr)
ind_amr = match(pos_nfe,pos_amr)
ind_eas = match(pos_nfe,pos_eas)
na_ind = is.na(ind_afr)|is.na(ind_amr)|is.na(ind_eas)
ind_nfe = which(!na_ind)
ind_afr = ind_afr[!na_ind]
ind_amr = ind_amr[!na_ind]
ind_eas = ind_eas[!na_ind]

chrXgnomad.nfe.SNP.FreqTable = chrXgnomad.nfe.SNP.FreqTable[ind_nfe,]
chrXgnomad.afr.SNP.FreqTable = chrXgnomad.afr.SNP.FreqTable[ind_afr,]
chrXgnomad.amr.SNP.FreqTable = chrXgnomad.amr.SNP.FreqTable[ind_amr,]
chrXgnomad.eas.SNP.FreqTable = chrXgnomad.eas.SNP.FreqTable[ind_eas,]

par = ispar38(chrXgnomad.nfe.SNP.FreqTable$POS)
target_ind1 = par=="NPR"|par=="PAR3"
target_ind2 = par=="PAR1"|par=="PAR2"


target_ind1 = target_ind1 & chrXgnomad.nfe.SNP.FreqTable$A_RAF>=0.05 & chrXgnomad.nfe.SNP.FreqTable$A_RAF<=0.95 &
  chrXgnomad.afr.SNP.FreqTable$A_RAF>=0.05 & chrXgnomad.afr.SNP.FreqTable$A_RAF<=0.95 &
  chrXgnomad.amr.SNP.FreqTable$A_RAF>=0.05 & chrXgnomad.amr.SNP.FreqTable$A_RAF<=0.95 &
  chrXgnomad.eas.SNP.FreqTable$A_RAF>=0.05 & chrXgnomad.eas.SNP.FreqTable$A_RAF<=0.95

target_ind2 = target_ind2 & chrXgnomad.nfe.SNP.FreqTable$A_RAF>=0.05 & chrXgnomad.nfe.SNP.FreqTable$A_RAF<=0.95 &
  chrXgnomad.afr.SNP.FreqTable$A_RAF>=0.05 & chrXgnomad.afr.SNP.FreqTable$A_RAF<=0.95 &
  chrXgnomad.amr.SNP.FreqTable$A_RAF>=0.05 & chrXgnomad.amr.SNP.FreqTable$A_RAF<=0.95 &
  chrXgnomad.eas.SNP.FreqTable$A_RAF>=0.05 & chrXgnomad.eas.SNP.FreqTable$A_RAF<=0.95


ind1length = length(which(target_ind1))
ind2length = length(which(target_ind2))
target_pos = c(chrXgnomad.nfe.SNP.FreqTable$POS[target_ind1],chrXgnomad.nfe.SNP.FreqTable$POS[target_ind2])

pop = c("nfe", "afr", "amr", "eas")
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

pop = c("afr", "amr", "eas")
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

new_pvalue = c(pval_afr, pval_amr, pval_eas)

write.csv(data.frame(POS=rep(target_pos,length(pop)),pvalue=new_pvalue),
          file="sdMAFcomparisonORIGIN_chrXgnomad.csv",row.names = F)



# Plots for application session (chrX sdMAF sum)
target_pos = c(chrXgnomad.nfe.SNP.FreqTable$POS[target_ind1],chrXgnomad.nfe.SNP.FreqTable$POS[target_ind2])
pop = c("nfe", "afr", "amr", "eas")
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

new_stat1 = stat1_nfe + stat1_afr + stat1_amr + stat1_eas
new_pvalue1 = -pchisq(new_stat1,df=length(pop),lower.tail = F,log.p=T)/log(10)
new_stat2 = stat2_nfe + stat2_afr + stat2_amr + stat2_eas
new_pvalue2 = -pchisq(new_stat2,df=length(pop),lower.tail = F,log.p=T)/log(10)

write.csv(data.frame(POS=target_pos,pvalue=c(new_pvalue1,new_pvalue2)),
          file="sdMAFsumORIGIN_chrXgnomad.csv",row.names = F)






# generate PLOT for analysis
sdMAFcomparisonORIGIN_chrXgnomad = read.csv("sdMAFcomparisonORIGIN_chrXgnomad.csv",stringsAsFactors = F)
sdMAFsumORIGIN_chrXgnomad = read.csv("sdMAFsumORIGIN_chrXgnomad.csv",stringsAsFactors = F)


tiff("chrX_gnomad.sdMAFcomparison_NFEvsothers_Manhattan(RAF_0.05).tiff",width=3000,height=4000,res=300)
outer = FALSE
line = -2
cex = 2
adj  = 0.025
par(mfrow=c(4,1))
nsnp = length(which(target_ind1)) + length(which(target_ind2))
char = c("A","B","C")
ind = 0
pop = c("afr", "amr", "eas")
ylim = max(10,ceiling(max(na.omit(sdMAFcomparisonORIGIN_chrXgnomad$pvalue)))*2)
#ylim = 400
for(p in pop)
{
  ind = ind + 1
  pos = sdMAFcomparisonORIGIN_chrXgnomad$POS[((ind-1)*nsnp+1):(ind*nsnp)]
  pval = sdMAFcomparisonORIGIN_chrXgnomad$pvalue[((ind-1)*nsnp+1):(ind*nsnp)]
  
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
m.data.temp2$LOGP[which(m.data.temp2$LOGP<1)] = 1

manhattan.semilog(m.data.temp2, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='X chromosome position (Mb)',ylab=expression(-log[10](italic(p))),
                  highlight = snphighlight2,xaxt="n",
                  ylim = c(1, 47982.36),
                  cex.axis=1,cex.lab=1.2,cex=0.5)
axis(1, at = seq(0, 160, by = 20), cex.axis=1,cex.lab=1.2,cex=0.5)
abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
title(outer=outer,adj=adj,main="D",cex.main=cex,col="black",font=2,line=line)
title(main="Multi-population sdMAF",outer=outer,cex.main=cex,col="black",font=2,line=line)
dev.off()








