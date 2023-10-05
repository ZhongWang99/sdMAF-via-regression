
source("../1000_genome/chr7/Semilog Manhattan Plot-copy.R")
source("../1000_genome/chr7/sdMAF comparison.R")
source("../sdMAF test.R")


########################################## replication part 1 ##########################################
sdMAFsumORIGIN_chr7_1kgenome = read.csv("../1000_genome/chr7/sdMAFsum_chr7.csv",stringsAsFactors = F)
pos_1kgenome = sdMAFsumORIGIN_chr7_1kgenome$POS

chr7_2_gnomad.nfe.SNP.FreqTable = read.csv("chr7_2_gnomad.nfe.SNP.FreqTable.csv",stringsAsFactors = F)
chr7_2_gnomad.afr.SNP.FreqTable = read.csv("chr7_2_gnomad.afr.SNP.FreqTable.csv",stringsAsFactors = F)
chr7_2_gnomad.amr.SNP.FreqTable = read.csv("chr7_2_gnomad.amr.SNP.FreqTable.csv",stringsAsFactors = F)
chr7_2_gnomad.eas.SNP.FreqTable = read.csv("chr7_2_gnomad.eas.SNP.FreqTable.csv",stringsAsFactors = F)
chr7_2_gnomad.sas.SNP.FreqTable = read.csv("chr7_2_gnomad.sas.SNP.FreqTable.csv",stringsAsFactors = F)

pos_nfe = chr7_2_gnomad.nfe.SNP.FreqTable$POS
pos_afr = chr7_2_gnomad.afr.SNP.FreqTable$POS
pos_amr = chr7_2_gnomad.amr.SNP.FreqTable$POS
pos_eas = chr7_2_gnomad.eas.SNP.FreqTable$POS
pos_sas = chr7_2_gnomad.sas.SNP.FreqTable$POS

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

chr7_2_gnomad.nfe.SNP.FreqTable = chr7_2_gnomad.nfe.SNP.FreqTable[ind_nfe,]
chr7_2_gnomad.afr.SNP.FreqTable = chr7_2_gnomad.afr.SNP.FreqTable[ind_afr,]
chr7_2_gnomad.amr.SNP.FreqTable = chr7_2_gnomad.amr.SNP.FreqTable[ind_amr,]
chr7_2_gnomad.eas.SNP.FreqTable = chr7_2_gnomad.eas.SNP.FreqTable[ind_eas,]
chr7_2_gnomad.sas.SNP.FreqTable = chr7_2_gnomad.sas.SNP.FreqTable[ind_sas,]


indlength = length(ind_target)
target_pos = chr7_2_gnomad.nfe.SNP.FreqTable$POS

pop = c("nfe", "eas", "amr", "afr", "sas")
for(p in pop) 
{
  chr7_2_gnomad.SNP.FreqTable = get(paste0("chr7_2_gnomad.",p,".SNP.FreqTable"))
  ind1 = grep("F_A1A1", colnames(chr7_2_gnomad.SNP.FreqTable))
  
  females = chr7_2_gnomad.SNP.FreqTable[,ind1:(ind1+2)]
  males = chr7_2_gnomad.SNP.FreqTable[,(ind1+3):(ind1+5)]
  
  assign(paste0('males_',p),males)
  assign(paste0('females_',p),females)
  print(paste0("Done with ",p," phase1..."))
}

pop = c("eas", "amr", "afr", "sas")
for(p in pop) 
{
  males = get(paste0('males_',p))
  females = get(paste0('females_',p))
  doc_p = numeric(indlength)
  for(i in 1:indlength)
  {
    doc_p[i] = sdMAFcomparison_A(list(c(as.numeric(females_nfe[i,]),as.numeric(males_nfe[i,])),
                                      c(as.numeric(females[i,]),as.numeric(males[i,]))))
  }
  assign(paste0('pval_',p),doc_p)
  print(paste0("Done with ",p," phase2..."))
}

new_pvalue = c()
for(p in pop) new_pvalue = c(new_pvalue, get(paste0("pval_",p)))

write.csv(data.frame(POS=rep(target_pos,length(pop)),pvalue=new_pvalue,pop=rep(pop,each=indlength)),
          file="sdMAFcomparisonORIGIN_chr7gnomad_match1kgenomeHC.csv",row.names = F)



# Plots for application session (chr7 sdMAF sum)
target_pos = chr7_2_gnomad.nfe.SNP.FreqTable$POS
pop = c("nfe", "eas", "amr", "afr", "sas")
for(p in pop) 
{
  chr7_2_gnomad.SNP.FreqTable = get(paste0("chr7_2_gnomad.",p,".SNP.FreqTable"))
  ind1 = grep("F_A1A1", colnames(chr7_2_gnomad.SNP.FreqTable))
  
  dt = chr7_2_gnomad.SNP.FreqTable[,ind1:(ind1+5)]
  
  pval = apply(dt,MARGIN = 1,FUN = wald.1df.hwd.auto)
  
  assign(paste0('stat_',p),qchisq(10^(-pval),df=1,lower.tail = F))
  print(paste0("Done with ",p,"..."))
}

new_stat = numeric(indlength)
for(p in pop) new_stat = new_stat + get(paste0("stat_",p))
new_pvalue = -pchisq(new_stat,df=length(pop),lower.tail = F,log.p=T)/log(10)

write.csv(data.frame(POS=target_pos,pvalue=new_pvalue),
          file="sdMAFsumORIGIN_chr7gnomad_match1kgenomeHC.csv",row.names = F)



# generate PLOT for analysis
sdMAFcomparisonORIGIN_chr7gnomad = read.csv("sdMAFcomparisonORIGIN_chr7gnomad_match1kgenomeHC.csv",stringsAsFactors = F)
sdMAFsumORIGIN_chr7gnomad = read.csv("sdMAFsumORIGIN_chr7gnomad_match1kgenomeHC.csv",stringsAsFactors = F)


tiff("chr7_gnomad.sdMAFcomparison_NFEvsothers_Manhattan(match1kgenomeHC).tiff",width=3000,height=4000,res=300)
outer = FALSE
line = -0.2
cex = 2
adj  = 0.025
nsnp = indlength
char = c("A","B","C","D")
ind = 0
pop = c("eas", "amr", "afr", "sas")
par(mfrow=c(length(pop)+1,1))
#ylim = 1000
ylim = max(10,ceiling(max(na.omit(sdMAFcomparisonORIGIN_chr7gnomad$pvalue)))*2)
#ylim = 400
for(p in pop)
{
  ind = ind + 1
  pos = sdMAFcomparisonORIGIN_chr7gnomad$POS[sdMAFcomparisonORIGIN_chr7gnomad$pop==p]
  pval = sdMAFcomparisonORIGIN_chr7gnomad$pvalue[sdMAFcomparisonORIGIN_chr7gnomad$pop==p]
  
  # sdMAF comparison test (Combined RAF>=0.05)
  m.data.temp = cbind.data.frame(CHR=7,BP=pos,
                                 SNP=NA,
                                 LOGP=pval)
  
  lgt = nrow(m.data.temp)
  snp.sp = as.character(sample(0:100000000,size=lgt,replace = F))
  m.data.temp$SNP = snp.sp
  m.data.temp$BP = m.data.temp$BP/1000000
  m.data.temp$LOGP[which(m.data.temp$LOGP<1)] = 1
  
  manhattan.semilog(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='Chromosome 7 position (Mb)',ylab=expression(-log[10](italic(p))),
                    xaxt="n",
                    ylim = c(1, ylim),
                    cex.axis=1,cex.lab=1.2,cex=0.5)
  axis(1, at = seq(0, 160, by = 20), cex.axis=1,cex.lab=1.2,cex=0.5)
  abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
  title(outer=outer,adj=adj,main=char[ind],cex.main=cex,col="black",font=2,line=line)
  title(main=paste(toupper(p),"vs","NFE sdMAF comparison"),outer=outer,cex.main=cex,col="black",font=2,line=line)
  
  print(paste0("Done with ",p,"..."))
}

m.data.temp2 = cbind.data.frame(CHR=7,BP=pos,
                                SNP=NA,
                                LOGP=sdMAFsumORIGIN_chr7gnomad$pvalue)

lgt = nrow(m.data.temp2)
snp.sp = as.character(sample(0:100000000,size=lgt,replace = F))
m.data.temp2$SNP = snp.sp
m.data.temp2$BP = m.data.temp2$BP/1000000
m.data.temp2$LOGP[m.data.temp2$LOGP<1] = 1
m.data.temp2$LOGP[m.data.temp2$LOGP>1000] = 1000

manhattan.semilog(m.data.temp2, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='Chromosome 7 position (Mb)',ylab=expression(-log[10](italic(p))),
                  xaxt="n",
                  ylim = c(1, 1000),
                  cex.axis=1,cex.lab=1.2,cex=0.5)
axis(1, at = seq(0, 160, by = 20), cex.axis=1,cex.lab=1.2,cex=0.5)
abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
title(outer=outer,adj=adj,main="E",cex.main=cex,col="black",font=2,line=line)
title(main="Multi-population sdMAF",outer=outer,cex.main=cex,col="black",font=2,line=line)
dev.off()


# Pval-Pval comparison plot (run on Mac)
require(scales)
sdMAFsumORIGIN_chr7_1kgenome = read.csv("../1000_genome/chr7/sdMAFsum_chr7.csv",stringsAsFactors = F)
sdMAFsumORIGIN_chr7gnomad = read.csv("sdMAFsumORIGIN_chr7gnomad_match1kgenomeHC.csv",stringsAsFactors = F)

save(ind_target,sdMAFsumORIGIN_chr7_1kgenome,sdMAFsumORIGIN_chr7gnomad,file="gnomad_chr7.RData")
load("gnomad_chr7.RData")

tiff("PPplot_chr7(sdMAFsum_1kgenomeHCvsgnomAD).tiff",width=2400,height=2000,res=300)
data1 = cbind.data.frame(p_1kgenome=sdMAFsumORIGIN_chr7_1kgenome$pvalue[ind_target], 
                         p_gnomAD=sdMAFsumORIGIN_chr7gnomad$pvalue)
data1$p_1kgenome[which(data1$p_1kgenome>1000)] = 1000
data1$p_gnomAD[which(data1$p_gnomAD>1000)] = 1000
data1$p_1kgenome[which(data1$p_1kgenome<1e-5)] = 1e-5
data1$p_gnomAD[which(data1$p_gnomAD<1e-5)] = 1e-5
ggplot(data1, aes(p_1kgenome, p_gnomAD)) + 
  xlab(expression(-log[10](italic(p['1k genome HC'])))) +
  ylab(expression(-log[10](italic(p['gnomAD'])))) +
  scale_x_continuous(trans='log10',breaks = c(1e-4,1e-2,1e0,1e2),labels = c('0.0001','0.01','1','100'), limits=c(0.00001,110)) +
  scale_y_continuous(trans='log10',breaks = c(1e-4,1e-2,1e0,1e2),labels = c('0.0001','0.01','1','100'), limits=c(0.00001,110)) +
  geom_hex(bins = 100) + 
  geom_abline(slope=1, intercept = 0, linetype="dashed", col='grey') +
  geom_hline(yintercept = -log10(5e-8), linetype="dotted", col='red') + 
  geom_vline(xintercept = -log10(5e-8), linetype="dotted", col='red')
dev.off()



########################################## replication part 2 ##########################################
# replicate the results for all autosomes

source("1000_genome/chr7/Semilog Manhattan Plot-copy.R")
source("1000_genome/chr7/sdMAF comparison.R")
source("sdMAF test.R")

chr = c(1:5)

for(c in chr)
{
  sdMAFsumORIGIN_chrA_1kgenome = read.csv(paste0("1000_genome/chr",c,"/sdMAFsum_chr",c,".csv"),stringsAsFactors = F)
  pos_1kgenome = sdMAFsumORIGIN_chrA_1kgenome$POS
  
  chrA_2_gnomad.nfe.SNP.FreqTable = read.csv(paste0("chr",c,"_2/chr",c,"_2_gnomad.nfe.SNP.FreqTable.csv"),stringsAsFactors = F)
  chrA_2_gnomad.afr.SNP.FreqTable = read.csv(paste0("chr",c,"_2/chr",c,"_2_gnomad.afr.SNP.FreqTable.csv"),stringsAsFactors = F)
  chrA_2_gnomad.amr.SNP.FreqTable = read.csv(paste0("chr",c,"_2/chr",c,"_2_gnomad.amr.SNP.FreqTable.csv"),stringsAsFactors = F)
  chrA_2_gnomad.eas.SNP.FreqTable = read.csv(paste0("chr",c,"_2/chr",c,"_2_gnomad.eas.SNP.FreqTable.csv"),stringsAsFactors = F)
  chrA_2_gnomad.sas.SNP.FreqTable = read.csv(paste0("chr",c,"_2/chr",c,"_2_gnomad.sas.SNP.FreqTable.csv"),stringsAsFactors = F)
  
  pos_nfe = chrA_2_gnomad.nfe.SNP.FreqTable$POS
  pos_afr = chrA_2_gnomad.afr.SNP.FreqTable$POS
  pos_amr = chrA_2_gnomad.amr.SNP.FreqTable$POS
  pos_eas = chrA_2_gnomad.eas.SNP.FreqTable$POS
  pos_sas = chrA_2_gnomad.sas.SNP.FreqTable$POS
  
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
  
  chrA_2_gnomad.nfe.SNP.FreqTable = chrA_2_gnomad.nfe.SNP.FreqTable[ind_nfe,]
  chrA_2_gnomad.afr.SNP.FreqTable = chrA_2_gnomad.afr.SNP.FreqTable[ind_afr,]
  chrA_2_gnomad.amr.SNP.FreqTable = chrA_2_gnomad.amr.SNP.FreqTable[ind_amr,]
  chrA_2_gnomad.eas.SNP.FreqTable = chrA_2_gnomad.eas.SNP.FreqTable[ind_eas,]
  chrA_2_gnomad.sas.SNP.FreqTable = chrA_2_gnomad.sas.SNP.FreqTable[ind_sas,]
  
  
  indlength = length(ind_target)
  target_pos = chrA_2_gnomad.nfe.SNP.FreqTable$POS
  
  pop = c("nfe", "eas", "amr", "afr", "sas")
  for(p in pop) 
  {
    chrA_2_gnomad.SNP.FreqTable = get(paste0("chrA_2_gnomad.",p,".SNP.FreqTable"))
    ind1 = grep("F_A1A1", colnames(chrA_2_gnomad.SNP.FreqTable))
    
    females = chrA_2_gnomad.SNP.FreqTable[,ind1:(ind1+2)]
    males = chrA_2_gnomad.SNP.FreqTable[,(ind1+3):(ind1+5)]
    
    assign(paste0('males_',p),males)
    assign(paste0('females_',p),females)
    print(paste0("Done with ",p," phase1..."))
  }
  
  pop = c("eas", "amr", "afr", "sas")
  for(p in pop) 
  {
    males = get(paste0('males_',p))
    females = get(paste0('females_',p))
    doc_p = numeric(indlength)
    for(i in 1:indlength)
    {
      doc_p[i] = sdMAFcomparison_A(list(c(as.numeric(females_nfe[i,]),as.numeric(males_nfe[i,])),
                                        c(as.numeric(females[i,]),as.numeric(males[i,]))))
    }
    assign(paste0('pval_',p),doc_p)
    print(paste0("Done with ",p," phase2..."))
  }
  
  new_pvalue = c()
  for(p in pop) new_pvalue = c(new_pvalue, get(paste0("pval_",p)))
  
  write.csv(data.frame(POS=rep(target_pos,length(pop)),pvalue=new_pvalue,pop=rep(pop,each=indlength)),
            file=paste0("chr",c,"_2/sdMAFcomparisonORIGIN_chr",c,"gnomad_match1kgenomeHC.csv"),row.names = F)
  
  
  
  # Plots for application session (chrA sdMAF sum)
  target_pos = chrA_2_gnomad.nfe.SNP.FreqTable$POS
  pop = c("nfe", "eas", "amr", "afr", "sas")
  for(p in pop) 
  {
    chrA_2_gnomad.SNP.FreqTable = get(paste0("chrA_2_gnomad.",p,".SNP.FreqTable"))
    ind1 = grep("F_A1A1", colnames(chrA_2_gnomad.SNP.FreqTable))
    
    dt = chrA_2_gnomad.SNP.FreqTable[,ind1:(ind1+5)]
    
    pval = apply(dt,MARGIN = 1,FUN = wald.1df.hwd.auto)
    
    assign(paste0('stat_',p),qchisq(10^(-pval),df=1,lower.tail = F))
    print(paste0("Done with ",p,"..."))
  }
  
  new_stat = numeric(indlength)
  for(p in pop) new_stat = new_stat + get(paste0("stat_",p))
  new_pvalue = -pchisq(new_stat,df=length(pop),lower.tail = F,log.p=T)/log(10)
  
  write.csv(data.frame(POS=target_pos,pvalue=new_pvalue),
            file=paste0("chr",c,"_2/sdMAFsumORIGIN_chr",c,"gnomad_match1kgenomeHC.csv"),row.names = F)
  
  
  
  # generate PLOT for analysis
  sdMAFcomparisonORIGIN_chrAgnomad = read.csv(paste0("chr",c,"_2/sdMAFcomparisonORIGIN_chr",c,"gnomad_match1kgenomeHC.csv"),stringsAsFactors = F)
  sdMAFsumORIGIN_chrAgnomad = read.csv(paste0("chr",c,"_2/sdMAFsumORIGIN_chr",c,"gnomad_match1kgenomeHC.csv"),stringsAsFactors = F)
  
  
  tiff(paste0("gnomad2results/chr",c,"_gnomad.sdMAFcomparison_NFEvsothers_Manhattan(match1kgenomeHC).tiff"),width=3000,height=4000,res=300)
  outer = FALSE
  line = -0.2
  cex = 2
  adj  = 0.025
  nsnp = indlength
  char = c("A","B","C","D")
  ind = 0
  pop = c("eas", "amr", "afr", "sas")
  par(mfrow=c(length(pop)+1,1))
  #ylim = 1000
  ylim = max(10,ceiling(max(na.omit(sdMAFcomparisonORIGIN_chrAgnomad$pvalue)))*2)
  #ylim = 400
  for(p in pop)
  {
    ind = ind + 1
    pos = sdMAFcomparisonORIGIN_chrAgnomad$POS[sdMAFcomparisonORIGIN_chrAgnomad$pop==p]
    pval = sdMAFcomparisonORIGIN_chrAgnomad$pvalue[sdMAFcomparisonORIGIN_chrAgnomad$pop==p]
    
    # sdMAF comparison test (Combined RAF>=0.05)
    m.data.temp = cbind.data.frame(CHR=c,BP=pos,
                                   SNP=NA,
                                   LOGP=pval)
    
    lgt = nrow(m.data.temp)
    snp.sp = as.character(sample(0:100000000,size=lgt,replace = F))
    m.data.temp$SNP = snp.sp
    m.data.temp$BP = m.data.temp$BP/1000000
    m.data.temp$LOGP[which(m.data.temp$LOGP<1)] = 1
    
    manhattan.semilog(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab=paste0('Chromosome ',c,' position (Mb)'),ylab=expression(-log[10](italic(p))),
                      ylim = c(1, ylim),
                      cex.axis=1,cex.lab=1.2,cex=0.5)
    abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
    title(outer=outer,adj=adj,main=char[ind],cex.main=cex,col="black",font=2,line=line)
    title(main=paste(toupper(p),"vs","NFE sdMAF comparison"),outer=outer,cex.main=cex,col="black",font=2,line=line)
    
    print(paste0("Done with ",p,"..."))
  }
  
  m.data.temp2 = cbind.data.frame(CHR=c,BP=pos,
                                  SNP=NA,
                                  LOGP=sdMAFsumORIGIN_chrAgnomad$pvalue)
  
  lgt = nrow(m.data.temp2)
  snp.sp = as.character(sample(0:100000000,size=lgt,replace = F))
  m.data.temp2$SNP = snp.sp
  m.data.temp2$BP = m.data.temp2$BP/1000000
  m.data.temp2$LOGP[m.data.temp2$LOGP<1] = 1
  m.data.temp2$LOGP[m.data.temp2$LOGP>1000] = 1000
  
  manhattan.semilog(m.data.temp2, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab=paste0('Chromosome ',c,' position (Mb)'),ylab=expression(-log[10](italic(p))),
                    ylim = c(1, 1000),
                    cex.axis=1,cex.lab=1.2,cex=0.5)
  abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
  title(outer=outer,adj=adj,main="E",cex.main=cex,col="black",font=2,line=line)
  title(main="Multi-population sdMAF",outer=outer,cex.main=cex,col="black",font=2,line=line)
  dev.off()


}


########################################## replication part 2 ##########################################

# Chromosome-wide plot


source("Semilog Manhattan Plot-copy.R")

chr = c(1:22)

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

tiff("AllAuto_gnomad.sdMAFsum_Manhattan(match1kgenomeHC).tiff",width=4000,height=2000,res=300)
outer = FALSE
line = -2
cex = 2
adj  = 0.025
manhattan.semilog(m.data.temp, chr="CHR", bp="BP", snp="SNP", p="LOGP", suggestiveline=FALSE, genomewideline=FALSE, logp=FALSE,xlab='Chromosome',ylab=expression(-log[10](italic(p))),
                  cex.axis=0.6,cex.lab=1.2,cex=0.5,ylim=c(1,401))
abline(h=-log(0.05/1000000,base=10),col='red',lty='dashed',lwd=1.5)
dev.off()











