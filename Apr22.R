source("1000_genome/chr7/Semilog Manhattan Plot-copy.R")
source("1000_genome/chr7/sdMAF comparison.R")
source("sdMAF test.R")

chr = c(2:5)

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
  
  
}




