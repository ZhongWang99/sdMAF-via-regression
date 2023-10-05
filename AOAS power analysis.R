
require(dplyr)

ispar38 <- function(pos)
{
  ifelse(pos<=2781479,"PAR1",ifelse(pos>=155701383,"PAR2",ifelse(pos>=89145000 & pos<=92745001,"PAR3","NPR")))
}

################################## Power analysis ##################################

# sampling z-values based on real data (from 3 by 2 plots)

snps_alt1 = c(2749353,2763977,2765220,2765769,2773342,2774420,74874310,155749677,155789452)
snps_alt2 = c(3065297,23111666,120067941,120389431,140013597,141047855,141578688,149768848,152781469)

basic = cbind.data.frame(pop = c("EUR", "EAS", "AMR", "AFR", "SAS"),
                         samplesize = c(503,504,347,661,489))
for(p in basic$pop) 
{
  chrX_p3hc.SNP.FreqTable = read.csv(file=paste0('/Users/MarkWang/StatGen/1000_genome/chrX_p3hc_',p,'.SNP.FreqTable.csv'),stringsAsFactors = F)
  assign(paste0(p,".FreqTable"), chrX_p3hc.SNP.FreqTable)
}

n_alt1.EUR = EUR.FreqTable %>% filter(POS%in%snps_alt1) %>% arrange(POS) %>% select(F_A1A1, F_A1A2, F_A2A2, M_A1A1.A1, M_A1A2, M_A2A2.A2)
n_alt1.EAS = EAS.FreqTable %>% filter(POS%in%snps_alt1) %>% arrange(POS) %>% select(F_A1A1, F_A1A2, F_A2A2, M_A1A1.A1, M_A1A2, M_A2A2.A2)
n_alt1.AMR = AMR.FreqTable %>% filter(POS%in%snps_alt1) %>% arrange(POS) %>% select(F_A1A1, F_A1A2, F_A2A2, M_A1A1.A1, M_A1A2, M_A2A2.A2)
n_alt1.AFR = AFR.FreqTable %>% filter(POS%in%snps_alt1) %>% arrange(POS) %>% select(F_A1A1, F_A1A2, F_A2A2, M_A1A1.A1, M_A1A2, M_A2A2.A2)
n_alt1.SAS = SAS.FreqTable %>% filter(POS%in%snps_alt1) %>% arrange(POS) %>% select(F_A1A1, F_A1A2, F_A2A2, M_A1A1.A1, M_A1A2, M_A2A2.A2)

n_alt2.EUR = EUR.FreqTable %>% filter(POS%in%snps_alt2) %>% arrange(POS) %>% select(F_A1A1, F_A1A2, F_A2A2, M_A1A1.A1, M_A1A2, M_A2A2.A2)
n_alt2.EAS = EAS.FreqTable %>% filter(POS%in%snps_alt2) %>% arrange(POS) %>% select(F_A1A1, F_A1A2, F_A2A2, M_A1A1.A1, M_A1A2, M_A2A2.A2)
n_alt2.AMR = AMR.FreqTable %>% filter(POS%in%snps_alt2) %>% arrange(POS) %>% select(F_A1A1, F_A1A2, F_A2A2, M_A1A1.A1, M_A1A2, M_A2A2.A2)
n_alt2.AFR = AFR.FreqTable %>% filter(POS%in%snps_alt2) %>% arrange(POS) %>% select(F_A1A1, F_A1A2, F_A2A2, M_A1A1.A1, M_A1A2, M_A2A2.A2)
n_alt2.SAS = SAS.FreqTable %>% filter(POS%in%snps_alt2) %>% arrange(POS) %>% select(F_A1A1, F_A1A2, F_A2A2, M_A1A1.A1, M_A1A2, M_A2A2.A2)


# statistical power (type-II error) (set1)
result_set1 = list()

M = 100000 # total times of simulation
for(i in 1:length(snps_alt1))
{
  par = ispar38(snps_alt1[i])
  result_set1[[i]] = list(POS=snps_alt1[i],PAR=par,sdMAFSUM=NULL,META=NULL,MEGA=NULL)
  for(p in basic$pop)
  {
    f.n_alt1 = get(paste0("n_alt1.",p))[i,1:3]
    assign(paste0("f.sample_alt1.",p), t(rmultinom(M, size=sum(f.n_alt1), prob=f.n_alt1/sum(f.n_alt1))))
    m.n_alt1 = get(paste0("n_alt1.",p))[i,4:6]
    assign(paste0("m.sample_alt1.",p), t(rmultinom(M, size=sum(m.n_alt1), prob=m.n_alt1/sum(m.n_alt1))))
  }
  
  meta_sdMAF = list()
  meta_wi = list()
  for(p in basic$pop) 
  {
    f.sample_alt1 = get(paste0("f.sample_alt1.",p))
    m.sample_alt1 = get(paste0("m.sample_alt1.",p))
    if(par=="PAR1"|par=="PAR2")
    {
      pval = apply(cbind(f.sample_alt1,m.sample_alt1),MARGIN = 1,FUN = wald.1df.hwd.auto)
      zval = apply(cbind(f.sample_alt1,m.sample_alt1),MARGIN = 1,FUN = wald.1df.hwd.auto_zval)
      sdMAF = apply(cbind(f.sample_alt1,m.sample_alt1),MARGIN = 1,FUN = sdMAF.auto)
    }
    else
    {
      pval = apply(cbind(f.sample_alt1,m.sample_alt1),MARGIN = 1,FUN = wald.1df.hwd.xchr)
      zval = apply(cbind(f.sample_alt1,m.sample_alt1),MARGIN = 1,FUN = wald.1df.hwd.xchr_zval)
      sdMAF = apply(cbind(f.sample_alt1,m.sample_alt1),MARGIN = 1,FUN = sdMAF.xchr)
    }
    meta_sdMAF[[p]] = sdMAF
    meta_se = sdMAF/zval
    meta_wi[[p]] = 1/meta_se^2
    assign(paste0('stat_',p),qchisq(-pval*log(10),df=1,lower.tail = F,log.p = T))
    print(paste0("Done with ",p,"..."))
  }
  
  f.sample_alt1 = f.sample_alt1.EUR + f.sample_alt1.EAS + f.sample_alt1.AMR + f.sample_alt1.AFR + f.sample_alt1.SAS
  m.sample_alt1 = m.sample_alt1.EUR + m.sample_alt1.EAS + m.sample_alt1.AMR + m.sample_alt1.AFR + m.sample_alt1.SAS
  
  if(par=="PAR1"|par=="PAR2") pval = apply(cbind(f.sample_alt1,m.sample_alt1),MARGIN = 1,FUN = wald.1df.hwd.auto)
  else pval = apply(cbind(f.sample_alt1,m.sample_alt1),MARGIN = 1,FUN = wald.1df.hwd.xchr)
  
  result_set1[[i]]$MEGA = pval
  
  new_stat = stat_EUR + stat_EAS + stat_AMR + stat_AFR + stat_SAS
  new_pvalue = -pchisq(new_stat,df=5,lower.tail = F,log.p=T)/log(10)
  result_set1[[i]]$sdMAFSUM = new_pvalue
  
  meta_wi_sum = numeric(M)
  meta_beta = numeric(M)
  for(p in basic$pop)
  {
    meta_wi_sum = meta_wi_sum + meta_wi[[p]]
    meta_beta = meta_beta + meta_wi[[p]]*meta_sdMAF[[p]]
  }
  meta_beta = meta_beta/meta_wi_sum
  meta_se = sqrt(1/meta_wi_sum)
  meta_z = meta_beta/meta_se
  meta_pvalue = -pchisq(meta_z^2,df=1,lower.tail = F,log.p=T)/log(10)
  result_set1[[i]]$META = meta_pvalue
}



# statistical power (type-II error) (set2)
result_set2 = list()

M = 100000 # total times of simulation
for(i in 1:length(snps_alt2))
{
  par = ispar38(snps_alt2[i])
  result_set2[[i]] = list(POS=snps_alt2[i],PAR=par,sdMAFSUM=NULL,META=NULL,MEGA=NULL)
  for(p in basic$pop)
  {
    f.n_alt2 = get(paste0("n_alt2.",p))[i,1:3]
    assign(paste0("f.sample_alt2.",p), t(rmultinom(M, size=sum(f.n_alt2), prob=f.n_alt2/sum(f.n_alt2))))
    m.n_alt2 = get(paste0("n_alt2.",p))[i,4:6]
    assign(paste0("m.sample_alt2.",p), t(rmultinom(M, size=sum(m.n_alt2), prob=m.n_alt2/sum(m.n_alt2))))
  }
  
  meta_sdMAF = list()
  meta_wi = list()
  for(p in basic$pop) 
  {
    f.sample_alt2 = get(paste0("f.sample_alt2.",p))
    m.sample_alt2 = get(paste0("m.sample_alt2.",p))
    if(par=="PAR1"|par=="PAR2")
    {
      pval = apply(cbind(f.sample_alt2,m.sample_alt2),MARGIN = 1,FUN = wald.1df.hwd.auto)
      zval = apply(cbind(f.sample_alt2,m.sample_alt2),MARGIN = 1,FUN = wald.1df.hwd.auto_zval)
      sdMAF = apply(cbind(f.sample_alt2,m.sample_alt2),MARGIN = 1,FUN = sdMAF.auto)
    }
    else
    {
      pval = apply(cbind(f.sample_alt2,m.sample_alt2),MARGIN = 1,FUN = wald.1df.hwd.xchr)
      zval = apply(cbind(f.sample_alt2,m.sample_alt2),MARGIN = 1,FUN = wald.1df.hwd.xchr_zval)
      sdMAF = apply(cbind(f.sample_alt2,m.sample_alt2),MARGIN = 1,FUN = sdMAF.xchr)
    }
    meta_sdMAF[[p]] = sdMAF
    meta_se = sdMAF/zval
    meta_wi[[p]] = 1/meta_se^2
    
    assign(paste0('stat_',p),qchisq(-pval*log(10),df=1,lower.tail = F,log.p = T))
    print(paste0("Done with ",p,"..."))
  }
  
  f.sample_alt2 = f.sample_alt2.EUR + f.sample_alt2.EAS + f.sample_alt2.AMR + f.sample_alt2.AFR + f.sample_alt2.SAS
  m.sample_alt2 = m.sample_alt2.EUR + m.sample_alt2.EAS + m.sample_alt2.AMR + m.sample_alt2.AFR + m.sample_alt2.SAS
  if(par=="PAR1"|par=="PAR2") pval = apply(cbind(f.sample_alt2,m.sample_alt2),MARGIN = 1,FUN = wald.1df.hwd.auto)
  else pval = apply(cbind(f.sample_alt2,m.sample_alt2),MARGIN = 1,FUN = wald.1df.hwd.xchr)
  
  result_set2[[i]]$MEGA = pval
  
  new_stat = stat_EUR + stat_EAS + stat_AMR + stat_AFR + stat_SAS
  new_pvalue = -pchisq(new_stat,df=5,lower.tail = F,log.p=T)/log(10)
  result_set2[[i]]$sdMAFSUM = new_pvalue
  
  meta_wi_sum = numeric(M)
  meta_beta = numeric(M)
  for(p in basic$pop)
  {
    meta_wi_sum = meta_wi_sum + meta_wi[[p]]
    meta_beta = meta_beta + meta_wi[[p]]*meta_sdMAF[[p]]
  }
  meta_beta = meta_beta/meta_wi_sum
  meta_se = sqrt(1/meta_wi_sum)
  meta_z = meta_beta/meta_se
  meta_pvalue = -pchisq(meta_z^2,df=1,lower.tail = F,log.p=T)/log(10)
  result_set2[[i]]$META = meta_pvalue
}


# output to xlsx
require(xlsx)

table_set1 = data.frame(POS=numeric(length(result_set1)),PAR=numeric(length(result_set1)),
                        sdMAFsum=numeric(length(result_set1)),META=numeric(length(result_set1)),MEGA=numeric(length(result_set1)))
table_set2 = data.frame(POS=numeric(length(result_set2)),PAR=numeric(length(result_set2)),
                        sdMAFsum=numeric(length(result_set2)),META=numeric(length(result_set2)),MEGA=numeric(length(result_set2)))

for(i in 1:length(result_set1))
{
  table_set1$POS[i] = result_set1[[i]]$POS
  table_set1$PAR[i] = result_set1[[i]]$PAR
  table_set1$sdMAFsum[i] = length(which(result_set1[[i]]$sdMAFSUM>-log10(5e-8)))/length(result_set1[[i]]$sdMAFSUM)
  table_set1$META[i] = length(which(result_set1[[i]]$META>-log10(5e-8)))/length(result_set1[[i]]$META)
  table_set1$MEGA[i] = length(which(result_set1[[i]]$MEGA>-log10(5e-8)))/length(result_set1[[i]]$MEGA)
}
for(i in 1:length(result_set2))
{
  table_set2$POS[i] = result_set2[[i]]$POS
  table_set2$PAR[i] = result_set2[[i]]$PAR
  table_set2$sdMAFsum[i] = length(which(result_set2[[i]]$sdMAFSUM>-log10(5e-8)))/length(result_set2[[i]]$sdMAFSUM)
  table_set2$META[i] = length(which(result_set2[[i]]$META>-log10(5e-8)))/length(result_set2[[i]]$META)
  table_set2$MEGA[i] = length(which(result_set2[[i]]$MEGA>-log10(5e-8)))/length(result_set2[[i]]$MEGA)
}

write.xlsx(table_set1,file="power_analysis_set1.xlsx")
write.xlsx(table_set2,file="power_analysis_set2.xlsx")





