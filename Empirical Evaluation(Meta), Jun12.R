# perform empirical evaluation for the tests appearing in the sdMAF regression paper
# BY PERMUTATION

gc_lambda <- function(p,df=5)
{
  if(df==1) med = 0.456
  else med = df*(1-2/(9*df))^3
  #round(qchisq(median(p),df=df,lower.tail=F, log.p=F)/med, digits=3)
  formatC(qchisq(median(p),df=df,lower.tail=F, log.p=F)/med, format="f", digits=3)
}

ispar38 <- function(pos)
{
  ifelse(pos<=2781479,"PAR1",ifelse(pos>=155701383,"PAR2",ifelse(pos>=89145000 & pos<=92745001,"PAR3","NPR")))
}

# 'Wald' type, 1 d.f. assuming HWD, Xchr 
wald.1df.hwd.xchr <- function(x)
  # 'Wald' type, 1 d.f. assuming HWD, Xchr 
  # x is a vector with elements F_A1A1,F_A1A2,F_A2A2,M_A1A1.A1,M_A1A2,M_A2A2.A2 respectively
{
  r0 = x[1]; r1 = x[2]; r2 = x[3]; s0 = x[4]; s1 = x[5]; s2 = x[6]
  s = s0+s2; r = r0+r1+r2
  pM = s2/s; pF = (0.5*r1+r2)/r 
  pAA.F = r2/r
  delta.F = pAA.F-pF^2
  stat = (pM-pF)^2/(1/s*(pM*(1-pM))+1/(2*r)*(pF*(1-pF)+delta.F))
  -pchisq(stat,df=1,lower.tail = F,log.p=T)/log(10)     # -log10
}

sdMAF.xchr <- function(x)
  # sex difference in MAF: F-M
  # x is a vector with elements F_A1A1,F_A1A2,F_A2A2,M_A1A1.A1,M_A1A2,M_A2A2.A2 respectively
{
  r0 = x[1]; r1 = x[2]; r2 = x[3]; s0 = x[4]; s1 = x[5]; s2 = x[6]
  s = s0+s2; r = r0+r1+r2
  pM = s2/s; pF = (0.5*r1+r2)/r 
  pF-pM
}

wald.1df.hwd.auto <- function(x)
  # 'Wald' type, 1 d.f. assuming HWD, Autosomal 
  # x is a vector with elements F_A1A1,F_A1A2,F_A2A2,M_A1A1.A1,M_A1A2,M_A2A2.A2 respectively
{
  r0 = x[1]; r1 = x[2]; r2 = x[3]; s0 = x[4]; s1 = x[5]; s2 = x[6]
  s = s0+s1+s2; r = r0+r1+r2
  pM = (0.5*s1+s2)/s; pF = (0.5*r1+r2)/r
  pAA.M = s2/s; pAA.F = r2/r
  delta.M = pAA.M-pM^2; delta.F = pAA.F-pF^2
  stat = (pM-pF)^2/(1/(2*s)*(pM*(1-pM)+delta.M)+1/(2*r)*(pF*(1-pF)+delta.F))
  -pchisq(stat,df=1,lower.tail = F,log.p=T)/log(10)     # -log10
}

sdMAF.auto <- function(x)
  # sex difference in MAF: F-M
  # x is a vector with elements F_A1A1,F_A1A2,F_A2A2,M_A1A1.A1,M_A1A2,M_A2A2.A2 respectively
{
  r0 = x[1]; r1 = x[2]; r2 = x[3]; s0 = x[4]; s1 = x[5]; s2 = x[6]
  s = s0+s1+s2; r = r0+r1+r2
  pM = (0.5*s1+s2)/s; pF = (0.5*r1+r2)/r
  pF-pM
}


meta_iv <- function(p_delta, n)
  # p_delta: c(p[i], delta[i])
  # p[i]: -log10
  # delta[i]: sdMAF
  # n: list of sample size
  # length(p_delta) = 2*length(n)
  # meta.p = apply(dt,MARGIN = 1,FUN=function(p_delta){meta(p_delta,samplesize)})
{
  l = length(n)
  p = p_delta[1:l]
  delta = p_delta[(l+1):(2*l)]
  z = qnorm(-p*log(10)-log(2),lower.tail = F,log.p = T)*sign(delta)
  se = delta/z
  w = 1/se^2
  delta_comb = sum(delta*w/sum(w))
  se_comb = sqrt(1/sum(w))
  Z = delta_comb/se_comb
  -pchisq(Z^2,df=1,lower.tail = F,log.p=T)/log(10)     # -log10
}

meta_ss <- function(p_delta, n)
  # p_delta: c(p[i], delta[i])
  # p[i]: -log10
  # delta[i]: sdMAF
  # n: list of sample size
  # length(p_delta) = 2*length(n)
  # meta.p = apply(dt,MARGIN = 1,FUN=function(p_delta){meta(p_delta,samplesize)})
{
  l = length(n)
  p = p_delta[1:l]
  delta = p_delta[(l+1):(2*l)]
  z = qnorm(-p*log(10)-log(2),lower.tail = F,log.p = T)*sign(delta)
  w = sqrt(n)
  Z = sum(z*w/sqrt(sum(n)))
  -pchisq(Z^2,df=1,lower.tail = F,log.p=T)/log(10)     # -log10
}

# 1. X chr, sdMAF summation (NPR SNPs), generate male sample and female sample
pop = c("EUR", "EAS", "AMR", "AFR", "SAS")
for(p in pop) 
{
  dt = read.csv(paste0("resimu/chrX_resimuFM_",p,".csv"))
  pval = apply(dt[,2:7],MARGIN = 1,FUN = wald.1df.hwd.xchr)
  sdMAF = apply(dt[,2:7],MARGIN = 1,FUN = sdMAF.xchr)
  assign(paste0('pval_',p),pval)
  assign(paste0('sdMAF_',p),sdMAF)
  
  print(paste0("Done with ",p,"..."))
}

new_dt = cbind.data.frame(pval_EUR, pval_EAS, pval_AMR, pval_AFR, pval_SAS, sdMAF_EUR, sdMAF_EAS, sdMAF_AMR, sdMAF_AFR, sdMAF_SAS)
new_pvalue = apply(new_dt,MARGIN = 1,FUN = function(x)meta_iv(x,c(503,504,347,661,489)))

write.csv(data.frame(POS=dt$POS,pvalue=new_pvalue),
          file="MetaIV_resimuFM_chrX.csv",row.names = F)


MetaIV_resimuFM_chrX = read.csv("MetaIV_resimuFM_chrX.csv")


pop = c("EUR", "EAS", "AMR", "AFR", "SAS")
for(p in pop) 
{
  dt = read.csv(paste0("resimu/chrX_resimuFM_",p,".csv"))
  pval = apply(dt[,2:7],MARGIN = 1,FUN = wald.1df.hwd.xchr)
  sdMAF = apply(dt[,2:7],MARGIN = 1,FUN = sdMAF.xchr)
  assign(paste0('pval_',p),pval)
  assign(paste0('sdMAF_',p),sdMAF)
  
  print(paste0("Done with ",p,"..."))
}

new_dt = cbind.data.frame(pval_EUR, pval_EAS, pval_AMR, pval_AFR, pval_SAS, sdMAF_EUR, sdMAF_EAS, sdMAF_AMR, sdMAF_AFR, sdMAF_SAS)
new_pvalue = apply(new_dt,MARGIN = 1,FUN = function(x)meta_ss(x,c(503,504,347,661,489)))

write.csv(data.frame(POS=dt$POS,pvalue=new_pvalue),
          file="MetaSS_resimuFM_chrX.csv",row.names = F)


MetaSS_resimuFM_chrX = read.csv("MetaSS_resimuFM_chrX.csv")

par(mfrow=c(1,2))
hist(10^(-MetaSS_resimuFM_chrX$pvalue),breaks=100,main="sample-size based")
hist(10^(-MetaIV_resimuFM_chrX$pvalue),breaks=100,main="inverse-variance based")


# chrX PAR
pop = c("EUR", "EAS", "AMR", "AFR", "SAS")
for(p in pop) 
{
  dt = read.csv(paste0("resimu/chrX.PAR_resimuFM_",p,".csv"))
  pval = apply(dt[,2:7],MARGIN = 1,FUN = wald.1df.hwd.auto)
  sdMAF = apply(dt[,2:7],MARGIN = 1,FUN = sdMAF.auto)
  assign(paste0('pval_',p),pval)
  assign(paste0('sdMAF_',p),sdMAF)
  
  print(paste0("Done with ",p,"..."))
}

new_dt = cbind.data.frame(pval_EUR, pval_EAS, pval_AMR, pval_AFR, pval_SAS, sdMAF_EUR, sdMAF_EAS, sdMAF_AMR, sdMAF_AFR, sdMAF_SAS)
new_pvalue = apply(new_dt,MARGIN = 1,FUN = function(x)meta_ss(x,c(503,504,347,661,489)))

write.csv(data.frame(POS=dt$POS,pvalue=new_pvalue),
          file="MetaSS_resimuFM_chrX.PAR.csv",row.names = F)


# chr7
pop = c("EUR", "EAS", "AMR", "AFR", "SAS")
for(p in pop) 
{
  dt = read.csv(paste0("../chr7/resimu/chr7_resimuFM_",p,".csv"))
  pval = apply(dt[,2:7],MARGIN = 1,FUN = wald.1df.hwd.auto)
  sdMAF = apply(dt[,2:7],MARGIN = 1,FUN = sdMAF.auto)
  assign(paste0('pval_',p),pval)
  assign(paste0('sdMAF_',p),sdMAF)
  
  print(paste0("Done with ",p,"..."))
}

new_dt = cbind.data.frame(pval_EUR, pval_EAS, pval_AMR, pval_AFR, pval_SAS, sdMAF_EUR, sdMAF_EAS, sdMAF_AMR, sdMAF_AFR, sdMAF_SAS)
new_pvalue = apply(new_dt,MARGIN = 1,FUN = function(x)meta_ss(x,c(503,504,347,661,489)))

write.csv(data.frame(POS=dt$POS,pvalue=new_pvalue),
          file="../chr7/MetaSS_resimuFM_chr7.csv",row.names = F)


MetaSS_resimuFM_chrX = read.csv("MetaSS_resimuFM_chrX.csv")
MetaSS_resimuFM_chrX.PAR = read.csv("MetaSS_resimuFM_chrX.PAR.csv")
MetaSS_resimuFM_chr7 = read.csv("../chr7/MetaSS_resimuFM_chr7.csv")

tiff("EmpEval_multipop.sdMAF_Meta.tiff",width=3000,height=3000,res=300)
par(mfrow=c(2,2))
outer = FALSE
line = 1
cex = 2
adj  = 0.025
hist(10^(-MetaSS_resimuFM_chr7$pvalue),breaks=50,xlab="p value",main="Chromosome 7",col="white")
title(outer=outer,adj=adj,main="A",cex.main=cex,col="black",font=2,line=line)
hist(10^(-MetaSS_resimuFM_chrX.PAR$pvalue),breaks=50,xlab="p value",main="X chromosome PAR",col="white")
title(outer=outer,adj=adj,main="B",cex.main=cex,col="black",font=2,line=line)
hist(10^(-MetaSS_resimuFM_chrX$pvalue),breaks=50,xlab="p value",main="X chromosome NPR",col="white")
title(outer=outer,adj=adj,main="C",cex.main=cex,col="black",font=2,line=line)
dev.off()











