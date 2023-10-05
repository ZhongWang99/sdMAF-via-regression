source("1000_genome/chr7/Semilog Manhattan Plot-copy.R")
source("1000_genome/chr7/sdMAF comparison.R")
source("sdMAF test.R")

chr = c(1:22)

for(c in chr)
{
  
  # generate PLOT for analysis
  sdMAFcomparisonORIGIN_chrAgnomad = read.csv(paste0("chr",c,"_2/sdMAFcomparisonORIGIN_chr",c,"gnomad_match1kgenomeHC.csv"),stringsAsFactors = F)
  sdMAFsumORIGIN_chrAgnomad = read.csv(paste0("chr",c,"_2/sdMAFsumORIGIN_chr",c,"gnomad_match1kgenomeHC.csv"),stringsAsFactors = F)
  
  
  tiff(paste0("gnomad2results/chr",c,"_gnomad.sdMAFcomparison_NFEvsothers_Manhattan(match1kgenomeHC).tiff"),width=3000,height=4000,res=300)
  outer = FALSE
  line = -0.2
  cex = 2
  adj  = 0.025
  nsnp = nrow(sdMAFsumORIGIN_chrAgnomad)
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
