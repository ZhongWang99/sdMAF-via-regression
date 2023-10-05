
####################################### prepare file for liftover #######################################

prep_liftover <- function(chr,pos,f="new.txt")
  # prepare liftover file for https://genome.ucsc.edu/cgi-bin/hgLiftOver
{
  c = as.character(chr)
  result = character(length(pos))
  for(i in 1:length(pos))
  {
    result[i] = paste0("chr",c," ",pos[i]," ",pos[i])
  }
  fileConn = file(f)
  writeLines(result, fileConn)
  close(fileConn)
}

decode_liftover <- function(f="new.txt",b1,b2)
  # process the output from https://genome.ucsc.edu/cgi-bin/hgLiftOver
  # b1: build before liftover
  # b2: build after liftover
  # return a dataframe
{
  data = read.table(f, sep = "" , header = F ,stringsAsFactors= F)
  result = data.frame(p1=data[,3],p2=unlist(lapply(str_split(data[,4],"-"), function(x) x[2])))
  colnames(result) = c(b1,b2)
  result
}


####################################### meta analysis #######################################


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













