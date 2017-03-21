
library(mgcv)

xx <- read.table("peptide.sequence",stringsAsFactors=F)
ss <- xx[1,1]

L <- list()

seqL <- nchar(ss)
for(ii in 1:seqL) if (substr(ss,ii,ii)=="K") L <- c(L,ii)

sites <- length(L)
L <- c(L,seqL+1)
N <- 2^sites

# B-ions

BM <- as.matrix(read.table("B_ions.matrix",sep=" "))
y <- as.matrix(read.table("B_ions.mz.intensities",sep="\t"))
w <- scan("weightsB.txt")

# normalize by total intensity
for(l in 1:seqL) for(i in 1:N) for(j in 1:N) BM[N*(l-1)+i,j] = BM[N*(l-1)+i,j] * w[l]  

BM2 <- matrix(0,N*sites,N)
y2 <- rep(0,N*sites)

for(site in 1:sites){
  for(l in L[[site]]:(L[[site+1]]-1)) for(i in 1:N) for (j in 1:N) BM2[(site-1)*N+i,j] = BM2[(site-1)*N+i,j] + BM[N*(l-1)+i,j] 
  for(l in L[[site]]:(L[[site+1]]-1)) for(i in 1:N)                 y2[(site-1)*N+i]   =  y2[(site-1)*N+i]   +  y[N*(l-1)+i,3] 
}

# Y-ions

YM <- as.matrix(read.table("Y_ions.matrix",sep=" "))
y3 <- as.matrix(read.table("Y_ions.mz.intensities",sep="\t"))
w <- scan("weightsY.txt")

# normalize by total intensity
for(l in 1:seqL) for(i in 1:N) for(j in 1:N) YM[N*(l-1)+i,j] = YM[N*(l-1)+i,j] * w[l]  

# aggregate Y-ions with the same pattern
YM2 <- matrix(0,N*sites,N)
y4 <- rep(0,N*sites)

for(l in 1:L[[1]]) for(i in 1:N) for (j in 1:N) YM2[i,j] = YM2[i,j] + YM[N*(l-1)+i,j]
for(l in 1:L[[1]]) for(i in 1:N)                 y4[i]   =  y4[i]   + y3[N*(l-1)+i,3]

if(sites>=2){
	for(site in 2:sites){
	  for(l in (L[[site-1]]+1):L[[site]]) 
	     for(i in 1:N) for (j in 1:N) 
		YM2[(site-1)*N+i,j] = YM2[(site-1)*N+i,j] + YM[N*(l-1)+i,j]
	  for(l in (L[[site-1]]+1):L[[site]]) 
	      for(i in 1:N)
                 y4[(site-1)*N+i]   =  y4[(site-1)*N+i]   + y3[N*(l-1)+i,3]
         }
}

MB <- list(y=y2, w=rep(1, length(y2)), X=BM2, C=matrix(1,1,N), p=rep(1/N, N), off=array(0,0), S=list(), sp=array(0,0), Ain=diag(N), bin=rep(0,N) )
MY <- list(y=y4, w=rep(1, length(y4)), X=YM2, C=matrix(1,1,N), p=rep(1/N, N), off=array(0,0), S=list(), sp=array(0,0), Ain=diag(N), bin=rep(0,N) )
M2 <- list(y=c(y4,y2), w=rep(1, 2*length(y4)), X=rbind(YM2,BM2), C=matrix(1,1,N), p=rep(1/N, N), off=array(0,0), S=list(), sp=array(0,0), Ain=diag(N), bin=rep(0,N) )

est2 <- tryCatch({pcls(M2)},error=function(err) { return(rep(NA,N)); })

#plot(est2)
#dev.new(); plot(rbind(YM2,BM2) %*% est2+0.001,c(y4,y2)+0.001,log="xy",xlim=c(1e-3,1e7),ylim=c(1e-3,1e7)); lines(c(1e-4,1e8),c(1e-4,1e8))
write.table(est2,file="result.txt",row.names=F, col.names=F)
