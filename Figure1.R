

two <- read.csv("two.csv")[,-1]
three <- read.csv("three.csv")[,-1]
hmm <- read.csv("hmm.csv")[,-1]

pdf("Figure1.pdf",height = 5,width=6.5)
par(mar=c(0,0,0,0),oma=c(0,0,0,0),fig=c(0,1,0,1))
plot(NA, NA, type="l", ylim=c(-20, 200), xlim=c(-48, 65 ),
     xlab="", ylab="",
     xaxs="i", yaxs="i",col="#EE9572",lwd=2,
     axes=FALSE, mgp=c(1.8,0.5,0), cex.lab=1.3)

segments(-20,0,-20,max(two[,2]),lwd=20,col="#5CACEE")
for(i in 1:dim(two)[1]){
  rect(-24,two[i,2]-0.2,-16,two[i,2]+0.2,border = NA,col="#912CEE")
}
text(-20,192,"Two-point",cex=1.5)
text(-20,-13,"181 cM",cex=1.5)

segments(15,0+32,15,max(hmm[,2])+32,lwd=20,col="#5CACEE")
for(i in 1:dim(hmm)[1]){
  rect(11,hmm[i,2]-0.2+32,19,hmm[i,2]+0.2+32,border = NA,col="#912CEE")
}

text(15,192,"HMM",cex=1.5)
text(15,-13,"127 cM",cex=1.5)
segments(50,0+50,50,max(three[,2])+50,lwd=20,col="#5CACEE")
for(i in 1:dim(three)[1]){
  rect(46,three[i,2]+50-0.2,54,three[i,2]+50+0.2,border = NA,col="#912CEE")
}

text(50,192,"Three-point",cex=1.5)
text(50,-13,"81 cM",cex=1.5)
text(-35,-13,"Length",cex=1.5)
for(i in 1:dim(two)[1]){
  segments(-14,two[i,2],9,hmm[i,2]+32,col="#EBEBEB")
}

for(i in 1:dim(three)[1]){
  segments(21,hmm[i,2]+32,44,three[i,2]+50,col="#EBEBEB")
}

rect(2-30,0,0.8-30,180,border = NA,col="grey")
for(i in 0:6){
  
  rect(-1.2-30,30*i,2-30,30*i-1.5,border = NA,col="grey")
  text(-7-30,30*i-0.75,i*30,cex=1.2)
}


dev.off()


