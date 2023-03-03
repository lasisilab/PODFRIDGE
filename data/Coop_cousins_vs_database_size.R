library("wesanderson")

path<-"/Users/tinalasisi-usc/GitHub/POPFORGE/data/"
savepath <- "/Users/tinalasisi-usc/GitHub/POPFORGE/output/"

p<-c(1:8)  ##cousin degree

N<-76e6

US_pop<-read.csv(paste(path,"US_popsize.csv",sep=""))
p_grandpar_gen<-1950-30*(p+1)
these.years<-match(p_grandpar_gen,US_pop$year)
 US_Ns<-US_pop$Population[these.years]
N<- US_Ns

N<- N* 0.5*0.9
N[US_Ns<1e6]<-1e6
DB.sizes<-c(1e6,5e6,10e6)


my.cols<-wes_palette("Darjeeling1")

################Figure 1
######################################################
#Probability you have a cousins in database
######################################################

layout(t(1:2))

plot(c(1,8),c(0,1),type="n",ylab="Probability of at least one p-th cousin in database",xlab="p (degree of cousin)")
sapply(1:length(DB.sizes),function(i){
	DB.size<-DB.sizes[i]
	prob.no.rellys<-exp(-2^(2*p-2)*DB.size/N)

 	points(p,1-prob.no.rellys,type="b",col=my.cols[i],pch=19,lwd=1.5)
 })



( 2^(p+1)-1)*2^p/N

legend(x="bottomright",legend=c("Database size (Millions)=",format(DB.sizes/1e6,dig=1)),col=c(NA,my.cols),pch=19)



###I choose 2^p great^(p-1) grandparent pairs, as do you there are N/2 couples in the population
num_cousins.hg<-function(p,N,DB.size){
 (1-dhyper(x=0, m=2^(p), n=((N/2)-2^(p)), k=2^(p), log = FALSE))*DB.size
}




######################################################
##expected number of pth cousins in sample
######################################################
plot(c(1,8),c(0,1000),type="n",ylab="Expected number of p-th cousins in database",xlab="p (degree of cousin)")
sapply(1:length(DB.sizes),function(i){

    num.cousins<-4^(p)*DB.sizes[i]/(N/2)
    points(p,num.cousins,type="b",col=my.cols[i],lwd=1.5,pch=19)
   # points(p,4^(p),type="b",col="black")
})

#E_p_cousins<-( 2^(p+1)-1)*2^p

dev.copy2pdf(file=paste(savepath,"Genealogical_cousins.pdf",sep=""))

######################################################
## Genetic overlap probs
######################################################
meiosis<-p+1
##expected number of blocks shared between cousins
E.num.blocks<-2*(33.8*(2*meiosis)+22)/(2^(2*meiosis-1))
##use Poisson assumption
Prob.genetic<-1-exp(-E.num.blocks)
prob.g.e.2.blocks<-1-sapply(E.num.blocks,function(expected.num){sum(dpois(0:1,expected.num))})
prob.g.e.3.blocks<-1-sapply(E.num.blocks,function(expected.num){sum(dpois(0:2,expected.num))})


################Figure 2
######################################################
##Prob. pth cousins share blocks
######################################################
layout(t(1))

my.cols2<-wes_palette("FantasticFox1")[3:5]


#png(file=paste(path,"Prob_cousin_detected.png",sep=""))
plot(c(1,8),c(0,1),type="n",ylab="Probability p-th cousin \"detectable\"",xlab="p (degree of cousin)")
points(p,Prob.genetic,col=my.cols2[1],pch=19,type="b",lwd=2)
points(p,prob.g.e.2.blocks,col=my.cols2[2],pch=19,type="b",lwd=2)
points(p,prob.g.e.3.blocks,col=my.cols2[3],pch=19,type="b",lwd=2)
legend(x="topright",legend=c("Cousins (w. >0 genomic blocks)","Cousins (w. >1 genomic blocks)","Cousins (w. >2 genomic blocks)"),col=my.cols2[1:3],lty=1)
dev.copy2pdf(file=paste(savepath,"Prob_cousin_detected.pdf",sep=""))

################Figure 3
######################################################
##expected number of pth GENETIC cousins in sample
######################################################
layout(t(1))



plot(c(1,8),c(0,350),type="n",ylab="Expected number of genetic p-th cousins in database",xlab="p (degree of cousin)")
sapply(1:length(DB.sizes),function(i){

    num.cousins<-4^(p)*DB.sizes[i]/(N/2)
  #  points(p,num.cousins,type="b",col=my.cols[i])
    #points(p,num.cousins*Prob.genetic,type="b",lty=1,col=my.cols[i])
	#points(p,num.cousins*prob.g.e.2.blocks,type="b",lty=2,col=my.cols[i])
	points(p,num.cousins*prob.g.e.3.blocks,type="b",lty=1,col=my.cols[i])

})
my.leg<-c(c("Genetic Cousins",">2 blocks"),
c("Database size (Millions)=",format(DB.sizes/1e6,dig=1)))
legend(x="topright",legend=my.leg,col=c(NA,rep("black",1),NA,my.cols),lty=c(NA,1:2,rep(NA,3)),pch=c(rep(NA,3),rep(19,3)))

dev.copy2pdf(file=paste(savepath,"E_genetic_cousin_detected.pdf",sep=""))

######################################################
##Egs of genomic overlap in cousins
######################################################
###This relies on functions from another git
meiosis<-2; old.relly<-"Grandmother";other.relly<-"1st cousin's"
layout(t(1:3))
plot.all.chr()
chr.chunks(my.families[,1],my.col="red",meiosis=meiosis,relly.pos=1);
text(0.7,20,paste("Your genome in\n your",old.relly),cex=1.5,col="red")
plot.all.chr()
chr.chunks(my.families[,2],my.col=adjustcolor("blue",.5),meiosis=meiosis,relly.pos=1);
text(0.7,20,paste("Your",other.relly, "\n genome in\n your",old.relly),cex=1.5,col="blue")
plot.all.chr()
chr.chunks(my.families[,1],my.col="red",meiosis=meiosis,relly.pos=1);
text(0.7,20,paste("Both your genomes \n in your \n",old.relly),cex=1.5,col="purple")
chr.chunks(my.families[,2],my.col=adjustcolor("blue",.5),meiosis=meiosis,relly.pos=1)
dev.copy2pdf(file=paste(savepath,"First_cousin_overlap.pdf",sep=""))


meiosis<-4; old.relly<-"Great,\n Great Grandmother";other.relly<-"3rd cousin's"
i=1;j=5;
layout(t(1:3))
plot.all.chr()
chr.chunks(my.families[,i],my.col="red",meiosis=meiosis,relly.pos=1);
text(0.7,20,paste("Your genome in\n your",old.relly),cex=1.5,col="red")
plot.all.chr()
chr.chunks(my.families[,j],my.col=adjustcolor("blue",.5),meiosis=meiosis,relly.pos=1);
text(0.7,20,paste("Your",other.relly, "\n genome in\n your",old.relly),cex=1.5,col="blue")
plot.all.chr()
chr.chunks(my.families[,i],my.col="red",meiosis=meiosis,relly.pos=1);
text(0.7,20,paste("Both your genomes \n in your \n",old.relly),cex=1.5,col="purple")
chr.chunks(my.families[,j],my.col=adjustcolor("blue",.5),meiosis=meiosis,relly.pos=1)
dev.copy2pdf(file=paste(savepath,"Third_cousin_overlap_1.pdf",sep=""))

meiosis<-4; old.relly<-"Great,\n Great Grandmother";other.relly<-"3rd cousin's"
i=1;j=21; layout(t(1:3))
plot.all.chr()
chr.chunks(my.families[,i],my.col="red",meiosis=meiosis,relly.pos=1);
text(0.7,20,paste("Your genome in\n your",old.relly),cex=1.5,col="red")
plot.all.chr()
chr.chunks(my.families[,j],my.col=adjustcolor("blue",.5),meiosis=meiosis,relly.pos=1);
text(0.7,20,paste("Your",other.relly, "\n genome in\n your",old.relly),cex=1.5,col="blue")
plot.all.chr()
chr.chunks(my.families[,i],my.col="red",meiosis=meiosis,relly.pos=1);
text(0.7,20,paste("Both your genomes \n in your \n",old.relly),cex=1.5,col="purple")
chr.chunks(my.families[,j],my.col=adjustcolor("blue",.5),meiosis=meiosis,relly.pos=1)
dev.copy2pdf(file=paste(savepath,"Third_cousin_overlap_2.pdf",sep=""))



