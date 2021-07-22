trial<-read.csv("Field trial results.csv", header=TRUE)
attach(trial)
standarde<-function(x,npar=TRUE,print=TRUE){
	result<-sd(x,na.rm=TRUE)/sqrt(length(x))
	return(result)
}
library(agricolae)

AUDPSresults<-NULL

dates<-c(71,91,104,119) # days

for(i in 1:120) {
	evaluation<-c(X16_6[i],X30_6[i],X13_07[i],X28_07[i])
	AUDPSresults<-c(AUDPSresults,audps(evaluation,dates))  
}

AUDPSresults<-as.numeric(gsub("evaluation","",AUDPSresults))

boxplot(AUDPSresults~Line)

#AUDPS


#dates<-c(71,91,104,119) # days
# example 1: evaluation - vector
#evaluation<-c(40,80,90)
#audps(evaluation,dates)
#audps(evaluation,dates,"relative")


#Subsets
susceptible=subset(trial,Genotype=="SxM_susceptible",drop=TRUE)
resistant=subset(trial,Genotype=="SxM_resistant",drop=TRUE)
#Problem - unwanted levels. Drop by converting genotype to character
resistant$Genotype<-as.factor(as.character(resistant$Genotype))
susceptible$Genotype<-as.factor(as.character(susceptible$Genotype))


bob<-tapply(AR1xAR1_model,Genotype,mean,na.rm=TRUE)
error1<-tapply(AR1xAR1_model,Genotype,standarde)
colours<-c("red","white")
colours<-c("white","red","white","red","white","red","white","red","white","white","white","red","red","white","white","white")



SxM=rbind(susceptible,resistant)

bob<-tapply(SxM$AUDPS,SxM$Genotype,mean,na.rm=TRUE)
error1<-tapply(SxM$AUDPS,SxM$Genotype,standarde)
colours<-c("grey","white")

barCenters<-barplot(height = bob,ylim=c(0,3000),ylab = "AUDPS",col=colours,las=2)
segments(barCenters, bob - error1 * 2, barCenters,bob + error1 * 2, lwd = 1.5)
arrows(barCenters, bob - error1 * 2, barCenters,bob + error1 * 2, lwd = 1.5, angle = 90,code = 3, length = 0.05)



#T.test

susceptible1=AUDPS[which(Genotype=="SxM_susceptible")]
resistant1=AUDPS[which(Genotype=="SxM_Rrs18")]
TTEST<-t.test(susceptible1,resistant1,var.equal = TRUE)


wd=500*0.03937
ht=450*0.03937
ps=14
postscript(file="Field_barchartcontrolsSxMzero3.eps",width=wd,height=ht,pointsize=ps)
par(mar=c(11,9,1,1))


boxplot(AUDPS~Genotype,ylab="AUDPS",las=2)

boxplot(AUDPS~Genotype,data=SxM,ylab="AUDPS")
#Add significance
#text(1,y=5,"hello",cex=2)
dev.off()

boxplot(c(1:10),ylim=c(0,12),axes=F)
text(11,"**",cex=2)



boxplot(Average_27_7_16~Genotype,xlab="Genotype",ylab="Score"

a1<-aov(AUDPS~Genotype)
a2<-aov(Average_27_7_16~Genotype)

#Post hoc tests

TukeyHSD(a1)
TukeyHSD(a2)

pairwise.t.test(AUDPS,Genotype,p.adj = "none")

pairwise.t.test(AUDPS,Genotype,p.adj = "hommel")

pairwise.t.test(Average_27_7_16,Genotype,p.adj = "none")



wd=270*0.03937
ht=150*0.03937
ps=12
postscript(file="Feild trial boxplot.eps",width=wd,height=ht,pointsize=ps)

boxplot(AUDPS~Genotype,xlab="Genotype",ylab="AUDPS")
dev.off()

