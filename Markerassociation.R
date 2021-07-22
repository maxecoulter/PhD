"""Script used for fine mapping of Rrs18 by REML. Results of analyses found in Coulter et al. 2019 (https://link.springer.com/article/10.1007/s00122-018-3262-8)
Author: Max Coulter"""




test<-read.csv("Association of all L73a data.csv",header=TRUE)
test<-read.csv("L73a for REML.csv",header=TRUE)
#test<-read.csv("L73a for REMLno7_06.csv",header=TRUE)

test<-read.csv("L73a for REMLno7_06no_parents.csv",header=TRUE)
test<-read.csv("L73a for REMLno7_06no_parents_all_predicted.csv",header=TRUE)


test2<-subset(test,Line!="01_03*",drop=TRUE)
test<-test2
#attach(test)
#Test if variance is equal
#var.test(Predicted.means..REML.~marker3)
#var.test(Predicted.means..REML.~marker30)
#Variance is equal


#ttest<-t.test(Predicted.means..REML.~test[,20],var.equal=T)
#str(ttest)
#ttest$p.value

#HIstogram
#hist(REML_1_2_and_DLA_19,breaks = 8, main = NULL,ylim = c(0,10),xlab = "Disease score",ylab = "Number of lines",col = "lightblue")

#susceptible=subset(test,X11264412=="b",drop=TRUE)
#resistant=subset(test,X11264412=="a",drop=TRUE)

#hist(susceptible$REML_1_2_and_DLA_19,breaks = 8, main = NULL,ylim = c(0,10),xlab = "Disease score",ylab = "Number of lines",col = "lightblue")
#hist(resistant$REML_1_2_and_DLA_19,breaks = 8, main = NULL,ylim = c(0,10),xlab = "Disease score",ylab = "Number of lines",col = "lightblue")


###Standard deviation for each line

listoflines<-NULL
listofsd<-NULL
test2<-subset(test,test$DLA==1.2,drop=TRUE)
test2<-subset(test,test$DLA==19,drop=TRUE)

for ( i in levels(test2$Line))
	{
	test3<-subset(test2,test2$Line==i,drop=TRUE)
	listoflines<-c(listoflines,i) 
	#listofsd<-c(listofsd,sd(test3$Greatest.lesion.size.0.removed,na.rm=TRUE))
	listofsd<-c(listofsd,sd(test3$Greatest.lesion.size,na.rm=TRUE))
	}

sdout<-cbind(listoflines,listofsd)
write(sdout,file="sdout_map_REML1_2_0s.csv")

write(sdout,file="sdout_map_REML19_0s.csv")








X9846178
colnames(test)
#Marker association analysis with t.test (REML data)
sig<-NULL
pos<-NULL
LOD<-NULL


#Phenotype Column names from test[,10-16]
#Map 1.2 is 10
#REML_predicted_means_17_08_17 is 11
#

#for (i in 2:9)
	#{
	#Box<-t.test(test[,17]~test[,i],var.equal=T)
	#sig<-c(sig,-log10(Box$p.value))
	#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
	#LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
	#position<-colnames(test)[i]
	#position<-as.numeric(gsub("X","",position))
	#pos<-c(pos,position)
	#}
#plot(pos,sig,type="l")
#plot(pos,LOD,type="l")
#############################REML######################
#Use lme4
#install.packages("lme4")
library(lme4)
testsubset<-test

#As<-subset(test,X11571800=="a",drop=TRUE)
#Bs<-subset(test,X11571800=="b",drop=TRUE)
#testsubset<-rbind(As,Bs)#Remove all parent genotypes from analysis

#Greatest.lesion.size.0.removed is y variable
test.model<- lmer(Greatest.lesion.size.0.removed ~ X11572955 + DLA + (1|Family) + (1|Line) + (1|Box), data=testsubset)
test.model.null<- lmer(Greatest.lesion.size.0.removed ~ DLA + (1|Family) + (1|Line) + (1|Box), data=testsubset)
sig<-anova(test.model,test.model.null)#This give significance of a particular marker


#Example below
##politeness.model  =  lmer(frequency  ~  attitude  + 
##+ (1|subject) + (1|scenario), data=politeness)


#to extract p value from anova:

#sig<-unlist(sig)
#log<-(-log10(as.numeric(sig["Pr(>Chisq)2"])))

#Loop to run this for all markers

sig<-NULL
pos<-NULL
LOD<-NULL
for (i in 1:11)
	{
	#create new dataframe with only one marker
	testsubset2<-data.frame(marker=testsubset[,i],DLA=testsubset[,18],Line=testsubset[,15],score=testsubset[,17],Box=testsubset[,12],Family=testsubset[,14])
	test.model<- lmer(score ~ marker + DLA + (1|Line) + (1|Box), data=testsubset2)
	test.model.null<- lmer(score ~ DLA + (1|Line) + (1|Box), data=testsubset2)
	sig1<-anova(test.model,test.model.null)#This give significance of a particular marker, refitted with ML
	sig1<-unlist(sig1)
	pvalue<-as.numeric(sig1["Pr(>Chisq)2"])
	sig<-c(sig,-log10(pvalue))
	
	#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
	LOD<-c(LOD,(qchisq(p=(pvalue), df=1, lower.tail=FALSE)/(2*log(10))))
	position<-colnames(testsubset)[i]
	position<-as.numeric(gsub("X","",position))
	pos<-c(pos,position)
	}
plot(pos,LOD,type="l")

#1.5 LOD support interval
high<-LOD[which.max(LOD)]

support<-(high-1.5)
 
supportl<-approx(x = LOD,y = pos, xout = support)
supportr<-approx(x = LOD[7:11],y = pos[7:11],xout = support) 
supportl<-as.numeric(supportl[2])
supportr<-as.numeric(supportr[2])

abline(v=supportl,lty="dashed")
abline(v=supportr,lty="dashed")

#l = 10974267 (04/01/18)
#r = 11578960 (04/01/18)

#l(18/5/18) = 10963646
#r (18/5/18) = 11578950



#Permutations
highestlods<-NULL
highestsigs<-NULL

for(i in 1:100){
	sig<-NULL
	pos<-NULL
	LOD<-NULL
	testsubset$random<-NULL
	random<-sample(testsubset[,17])
	testsubset$random<-c(random)
	for(i in 1:11){
		#create new dataframe with only one marker
		testsubset2<-data.frame(marker=testsubset[,i],DLA=testsubset[,18],Line=testsubset[,15],score=testsubset[,17],Box=testsubset[,12],Family=testsubset[,14],random=testsubset[,19])
		test.model<- lmer(random ~ marker + DLA + (1|Line) + (1|Box), data=testsubset2)
		test.model.null<- lmer(random ~ DLA + (1|Line) + (1|Box), data=testsubset2)
		sig1<-anova(test.model,test.model.null)#This give significance of a particular marker
		sig1<-unlist(sig1)
		pvalue<-as.numeric(sig1["Pr(>Chisq)2"])
		sig<-c(sig,-log10(pvalue))
	
		#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
		LOD<-c(LOD,(qchisq(p=(pvalue), df=1, lower.tail=FALSE)/(2*log(10))))
		position<-colnames(testsubset)[i]
		position<-as.numeric(gsub("X","",position))
		pos<-c(pos,position)
	}

	
	high<-LOD[which.max(LOD)]
	highestlods<-c(highestlods,high)
	highsigs<-sig[which.max(sig)]
	highestsigs<-c(highestsigs,highsigs)
}


highestlods<-sort(highestlods, decreasing=TRUE)
fifthpercentileLOD<-highestlods[5]

highestsigs<-sort(highestsigs, decreasing=TRUE)
fifthpercentilesig<-highestsigs[5]
#rm(abline)
abline(h=fifthpercentilesig,lty="dashed")
abline(h=fifthpercentileLOD,lty="longdash")

#Make smarter ggplot
install.packages("ggplot2")
library(ggplot2)
#Change positions to megabases
pos<-pos/1000000
graph<-cbind(pos,LOD)
graph<-as.data.frame(graph)

P<-ggplot(data=graph,aes(x=pos, y=LOD, group=1)) + geom_line(size=1)+ xlab("Physical map position (Mbp)")+ylab("LOD") +
theme_bw()+ scale_x_continuous(labels = comma)+theme(axis.text=element_text(size=16),axis.title=element_text(size=16))

#theme(panel.background = element_rect(fill = 'white', colour = 'black'))
#Nt<-4
#v1<-c(0:Nt)*1000000
#axis(side = 1, at= v1, labels =(v1/1000000))
#lab=v1/1000000
require(scales)
#p + scale_x_continuous(labels = comma)
#xlab("Physical map position (bp)")+ylab("LOD")

P + geom_hline(yintercept=fifthpercentileLOD,lty="dotted")+
geom_vline(xintercept=(supportl/1000000),lty="dashed")+
geom_vline(xintercept=(supportr/1000000),lty="dashed")

wd=129*0.03937
ht=129*0.03937
ps=12
postscript(file="L73a_fine_map_REML_18_5.eps",width=wd,height=ht,pointsize=ps)

P + geom_hline(yintercept=fifthpercentileLOD,lty="dotted")+
geom_vline(xintercept=supportl/1000000,lty="dashed")+
geom_vline(xintercept=supportr/1000000,lty="dashed")


dev.off()


##################################

#DLA19 on its own

pos<-NULL
LOD<-NULL
sig<-NULL

for (i in 2:9)
	{
	Box<-t.test(test[,12]~test[,i],var.equal=T)
	sig<-c(sig,-log10(Box$p.value))
	#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
	LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
	position<-colnames(test)[i]
	position<-as.numeric(gsub("X","",position))
	pos<-c(pos,position)
	}



plot(pos,LOD,type="l")
wd=129*0.03937
ht=129*0.03937
ps=12
postscript(file="L73a_fine_map_DLA19_no7_06big.eps",width=wd,height=ht,pointsize=ps)
dev.off()
P + geom_hline(yintercept=fifthpercentileLOD,lty="dashed")
P
plot(pos,LOD,type="l",xlab="Physical position (Mbp) on 6HS",lwd=2,ylab="LOD",xaxt="n")
#number of labels
Nt<-5
v1<-c(0:Nt)*(7e+07/Nt)
axis(side = 1, at= v1, labels =(v1/1000000))  

highestlods<-NULL
highestsigs<-NULL

for(i in 1:100){
	sig<-NULL
	pos<-NULL
	LOD<-NULL
	test$random<-NULL
	random<-sample(test[,12])
	test$random<-c(random)
	for(k in 2:9){
		Box<-t.test(random~test[,k],var.equal=T)
		sig<-c(sig,-log10(Box$p.value))
		#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
		LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
		position<-colnames(test)[k]
		position<-as.numeric(gsub("X","",position))
		pos<-c(pos,position)
	}
	high<-LOD[which.max(LOD)]
	highestlods<-c(highestlods,high)
	highsigs<-sig[which.max(sig)]
	highestsigs<-c(highestsigs,highsigs)
}


highestlods<-sort(highestlods, decreasing=TRUE)
fifthpercentileLOD<-highestlods[5]





##############################################
#Add 1.2 and DLA 19 together
test2<-test[,1:9]


DLA1_2<-test[,10]
DLA_19<-test[,12]
DLA1_19<-rbind(test2,test2)


both<-c(DLA1_2,DLA_19)

DLA1_19<-cbind(DLA1_19,both)


sig<-NULL
pos<-NULL
LOD<-NULL
for (i in 2:9)
	{
	Box<-t.test(DLA1_19[,10]~DLA1_19[,i],var.equal=T)
	sig<-c(sig,-log10(Box$p.value))
	#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
	LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
	position<-colnames(test)[i]
	position<-as.numeric(gsub("X","",position))
	pos<-c(pos,position)
	}
plot(pos,LOD,type="l")

#REML
pos<-NULL
LOD<-NULL
sig<-NULL

for (i in 2:9)
	{
	Box<-t.test(test[,11]~test[,i],var.equal=T)
	sig<-c(sig,-log10(Box$p.value))
	#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
	LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
	position<-colnames(test)[i]
	position<-as.numeric(gsub("X","",position))
	pos<-c(pos,position)
	}

plot(pos,LOD,type="l")



#This gives high -Log10 and good results

high<-LOD[which.max(LOD)]
abline(h=high-2,lty="dashed")

postscript(file="L73a_fine_map_DLA19_Map_1_2_together_no7_06big.eps",width=wd,height=ht,pointsize=ps)
dev.off()

#delete rows

#Good rows:1-34, 53-86

DLA1_9_1<-DLA1_19[1:34,]
DLA1_9_2<-DLA1_19[53:86,]
DLA1_9<-rbind(DLA1_19[1:34,],DLA1_19[53:86,])
DLA1_19<-DLA1_9

#Permutations
highestlods<-NULL
highestsigs<-NULL

for(i in 1:100){
	sig<-NULL
	pos<-NULL
	LOD<-NULL
	DLA1_19$random<-NULL
	random<-sample(DLA1_19[,10])
	DLA1_19$random<-c(random)
	for(k in 2:9){
		Box<-t.test(random~DLA1_19[,k],var.equal=T)
		sig<-c(sig,-log10(Box$p.value))
		#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
		LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
		position<-colnames(test)[k]
		position<-as.numeric(gsub("X","",position))
		pos<-c(pos,position)
	}
	high<-LOD[which.max(LOD)]
	highestlods<-c(highestlods,high)
	highsigs<-sig[which.max(sig)]
	highestsigs<-c(highestsigs,highsigs)
}


highestlods<-sort(highestlods, decreasing=TRUE)
fifthpercentileLOD<-highestlods[5]

highestsigs<-sort(highestsigs, decreasing=TRUE)
fifthpercentilesig<-highestsigs[5]

abline(h=fifthpercentilesig,lty="dashed")
abline(h=fifthpercentileLOD,lty="dashed")




DLA1_2<-cbind(test2,DLA1_2)
DLA_19<-cbind(test2,DLA_19)
DLA1_19<-rbind(DLA1_2,DLA_19)













#0 lesions removed
sig<-NULL
pos<-NULL
LOD<-NULL

for (i in 9:54)
	{
	Box<-t.test(Average.phenotype..0.leaves.removed.~test[,i],var.equal=T)
	sig<-c(sig,-log10(Box$p.value))
	#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
	LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
	position<-colnames(test)[i]
	position<-as.numeric(gsub("X","",position))
	pos<-c(pos,position)
	}
#Permutations
highestlods<-NULL

for(i in 1:100){
	sig<-NULL
	pos<-NULL
	LOD<-NULL
	test$random<-NULL
	random<-sample(Average.phenotype..0.leaves.removed.)
	test$random<-c(random)
	for(k in 9:54){
		Box<-t.test(random~test[,k],var.equal=T)
		sig<-c(sig,-log10(Box$p.value))
		#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
		LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
		position<-colnames(test)[k]
		position<-as.numeric(gsub("X","",position))
		pos<-c(pos,position)
	}
	high<-LOD[which.max(LOD)]
	highestlods<-c(highestlods,high)
}


wd=270*0.03937
ht=150*0.03937
ps=12
postscript(file="L73a0leavesremoved_5percentpermute.eps",width=wd,height=ht,pointsize=ps)





highestlods<-sort(highestlods, decreasing=TRUE)
fifthpercentileLOD<-highestlods[5]

highestlods<-c(highestlods,(LOD[which.max(LOD)])


plot(pos,sig,type="l",xlab="Physical position (Mbp) on 6HS",lwd=2,ylab="-Log10(P)",xaxt="n")
plot(pos,LOD,type="l",xlab="Physical position (Mb) on 6HS",lwd=2,ylab="LOD",xaxt="n")
abline(h=fifthpercentileLOD,lty="dashed")
#95% confidence intervals
abline(v=9300000)
abline(v=11900000)
#number of labels
Nt<-5
v1<-c(0:Nt)*(7e+07/Nt)
axis(side = 1, at= v1, labels =(v1/1000000))

#Calulate percentage of variance
str(aov(Average.phenotype..0.leaves.removed.~X9846099))
sumsquares<-summary(aov(Average.phenotype..0.leaves.removed.~X9846099))[[1]]$'Sum Sq'
variance_L73a<-((sumsquares[1]/(sumsquares[1]+sumsquares[2]))*100)
#52% variance
plot(pos,sig,type="l",xlab="Physical position (Mbp) on 6HS",lwd=2,ylab="-Log10(P)",xaxt="n")
#number of labels
Nt<-5
v1<-c(0:Nt)*(7e+07/Nt)
axis(side = 1, at= v1, labels =(v1/1000000))  




#Permutations for 271
highestlods<-NULL

for(i in 1:100){
	sig<-NULL
	pos<-NULL
	LOD<-NULL
	test$random<-NULL
	random<-sample(X271_Bianca)
	test$random<-c(random)
	for(k in 9:54){
		Box<-t.test(random~test[,k],var.equal=T)
		sig<-c(sig,-log10(Box$p.value))
		#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
		LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
		position<-colnames(test)[k]
		position<-as.numeric(gsub("X","",position))
		pos<-c(pos,position)
	}
	high<-LOD[which.max(LOD)]
	highestlods<-c(highestlods,high)
}

highestlods<-sort(highestlods, decreasing=TRUE)
fifthpercentileLOD<-highestlods[5]




sig<-NULL
pos<-NULL
LOD<-NULL
#271_Bianca
for (i in 9:54)
	{
	Box<-t.test(X271_Bianca~test[,i],var.equal=T)
	sig<-c(sig,-log10(Box$p.value))
	#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
	LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
	position<-colnames(test)[i]
	position<-as.numeric(gsub("X","",position))
	pos<-c(pos,position)
	}
#Calulate percentage of variance
str(aov(X271_Bianca~X9846099))
sumsquares<-summary(aov(X271_Bianca~X9846099))[[1]]$'Sum Sq'
variance_271<-((sumsquares[1]/(sumsquares[1]+sumsquares[2]))*100)
#75% variance

#Plot graphs

wd=270*0.03937
ht=150*0.03937
ps=12
postscript(file="Bianca_271_5percentpermute.eps",width=wd,height=ht,pointsize=ps)


plot(pos,sig,type="l",xlab="Physical position (Mbp) on 6HS",lwd=2,ylab="-Log10(P)",xaxt="n")
plot(pos,LOD,type="l",xlab="Physical position (Mbp) on 6HS",lwd=2,ylab="LOD",xaxt="n")
abline(h=fifthpercentileLOD,lty="dashed")

###LOD 95% flanking markers: 8757294,12052797
#number of labels
Nt<-14
v1<-c(0:Nt)*(7e+07/Nt)
axis(side = 1, at= v1, labels =(v1/1000000))
#Add confidence line
#Add lines for 95% confidence intervals
abline(h=4.5)
#95% confidence intervals
abline(v=9300000)
abline(v=11900000)
dev.off()



#For Lf_Bianca
sig<-NULL
pos<-NULL
LOD<-NULL

#LFL_Bianca
for (i in 9:54)
	{
	Box<-t.test(LfL.12.F_Bianca~test[,i],var.equal=T)
	sig<-c(sig,-log10(Box$p.value))
	#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
	LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
	position<-colnames(test)[i]
	position<-as.numeric(gsub("X","",position))
	pos<-c(pos,position)
	}

#Calulate percentage of variance

sumsquares<-summary(aov(LfL.12.F_Bianca~X9846099))[[1]]$'Sum Sq'
variance_LFL<-((sumsquares[1]/(sumsquares[1]+sumsquares[2]))*100)

#Permutations
highestlods<-NULL

for(i in 1:100){
	sig<-NULL
	pos<-NULL
	LOD<-NULL
	test$random<-NULL
	random<-sample(LfL.12.F_Bianca)
	test$random<-c(random)
	for(k in 9:54){
		Box<-t.test(random~test[,k],var.equal=T)
		sig<-c(sig,-log10(Box$p.value))
		#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
		LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
		position<-colnames(test)[k]
		position<-as.numeric(gsub("X","",position))
		pos<-c(pos,position)
	}
	high<-LOD[which.max(LOD)]
	highestlods<-c(highestlods,high)
}






highestlods<-sort(highestlods, decreasing=TRUE)
fifthpercentileLOD<-highestlods[5]

wd=240*0.03937
ht=100*0.03937
ps=12
postscript(file="LfL_new.eps",width=wd,height=ht,pointsize=ps)



plot(pos,LOD,type="l",xlab="Physical position (Mb) on 6HS",lwd=2,ylab="LOD",xaxt="n")
abline(h=fifthpercentileLOD,lty="dashed")

plot(pos,sig,type="l",xlab="Physical position (Mbp) on 6HS",lwd=2,ylab="-Log10(P)",xaxt="n")
#number of labels
Nt<-14
v1<-c(0:Nt)*(7e+07/Nt)
axis(side = 1, at= v1, labels =(v1/1000000))
abline(v=9300000)
abline(v=11900000)
dev.off()

for (i in 1:7)
	{
	LOD<-c(LOD,(qchisq((10^-lug[i,]),df=1,lower.tail=FALSE)/(2*log(10))))
	}



