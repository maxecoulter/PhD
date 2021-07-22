test<-read.csv("Hets_data_for_R.csv")

library(lme4)

test2<-subset(test,Experiment="hetrozygote",drop=TRUE)

test3<-subset(test2,Genotype=="Steptoe",drop=TRUE)
test4<-subset(test2,Genotype=="Hetrozygote",drop=TRUE)
test5<-subset(test2,Genotype=="Morex",drop=TRUE)

test2<-rbind(test3,test4,test5)


#test2<-subset(subset(subset(test2,Genotype=="Steptoe"),Genotype=="Hetrozygote"),Genotype=="Morex")
#test2<-subset(test2,subset = Genotype %in% c(Steptoe,Hetrozygote,Morex))
testsubset<-test


#Columns: Family, Sub_family, Genotype, Experiment, Box, Position, Day_**, Greatest, Greatest_0_removed 
test.model<- lmer(Greatest_0_removed ~ Genotype + (1|Family/Sub_family/Line), data=test2)
test.model.null<- lmer(Greatest_0_removed ~ (1|Family/Sub_family/Line), data=test2)
sig<-anova(test.model,test.model.null)#This give significance

##Results from Greatest_0_removed, random model = f/sf/l: Steptoe/Morex pvalue = 0.0054,het/Morex pvalue = 0.017,Steptoe/het pvalue = 0.4944
##Results from Greatest_0_removed, random model = f/sf: Steptoe/Morex pvalue = 0.0002107,het/Morex pvalue = 0.002598,Steptoe/het pvalue = 0.4626

#Line has a significant effect in model - p = 0.0005208
#This is more significant than genotype - p = 0.006065


##LSD

library(predictmeans)
predictmeans(test.model,"Genotype",level=0.005)
#Result predicted means:
#Hetrozygote  8.4824
#Morex 10.5214
#Steptoe 7.9610
#Differences:Steptoe/Morex: 2.5604, Steptoe/Hetrozygote: 0.5214, Morex/Hetrozygote:2.039, Aveg LSD (0.05): 1.69934,Aveg LSD (0.01):2.272,Aveg LSD (0.005):2.49511


###############map lines###########################

test<-read.csv("map_forREML.csv")

#drop 16_4_5
#test<-subset(test,Line!="16_4 5",drop=TRUE)


X9846178
colnames(test)
#Marker association analysis with t.test (REML data)
sig<-NULL
pos<-NULL
LOD<-NULL


#Phenotype Column names from test[,10-16]

#

for (i in 1:11)
	{
	as<-subset(test,test[,i]=="a",drop=TRUE)
	bs<-subset(test,test[,i]=="b",drop=TRUE)
	test2<-rbind(as,bs)
	Box<-t.test(test2[,17]~test2[,i],var.equal=T)
	sig<-c(sig,-log10(Box$p.value))
	#Convert P value to lod score (qchisq returns inverse of probability of chisq distribution, the LRT))
	LOD<-c(LOD,(qchisq(p=(Box$p.value), df=1, lower.tail=FALSE)/(2*log(10))))
	position<-colnames(test2)[i]
	position<-as.numeric(gsub("X","",position))
	pos<-c(pos,position)
	}
plot(pos,sig,type="l")
plot(pos,LOD,type="l")
plot(pos,LOD,type="l",xlab="Physical position (Mbp) on 6HS",lwd=2,ylab="LOD",xaxt="n")
Nt<-5
v1<-c(0:Nt)*(7e+07/Nt)
axis(side = 1, at= v1, labels =(v1/1000000)) 



#as<-subset(test,test[,2]=="a",drop=TRUE)
#bs<-subset(test,test[,2]=="b",drop=TRUE)
#test2<-rbind(as,bs)

#Box<-t.test(test2[,22]~test2[,2],var.equal=T)

#ggplot
install.packages("ggplot2")
library(ggplot2)
#Change positions to megabases
pos<-pos/1000000
graph<-cbind(pos,LOD)
graph<-as.data.frame(graph)
require(scales)
P<-ggplot(data=graph,aes(x=pos, y=LOD, group=1)) + geom_line(size=1)+ xlab("Physical map position (Mbp)")+ylab("LOD") +
theme_bw()+ scale_x_continuous(labels = comma)+theme(axis.text=element_text(size=16),axis.title=element_text(size=16))




