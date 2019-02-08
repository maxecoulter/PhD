 -*- coding: utf-8 -*-
#"""
#Created on Mon Jan 23 14:25:07 2017

#@author: MC42302
#"""

#Open multiple Fastq files of sequence

#20 samples, 2 fastq files per sample
#sample name:MC_01_1P.fastq - first number 01-20, second number 1P or 2P

import os

os.chdir('/mnt/shared/scratch/mb40521/201608_RhyncResCIRC/RNAseq/trimmedReads')

samplelist=[]
readpairs=[1,2]
for i in range(1,21):
    if len(str(i))==1:
        p="0"+str(i)
        samplelist.append(p)
    else:
        samplelist.append(i)



#Histogram has 76 classes
index=[]
for i in range(2,79):
    index.append(i)

listofGstrings=[]
for i in range(2,79):
    listofGstrings.append(0)
Gdic=dict(zip(index,listofGstrings))





        
#for each sequence of 4 lines count number of Gs
samplelist=["08"] #test
readpairs=[1,2]  #test
readcount=0
for sample in samplelist:
    for read in readpairs:
        Gcount=0
        print("C count")
        print("MC_"+str(sample)+"_"+str(read)+"P.fastq")
        linecount=0
        for line in open("MC_"+str(sample)+"_"+str(read)+"P.fastq"):
        #for line in open("MC_01_1Phead2.fastq"):  #test
            if line.startswith("@"):
                linecount=2
            elif linecount==2:
                Gcount=0
                #G string counting algorithm
                for li in line:
                    if li=="C":
                        Gcount+=1
                    elif Gcount>=2:
                        if li=="C":
                            Gcount+=1
                        else:
                            Gdic[Gcount]+=1
                            Gcount=0
                    else:
                        Gcount=0
                if Gcount!=0:
                    Gdic[Gcount]+=1
                    Gcount=0
                linecount+=1
                        
            else:
                linecount+=1










print("Histogram output:")

for ind in index:
    print("Number of G strings with "+str(ind)+"Gs: "+str(Gdic[ind]))






            
histoutput=open("Histogram_of_G08s","w")
histoutput.write("Histogram output:\n")
for ind in index:
    histoutput.write("Number of G strings with "+str(ind)+"Gs: "+str(Gdic[ind])+"\n")

histoutput.close()
    



                
#Once all files done, make histogram of list of numbers of Gs
