# -*- coding: utf-8 -*-
"""


@author: MC42302
"""

#Get rid of reads with G strings

import os
threshold=21

os.chdir('/mnt/shared/scratch/mb40521/201608_RhyncResCIRC/RNAseq/mappingPseudomols')

#line by line
samplelist=[]
readpairs=[1,2]
for i in range(1,21):
    if len(str(i))==1:
        p="0"+str(i)
        samplelist.append(p)
    else:
        samplelist.append(i)

#samplelist=["01"]
#readpairs=[1,2]

readcounts=0
badcounts=0

##MC_01/MC_01Unmapped.out.mate1
##C_01Unmapped.out.mate2

for sample in samplelist:
    for read in readpairs:
        Gcount=0
        print("MC_"+str(sample)+"_"+str(read)+"P.fastq")
        #outfile=open("cleanedUnmappedReads/MC_"+str(sample)+"_"+str(read)+"P_Gsremoved.fastq","w")
        #outfile=open("MC_01_1Phead2Gs_removed.fastq","w") #test
        outputfile2=open("cleanedUnmappedReads/G_string_removerpythonout.txt","a")
        #outfilebad=open("cleanedUnmappedReads/removed_G_reads.fastq","w")
        readings=""
        stop=0
        linecount=0
        badcount=0
        readcount=0
        for line in open("MC_"+str(sample)+"/"+"MC_"+str(sample)+"Unmapped.out.mate"+str(read)):
        #for line in open("MC_01_1Phead2.fastq"): #test
            #print repr(str(line))
            linecount+=1
            readings+=line
            #print(str(linecount))
            if linecount==2:
                Gcounts=[]
                Gcount=0
                #G string counting algorithm
                for li in line:
                    if li=="G":
                        Gcount+=1
                    elif Gcount>=2:
                        if li=="G":
                            Gcount+=1
                        else:
                            Gcounts.append(Gcount)
                            Gcount=0
                    else:
                        Gcount=0
                if Gcount!=0:
                    Gcounts.append(Gcount)
                    Gcount=0                
            elif linecount==4:
                readcount+=1
                readcounts+=1
                Gcounts.append(0)
                if max(Gcounts)>=threshold:
                    #outfilebad.write(str(readings))
                    badcount+=1
                    badcounts+=1
                    readings=""
                else:
                    #outfile.write(str(readings))
                    readings=""
                linecount=0
            else:
                pass
        
        percentage=(float(badcount)/float(readcount))*100
        outputfile2.write("MC_"+str(sample)+"_"+str(read)+"P.fastq percentage reads removed\n"+str(percentage)+"\n")
        #outfile.close()

outputfile2.write("Total reads: "+str(readcounts)+"\n"+"Total reads removed: "+str(badcounts)+"\n"+"Total percentage of reads removed: "+str((float(badcounts)/float(readcounts))*100))
outputfile2.close()
#outfilebad.close()
