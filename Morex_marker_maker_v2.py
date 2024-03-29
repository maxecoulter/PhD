"""Program for creating kasp primers from vcf.
Authors Max Coulter, Dr Matthew Moscou"""


#For creating KASP markers from a .vcf file, with Morex barley pseudomolecules reference genome (Mascher 2017). Primers are then tested for specificity using BLASTn.
#Program written for use on the linux cluster at the James Hutton Institute.
# In my project, the reference allele was from cv. Morex, while the snp allele was from cv. Steptoe.
#Program creates kasp primers, and is based on a previous script developed by Matthew Mouscou. 



#Program parameters
#python 3 path linux: /mnt/apps/python/3.5/bin/python3.5

#To run on linux:
#Usage:
#python morex_marker_makerv2.py <inputfilename>

#Input file: VCF column order (tab delimited): Chromosome,physical position,reference allele, snp allele. Gene reference must be in column 8
filt=2 #Filtering stingency (bases)


#Blast result filtering parameters



#Modules to import

import subprocess
import os
import optparse
from optparse import OptionParser 
import re
import math

myworkingdirectory=os.getcwd()

# import arguments and options
parser = OptionParser()
(options, args) = parser.parse_args()


# global variables
cmp_DNA = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'M':'K', 'R':'Y', 'S':'S', 'W':'W', 'Y':'R', 'K':'M', 'V':'B', 'H':'D', 'D':'H', 'B':'V', 'X':'X', 'N':'N'}
IUPAC = {'AG':'R', 'GA':'R', 'CT':'Y', 'TC':'Y', 'CA':'M', 'AC':'M', 'TG':'K', 'GT':'K', 'TA':'W', 'AT':'W', 'CG':'S', 'GC':'S'}
DNA = ['A', 'C', 'G', 'T']
IUPAC_dinucleotides = ['R', 'Y', 'M', 'K', 'W', 'S']
IUPAC_SNP = {'R':['A', 'G'], 'Y':['C', 'T'], 'M':['A', 'C'], 'K':['G', 'T'], 'W':['A', 'T'], 'S':['C', 'G']}
VIC = 'GAAGGTCGGAGTCAACGGATT'
FAM = 'GAAGGTGACCAAGTTCATGCT'

#Functions


def manifestmaker(index,indexlist,chromdic,genomedic,refdic,genedic,IUPAC,listofmanifests,referencesnps,output):
    print("processing snp "+str(index))
    chromnumber=chromdic[index]
    chrom=genomedic[chromnumber]
    index2=index-1
    assert(chrom[index2]==refdic[index])
    sequence=chrom[index2-70:index2+70]
    sequence=sequence.upper()
    #Add in other SNPs into the manifest sequence as 'N'
    snpindex=index
    listofclosesnps=[]
    for ind in indexlist:
        if ind>=(snpindex-70) and ind<=(snpindex+70):
            if ind!=snpindex:
                listofclosesnps.append(ind)
    for inde in listofclosesnps:
        manpos=(inde-snpindex)+70
        sequencelist=list(sequence)
        try:
            sequencelist[manpos]="N"
        except IndexError:
            pass
        sequence="".join(sequencelist)    
    sequencelist=list(sequence)
    IUPACcodeindex=refdic[index]+snpdic[index]
    sequencelist[70]=IUPAC[IUPACcodeindex]
    #sequencelist[70]="["+refdic[index]+"/"+snpdic[index]+"]"
    manifest="".join(sequencelist)
    listofmanifests.append(manifest)
    referencesnps.append(genedic[index]+"_"+str(index))
    print(genedic[index]+"_"+str(index))
    output.write(str(genedic[index])+"\t"+str(index)+"\t"+str(manifest)+"\n")
    manifests_references=[listofmanifests,referencesnps]
    return manifests_references


    



# Kasp marker maker FUNCTIONS
# reverse
# input : DNA sequence
# output : reverse of said DNA sequence
def reverse(orig_DNA):
	rev_DNA = ''
	for index in range(len(orig_DNA)):
		rev_DNA += orig_DNA[len(orig_DNA)-index-1]
	return rev_DNA

# reverse complement
# input : DNA sequence
# output : reverse complement of said DNA sequence
def reverse_complement(orig_DNA):
	rev_comp_DNA = ''
	for index in range(len(orig_DNA)):
		rev_comp_DNA += cmp_DNA[orig_DNA[len(orig_DNA)-index-1]]
	return rev_comp_DNA

# primer3 output
def primer3_output(output_stream, markerID, template, right_or_force):
    output_stream.write('SEQUENCE_ID=' + markerID + '\n')
    output_stream.write('SEQUENCE_TEMPLATE=' + template + '\n')
	
    if right_or_force == 'right':
        output_stream.write('PRIMER_PICK_LEFT_PRIMER=0' + '\n')
    elif right_or_force == 'force':
        output_stream.write('PRIMER_PICK_LEFT_PRIMER=1' + '\n')
        output_stream.write('SEQUENCE_FORCE_LEFT_START=0' + '\n')
        output_stream.write('SEQUENCE_FORCE_LEFT_END=19' + '\n')

    output_stream.write('PRIMER_PICK_INTERNAL_OLIGO=0' + '\n')
    output_stream.write('PRIMER_PICK_RIGHT_PRIMER=1' + '\n')
    output_stream.write('PRIMER_OPT_SIZE=18' + '\n')
    output_stream.write('PRIMER_MIN_SIZE=15' + '\n')
    output_stream.write('PRIMER_MAX_SIZE=21' + '\n')
    output_stream.write('PRIMER_MAX_NS_ACCEPTED=0' + '\n')
    output_stream.write('PRIMER_LIBERAL_BASE=1' + '\n')

    if right_or_force == 'right':
        output_stream.write('PRIMER_PRODUCT_SIZE_RANGE=21-100' + '\n')
    elif right_or_force == 'force':
        output_stream.write('PRIMER_PRODUCT_SIZE_RANGE=41-100' + '\n')
	
    output_stream.write('P3_FILE_FLAG=0' + '\n')
    output_stream.write('PRIMER_EXPLAIN_FLAG=0' + '\n')
    output_stream.write('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/mnt/shared/scratch/mb40521/software/primer3-2.3.7/src/primer3_config/' + '\n')
    #output_stream.write('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=C:/release-2.3.6/primer3_config/' + '\n')
    #output_stream.write('PRIMER_THERMODYNAMIC_PARAMETERS_PATH='+config_dir + '\n')
    output_stream.write('=' + '\n')

    return 1

##########################################################################################################################################################
#Upload morex chromosomes




#Open genome file
#Creates dictionary of chromosmes
#Memory hungry - linux cluster only!
print("Opening genome files...")
genome=[]
for i in range(1,8):
    chrom=open("chr"+str(i)+"H.fasta")
    chrom=chrom.read()
    chrom=chrom.replace("\n","")
    genome.append(chrom)

chrom=open("chrUn.fasta")
chrom=chrom.read()
chrom=chrom.replace("\n","")
genome.append(chrom)


genomeindex=[]
for i in range(1,8):
    name="chr"+str(i)+"H"
    genomeindex.append(name)

genomeindex.append("chrUn")

genomedic=dict(zip(genomeindex,genome))

#VCF input or FASTA input?

os.chdir(myworkingdirectory)

linecount=0
while linecount==0:
    for line in open(args[0]):
        if line.startswith(">"):
            filetype="FASTA"
        else:
            filetype="VCF"
        linecount=1

        


if filetype=="VCF":
    
    #Dictionaries
    indexlist=[]
    chromlist=[]
    referencelist=[]
    snplist=[]
    genelist=[]
    count=0
    #Parse vcf input
    print("Parsing vcf input")
    for line in open(args[0]):
        if count==0:
            count+=1
        else:
            line=line.replace("\n","")
            splits=line.split("\t")
            chromlist.append(splits[0])
            indexlist.append(int(splits[1]))
            referencelist.append(splits[2])
            snplist.append(splits[3])
            genelist.append(splits[7])
    print("Creating dictionaries")
    refdic=dict(zip(indexlist,referencelist))
    snpdic=dict(zip(indexlist,snplist))
    genedic=dict(zip(indexlist,genelist))
    chromdic=dict(zip(indexlist,chromlist))

    #VCF input:Create manifests with UIPAC code as SNP. Generate file

    output=open("List_of_manifests_and_positions.txt","w")
    output.write("Gene\tSNP position\tSNP manifest\n")
    print("Making manifests...")
    listofmanifests=[]
    referencesnps=[]
    for index in indexlist:
        manifests_references=manifestmaker(index,indexlist,chromdic,genomedic,refdic,genedic,IUPAC,listofmanifests,referencesnps,output)
        listofmanifests=manifests_references[0]
        referencesnps=manifests_references[1]
    output.close()

    manifestdic=dict(zip(indexlist,listofmanifests))

    #Output to fasta file
    output=open("Manifests.fasta","w")
    print("Outputting manifests to fasta...")
    for index in indexlist:
        print("processing snp "+str(index))
        #chromnumber=chromdic[index]
        #chrom=genomedic[chromnumber]
        #index2=index-1
        #assert(chrom[index2]==refdic[index])
        #sequence=chrom[index2-70:index2+70]
        output.write(">"+str(genedic[index])+"_"+str(index)+"\n"+str(manifestdic[index])+"\n")
    output.close()
    
    #dictionary of manifests
    manifestdic=dict(zip(referencesnps,listofmanifests))
    inputfile="Manifests.fasta"
    
else:
    inputfile="Inputfile.fasta"


#Generate KASP primers from either FASTA file or newly generated output from VCF input

# import arguments and options
#parser = OptionParser()
#(options, args) = parser.parse_args()


#Create string containing Primer3 config directory
#proc=subprocess.Popen(['where', 'primer3_core'],stdout=subprocess.PIPE)
#config_dir = str (proc.stdout.read())
#config_dir = config_dir[2:len(config_dir)-21] + 'primer3_config\\'

# dictionaries
marker_sequence = {}		# marker -> sequence
marker_SNP_positions = {}	# marker -> [SNP positions, ...]
SNP_marker_information = {}	# SNP_ID -> [marker, SNP position, {F:[s1, s2, r1], R:[s1, s2, r1]}]

print("Creating kasp primers...")
# STEP 1. Import IUPAC sequences
gene_sequence = open(inputfile, 'r')
#gene_sequence = open(args[0], 'r')
marker_order = []


line = gene_sequence.readline()

while line:
        
	line = ''.join(str.split(line))
	if len(line) > 0:
		if line[0] == '>':
			marker = line[1:]
			marker_sequence[marker] = ''
			marker_order.append(marker)
		else:
			#marker_sequence[marker] += str.replace(line, '-', 'N')
			marker_sequence[marker] = line.replace('-', 'N')
                        

	line = gene_sequence.readline()
	

gene_sequence.close()

# STEP 2. Identify positions of every SNP or ambiguous sequence (if they are more than 20 bases from the start or end of the sequence)
for marker in marker_sequence.keys():
	marker_SNP_positions[marker] = []

	for SNP in IUPAC_dinucleotides:
		for position in re.finditer(SNP, marker_sequence[marker]):
			if position.start() >= 20 and position.start() <= (len(marker_sequence[marker]) - 20):
				marker_SNP_positions[marker].append(position.start())

# STEP 3. Output sequencing primer to primer3
print("Creating primer3 input files...")
right_primers = open('right_primer_primer3_input.txt', 'w')
force_left_primers = open('force_left_primer_primer3_input.txt', 'w')

for marker in marker_SNP_positions.keys():
	sequence_RC = reverse_complement(marker_sequence[marker])
	
	for SNP in marker_SNP_positions[marker]:
		SNP_marker_information[marker + '_' + str(SNP)] = [marker, SNP, {'F':['', '', '', False], 'R':['', '', '', False]}, 0]
		SNP_RC = len(marker_sequence[marker]) - (SNP + 1)
		
		# forward sequence
		SNP_marker_information[marker + '_' + str(SNP)][2]['F'][0] = marker_sequence[marker][SNP - 19:SNP] + IUPAC_SNP[marker_sequence[marker][SNP]][0]
		SNP_marker_information[marker + '_' + str(SNP)][2]['F'][1] = marker_sequence[marker][SNP - 19:SNP] + IUPAC_SNP[marker_sequence[marker][SNP]][1]

		primer3_output(right_primers, marker + '_' + str(SNP) + '_F_r1', marker_sequence[marker][SNP + 1:SNP + 101], 'right')
		primer3_output(force_left_primers, marker + '_' + str(SNP) + '_F_r1', marker_sequence[marker][SNP - 19:SNP] + IUPAC_SNP[marker_sequence[marker][SNP]][0] + marker_sequence[marker][SNP + 1:SNP + 101], 'force')

		# reverse sequence
		SNP_marker_information[marker + '_' + str(SNP)][2]['R'][0] = sequence_RC[SNP_RC - 19:SNP_RC] + IUPAC_SNP[sequence_RC[SNP_RC]][0]
		SNP_marker_information[marker + '_' + str(SNP)][2]['R'][1] = sequence_RC[SNP_RC - 19:SNP_RC] + IUPAC_SNP[sequence_RC[SNP_RC]][1]

		primer3_output(right_primers, marker + '_' + str(SNP) + '_R_r1', sequence_RC[SNP_RC + 1:SNP_RC + 101], 'right')
		primer3_output(force_left_primers, marker + '_' + str(SNP) + '_R_r1', sequence_RC[SNP_RC - 19:SNP_RC] + IUPAC_SNP[sequence_RC[SNP_RC]][0] + sequence_RC[SNP_RC + 1:SNP_RC + 101], 'force')

right_primers.close()
force_left_primers.close()


# STEP 4. Run primer3
print("Running primer3...")
subprocess.getstatusoutput('/mnt/shared/scratch/mb40521/software/primer3-2.3.7/src/primer3_core < right_primer_primer3_input.txt > right_primer_primer3_output.txt')

subprocess.getstatusoutput('/mnt/shared/scratch/mb40521/software/primer3-2.3.7/src/primer3_core < force_left_primer_primer3_input.txt > force_left_primer_primer3_output.txt')

#subprocess.getstatusoutput('primer3_core < right_primer_primer3_input.txt > right_primer_primer3_output.txt')
#subprocess.getstatusoutput('primer3_core < force_left_primer_primer3_input.txt > force_left_primer_primer3_output.txt')

# STEP 5. Parse primer3 output
print("Parsing primer3 output")
force_left_primer_output = open('force_left_primer_primer3_output.txt', 'r')
right_primer_output = open('right_primer_primer3_output.txt', 'r')

# PART 5A. Parse primer3 output for forced left primer
line = force_left_primer_output.readline()

while line:
	line = line.replace('\n', '')
	sline = line.split('=')


	if sline[0] == 'SEQUENCE_ID':
		SNP = sline[1][:len(sline[1]) - 5]
		strand = sline[1][len(sline[1]) - 4]
		primer = sline[1][(len(sline[1]) - 2):len(sline[1])]

	if sline[0] == 'PRIMER_RIGHT_0_SEQUENCE':
		SNP_marker_information[SNP][2][strand][2] = sline[1]
		SNP_marker_information[SNP][2][strand][3] = True

	if sline[0] == 'PRIMER_PAIR_0_PRODUCT_SIZE':
		SNP_marker_information[SNP][3] = int(sline[1])

	line = force_left_primer_output.readline()

force_left_primer_output.close()

# PART 5B. Parse primer3 output for right primer design only
line = right_primer_output.readline()

while line:
	line = line.replace('\n', '')
	sline = line.split('=')

	if sline[0] == 'SEQUENCE_ID':
		SNP = sline[1][:len(sline[1]) - 5]
		strand = sline[1][len(sline[1]) - 4]
		primer = sline[1][(len(sline[1]) - 2):len(sline[1])]

	if not SNP_marker_information[SNP][2][strand][3]: 
		if sline[0] == 'PRIMER_RIGHT_0_SEQUENCE':
			SNP_marker_information[SNP][2][strand][2] = sline[1]
			SNP_marker_information[SNP][2][strand][3] = False

		if sline[0] == 'PRIMER_RIGHT_0':
			primer_start, primer_length = sline[1].split(',')
			SNP_marker_information[SNP][3] = int(primer_start) + int(primer_length) + 20

	line = right_primer_output.readline()

right_primer_output.close()

# STEP 6. Export primer designs
#	Overall: In lieu of designing logic based selection criteria, I will output a single file to allow human curation
#		A. 1 primer design per gene
#		B. 1 primer design per SNP
#		C. All primer designs (including forward and reverse)
#	Main function: Remove all SNPs with any IUPAC sequence

#primer_design_file = open(args[1], 'w')
print("Outputting primer designs...")
primer_design_file=open("Designed_primers.txt",'w')

primer_design_file.write('marker')
primer_design_file.write('\t' + 'template')
primer_design_file.write('\t' + 'SNP_position')
primer_design_file.write('\t' + 'strand')
primer_design_file.write('\t' + 'product_size')
primer_design_file.write('\t' + 'quality')
primer_design_file.write('\t' + 'SNP_primer_1_name')
primer_design_file.write('\t' + 'SNP_primer_1_sequence')
primer_design_file.write('\t' + 'SNP_primer_2_name')
primer_design_file.write('\t' + 'SNP_primer_2_sequence')
primer_design_file.write('\t' + 'reverse_primer_1_name')
primer_design_file.write('\t' + 'reverse_primer_1_sequence')
primer_design_file.write('\n')

for marker in marker_order:
	for SNP in marker_SNP_positions[marker]:
		for strand in ['F', 'R']:
			primer_QC = True

			for primer_index in range(3):
				if len(set(SNP_marker_information[marker + '_' + str(SNP)][2][strand][primer_index]) - set(DNA)) > 0:
					primer_QC = False
				if len(SNP_marker_information[marker + '_' + str(SNP)][2][strand][primer_index]) == 0:
					primer_QC = False

			if primer_QC:
				marker_ID = marker + '_' + str(SNP) + '_' + strand
				primer_design_file.write(marker_ID)
				primer_design_file.write('\t' + marker)
				primer_design_file.write('\t' + str(SNP))
				primer_design_file.write('\t' + strand)
				primer_design_file.write('\t' + str(SNP_marker_information[marker + '_' + str(SNP)][3]))
				primer_design_file.write('\t' + str(int(SNP_marker_information[marker + '_' + str(SNP)][2][strand][3])))
				primer_design_file.write('\t' + marker_ID + '_' + 's1')
				primer_design_file.write('\t' + VIC + SNP_marker_information[marker + '_' + str(SNP)][2][strand][0])
				primer_design_file.write('\t' + marker_ID + '_' + 's2')
				primer_design_file.write('\t' + FAM + SNP_marker_information[marker + '_' + str(SNP)][2][strand][1])
				primer_design_file.write('\t' + marker_ID + '_' + 'r1')
				primer_design_file.write('\t' + SNP_marker_information[marker + '_' + str(SNP)][2][strand][2])
				primer_design_file.write('\n')

primer_design_file.close()

###########################################################################################################################


#Create fasta file of primer output

print("Making fasta file of primers...")

output=open("Primersforblast.fasta","w")
count=0
primersequences=[]
shortprimersequences=[] #vic/fam sequence removed
for line in open("Designed_primers.txt"):
    if count==0:
        count+=1
    else:
        linesplit=line.split("\t")
        output.write(">"+linesplit[6]+"\n")
        output.write(linesplit[7][21:]+"\n") #Remove vic/fam sequence
        primersequences.append(linesplit[7])
        shortprimersequences.append(linesplit[7][21:])
        output.write(">"+linesplit[8]+"\n")
        output.write(linesplit[9][21:]+"\n")
        primersequences.append(linesplit[9])
        shortprimersequences.append(linesplit[9][21:])
        output.write(">"+linesplit[10]+"\n")
        output.write(linesplit[11])
        primersequences.append(linesplit[11].replace("\n",""))
        shortprimersequences.append(linesplit[11].replace("\n",""))
output.close()

##########################################################################################################################

#Primer output BLAST against v4 Chromosomes
print("BLASTing primers...")

os.system('blastn -query Primersforblast.fasta -outfmt 10 -subject 150831_barley_pseudomolecules.fasta -task blastn-short -out BLASTresults.txt')


print("Filtering blast results...")
   

             
#Parse BLAST results



count=0
indexlist=[]
chromlist=[]
percenthitlist=[]
noofhitslist=[]
snpreferencelist=[]
querystart=[]
queryend=[]
for line in open("BLASTresults.txt"):
    count=count+1
    splits=line.split(",")
    indexlist.append(count)
    snpreferencelist.append(splits[0])
    chromlist.append(splits[1])
    percenthitlist.append(float(splits[2]))
    noofhitslist.append(int(splits[3]))
    querystart.append(int(splits[6]))
    queryend.append(int(splits[7]))
                    


referencedic=dict(zip(indexlist,snpreferencelist))
chromdic=dict(zip(indexlist,chromlist))
percenthitdic=dict(zip(indexlist,percenthitlist))
noofhitsdic=dict(zip(indexlist,noofhitslist))
querystartdic=dict(zip(indexlist,querystart))
queryenddic=dict(zip(indexlist,queryend))                    

#Create list of snp references without duplicates
references=[]
for snp in snpreferencelist:
    if snp in references:
        pass
    else:
        references.append(snp)

primerdic=dict(zip(references,primersequences))
shortprimerdic=dict(zip(references,shortprimersequences))

#primers either 41 or 18 long
#Filter Primer sets
output=open("Marker_maker_primers_filtered.txt","w")
output.write("Primer set\tS1 primer\tNumber of hits from best alternative blast hit\tFiltering test\tS2 primer\tNumber of hits from best alternative blast hit\tFiltering test\tr1 primer\tNumber of hits from best alternative blast hit\tFiltering test\n")
badsnp=0
listofsnpfinfo=[]
count=0
filt-=1
for snp in references:
    count+=1
    primerlength=len(shortprimerdic[snp])
    print(str(snp)+" length: "+str(primerlength))
    primerpercentthresh=((primerlength-filt)/primerlength)*100
    primernohitsthreshold=float(primerlength-filt)
    qualityscores=[]
    print("processing snp "+str(snp))


    for index in indexlist:
        if referencedic[index]==snp:
            totalhits=(noofhitsdic[index]*percenthitdic[index])/100
            totalhitsi=int(totalhits)
            difference=totalhits-totalhitsi
            if difference>=0.5:
                totalhits=math.ceil(totalhits)
            else:
                totalhits=totalhitsi
            qualityscores.append(totalhits)
    qualityscores.sort()
    qualityscore=qualityscores[-2]
    if qualityscore>=primernohitsthreshold:
        badsnp=1

    
        
        

                            
    if badsnp==0:
        snpinfo="Primer ok"
    else:
        snpinfo="Primer failed filtering test"
        badsnp=0


    if count==1:
        output.write(str(snp.replace(str(snp[-2]+snp[-1]),""))+"\t"+str(primerdic[snp])+"\t"+str(qualityscore)+"/"+str(primerlength)+"\t"+snpinfo+"\t")
    elif count<=2:
        output.write(str(primerdic[snp])+"\t"+str(qualityscore)+"/"+str(primerlength)+"\t"+snpinfo+"\t")
    else:
        output.write(str(primerdic[snp])+"\t"+str(qualityscore)+"/"+str(primerlength)+"\t"+snpinfo+"\n")
        count=0
output.close()       
