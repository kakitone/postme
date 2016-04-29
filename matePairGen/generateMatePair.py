from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

folderName = ""
filename = ""


def convertLR2MP(folderName, longReadFilename, pairEndOutputHeader, readLength, gapLength): 
	'''
	Goal : convert Long reads into matepair library. 
	Input :  longread.fasta
	Output : pair_1.fasta, pair_2.fasta
	'''
	pairEnd1List, pairEnd2List = [], []

	for record in SeqIO.parse(folderName + longReadFilename, "fasta"):
		if len(record.seq) >= 2*readLength + gapLength:
			record.seq = record.seq[0:2*readLength + gapLength] 
			beginSeq = SeqRecord(Seq(str(record.seq[0:readLength]), generic_dna), description="", id="read_" +  "%016d" % (len(pairEnd1List) + 1) + "/1" )
			beginSeq.letter_annotations["phred_quality"] = [40 for i in range(len(beginSeq.seq))]

			endSeq = SeqRecord(Seq(str(record.reverse_complement().seq[0:readLength]), generic_dna), description="", id="read_" +  "%016d" % (len(pairEnd2List) + 1) + "/2")
			endSeq.letter_annotations["phred_quality"] = [40 for i in range(len(endSeq.seq))]

			pairEnd1List.append(beginSeq)
			pairEnd2List.append(endSeq)


	SeqIO.write(pairEnd1List, folderName + pairEndOutputHeader + "_1.fastq" , "fastq")
	SeqIO.write(pairEnd2List, folderName + pairEndOutputHeader + "_2.fastq", "fastq")


def unitTest():
	folderName = "/tmp/"
	longReadFilename = "longread.fasta"
	pairEndOutputHeader = "pairend"
	readLength =  3

	SeqIO.write([SeqRecord(Seq("ACCTTTGGA", generic_dna), description="", id="segkk_0")], folderName + longReadFilename  , "fasta")
	
	### Test 1 
	gapLength = 3
	convertLR2MP(folderName, longReadFilename, pairEndOutputHeader, readLength, gapLength)
	
	mySeq = SeqIO.parse(folderName + pairEndOutputHeader + "_1.fastq", "fastq").next()
	str1 = str(mySeq.seq)
	qual1 = mySeq.letter_annotations["phred_quality"]
	

	mySeq = SeqIO.parse(folderName + pairEndOutputHeader + "_2.fastq", "fastq").next()
	str2 = str(mySeq.seq)
	qual2 = mySeq.letter_annotations["phred_quality"]

	assert(str1 == "ACC")
	assert(qual1 == [40 , 40, 40])
	assert(str2 == "TCC")
	assert(qual2 == [40 , 40, 40])

	### Test 2 
	gapLength = 4
	convertLR2MP(folderName, longReadFilename, pairEndOutputHeader, readLength, gapLength)

	assert(len(list(SeqIO.parse(folderName + pairEndOutputHeader + "_1.fastq", "fastq"))) == 0) 
	assert(len(list(SeqIO.parse(folderName + pairEndOutputHeader + "_2.fastq", "fastq"))) == 0) 


	### Test 3 
	gapLength = 2
	convertLR2MP(folderName, longReadFilename, pairEndOutputHeader, readLength, gapLength)
	
	mySeq = SeqIO.parse(folderName + pairEndOutputHeader + "_1.fastq", "fastq").next()
	str1 = str(mySeq.seq)
	qual1 = mySeq.letter_annotations["phred_quality"]
	

	mySeq = SeqIO.parse(folderName + pairEndOutputHeader + "_2.fastq", "fastq").next()
	str2 = str(mySeq.seq)
	qual2 = mySeq.letter_annotations["phred_quality"]
	
	assert(str1 == "ACC")
	assert(qual1 == [40 , 40, 40])
	assert(str2 == "CCA")
	assert(qual2 == [40 , 40, 40])
	
	print "Passing all test cases"

unitTest()



folderName = "/home/kakitone/Apr29-2016/datafolder3/"
longReadFilename = "EC.fasta"
readLength = 100

for gapLength in [600, 1200, 1800, 2400, 3000]:
	pairEndOutputHeader = "pairGap" + str(gapLength)
	convertLR2MP(folderName, longReadFilename, pairEndOutputHeader, readLength, gapLength)







