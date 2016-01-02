'''
Library functions for convenience
'''
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def transformFileHeaders(folderName, inputFastaName, outputFastaName):
    targetToSourceNameDic = {}
    targetList = list(SeqIO.parse(folderName + filename, "fasta"))
    for eachRecord in targetList:
    	tmpName = "Seg" + str(len(targetToSourceNameDic))
        targetToSourceNameDic[tmpName]  = eachRecord.id
        eachRecord.id = tmpName

    SeqIO.write(targetList, folderName + outputFastaName , "fasta")
    return targetToSourceNameDic

def trailingFolderCorrection(folderName):
    return folderName if folderName[-1] == "/" else folderName + "/"
  