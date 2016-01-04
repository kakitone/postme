'''
Library functions for convenience
'''
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json

def dumpDataToJson(folderName , filename, dataList) :
    with open(folderName + filename, 'w') as outfile:
        json.dump(dataList, outfile)

def readInJSON(folderName, filename):
    return json.load(open(folderName + filename, 'r'))
    
def transformFileHeaders(folderName, inputFastaName, outputFastaName, noAlignment):
    if noAlignment == False:
        targetToSourceNameDic = {}
        targetList = list(SeqIO.parse(folderName + inputFastaName, "fasta"))
        for eachRecord in targetList:
            tmpName = "Seg" + str(len(targetToSourceNameDic))
            targetToSourceNameDic[tmpName]  = eachRecord.id
            eachRecord.id = tmpName

        SeqIO.write(targetList, folderName + outputFastaName , "fasta")
        dumpDataToJson(folderName, "targetToSourceNameDic_" + inputFastaName + "_" + outputFastaName + ".json", targetToSourceNameDic)
        return targetToSourceNameDic
    else:
        return readInJSON(folderName, "targetToSourceNameDic_" + inputFastaName + "_" + outputFastaName + ".json")

def trailingFolderCorrection(folderName):
    return folderName if folderName[-1] == "/" else folderName + "/"
  