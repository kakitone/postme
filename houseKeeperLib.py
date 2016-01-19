'''
Library functions for convenience
'''
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import groupby
import json
from operator import itemgetter


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
  
def extractMergeCandidStructToList(candidatesStructList):
    returnList = [] 
    for eachRecord in candidatesStructList:
        returnList += eachRecord[0]
    return returnList

def removeRedundantList(sourceList):
    targetList = []
    sourceList.sort()
    for key, items in groupby(sourceList):
        targetList.append(key)
    return targetList 
