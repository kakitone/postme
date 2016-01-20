'''
Input : condenseCandidatesList, coverage/spanReadNum from contigGraph
Output: scoreList = [pairs, c_score, N_score]
Algorithm:
    1. Extract the information
    2. Perform CTest
    3. Output the score
'''
from Bio import SeqIO
from itertools import groupby
from operator import itemgetter
import scipy.stats


def assignLengthToContigNodes(G, folderName, contigFileName):
    for record in SeqIO.parse(folderName + contigFileName, "fasta"): 
        G.dicOfContigNodes[record.id].contigLength = len(record.seq)

def assignCoverageFromDataList(G, dataList,  folderName, contigFileName):
    contigReadCountDic = {}
    assignLengthToContigNodes(G, folderName, contigFileName)
    dataList.sort(key = itemgetter(-1, -3))
    
    for readName, mummerRecordList in groupby(dataList, itemgetter(-1)):
        topRecord = mummerRecordList.next()
        contigName = topRecord[-2]
        G.dicOfContigNodes[contigName].readToContigCount +=  1

def assignCoverageFromHeader(G, folderName, contigFileName, targetToSourceContigsNamesDic):
    assignLengthToContigNodes(G, folderName, contigFileName)
    for targetName in targetToSourceContigsNamesDic: 
        G.dicOfContigNodes[targetName].readToContigCount = float(targetToSourceContigsNamesDic[targetName].split("_")[5]) * G.dicOfContigNodes[targetName].contigLength/200

def calculateConfidenceScore(G, condenseCandidatesStructList):
    scoreList = []
    
    for eachcandidatestruct in condenseCandidatesStructList:
        tmpStruct = [[], eachcandidatestruct[1]]
        for eachcandidate in eachcandidatestruct[0]:
            infoList = eachcandidate.split("~") 
            contig1Name, contig2Name, mCount = infoList[0][0:-2], infoList[1][0:-2], int(infoList[2])
                    
            xTmp = G.dicOfContigNodes[contig1Name].readToContigCount
            nTmp = G.dicOfContigNodes[contig1Name].readToContigCount + G.dicOfContigNodes[contig2Name].readToContigCount
            pTmp = G.dicOfContigNodes[contig1Name].contigLength *1.0/ (G.dicOfContigNodes[contig1Name].contigLength + G.dicOfContigNodes[contig2Name].contigLength)
            
            mScore = mCount
            pvalue = scipy.stats.binom_test(xTmp, n=nTmp, p=pTmp)

            tmpStruct[0].append([eachcandidate, 1-pvalue, mScore])

        scoreList.append(tmpStruct)
        
    return scoreList

