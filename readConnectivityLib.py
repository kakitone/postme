'''
Input : mummerDataList
Output : spanReadsList, gapLookUpDic
Algorithms: 
    1. Group by readName
    2. Identify L/R attaching contigs
    3. Decide on whether it is confirming reads
    4. Form spanReadsList, contigGapReadLookUpDic
'''

from itertools import groupby
from operator import itemgetter

def decideMatch(recordList, overlapThres, matchThres):
    #381     6380  |        1     5994  |     6000     5994  |    97.92  |  5000000     5994  | Segkk0  Segkk6408
    isMatch, associatedRecord =  "N", []
    for record in recordList:
        contigStart, contigEnd, readStart, readEnd  = min(record[0:2]), max(record[0:2]), min(record[2:4]), max(record[2:4])
        contigLen, readLen = record[7], record[8]
        matchLen = (record[4] + record[5])/2
        if readStart < overlapThres and contigEnd > contigLen - overlapThres and matchLen > matchThres:
            isMatch, associatedRecord = "L", record
            break

        if contigStart < overlapThres and readEnd > readLen - overlapThres and matchLen > matchThres:
            isMatch, associatedRecord = "R", record
            break

    return isMatch, associatedRecord

def findConnectingReadsList(mummerDataList):
    mummerDataList.sort(key=itemgetter(-2, -1))
    overlapThres, matchThres  = 30, 100

    connectingReadsList = []
    for key, overlapList in groupby(mummerDataList, itemgetter(-2, -1)):
        contigName, readName = key[0], key[1]
        isMatch, associatedRecord = decideMatch(overlapList, overlapThres, matchThres)
        if not isMatch == "N":
            connectingReadsList.append([readName, isMatch, contigName, associatedRecord])

    return connectingReadsList

def isReverse(record):
    checkIsReverse = True if record[2] > record[3] else False
    contigName = record[-2] + "_d" if checkIsReverse else record[-2] + "_p"
    return contigName

def formatSpanReadReturn(recordList):
    spanReadRecord, gapRecord = [], []
    
    readName = recordList[0][0]
    contig0, side0,  record0 = recordList[0][2], recordList[0][1], recordList[0][3] 
    contig1, side1,  record1 = recordList[1][2], recordList[1][1], recordList[1][3]
    name0, name1 = isReverse(record0), isReverse(record1)
    
    spanReadRecord = [name0, name1, readName] if side0 == "L" else [name1, name0, readName]
    gapRecordKey = name0 + "-" + name1 if side0 == "L" else name1 + "-" + name0
    gapRecordValue = [record0, record1] if side0 == "L" else [record1, record0]

    return spanReadRecord, gapRecordKey, gapRecordValue

def findSpanReadsList(connectingReadsList):
    spanReadsList, contigGapReadLookUpDic = [], {}
    # ['Segkk10070', 'L', 'Segkk0', [2496339, 2500001, 1, 3667, 3663, 3667, 98.11, 2500001, 6002, 'Segkk0', 'Segkk10070']] 265
    connectingReadsList.sort(key = itemgetter(0))
    for readName, overlappingContigsRecord in groupby(connectingReadsList, itemgetter(0)):
        recordList = list(overlappingContigsRecord)
        if len(recordList) == 2: 
            if recordList[0][1] != recordList[1][1]:
                spanReadRecord, gapRecordKey, gapRecordValue = formatSpanReadReturn(recordList)
                
                spanReadsList.append(spanReadRecord)

                contigGapReadLookUpDic[gapRecordKey] = \
                    contigGapReadLookUpDic[gapRecordKey] + [gapRecordValue]  if gapRecordKey in contigGapReadLookUpDic else  [gapRecordValue] 

    return spanReadsList, contigGapReadLookUpDic
