from itertools import groupby
from operator import itemgetter

import houseKeeperLib
import readConnectivityLib


class dummyNodeController(object):
    def __init__(self):
        self.name = ""
        self.realToDummyDic = {}
        self.dummyToRealDic = {}
        self.header = "Dummy"
        self.ctr = -1

    def addDummy(self, realName, toAddDummy):
        if toAddDummy or not realName in self.realToDummyDic:
            self.ctr += 1
            dummyName = self.header  + str(self.ctr)
            self.dummyToRealDic[dummyName] = realName
            self.realToDummyDic[realName] = self.realToDummyDic[realName] + [dummyName] \
                                            if realName in self.realToDummyDic else [ dummyName ]
            
            return dummyName 
        else:
            return self.realToDummyDic[realName] 

    def addAllContigs(self, contigsNamesList):
        for contigName in contigsNamesList:
            if not contigName in self.realToDummyDic:
                self.addDummy(contigName)

    def R2DLookUp(self, realName):
        return self.realToDummyDic[realName]

    def D2RLookUp(self, dummyName):
        return self.dummyToRealDic[dummyName]

    def getNameList(self):
        return [eachDummyName for eachDummyName in self.dummyToRealDic]

    def mapReadingList(self, readingListWithDummy):
        readingList = []
        for eachPath in readingListWithDummy:
            tmpPath = []
            for eachContigName in eachPath:
                tmpPath.append(self.D2RLookUp(eachContigName[0:-2]) + eachContigName[-2:]  )

            readingList.append(tmpPath)
        return readingList

    def createLocalR2DMapping(self, connectScoreList, toAddDummy):
        tmpR2DDic = {}
        tmpContigsList = [] 

        for eachRecord in connectScoreList:
            print connectScoreList
            infoList = eachRecord[0].split("~")
            tmpContigsList.append(infoList[0][0:-2])

        tmpContigsList.append(infoList[1][0:-2])
        
        for eachContigName in tmpContigsList:
            dummyName = self.addDummy(eachContigName, toAddDummy)
            tmpR2DDic[eachContigName] = dummyName

        return tmpR2DDic

def extendConnectivityFromReads(candidatesStructList, connectingReadsList, contigsNamesList):
    condenseCandidatesList = houseKeeperLib.extractMergeCandidStructToList(candidatesStructList)
    unUsedContigsDic = findUnUsedContigs(condenseCandidatesList, contigsNamesList)
    setCoverStructList = findSetCoverBaseLine(unUsedContigsDic, connectingReadsList)
    return setCoverStructList + candidatesStructList 

def findUnUsedContigs(condenseCandidatesList, contigsNamesList):
    # condenseCandidatesList == ['ContigDummyL_R~ContigDummyR_L~1']
    unUsedContigsDic = { contigName : True for contigName in contigsNamesList }

    for eachCondensePair in condenseCandidatesList:
        parsedList = eachCondensePair.split("~")
        for i in range(2):
            unUsedContigsDic[parsedList[i][0:-2]] = False 
    
    return unUsedContigsDic

def findSetCoverBaseLine(unUsedContigsDic, connectingReadsList):
    setCoverStructList, setRecords = [], transformConnectingReadsToSetStructure(connectingReadsList)

    setRecords.sort(key = itemgetter(2), reverse = True)
    for eachSetRecord in setRecords:
        allNotUsed = True
        for eachContigName in eachSetRecord[0]:
            if unUsedContigsDic[eachContigName] == False:
                allNotUsed = False

        if allNotUsed:
            setCoverStructList += [eachSetRecord[1], True]
            for eachContigName in eachSetRecord[0]:
                unUsedContigsDic[eachContigName] = False

    return setCoverStructList

def findSetCoverGreedy(unUsedContigsDic, connectingReadsList):
    setCoverList, setRecords = [], transformConnectingReadsToSetStructure(connectingReadsList)
    toCoverList, setsToSelect, currentCoveredList = [], [], []
    
    for eachContig in unUsedContigsDic:
        if unUsedContigsDic[eachContig]:
            toCoverList.append(eachContig)
            setsToSelect.append([ [eachContig] , [], -1 ])

    setsToSelect += setRecords

    toCoverSet = set(toCoverList)
    while len(toCoverSet) > 0:
        recordItem = findMostCostEffectiveSet(toCoverSet, setsToSelect) 
        setCoverList += recordItem[1]
        toCoverSet = toCoverSet.difference(set(recordItem[0])) 

    return setCoverList

def findMostCostEffectiveSet(toCoverSet, setsToSelect):
    recordItem, maxLength= [], -1 

    for eachRecord in setsToSelect:
        if len(set(eachRecord[0]).intersection(toCoverSet)) + penalizeSingleton(eachRecord) > maxLength:
            recordItem = eachRecord
            maxLength = len(set(eachRecord[0]).intersection(toCoverSet)) + penalizeSingleton(eachRecord)

    return recordItem 

def penalizeSingleton(record):
    return -0.5 if record[-1] < 0 else 0

def transformConnectingReadsToSetStructure(connectingReadsList):
    # Input :  connectingReadsList.append(['ReadDummy', 'B', 'ContigDummyB1', contigDummyBRecord1])
    # Output : [['Contig1', 'Contig2', 'Contig3'], [['Contig1_L', 'Contig2_R'], ['Contig2_L', 'Contig3_R']] ]
    spanReadsInfoList, linkedContigsInfoList = readConnectivityLib.extractInfoFromSpanRead(connectingReadsList)
    setStructures = []
    for eachlinkedInfo in linkedContigsInfoList:
        linkList = [[], [], -1]
        for eachContig in eachlinkedInfo[0]:
            linkList[0].append(eachContig[0:-2])

        for i in range(len(eachlinkedInfo[0]) - 1):
            linkList[1].append(convertPD2LR(eachlinkedInfo[0][i], "R") + "~" + convertPD2LR(eachlinkedInfo[0][i+1], "L") + "~" + "3")

        linkList[2] = len(eachlinkedInfo[0])
        setStructures.append(linkList)

    return setStructures

def convertPD2LR(pdContigName, side):
    if pdContigName[-1] == "p":
        return pdContigName[0:-1] + "R" if side == "R" else pdContigName[0:-1] + "L" 
    elif pdContigName[-1] == "d":
        return pdContigName[0:-1] + "L" if side == "R" else pdContigName[0:-1] + "R"
    else:
        print "Error : not pd"
        assert(False)

def assignRepeatedNodesToDummy(scoreStructList):
    scoreListWithDummy, dummyNodeDataRobot  = [], dummyNodeController()
    
    for eachScoreStruct in scoreStructList:
        connectScoreList = eachScoreStruct[0]
        toAddDummy = eachScoreStruct[1]
        tmpDummyDic = dummyNodeDataRobot.createLocalR2DMapping(connectScoreList, toAddDummy)
        
        for eachScoreItem in connectScoreList:
            infoList = eachScoreItem[0].split("~")
            name1 = tmpDummyDic[infoList[0][0:-2]] +  infoList[0][-2:]
            name2 = tmpDummyDic[infoList[1][0:-2]] +  infoList[1][-2:]
            mScore = infoList[2]
            newName = name1  + "~" + name2 + "~" + mScore 
            scoreListWithDummy.append([newName, eachScoreItem[1], eachScoreItem[2]])

    return scoreListWithDummy, dummyNodeDataRobot

