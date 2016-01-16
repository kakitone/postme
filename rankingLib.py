'''
Input : contigsFilename, readsFilename, scoreList, contigGapReadLookUpDic
Output : improved.fasta, scoreList.json
Algorithm : 
    1. Sort the score, cut off and form mergeList
    2. Form a contigGraph with mergeList
    3. Find the starter in graph
    4. Read the contigNodes concatenation, transform into readingList [1_p -> 2_d -> 3_p]
    5. Fill in the gap with the contigGapReadLookUpDic and readingList
'''
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from itertools import groupby
import json
from operator import itemgetter

import houseKeeperLib
import graphLib
import setCoverLib

class condenseEdgesOnlyGraph(graphLib.contigGraph):
    def formGraph(self, mergeList):
        for eachMergeItem in mergeList:
            infoList = eachMergeItem[0].split("~")
            bContigName, bSide, eContigName, eSide \
                = infoList[0][0:-2], infoList[0][-1], infoList[1][0:-2], infoList[1][-1]

            bContainer = self.returnContainer(bContigName, bSide)
            eContainer = self.returnContainer(eContigName, eSide)

            self.attachContainers(bContainer, eContainer, "dummyRead")

    def findStarterList(self):
        self.starterList = []
        for contigName in self.dicOfContigNodes:
            lenL = len(self.returnContainer(contigName, "L").connectedContigsDic)
            lenR = len(self.returnContainer(contigName, "R").connectedContigsDic)
            if lenL + lenR == 1 :
                self.starterList = \
                 self.starterList + [contigName + "_R"] if lenL == 1 else self.starterList + [contigName + "_L"]

    def formContigReadingList(self):
        readingList, checkDic = [],  {}

        for eachContigName in self.dicOfContigNodes:
            checkDic[eachContigName] = False

        for eachStarter in self.starterList:
            if checkDic[eachStarter[0:-2]] == False:
                path = self.returnCondensePath(eachStarter[0:-2], eachStarter[-1])
                readingList.append(path)
                for eachContigName in path:
                    checkDic[eachContigName[0:-2]] = True

        for eachContigName in checkDic:
            if checkDic[eachContigName] == False:
                readingList.append([eachContigName + "_p"])
                checkDic[eachContigName] = True

        return readingList

    def returnCondensePath(self, nextName, nextInSide):
        currentName = nextName
        currentSide = "R" if nextInSide == "L" else "L"
        returnName = currentName + "_p" if nextInSide == "L" else currentName + "_d" 

        nextDic = self.returnContainer(currentName, currentSide).connectedContigsDic
        
        if len(nextDic) == 0:
            return [returnName]
        else:
            containerName = nextDic.keys()[0]
            nextName, nextInSide = containerName[0:-2], containerName[-1]
            return [returnName] + self.returnCondensePath(nextName, nextInSide)

    def mergeListToReadingList(self, mergeList):
        self.formGraph(mergeList)
        self.findStarterList()
        return self.formContigReadingList()

def negateSide(side):
    return "d" if side == "p" else "p"

def readContigGraphToContigs(folderName,  contigsFilename, readsFilename, contigGapReadLookUpDic, readingList, outputContigsFilename):
    contigsDic, readsDic = {}, {}
    for record in SeqIO.parse(folderName + contigsFilename, "fasta"):
        contigsDic[record.id] = record

    for record in SeqIO.parse(folderName + readsFilename, "fasta"):
        readsDic[record.id] = record

    improvedList = []
    for eachPath in readingList:
        contigName, orientation = eachPath[0][0:-2], eachPath[0][-1]
        newContig = contigsDic[contigName].seq if orientation == "p" else contigsDic[contigName].reverse_complement().seq
        
        for i in range(1, len(eachPath)):
        
            prevContigName, prevSide, nextContigName, nextContigSide \
                = eachPath[i-1][0:-2], eachPath[i-1][-1], eachPath[i][0:-2], eachPath[i][-1]

            keyTrial = prevContigName + "_" + prevSide + "-" + nextContigName + "_" + nextContigSide
            keyTrialRev = nextContigName + "_" + negateSide(nextContigSide) + "-" + prevContigName + "_" + negateSide(prevSide) 

            if keyTrial in contigGapReadLookUpDic:
                readOrientation = "p"
                gapRecord = contigGapReadLookUpDic[keyTrial][0]
                sInfo = gapRecord[0][2:4]
                eInfo = gapRecord[1][2:4]
                s, e = max(sInfo), min(eInfo)
            elif keyTrialRev in contigGapReadLookUpDic:
                readOrientation = "d"
                gapRecord = contigGapReadLookUpDic[keyTrialRev][0]
                eInfo = [gapRecord[0][8] - gapRecord[0][2], gapRecord[0][8] - gapRecord[0][3]]
                sInfo  = [gapRecord[1][8] - gapRecord[1][2], gapRecord[1][8] - gapRecord[1][3]]
                s, e = max(sInfo), min(eInfo)
            else:
                print "Error: Gap not found "
                assert(False)

            readName = gapRecord[0][-1]
            
            if s < e:
                gapContent = readsDic[readName].seq[s:e-1] if readOrientation == "p" else readsDic[readName].reverse_complement().seq[s:e-1]
                nextContigResidual = contigsDic[nextContigName].seq if nextContigSide == "p" else contigsDic[nextContigName].reverse_complement().seq
            else:
                gapContent = ""
                nextContigResidual = contigsDic[nextContigName].seq[s-e + 1:] if nextContigSide == "p" else contigsDic[nextContigName].reverse_complement().seq[s-e+ 1:]

            newContig = newContig + gapContent + nextContigResidual
            
        improvedList.append(SeqRecord(Seq(str(newContig), generic_dna), description="", id="Segkk" + str(len(improvedList))))

    SeqIO.write(improvedList, folderName + outputContigsFilename , "fasta")

def cutOffToFormMergeList(scoreList, mScoreThres, conScoreThres ):
    mergeList = []
    for eachPotentialMerge in scoreList: 
        if eachPotentialMerge[-1] > mScoreThres or (eachPotentialMerge[-1] == mScoreThres and eachPotentialMerge[-2] >= conScoreThres ) :
            mergeList.append(eachPotentialMerge)
    return mergeList

def rankAndMerge(folderName, contigsNamesList, contigsFilename, readsFilename, scoreList, contigGapReadLookUpDic, mScoreThres, conScoreThres, scoreListOutputName, outputContigsFilename, dummyNodeDataRobot):
    
    scoreList.sort(key = itemgetter(-1, -2), reverse = True)
    houseKeeperLib.dumpDataToJson(folderName , scoreListOutputName, scoreList)
    mergeList = cutOffToFormMergeList(scoreList, mScoreThres, conScoreThres )

    nameList = dummyNodeDataRobot.getNameList()

    GCondensed = condenseEdgesOnlyGraph(nameList)

    readingListWithDummy = GCondensed.mergeListToReadingList(mergeList)

    readingList = dummyNodeDataRobot.mapReadingList(readingListWithDummy)
    
    readContigGraphToContigs(folderName,  contigsFilename, readsFilename, contigGapReadLookUpDic, readingList, outputContigsFilename)










    
