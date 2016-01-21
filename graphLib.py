'''
Input :  spanReadsList
Output: contigGraph with contigNode
Algorithm: 
    1. Form the basic data structure
    2. Establish connectivity
'''
import houseKeeperLib

class contigGraph(object):
    def __init__(self, contigsNamesList):
        self.dicOfContigNodes = {}
        for contigName in contigsNamesList:
            self.dicOfContigNodes[contigName] = contigNode(contigName) 

    def returnContainer(self, contigName, side):
        assert(side in ["L", "R"])
        if side == "L":
            return self.dicOfContigNodes[contigName].leftEndContainer
        elif side == "R":
            return self.dicOfContigNodes[contigName].rightEndContainer

    def attachContainers(self, bContainer, eContainer, readName):
        bContainer.addNameToContainer(eContainer.myName, readName)
        eContainer.addNameToContainer(bContainer.myName, readName)

    def spanReadsToEdges(self, spanReadsList):
        for eachSpanReadRecord in spanReadsList:
            self.addSpanReadRecordToEdge(eachSpanReadRecord)

    def addSpanReadRecordToEdge(self, eachSpanReadRecord):
        readName = eachSpanReadRecord[-1]

        bContigName, bParity, eContigName, eParity \
            = eachSpanReadRecord[0][0:-2], eachSpanReadRecord[0][-1], eachSpanReadRecord[1][0:-2], eachSpanReadRecord[1][-1]

        bContainer = self.returnContainer(bContigName, "L")  if bParity == "d" else  self.returnContainer(bContigName, "R")
        eContainer = self.returnContainer(eContigName, "R")  if eParity == "d" else  self.returnContainer(eContigName, "L")
        
        self.attachContainers(bContainer, eContainer, readName)
        
    def findCondenseCandidatesList(self):
        condenseCandidatesList = []
        for contigName in self.dicOfContigNodes:
            potentialCondenseCandidates = self.dicOfContigNodes[contigName].findPotentialCondenseCandidateAtNode()
            condenseCandidatesList += self.confirmPotentialCondenseCandidates(potentialCondenseCandidates)

        return houseKeeperLib.removeRedundantList(condenseCandidatesList)
        
    def confirmPotentialCondenseCandidates(self, potentialCondenseCandidates):
        myList = []

        for eachPotential in potentialCondenseCandidates:
            myName, toCheckName = eachPotential[0], eachPotential[1]
    
            nameList = [myName, toCheckName]
            nameList.sort()
            returnName = nameList[0] + "~" + nameList[1] 
            
            contigName, side = toCheckName[0:-2], toCheckName[-1]
            
            if side == "L":
                if len(self.dicOfContigNodes[contigName].leftEndContainer.connectedContigsDic) == 1:
                    mNum = len(self.dicOfContigNodes[contigName].leftEndContainer.connectedContigsDic.itervalues().next())
                    myList.append([[returnName + "~" + str(mNum)], False]) 
            else:
                if len(self.dicOfContigNodes[contigName].rightEndContainer.connectedContigsDic) == 1:
                    mNum = len(self.dicOfContigNodes[contigName].rightEndContainer.connectedContigsDic.itervalues().next())
                    myList.append([[returnName + "~" + str(mNum)], False])   
                    
        return myList
        
    def findEdgeMultiplicity(self):
        multiplicityDic =  {}
        for contigName in self.dicOfContigNodes:
            myContainer = self.dicOfContigNodes[contigName].leftEndContainer
            for neighborName in myContainer.connectedContigsDic:
                multiplicityDic[myContainer.myName + "," + neighborName] = len(myContainer.connectedContigsDic[neighborName])

            myContainer = self.dicOfContigNodes[contigName].rightEndContainer
            for neighborName in myContainer.connectedContigsDic:
                multiplicityDic[myContainer.myName + "," + neighborName] = len(myContainer.connectedContigsDic[neighborName])

        return multiplicityDic

class contigNode(object):
    def __init__(self, contigName):
        self.name = contigName
        self.contigLength = 0
        self.readToContigCount = 0 
        self.leftEndContainer = endPtContainer(contigName, "L")
        self.rightEndContainer = endPtContainer(contigName, "R")

    def findPotentialCondenseCandidateAtNode(self):
        return self.findPotentialCondenseCandidateAtContainer(self.leftEndContainer) \
                + self.findPotentialCondenseCandidateAtContainer(self.rightEndContainer)

    def findPotentialCondenseCandidateAtContainer(self, targetEndPtContainer):
        return targetEndPtContainer.getCondenseCandidateList()

class endPtContainer(object):
    def __init__(self, contigName, side):
        self.myName = contigName + "_" + side
        self.connectedContigsDic = {}

    def addNameToContainer(self, containerName, readName):
        self.connectedContigsDic[containerName] = \
            self.connectedContigsDic[containerName] + [readName] if containerName in self.connectedContigsDic else [readName] 

    def getCondenseCandidateList(self):
        if len(self.connectedContigsDic) != 1 :
            return []
        else:
            tmpList = [self.myName, self.connectedContigsDic.keys()[0]]
            return [ tmpList ]

def formContigGraph(spanReadsList, contigsNamesList):
    G = contigGraph(contigsNamesList)
    G.spanReadsToEdges(spanReadsList)
    return G