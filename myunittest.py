# Example usage : python myunittest.py -a  /Users/kakitlam/Desktop/allinone/experimentBench/MUMmer3.23/ -f /Users/kakitlam/Desktop/baseline/tmp/

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import os
import unittest

import alignmentLib
import cTestLib
import graphLib
import rankingLib
import readConnectivityLib
import setCoverLib

class baselineAlgoTest(unittest.TestCase):
    
    def setUp(self):
        print "Set up : Started"
        
        self.folderName = "/tmp/testdir/"
        self.mummerLink = "/tmp/MUMmer3.23/"
        self.readsFilename = "LR.fasta"
        self.contigsFilename = "mFixed2.fasta"
        self.outputHeader = "outputHeader"
        self.splitNum = 1
        self.parallelNum = 1
        self.debug = False
        
        os.system("mkdir "+ self.folderName)
        
        print "Set up  : Done "
    
    def createSimpleFasta(self):
        contigContent = "AGAATCAATAAAATGTTTGCAGACATCATAAGAGATCTGCCGCGCCTATACACGATGCCTTGATGATAATGTGCAGTTTGATAACTGTGCATTAGACTTC"
        readContent = contigContent[10:80]
        SeqIO.write([SeqRecord(Seq(contigContent, generic_dna), id="ContigDummy", description="")], self.folderName + self.contigsFilename , "fasta")
        SeqIO.write([SeqRecord(Seq(readContent, generic_dna), id="ReadDummy", description="")], self.folderName + self.readsFilename , "fasta")
    
    def test_extractRead2Contig(self):
        self.createSimpleFasta()   
        dataList = alignmentLib.extractRead2Contig(self.folderName, self.mummerLink, self.readsFilename, \
                    self.contigsFilename, self.splitNum, self.outputHeader, self.parallelNum, self.debug)
        assert(dataList == [[11, 80, 1, 70, 70, 70, 100.0, 100, 70, 'ContigDummy', 'ReadDummy']])
    
    def test_findConnectingReadsList(self):
        dataList = [[601, 800, 1, 200, 200, 200, 100.0, 800, 400, 'ContigDummy', 'ReadDummy']]
        connectingReadsList = readConnectivityLib.findConnectingReadsList(dataList)
        assert(connectingReadsList == [['ReadDummy', 'L', 'ContigDummy', [601, 800, 1, 200, 200, 200, 100.0, 800, 400, 'ContigDummy', 'ReadDummy']]])

    def test_findConnectingReadsListEmbed(self):
        dataList = [[1, 200, 201, 400, 200, 200, 100.0, 200, 800, 'ContigDummy', 'ReadDummy']]
        connectingReadsList = readConnectivityLib.findConnectingReadsList(dataList)
        assert(connectingReadsList == [['ReadDummy', 'B', 'ContigDummy',[1, 200, 201, 400, 200, 200, 100.0, 200, 800, 'ContigDummy', 'ReadDummy']]])

    def test_findSpanReadsList(self):
        connectingReadsList = []
        contigDummyLRecord, contigDummyRRecord = [601, 800, 1, 200, 200, 200, 100.0, 800, 400, 'ContigDummyL', 'ReadDummy'],  [1, 200, 201, 400, 200, 200, 100.0, 800, 400, 'ContigDummyR', 'ReadDummy']
        connectingReadsList.append(['ReadDummy', 'L', 'ContigDummyL', contigDummyLRecord])
        connectingReadsList.append(['ReadDummy', 'R', 'ContigDummyR', contigDummyRRecord])
        spanReadsList, contigGapReadLookUpDic = readConnectivityLib.findSpanReadsList(connectingReadsList)
        
        assert(spanReadsList == [  ['ContigDummyL_p', 'ContigDummyR_p', 'ReadDummy'] ] )
        assert(len(contigGapReadLookUpDic) == 1)
        assert(contigGapReadLookUpDic['ContigDummyL_p-ContigDummyR_p'].sort() == [[contigDummyLRecord,contigDummyRRecord]].sort())

    def test_findSpanReadsListEmbed(self):
        connectingReadsList = []
        contigDummyLRecord = [601, 800, 1, 200, 200, 200, 100.0, 800, 400, 'ContigDummyL', 'ReadDummy']  
        contigDummyBRecord1 = [1, 200, 101, 300, 200, 200, 100.0, 200, 400, 'ContigDummyB1', 'ReadDummy']
        contigDummyBRecord2 = [1, 200, 350, 151, 200, 200, 100.0, 200, 400, 'ContigDummyB2', 'ReadDummy']
        contigDummyRRecord = [1, 200, 201, 400, 200, 200, 100.0, 800, 400, 'ContigDummyR', 'ReadDummy']
        
        connectingReadsList.append(['ReadDummy', 'L', 'ContigDummyL', contigDummyLRecord])
        connectingReadsList.append(['ReadDummy', 'B', 'ContigDummyB1', contigDummyBRecord1])
        connectingReadsList.append(['ReadDummy', 'B', 'ContigDummyB2', contigDummyBRecord2])
        connectingReadsList.append(['ReadDummy', 'R', 'ContigDummyR', contigDummyRRecord])
        spanReadsList, contigGapReadLookUpDic = readConnectivityLib.findSpanReadsList(connectingReadsList)
        
        expectedSpanReadsList =  [['ContigDummyL_p', 'ContigDummyB1_p', 'ReadDummy'], \
                                  ['ContigDummyB1_p', 'ContigDummyB2_d', 'ReadDummy'], \
                                  ['ContigDummyB2_d', 'ContigDummyR_p', 'ReadDummy']]
        
        assert(spanReadsList.sort() == expectedSpanReadsList.sort())
        assert(len(contigGapReadLookUpDic) == 3)

        assert(contigGapReadLookUpDic['ContigDummyL_p-ContigDummyB1_p'].sort() == [[contigDummyLRecord,contigDummyBRecord1]].sort())
        assert(contigGapReadLookUpDic['ContigDummyB1_p-ContigDummyB2_d'].sort() == [[contigDummyBRecord1,contigDummyBRecord2]].sort())
        assert(contigGapReadLookUpDic['ContigDummyB2_d-ContigDummyR_p'].sort() == [[contigDummyBRecord2,contigDummyRRecord]].sort())
    
    def test_findContigsNames(self):
        self.createSimpleFasta()   
        contigsNamesList = alignmentLib.findContigsNames(self.folderName, self.contigsFilename)
        assert(contigsNamesList == ['ContigDummy'] )

    def test_formContigGraph(self):
        spanReadsList, contigsNamesList = [['ContigDummyL_p', 'ContigDummyR_p', 'ReadDummy']], ['ContigDummyL', 'ContigDummyR']
        G = graphLib.formContigGraph(spanReadsList, contigsNamesList)
        tmpObjL = G.dicOfContigNodes['ContigDummyL']
        tmpObjR = G.dicOfContigNodes['ContigDummyR']
        
        assert(len(tmpObjL.leftEndContainer.connectedContigsDic) == 0)
        assert(len(tmpObjL.rightEndContainer.connectedContigsDic) == 1 )
        assert(len(tmpObjR.leftEndContainer.connectedContigsDic) == 1)
        assert(len(tmpObjR.rightEndContainer.connectedContigsDic) == 0 )

    def test_findCondenseCandidatesList(self):
        spanReadsList, contigsNamesList = [['ContigDummyL_p', 'ContigDummyR_p', 'ReadDummy']], ['ContigDummyL', 'ContigDummyR']
        G = graphLib.formContigGraph(spanReadsList, contigsNamesList)
        condenseCandidatesList = G.findCondenseCandidatesList()  

        assert(condenseCandidatesList == [[['ContigDummyL_R~ContigDummyR_L~1'], False]])
    
    def test_transformConnectingReadsToSetStructure(self):
        connectingReadsList = []
        contigDummyLRecord = [601, 800, 1, 200, 200, 200, 100.0, 800, 400, 'ContigDummyL', 'ReadDummy']  
        contigDummyBRecord1 = [1, 200, 101, 300, 200, 200, 100.0, 200, 400, 'ContigDummyB1', 'ReadDummy']
        contigDummyBRecord2 = [1, 200, 350, 151, 200, 200, 100.0, 200, 400, 'ContigDummyB2', 'ReadDummy']
        contigDummyRRecord = [1, 200, 201, 400, 200, 200, 100.0, 800, 400, 'ContigDummyR', 'ReadDummy']
        
        connectingReadsList.append(['ReadDummy', 'L', 'ContigDummyL', contigDummyLRecord])
        connectingReadsList.append(['ReadDummy', 'B', 'ContigDummyB1', contigDummyBRecord1])
        connectingReadsList.append(['ReadDummy', 'B', 'ContigDummyB2', contigDummyBRecord2])
        connectingReadsList.append(['ReadDummy', 'R', 'ContigDummyR', contigDummyRRecord])

        setInfo =  setCoverLib.transformConnectingReadsToSetStructure(connectingReadsList)

        setOfElements = ['ContigDummyL', 'ContigDummyB1', 'ContigDummyB2', 'ContigDummyR']
        linkageAlongList = ['ContigDummyL_R~ContigDummyB1_L~3', 'ContigDummyB1_R~ContigDummyB2_R~3', 'ContigDummyB2_L~ContigDummyR_L~3']
        recordItem = [setOfElements, linkageAlongList, len(setOfElements)]
        expectedSetInfo = [recordItem]

        assert(expectedSetInfo == setInfo)

    def test_findUnUsedContigs(self):
        condenseCandidatesList, contigsNamesList = ['ContigDummyL_R~ContigDummyR_L~1'], ['ContigDummyL', 'ContigDummyR', 'ContigDummyB']
        unUsedContigsDic = setCoverLib.findUnUsedContigs(condenseCandidatesList, contigsNamesList)
        expectedUnUsedContigsDic = {'ContigDummyL': False, 'ContigDummyB': True, 'ContigDummyR': False}
        assert(unUsedContigsDic == expectedUnUsedContigsDic)

    def test_findSetCoverBaseLine(self):
        unUsedContigsDic, connectingReadsList = {'ContigDummyL': True, 'ContigDummyB1': True, 'ContigDummyB2': True, 'ContigDummyR': True} , []
        contigDummyLRecord = [601, 800, 1, 200, 200, 200, 100.0, 800, 400, 'ContigDummyL', 'ReadDummy']  
        contigDummyBRecord1 = [1, 200, 101, 300, 200, 200, 100.0, 200, 400, 'ContigDummyB1', 'ReadDummy']
        contigDummyBRecord2 = [1, 200, 350, 151, 200, 200, 100.0, 200, 400, 'ContigDummyB2', 'ReadDummy']
        contigDummyRRecord = [1, 200, 201, 400, 200, 200, 100.0, 800, 400, 'ContigDummyR', 'ReadDummy']

        connectingReadsList.append(['ReadDummy', 'L', 'ContigDummyL', contigDummyLRecord])
        connectingReadsList.append(['ReadDummy', 'B', 'ContigDummyB1', contigDummyBRecord1])
        connectingReadsList.append(['ReadDummy', 'B', 'ContigDummyB2', contigDummyBRecord2])
        connectingReadsList.append(['ReadDummy', 'R', 'ContigDummyR', contigDummyRRecord])
        
        setCoverList = setCoverLib.findSetCoverBaseLine(unUsedContigsDic, connectingReadsList)

        expectedSetCoverList = [[['ContigDummyL_R~ContigDummyB1_L~3', 'ContigDummyB1_R~ContigDummyB2_R~3', 'ContigDummyB2_L~ContigDummyR_L~3'], True]]

        assert(expectedSetCoverList == setCoverList)

    def test_findSetCoverGreedy(self):

        unUsedContigsDic, connectingReadsList = {'ContigDummyL': True, 'ContigDummyB1': True, 'ContigDummyB2': True, 'ContigDummyR': True} , []
        contigDummyLRecord = [601, 800, 1, 200, 200, 200, 100.0, 800, 400, 'ContigDummyL', 'ReadDummy']  
        contigDummyBRecord1 = [1, 200, 101, 300, 200, 200, 100.0, 200, 400, 'ContigDummyB1', 'ReadDummy']
        contigDummyRRecord = [1, 200, 201, 400, 200, 200, 100.0, 800, 400, 'ContigDummyR', 'ReadDummy']

        contigDummyBRecord2 = [1, 200, 350, 151, 200, 200, 100.0, 200, 400, 'ContigDummyB2', 'ReadDummy2']
        contigDummyBRecord3 = [1, 200, 101, 300, 200, 200, 100.0, 200, 400, 'ContigDummyB1', 'ReadDummy2']

        connectingReadsList.append(['ReadDummy', 'L', 'ContigDummyL', contigDummyLRecord])
        connectingReadsList.append(['ReadDummy', 'B', 'ContigDummyB1', contigDummyBRecord1])
        connectingReadsList.append(['ReadDummy', 'R', 'ContigDummyR', contigDummyRRecord])
        
        connectingReadsList.append(['ReadDummy2', 'B', 'ContigDummyB2', contigDummyBRecord2])
        connectingReadsList.append(['ReadDummy2', 'B', 'ContigDummyB1', contigDummyBRecord3])

        setCoverList = setCoverLib.findSetCoverGreedy(unUsedContigsDic, connectingReadsList)
        
        expectedSetCoverList =  ['ContigDummyL_R~ContigDummyB1_L~3', 'ContigDummyB1_R~ContigDummyR_L~3', 'ContigDummyB1_R~ContigDummyB2_R~3']       
        
        assert(expectedSetCoverList.sort() == setCoverList.sort())

    def test_assignCoverageFromDataList(self):
        dataList, contigList = [ [1, 6, 1, 6, 6, 6, 100.0, 6, 6, 'ContigDummyL', 'ReadDummy'] ], []

        contigList.append(SeqRecord(Seq("AAACCC", generic_dna), id="ContigDummyL", description=""))
        contigList.append(SeqRecord(Seq("CCCTTTT", generic_dna), id="ContigDummyR", description=""))
        SeqIO.write(contigList, self.folderName + self.contigsFilename , "fasta")
        
        spanReadsList, contigsNamesList = [['ContigDummyL_p', 'ContigDummyR_p', 'ReadDummy']], ['ContigDummyL', 'ContigDummyR']
        G = graphLib.formContigGraph(spanReadsList, contigsNamesList)
        cTestLib.assignCoverageFromDataList(G, dataList, self.folderName, self.contigsFilename)

        assert(G.dicOfContigNodes['ContigDummyL'].contigLength == 6)
        assert(G.dicOfContigNodes['ContigDummyR'].contigLength == 7)
        assert(G.dicOfContigNodes['ContigDummyL'].readToContigCount == 1)
        assert(G.dicOfContigNodes['ContigDummyR'].readToContigCount == 0)
       
    def test_calculateConfidenceScore(self):
        condenseCandidatesList = [[['ContigDummyL_R~ContigDummyR_L~1'], False]]
        spanReadsList, contigsNamesList = [['ContigDummyL_p', 'ContigDummyR_p', 'ReadDummy']], ['ContigDummyL', 'ContigDummyR']
        
        dataList, contigList = [ [1, 6, 1, 6, 6, 6, 100.0, 6, 6, 'ContigDummyL', 'ReadDummy'] ], []
        contigList.append(SeqRecord(Seq("AAACCC", generic_dna), id="ContigDummyL", description=""))
        contigList.append(SeqRecord(Seq("CCCTTTT", generic_dna), id="ContigDummyR", description=""))
        SeqIO.write(contigList, self.folderName + self.contigsFilename , "fasta")

        G = graphLib.formContigGraph(spanReadsList, contigsNamesList)        
        cTestLib.assignCoverageFromDataList(G, dataList, self.folderName, self.contigsFilename)
        scoreStructList = cTestLib.calculateConfidenceScore(G, condenseCandidatesList)

        assert(scoreStructList[0][0][0][0] == 'ContigDummyL_R~ContigDummyR_L~1' )
        assert(abs(scoreStructList[0][0][0][1] -  0.53846153846153844) < 0.01)
        assert(scoreStructList[0][0][0][2] == 1)
        assert(scoreStructList[0][1] == False)

    def test_assignRepeatedNodesToDummy(self):
        scoreList = [[[['ContigDummyL_R~ContigDummyR_L~1' , 1 , 1]], False ] ] 
        scoreListWithDummy, dummyNodeDataRobot = setCoverLib.assignRepeatedNodesToDummy(scoreList)
        expectedScoreListWithDummy = [['Dummy0_R~Dummy1_L~1', 1, 1]]

        assert( "ContigDummyL" == dummyNodeDataRobot.D2RLookUp("Dummy0"))
        assert( "ContigDummyR" == dummyNodeDataRobot.D2RLookUp("Dummy1"))
        
        assert( ["Dummy0"] == dummyNodeDataRobot.R2DLookUp("ContigDummyL"))
        assert( ["Dummy1"] == dummyNodeDataRobot.R2DLookUp("ContigDummyR"))
        
        assert(expectedScoreListWithDummy == scoreListWithDummy)

    def test_rankAndMerge(self):
        contigList = []
        contigList.append(SeqRecord(Seq("AAACCC", generic_dna), id="ContigDummyL", description=""))
        contigList.append(SeqRecord(Seq("CCCTTTT", generic_dna), id="ContigDummyR", description=""))
        SeqIO.write(contigList, self.folderName + self.contigsFilename , "fasta")
        
        SeqIO.write([SeqRecord(Seq("CCCGGGCCC", generic_dna), id="ReadDummy", description="")], self.folderName + self.readsFilename , "fasta")
        
        scoreList = [ ['ContigDummyL_R~ContigDummyR_L~1' , 1 , 1] ] 

        contigGapReadLookUpDic = {}
        contigDummyLRecord, contigDummyRRecord = [4, 6, 1, 3, 3, 3, 100.0, 6, 9, 'ContigDummyL', 'ReadDummy'],  [1, 3, 7, 9, 3, 3, 100.0, 7, 9, 'ContigDummyR', 'ReadDummy']
        contigGapReadLookUpDic['ContigDummyL_p-ContigDummyR_p'] = [[contigDummyLRecord, contigDummyRRecord]]

        contigsNamesList = alignmentLib.findContigsNames(self.folderName, self.contigsFilename)
        
        dummyNodeDataRobot = setCoverLib.dummyNodeController()

        dummyNodeDataRobot.realToDummyDic = {'ContigDummyL': 'ContigDummyL','ContigDummyR': 'ContigDummyR'}
        dummyNodeDataRobot.dummyToRealDic = {'ContigDummyL': 'ContigDummyL','ContigDummyR': 'ContigDummyR'}

        rankingLib.rankAndMerge(self.folderName,contigsNamesList, self.contigsFilename, self.readsFilename, scoreList, contigGapReadLookUpDic, 1, 0.95, "scoreList.json", "improved.fasta", dummyNodeDataRobot)
         
        expectedContig= "AAACCC" + "GGG" + "CCCTTTT"
        records = list(SeqIO.parse(self.folderName + "improved.fasta", "fasta"))
        assert(expectedContig == str(records[0].seq))

    def tearDown(self):
        print "Teardown : Started"
        os.system("rm -rf " + self.folderName)
        print "Teardown : Done"
        
def main():
    unittest.main()
    
if __name__ == '__main__':
    main()





