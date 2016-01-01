import alignmentLib
import cTestLib
import graphLib
import rankingLib
import readConnectivityLib

def reportGraph(GInput):
    for eachitem in GInput.dicOfContigNodes:
        tmpObj = GInput.dicOfContigNodes[eachitem]
        print tmpObj.name, len(tmpObj.leftEndContainer.connectedContigsDic), len(tmpObj.rightEndContainer.connectedContigsDic)
        for item in tmpObj.leftEndContainer.connectedContigsDic:
            print "\t Left",  len(tmpObj.leftEndContainer.connectedContigsDic[item])
        for item in tmpObj.rightEndContainer.connectedContigsDic:
            print "\t Right",  len(tmpObj.rightEndContainer.connectedContigsDic[item])

def tempDebug():
    folderName = "/Users/kakitlam/Desktop/baseline/datafolder/"
    mummerLink = "/Users/kakitlam/Desktop/allinone/experimentBench/MUMmer3.23/"
    readsFilename = "LR.fasta"
    contigsFilename = "mFixed2.fasta"
    outputHeader = "outputHeader"
    splitNum = 20
    parallelNum = 4
    debug = True
    
    dataList = alignmentLib.extractRead2Contig(folderName, mummerLink, readsFilename, \
                contigsFilename, splitNum, outputHeader, parallelNum, debug )

    print dataList[0]
    
    connectingReadsList = readConnectivityLib.findConnectingReadsList(dataList)
    print connectingReadsList[0], len(connectingReadsList)

    spanReadsList, contigGapReadLookUpDic = readConnectivityLib.findSpanReadsList(connectingReadsList)
    print len(spanReadsList),  spanReadsList[0], len(contigGapReadLookUpDic)

    contigsNamesList = alignmentLib.findContigsNames(folderName, contigsFilename)

    G = graphLib.formContigGraph(spanReadsList, contigsNamesList)
    reportGraph(G)

    condenseCandidatesList = G.findCondenseCandidatesList()

    cTestLib.assignCoverageFromDataList(G, dataList,folderName, contigsFilename)

    scoreList = cTestLib.calculateConfidenceScore(G, condenseCandidatesList)
    print scoreList 

    rankingLib.rankAndMerge(folderName,contigsNamesList, contigsFilename, readsFilename, scoreList, contigGapReadLookUpDic)


tempDebug()














