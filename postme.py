'''
Typical usage : 
python postme.py \
    -o datafolder/ \
    -a ~/Desktop/allinone/experimentBench/MUMmer3.23/ \
    -c mFixed2.fasta \
    -r LR.fasta
'''
### Standard librariess
import argparse
import json
import time
import os

### Local libraries
import alignmentLib
import cTestLib
import graphLib
import houseKeeperLib
import rankingLib
import readConnectivityLib
import setCoverLib

def mainFlow(folderName, mummerLink, inputContigsFilename, inputReadsFilename, useSpades, noAlignment, scoreListOutputName, outputContigsFilename, mScoreThres, conScoreThres):
    outputHeader, splitNum, parallelNum = "readToContigHeader",  20, 20  
    contigsFilename, readsFilename= "tmp" + inputContigsFilename , "tmp" + inputReadsFilename

    targetToSourceContigsNamesDic = houseKeeperLib.transformFileHeaders(folderName, inputContigsFilename, contigsFilename, noAlignment)
    targetToSourceReadsNamesDic = houseKeeperLib.transformFileHeaders(folderName, inputReadsFilename, readsFilename, noAlignment)

    dataList = alignmentLib.extractRead2Contig(folderName, mummerLink, readsFilename, contigsFilename, splitNum, outputHeader, parallelNum, noAlignment )
    
    connectingReadsList = readConnectivityLib.findConnectingReadsList(dataList)
    
    spanReadsList, contigGapReadLookUpDic = readConnectivityLib.findSpanReadsList(connectingReadsList)
    
    contigsNamesList = alignmentLib.findContigsNames(folderName, contigsFilename)
    
    G = graphLib.formContigGraph(spanReadsList, contigsNamesList)
    
    condenseCandidatesList = G.findCondenseCandidatesList()

    potentialMergesList = setCoverLib.extendConnectivityFromReads(condenseCandidatesList, connectingReadsList, contigsNamesList)
    
    if useSpades == True:
        cTestLib.assignCoverageFromHeader(G, folderName, contigsFilename, targetToSourceContigsNamesDic)
    else:
        cTestLib.assignCoverageFromDataList(G, dataList,folderName, contigsFilename)
    
    scoreList = cTestLib.calculateConfidenceScore(G, potentialMergesList)
    
    scoreListWithDummy, dummyNodeDataRobot = setCoverLib.assignRepeatedNodesToDummy(scoreList)

    rankingLib.rankAndMerge(folderName,contigsNamesList, contigsFilename, readsFilename, scoreListWithDummy, contigGapReadLookUpDic, mScoreThres, conScoreThres, scoreListOutputName, outputContigsFilename, dummyNodeDataRobot)

parser = argparse.ArgumentParser(description='PostMe')
parser.add_argument('-o', '--outputFolder', help= 'Output folder path', required=False)
parser.add_argument('-a', '--alignmentPath', help= 'MUMmer path', required=False)
parser.add_argument('-r', '--readsFilename', help= 'Filename for reads file', required=False)
parser.add_argument('-c', '--contigsFilename', help= 'Filename for contigs file', required=False)
parser.add_argument('-s', '--inputCoverageViaSpadesHeader', help= 'Input coverage from spades header files (T/F)', required=False)
parser.add_argument('-na', '--noAlignment', help= 'Use existing alignment file (T/F)', required=False)
parser.add_argument('-sn', '--scoreListOutputName', help= 'Output filename for scoreList', required=False)
parser.add_argument('-on', '--outputContigsFilename', help= 'Output filename for the improved contigs in FASTA format', required=False)
parser.add_argument('-ms', '--mScoreCutOff', help= 'Number of spanning reads cutoff', required=False)
parser.add_argument('-cs', '--cScoreCutOff', help= 'Abundance confidence score cutoff', required=False)

args = vars(parser.parse_args())
folderName, mummerLink = houseKeeperLib.trailingFolderCorrection(args['outputFolder']) , houseKeeperLib.trailingFolderCorrection(args['alignmentPath'])
contigsFilename, readsFilename =  args['contigsFilename'], args['readsFilename']

useSpades = True if args['inputCoverageViaSpadesHeader'] == 'T' else False
noAlignment = True if args['noAlignment'] == 'T' else False
scoreListOutputName = args['scoreListOutputName'] if args['scoreListOutputName'] != None else "scoreList.json"
outputContigsFilename = args['outputContigsFilename'] if args['outputContigsFilename'] != None else "improved.fasta"
mScoreThres = int(args['mScoreCutOff']) if args['mScoreCutOff'] != None  else  2 
conScoreThres = float(args['cScoreCutOff']) if args['cScoreCutOff'] != None  else  0.95

t0 = time.time()
mainFlow(folderName, mummerLink, contigsFilename, readsFilename, useSpades, noAlignment, scoreListOutputName, outputContigsFilename, mScoreThres, conScoreThres)
print  "Time", time.time() - t0
