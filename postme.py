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

def mainFlow(folderName, mummerLink, inputContigsFilename, inputReadsFilename, useSpades):
    outputHeader, splitNum, parallelNum, debug = "readToContigHeader",  20, 4, False    
    contigsFilename, readsFilename= "tmp" + inputContigsFilename , "tmp" + inputReadsFilename

    targetToSourceContigsNamesDic = houseKeeperLib.transformFileHeaders(folderName, inputContigsFilename, contigsFilename)
    targetToSourceReadsNamesDic = houseKeeperLib.transformFileHeaders(folderName, inputReadsFilename, readsFilename)

    dataList = alignmentLib.extractRead2Contig(folderName, mummerLink, readsFilename, contigsFilename, splitNum, outputHeader, parallelNum, debug )
    
    connectingReadsList = readConnectivityLib.findConnectingReadsList(dataList)
    
    spanReadsList, contigGapReadLookUpDic = readConnectivityLib.findSpanReadsList(connectingReadsList)
    
    contigsNamesList = alignmentLib.findContigsNames(folderName, contigsFilename)
    
    G = graphLib.formContigGraph(spanReadsList, contigsNamesList)
    
    condenseCandidatesList = G.findCondenseCandidatesList()

    if useSpades == "T":
        cTestLib.assignCoverageFromHeader(G, folderName, contigsFilename, targetToSourceContigsNamesDic)
    else:
        cTestLib.assignCoverageFromDataList(G, dataList,folderName, contigsFilename)
    
    scoreList = cTestLib.calculateConfidenceScore(G, condenseCandidatesList)
    
    rankingLib.rankAndMerge(folderName,contigsNamesList, contigsFilename, readsFilename, scoreList, contigGapReadLookUpDic)

parser = argparse.ArgumentParser(description='PostMe')
parser.add_argument('-o', '--outputFolder', help= 'Output folder path', required=False)
parser.add_argument('-a', '--alignmentPath', help= 'MUMmer path', required=False)
parser.add_argument('-r', '--readsFilename', help= 'Filename for reads file', required=False)
parser.add_argument('-c', '--contigsFilename', help= 'Filename for contigs file', required=False)
parser.add_argument('-s', '--inputCoverageViaSpadesHeader', help= 'Input coverage from spades header files (T/F)', required=False)

args = vars(parser.parse_args())
folderName, mummerLink = houseKeeperLib.trailingFolderCorrection(args['outputFolder']) , houseKeeperLib.trailingFolderCorrection(args['alignmentPath'])
contigsFilename, readsFilename =  args['contigsFilename'], args['readsFilename']
useSpades = args['inputCoverageViaSpadesHeader']

t0 = time.time()
mainFlow(folderName, mummerLink, contigsFilename, readsFilename, useSpades)
print  "Time", time.time() - t0
