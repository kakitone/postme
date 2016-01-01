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
import rankingLib
import readConnectivityLib

def mainFlow(folderName, mummerLink, contigsFilename, readsFilename):
    outputHeader, splitNum, parallelNum, debug= "readToContigHeader",  20, 4, False    
    dataList = alignmentLib.extractRead2Contig(folderName, mummerLink, readsFilename, contigsFilename, splitNum, outputHeader, parallelNum, debug )
    connectingReadsList = readConnectivityLib.findConnectingReadsList(dataList)
    spanReadsList, contigGapReadLookUpDic = readConnectivityLib.findSpanReadsList(connectingReadsList)
    contigsNamesList = alignmentLib.findContigsNames(folderName, contigsFilename)
    G = graphLib.formContigGraph(spanReadsList, contigsNamesList)
    condenseCandidatesList = G.findCondenseCandidatesList()
    cTestLib.assignCoverageFromDataList(G, dataList,folderName, contigsFilename)
    scoreList = cTestLib.calculateConfidenceScore(G, condenseCandidatesList)
    rankingLib.rankAndMerge(folderName,contigsNamesList, contigsFilename, readsFilename, scoreList, contigGapReadLookUpDic)


parser = argparse.ArgumentParser(description='PostMe')
parser.add_argument('-o', '--outputFolder', help= 'Output folder path', required=False)
parser.add_argument('-a', '--alignmentPath', help= 'MUMmer path', required=False)
parser.add_argument('-r', '--readsFilename', help= 'Filename for reads file', required=False)
parser.add_argument('-c', '--contigsFilename', help= 'Filename for contigs file', required=False)

args = vars(parser.parse_args())
folderName, mummerLink = args['outputFolder'], args['alignmentPath']
contigsFilename, readsFilename =  args['contigsFilename'], args['readsFilename']

t0 = time.time()
mainFlow(folderName, mummerLink, contigsFilename, readsFilename)
print  "Time", time.time() - t0
