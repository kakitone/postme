
'''
Typical usage : 

python dataAnalysisLib.py -o /home/kakitone/tmpkk/datafolder/  -a /usr/bin/  -r raw_reads.fasta -c contigs.fasta -q /home/kakitone/quast-2.3/ -f reference.fasta -na T

Input : contigs.fasta, reference.fasta
Output : settingsList.json, improved[0..N].fasta, performance.json, performance.png
Algorithm : 
    1. Run PostMe in batch
    2. Run QUAST analysis in a batch
    3. Extract QUAST analysis and generate plot

'''

import argparse
import json
import matplotlib.pyplot as plt
import os
import re
import time

import houseKeeperLib

def batchPostMe(folderName, mummerLink, contigsFilename, readsFilename, noAlignment):
    commandsList, settingsList = [],  []

    if not noAlignment:
        commandsList.append("python postme.py " \
            + " -o "+ folderName                \
            + " -a " + mummerLink               \
            + " -r " +  readsFilename           \
            + " -c " + contigsFilename          \
            + " -s T -na F -sn scoreList.json -on improved.fasta -ms 2 -cs 0.95")
        

    for i in range(1, 3):
        for j in range(5):
            mScoreThres, cScoreThres,  =  i, 0.95 - j*0.05
            settingsList.append(["score_" +  str(len(settingsList)) + ".fasta", \
                                 "imp_" + str(len(settingsList)) + ".fasta",    \
                                 mScoreThres, cScoreThres])

    commonHeader = "python postme.py " \
        + " -o "+ folderName           \
        + " -a " + mummerLink          \
        + " -r " +  readsFilename      \
        + " -c " + contigsFilename     \
        + " -s T -na T "

    for eachSetting in settingsList:
        command = commonHeader                       \
                    + " -sn " + str(eachSetting[0])  \
                    + " -on " + str(eachSetting[1])  \
                    + " -ms " + str(eachSetting[2])  \
                    + " -cs " + str(eachSetting[3]) 
        
        commandsList.append(command) 

    for eachCommand in commandsList: 
        os.system(eachCommand)

    houseKeeperLib.dumpDataToJson(folderName, "settingsList.json" ,settingsList)

def batchQuastAnalysis(quastPath, folderName, contigsFilename, referenceFilename):
    # Run Quast in a batch
    settingsList = houseKeeperLib.readInJSON(folderName, "settingsList.json")
    fileNames = ""
    for eachSetting in settingsList:
        fileNames = fileNames + " " + folderName + eachSetting[1] + " " 


    os.system("python " + quastPath + "quast.py " \
               + " -R " + folderName + referenceFilename \
               + " " + folderName + contigsFilename \
               + " " + fileNames)    
    
    # Parse Quast and review results 
    quastReport =  "quast_results/latest/report.txt"
    misassemblyList = []
    ncontigList = [] 
    
    with open(quastReport) as f2:
        for line in f2:
            result = re.match("# misassemblies (.*?) .*", line)
            if result != None: 
                misassemblyList = result.group(0).split()[2:]
                
            result = re.match("# contigs \(>= 0 bp\) (.*?) .*", line)
            if result != None: 
                ncontigList = result.group(0).split()[5:]

    ncontigList = [int(eachStr) for eachStr in ncontigList]
    misassemblyList = [int(eachStr) for eachStr in misassemblyList]


    # Plot graphs
    plt.scatter(misassemblyList, ncontigList, color = 'r')
    plt.scatter([152], [4439], color = 'b')

    plt.xlabel("Misassemblies Number")
    plt.ylabel("NContig")

    print ncontigList, misassemblyList
    plt.savefig("results")

parser = argparse.ArgumentParser(description='PostMe')
parser.add_argument('-o', '--outputFolder', help= 'Output folder path', required=False)
parser.add_argument('-a', '--alignmentPath', help= 'MUMmer path', required=False)
parser.add_argument('-r', '--readsFilename', help= 'Filename for reads file', required=False)
parser.add_argument('-c', '--contigsFilename', help= 'Filename for contigs file', required=False)
parser.add_argument('-na', '--noAlignment', help= 'Use existing alignment file (T/F)', required=False)
parser.add_argument('-q', '--quastPath', help= 'Path to QUAST', required=False)
parser.add_argument('-f', '--referenceFilename', help= 'Filename for reference file', required=False)

args = vars(parser.parse_args())
folderName, mummerLink = houseKeeperLib.trailingFolderCorrection(args['outputFolder']) , houseKeeperLib.trailingFolderCorrection(args['alignmentPath'])
contigsFilename, readsFilename =  args['contigsFilename'], args['readsFilename']
quastPath = houseKeeperLib.trailingFolderCorrection(args['quastPath'])
referenceFilename = args['referenceFilename']
noAlignment = True if args['noAlignment'] == 'T' else False

t0 = time.time()
batchPostMe(folderName, mummerLink, contigsFilename, readsFilename, noAlignment)
batchQuastAnalysis(quastPath, folderName, contigsFilename, referenceFilename)
print  "Time", time.time() - t0

