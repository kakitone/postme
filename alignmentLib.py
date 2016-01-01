'''
Input : raw_reads.fasta, contigs.fasta
Output : dataList of mummerData
Algorithm:
    1. Split raw_reads into multiple files
    2. Align split_raw_reads against the contigs
    3. Parse the combined data
    4. Output the dataList
'''
from Bio import SeqIO
from multiprocessing import Pool
import os
import re

def splitFasta(folderName, filename, splitNum, isDebug):
    header, listOfFilenames = "split_", []

    if not isDebug:
        records = list(SeqIO.parse(folderName + filename, "fasta"))
        chunkSize = len(records)/splitNum

        indicesList = [[i*chunkSize, (i+1)*chunkSize] for i in range(splitNum)]
        indicesList[-1][1] = len(records) 
                
        for i in range(splitNum):
            begin, end = indicesList[i][0], indicesList[i][1]
            splitFilename = header + "%02d" % i + "_" + filename
            SeqIO.write(records[begin:end], folderName + splitFilename , "fasta")
            listOfFilenames.append(splitFilename)
    else:
        listOfFilenames = [ header + "%02d" % i + "_" + filename for i in range(splitNum) ]

    return listOfFilenames

def runSingleAlign(folderName, mummerLink,  referenceName, queryName, outputHeader):
    command = mummerLink + "nucmer --maxmatch  -p " + folderName + outputHeader \
             + " " + folderName + referenceName + " " + folderName + queryName
    os.system(command)
    
    command = mummerLink + "show-coords -l -r " + folderName + outputHeader + \
                ".delta > " + folderName + outputHeader + "Out"
    os.system(command)

def executeJobs(args):
    func, parameters = args
    func(*parameters)

def parallelAlign(folderName,mummerLink, contigsFilename, listOfReadsFilenames, outputHeader, nProc, isDebug):
    mummerOutFilenameList, workerList = [] , []

    for i in range(len(listOfReadsFilenames)):
        workerList.append([folderName, mummerLink,  contigsFilename, \
            listOfReadsFilenames[i], outputHeader + "%02d" % i ])
        mummerOutFilenameList.append(outputHeader + "%02d" % i + "Out")

    if not isDebug:
        p = Pool(processes=nProc)
        jobList = []
        for eachjob in workerList:
            folderName, mummerLink,  referenceName, queryName, outputHeader = eachjob
            jobList.append((runSingleAlign, (folderName, mummerLink,  referenceName, queryName, outputHeader)))

        p.map_async(executeJobs, jobList, chunksize=max(1,len(jobList)/nProc))
        p.close()
        p.join()

    return mummerOutFilenameList

def convertMUMmerRecord(record):
    rawRecordList = re.split('[\t | \n ]',record)
    recordList = filter(None, rawRecordList)
    formattedrecord = [ int(recordList[i]) for i in range(6)]
    formattedrecord.append(float(recordList[6]))
    for i in range(7,9):
        formattedrecord.append(int(recordList[i]))
    formattedrecord += recordList[9:]
    return formattedrecord

def combineAlignmentData(folderName, mummerOutFilenameList):
    dataList = []

    for eachfile in mummerOutFilenameList:
        f = open(folderName + eachfile, 'r')
        
        for i in range(6):
            tmp = f.readline()

        while len(tmp) > 0:
            dataList.append(convertMUMmerRecord(tmp))
            tmp = f.readline().rstrip()
        
        f.close()

    return dataList

def extractRead2Contig(folderName, mummerLink, readsFilename, contigsFilename, splitNum, outputHeader, parallelNum, isDebug):
    splitReadsFilenameList = splitFasta(folderName, readsFilename, splitNum, isDebug)        
    mummerOutFilenameList = parallelAlign(folderName, mummerLink, contigsFilename, splitReadsFilenameList, outputHeader, parallelNum, isDebug)
    dataList = combineAlignmentData(folderName, mummerOutFilenameList)
    return dataList

def findContigsNames(folderName, filename):
    return [record.id for record in SeqIO.parse(folderName + filename, "fasta")]
    