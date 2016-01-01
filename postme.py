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

def mainFlow(folderName, mummerLink):
	print "Welcome to PostMe !!!"

parser = argparse.ArgumentParser(description='PostMe')
parser.add_argument('-o', '--outputFolder', help= 'Output folder path', required=False)
parser.add_argument('-a', '--alignmentPath', help= 'MUMmer path', required=False)

args = vars(parser.parse_args())
folderName, mummerLink = args['outputFolder'], args['alignmentPath']

t0 = time.time()
mainFlow(folderName, mummerLink)
print  "Time", time.time() - t0
