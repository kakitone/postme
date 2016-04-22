import math
import numpy as np  
from operator import itemgetter
import matplotlib.pyplot as plt

def findY(chat, l, L):
	return chat * l * 1.0 /L

def computeLogLikelihood(c1, l1, c2, l2, L, c1hat, c2hat):
	y1hat, y2hat = findY(c1hat, l1, L), findY(c2hat, l2, L)
	y1, y2 = c1*l1*1.0/L, c2*l2*1.0/L
	return -y1 - y2 + y1hat*math.log(y1) + y2hat*math.log(y2)

def findMLEFit(c1hat, c2hat, l1, l2, L, eps, rangeList):
	ck = False
	for eachitem in rangeList:
		if eachitem[0] < c1hat*1.0/c2hat < eachitem[1]:
			ck = True

	if ck :
		c1, c2 = c1hat , c2hat 
		return computeLogLikelihood(c1, l1, c2, l2, L, c1hat, c2hat)

	else:
		c2 = (c1hat*l1 + c2hat*l2)*1.0/((1-eps)*l1 + l2)
		c1 = c2*(1- eps)
		ml1 = computeLogLikelihood(c1, l1, c2, l2, L, c1hat, c2hat)

		c2 = (c1hat*l1 + c2hat*l2)*1.0/((1+eps)*l1 + l2)
		c1 = c2*(1+ eps)
		ml2 = computeLogLikelihood(c1, l1, c2, l2, L, c1hat, c2hat)

		return max(ml1, ml2)

def computeLogRatio(l1, l2, c1hat, c2hat, L):
	eps = 0.05

	omega = [ [1- eps , 1+ eps] ] 
	notomega = [ [0, 1-eps], [1+ eps , 10**9]]

	return findMLEFit(c1hat, c2hat, l1, l2, L, eps, omega) - findMLEFit(c1hat, c2hat, l1, l2, L, eps, notomega)

def unitTest():
	L = 100
	for j in range(5):
		for i in range(5):
			cov1 = 50.0  
			cov2 = 50.0 + 2.0*i
			l1 = 30*(10**j)
			l2 = 30*(10**j)

			logPdiff = computeLogRatio(l1, l2, cov1, cov2, L)
			print "cov1, cov2, l1, l2, logPdiff: ", cov1, cov2, l1, l2, logPdiff 

		print ""


def simulation(l1List, l2List , L):
	numberOfClusters = 5
	cList = [5*i + 5 for i in range(numberOfClusters)]
	resultList = []
	numberOfRounds = 100000
	
	TP = 0
	for i in range(numberOfRounds):	
		index1, index2 = np.random.randint(numberOfClusters), np.random.randint(numberOfClusters) 
		ellIndex1, ellIndex2 = np.random.randint(len(l1List)), np.random.randint(len(l2List))

		l1, l2 = l1List[ellIndex1], l2List[ellIndex2]

		y1, y2 = np.random.poisson(cList[index1] * l1/L), np.random.poisson(cList[index2]* l2/L)

		c1, c2 = y1*1.0*L/l1, y2*1.0*L/l2

		logPdiff = computeLogRatio(l1, l2, c1, c2, L)
		diffDivMax = abs(c1 - c2)/ max(c1, c2) 
		diffDivMin = abs(c1 - c2)/ min(c1, c2)

		if index1 == index2 : 
			TP += 1

		#print "logPdiff, index1 == index2, index1, index2 , y1, y2, c1, c2 : ", logPdiff, index1 == index2, index1, index2 , y1, y2, c1, c2
		resultList.append([l1, l2, c1, c2, index1, index2, logPdiff, diffDivMax, diffDivMin, 1 if index1 == index2 else 0])

	#resultList.sort(key = itemgetter(6), reverse= True)
	ROCLogDiffx1, ROCLogDiffy1 = findResultByMetric(resultList, "log", TP)
	ROCLogDiffx2, ROCLogDiffy2 = findResultByMetric(resultList, "max", TP)
	ROCLogDiffx3, ROCLogDiffy3 = findResultByMetric(resultList, "min", TP)

	if True:
		resultList.sort(key=itemgetter(6), reverse = True)
		for i in range(10):
			print i*numberOfRounds, resultList[i*numberOfRounds/10],  ROCLogDiffx1[i*numberOfRounds/10], ROCLogDiffy1[i*numberOfRounds/10]

		plt.scatter(ROCLogDiffx1, ROCLogDiffy1, color = 'r')
		plt.scatter(ROCLogDiffx2, ROCLogDiffy2, color = 'g')
		plt.scatter(ROCLogDiffx3, ROCLogDiffy3, color = 'b', marker='x')

		plt.xlabel("Recall")
		plt.ylabel("Precision")
		plt.show()


def findResultByMetric(resultList, mytype, TP):

	if mytype == "log":
		resultList.sort(key = itemgetter(6), reverse = True)
	elif mytype == "max":
		resultList.sort(key = itemgetter(7))	
	elif mytype == "min":
		resultList.sort(key = itemgetter(8))		

	ROCLogDiffx, ROCLogDiffy = [], []
	
	PTP = 0

	for i in range(len(resultList)):
		eachitem = resultList[i]
		if eachitem[-1] == 1:
			PTP += 1

		precision, recall = PTP*1.0/(i+1)  , PTP*1.0/(TP) 

		ROCLogDiffx.append(recall)
		ROCLogDiffy.append(precision)

	return ROCLogDiffx, ROCLogDiffy

unitTest()	 
simulation([500 + i*500 for i in range(10)], [500+ i*500 for i in range(10)] , 100)

