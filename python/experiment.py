#!/usr/bin/python

import sys
import math
from bandit import Bandit, RandomBandit
import algorithms
import logging
import getopt
import random
import ast
import os

class prettyfloat(float):
    def __repr__(self):
        return "%0.2f" % self

def usage():
    """prints the usage message"""
    print "Usage:\n\tpython <exp.py> --min <min-pulls> --max <max-pulls> --step <pull-step> --repeat <repetitions> --plot <output-file.png> <input-file.txt> <algorithm-name> [algorithm-name]...\nor\n\tbanditExperiment([--min, <min-pulls>, --max, <max-pulls>, --step, <pull-step>, --repeat, <repetitions>, --plot, <output-file.png>, <input-file.txt>, <algorithm-name>, [algorithm-name],...])"

 
def findBestArm(meanList):
    return algorithms.maxIndex(meanList)


def experimentFunc(pullsNum, algoname, bandit):
    """excecute the experiment with pullsNum number of arm pulls according to the given algorithm.
    returns the index of the best arm"""
    
    algorithm = getattr(algorithms, algoname) 
     
    counts = [0 for _i in range(bandit.getArmsNum())]
    means = [0.0 for _i in range(bandit.getArmsNum())]
    cumulative_regret = 0.0
     
    for _pull in range (0, int(pullsNum)):        
        stats = zip(counts, means)
        arm = algorithm(stats)
        cumulative_regret += bandit.calcRegret(arm)
        
        pullResult = bandit.pullArm(arm)
        
        means[arm] = (counts[arm] * means[arm] + pullResult) / (counts[arm] + 1)
        counts[arm] = counts[arm] + 1 
    
    simple_regret = bandit.calcRegret(findBestArm(means))

    return simple_regret, cumulative_regret
    

def calcAverageRegret(pullsNum, repetitions, algorithms, bandit, randomBandit):
    """ execute experimentFunc for a given repetitions num and return the average regret"""
    # initialize sum of the regrets in all repetitions
    algoNum = len(algorithms)
    
    avgSimpleRegrets = [0.0 for _i in range(algoNum)]    
    avgCumulativeRegrets = [0.0 for _i in range(algoNum)]    
    
    # execute the experiment Function repetitions times
    for algo in range(algoNum):        
        simple_regret = 0.0
        cumulative_regret = 0.0
        
        for _repeat in range(repetitions):                    
            if randomBandit:
                bandit = RandomBandit(bandit.getArmsNum())
            
            sr, cr = experimentFunc(pullsNum, algorithms[algo], bandit)
            simple_regret += sr
            cumulative_regret += cr
            
        avgSimpleRegrets[algo] =  simple_regret / repetitions
        avgCumulativeRegrets[algo] = cumulative_regret / repetitions
                
    return avgSimpleRegrets, avgCumulativeRegrets
     
    
def readAvgFromFile(fileName):
    """returns a list of averages that was read from a given file"""
    banditFile = open(fileName, 'r')
    # string = banditFile.read()
    avgList = []
    
    logging.info('Reading from file.')
    
    for line in banditFile:
        avgList.extend([float(num) for num in line.split()])
    banditFile.close()
    
    logging.info('Done reading from file.')
    
    return avgList 
      
    
def experimentMainLoop(minPulls,maxPulls, pullStep, repetitions, bandit, algorithms, randomBandit):
    """the outer loop of the experiment"""
    ####print first row in results table:
    printRow = "%-10s " % "#Pulls"
    for algoname in algorithms:
        printRow += "%-10s %-10s" % (algoname + "[SR]", algoname + "[CR]")
    print "\n" + printRow
    
    pullsNum = minPulls
    
    logging.info('Start experiment.') 
    
    while pullsNum <= maxPulls:
        
        printRow = "%-10d " % pullsNum
        
        logging.info('Number of pulls: %f', pullsNum)
        avgSimpleRegrets, avgCumulativeRegrets = calcAverageRegret(pullsNum, repetitions, algorithms, bandit, randomBandit)
        
        for algo in range(len(algorithms)):
            printRow += "%-10f %-10f" % (avgSimpleRegrets[algo], avgCumulativeRegrets[algo])
            
        print printRow
        
        pullsNum = math.ceil(pullsNum * pullStep) 
        
    logging.info('Done experimant')
    
        
def banditExperiment(argsList):
    """receives list of arguments and perform the experiment"""
     
    if len(argsList) < 2:
        usage()
        sys.exit(2)
        
    logfile = 'experiment-' + str(os.getpid()) + '.log'
    logging.basicConfig(filename = logfile, level = logging.INFO)
    
    #defult values: 
    
    minPulls = 10
    maxPulls = 10000
    pullStep = 2
    repeatitions = 1000
    randomBandit = False
    numBandits = 4
    randomSeed = -1

    try:
        opts, args = getopt.getopt(argsList, "", ["min=", "max=", "step=", "repeat=", "random=", "bandits=", "seed="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
    for opt, arg in opts:
        if opt == '--min':
            minPulls = int(arg)
        elif opt == '--max':
            maxPulls = int(arg)
        elif opt == '--step':
            pullStep = float(arg)
        elif opt == '--repeat':
            repeatitions = int(arg)
        elif opt == '--random':
            randomBandit = ast.literal_eval(arg)
        elif opt == '--bandits':
            numBandits = int(arg)
        elif opt == '--seed':
            randomSeed = int(arg)

    avgFileName = args[0]
    algoNames = args[1:] # may be few algorithms

    if randomSeed != -1:
        random.seed(randomSeed)  # use a constant seed

    if randomBandit:
        bandit = RandomBandit(numBandits)
    else:
        avgList = readAvgFromFile(avgFileName)
        bandit = Bandit(avgList)

    experimentMainLoop(minPulls, maxPulls, pullStep, repeatitions, bandit, algoNames, randomBandit)
    
    
if __name__ == "__main__":
    banditExperiment(sys.argv[1:])
