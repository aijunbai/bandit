### Bandit algorithms

"""Bandit algorithms

An algorithm receives a list of lists of pull outcomes
and returns the index of the next arm to pull."""

from math import sqrt
from random import random, randrange, betavariate

# def mean(a):
#     return (float(a.count(1)) / len(a)) if len(a) > 0 else 0.0

def maxIndex(a):
    return max(enumerate(a), key=lambda x: x[1])[0]

def ThompsonSampling(stats):
    alphabeta = [(1 + s[0] * s[1], 1 + s[0] - s[0] * s[1]) for s in stats]
    samples = [(betavariate(ab[0], ab[1]) if (ab[0] + ab[1] > 2) else float("inf")) for ab in alphabeta]
    
    return maxIndex(samples)

def RoundRobin(stats):
    """ receives a lists of arm pulling results and return a number of arm to pull
    just for testing"""
    k = len(stats)
    n = sum(s[0] for s in stats) # total number of pulls

    return n % k

def Randomized(stats):
    """ receives a lists of arm pulling results and return a number of arm to pull
    just for testing"""
     
    k = len(stats)
    return randrange(k)

def UCB(stats):
    """Upper Confidence Bounds 1"""
    k = len(stats)
    n = sum(s[0] for s in stats) # total number of pulls
    ibest = -1
    xbest = -1
    
    for i in range(k):        
        ni = stats[i][0]
        if ni == 0:     # if there is a completely unexplored arm,
            return i  # return the arm
        x = stats[i][1] + sqrt(2.0 * sqrt(n) / ni)
        if x > xbest: # otherwise, maximize the upper confidence bound
            xbest = x
            ibest = i
    return ibest

def Greedy(stats):
    """0.5-greedy: pull best arm with probability 0.5,
    rest of the arms with equal probability"""
    k = len(stats)
    imax = -1
    xmax = -1
    
    if random() > 0.5:
        for i in range(k):
            if stats[i][0] == 0:
                return i
            x = stats[i][1]
            if x > xmax:
                imax = i
                xmax = x
        return imax
    else:
        return randrange(k)
    
