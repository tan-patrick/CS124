import numpy as np
import scipy.stats as stat 
from copy import copy, deepcopy
import math
import time
import timeit

def createHaplotypes(numSNPs, probHet):
    separation = np.random.geometric(probHet, numSNPs)
    haplotype1 = []
    haplotype2 = []
    curLoc = 0
    for i in separation:
        curLoc = curLoc + i
        if(np.random.randint(0,2) == 0):
            haplotype1.append([0, curLoc])
            haplotype2.append([1, curLoc])
        else:
            haplotype1.append([1, curLoc])
            haplotype2.append([0, curLoc])
    return [haplotype1, haplotype2, curLoc]

def getReads(haplotypes, numReads, readLength):
    maxLength = haplotypes[2]
    reads = []
    loc = []
    for i in range(0, numReads):
        startAt = np.random.randint(0, maxLength) - readLength / 2
        loc.append(startAt)
        endAt = startAt + readLength
        curSNP = 0
        read = []
        if(np.random.randint(0,2) == 0):
            for i in haplotypes[0]:
                if(i[1] >= startAt):
                    read.append(curSNP)
                    while(curSNP < len(haplotypes[0]) and haplotypes[0][curSNP][1] < endAt):
                        read.append(haplotypes[0][curSNP][0])
                        curSNP = curSNP + 1
                    break
                curSNP = curSNP + 1
        else:
            for i in haplotypes[1]:
                if(i[1] >= startAt):
                    read.append(curSNP)
                    while(curSNP < len(haplotypes[1]) and haplotypes[1][curSNP][1] < endAt):
                        read.append(haplotypes[1][curSNP][0])
                        curSNP = curSNP + 1
                    break
                curSNP = curSNP + 1 
        reads.append(read)
    return [reads, loc]

def bruteForce(reads, numSNPs):
    for num in range(0, pow(2, numSNPs)):
        binary = format(num, '0' + str(numSNPs) + 'b')
        haplotype1 = list(binary)
        haplotype1 = [int(i) for i in haplotype1]
        haplotype2 = []
        nextBin = False
        for i in haplotype1:
            if(i == 0):
                haplotype2.append(1)
            else:
                haplotype2.append(0)
        #print "-----------"
        #print haplotype1
        #print haplotype2
        for read in reads:
            #print read
            success = True
            start = read[0]
            if(len(read) >= 2):
                if(haplotype1[start] == read[1]):
                    for i in range(1, len(read)):
                        if(haplotype1[start + i - 1] != read[i]):
                            success = False
                            break
                else:
                    for i in range(1, len(read)):
                        if(haplotype2[start + i - 1] != read[i]):
                            success = False
                            break
                if(success == False):
                    nextBin = True
                    break
            #print "Round Done"
        if(nextBin == False):
            return [haplotype1, haplotype2]
            
def greedyAssembly(reads, numSNPs):
    haplotype1 = [-1] * numSNPs
    haplotype2 = [-1] * numSNPs
    currentLoc = reads[0][0]
    for i in reads[0][1:]:
        if(i == 0):
            haplotype1[currentLoc] = 0
            haplotype2[currentLoc] = 1
        else:
            haplotype1[currentLoc] = 1
            haplotype2[currentLoc] = 0   
        currentLoc = currentLoc + 1     
        
    #print reads[0]
    #print haplotype1
    #print haplotype2
    #print "Round Done" 
    
    for read in reads[1:]:
        #print read
        
        if(len(read) >= 2):
            currentLoc = read[0]
            if(haplotype1[currentLoc] == -1):
                for i in read[1:]:
                    if(i == 0):
                        haplotype1[currentLoc] = 0
                        haplotype2[currentLoc] = 1
                    else:
                        haplotype1[currentLoc] = 1
                        haplotype2[currentLoc] = 0   
                    currentLoc = currentLoc + 1 
            elif(haplotype1[currentLoc] == read[1]):
                for i in read[1:]:
                    if(i == 0):
                        haplotype1[currentLoc] = 0
                        haplotype2[currentLoc] = 1
                    else:
                        haplotype1[currentLoc] = 1
                        haplotype2[currentLoc] = 0   
                    currentLoc = currentLoc + 1 
            else:
                for i in read[1:]:
                    if(i == 0):
                        haplotype1[currentLoc] = 1
                        haplotype2[currentLoc] = 0
                    else:
                        haplotype1[currentLoc] = 0
                        haplotype2[currentLoc] = 1   
                    currentLoc = currentLoc + 1 
            #print haplotype1
            #print haplotype2
            #print "Round Done"
    for i in range(0,len(haplotype1)):
        if(haplotype1[i] == -1):
            haplotype1[i] = 0
            haplotype2[i] = 1
    return [haplotype1, haplotype2]
                
def numErrors(haplotypes, solution):
    numErrors = 0
    haplotype1 = haplotypes[0]
    haplotype2 = haplotypes[1]
    sol1 = solution[0]
    sol2 = solution[1]
    currentLoc = 0
    currentCheck = 0
    if(sol1[0] == haplotype1[0][0]):
        currentCheck = 0
    else:
        currentCheck = 1
        
    for i in sol1:
        if(currentCheck == 0):
            if(i != haplotype1[currentLoc][0]):
                currentCheck = 1
                numErrors = numErrors + 1
        else:
            if(i != haplotype2[currentLoc][0]):
                currentCheck = 0
                numErrors = numErrors + 1
        currentLoc = currentLoc + 1
    return numErrors

def haplotypeAssembly(numSNPs, percentHet, numReads, readLength):
    haplotypes = createHaplotypes(numSNPs, percentHet)
    readRet = getReads(haplotypes, numReads, readLength)
    loc = readRet[len(readRet) - 1]
    reads = readRet[:-1][0]
    print "Presorted Reads:", reads
    reads = sorted(reads, key=lambda row:row[0])
    h1 = [i[0] for i in haplotypes[0]]
    h2 = [i[0] for i in haplotypes[1]]
    haplotypeOnly = [h1,h2]
    print "Haplotypes:", haplotypes
    print "Reads:", reads
    #print "Read Locs:", loc
    print "Only Haplotypes:", haplotypeOnly
    t = time.clock()
    bruteSol = bruteForce(reads, numSNPs)
    print "Brute Force took", time.clock() - t, "to run"
    print "BruteSol:", bruteSol
    print "Brute Force Errors", numErrors(haplotypes, bruteSol)
    t = time.clock()
    greedySol = greedyAssembly(reads, numSNPs)
    print "Greedy Algorithm took", time.clock() - t, "to run"
    print "GreedySol:", greedySol
    print "Greedy Solution Errors", numErrors(haplotypes, greedySol)

haplotypeAssembly(20, .05, 50, 2)