import numpy as np
import scipy.stats as stat 
from copy import copy, deepcopy
import math

#1A
disease_status = np.random.binomial(1,0.1,1000) #only one trial since one chromosome per individual
disease_status.sort() # first 900 people has 0, so they don't have disease 
print "1A"
print "There are", list(disease_status).count(1), "individuals with the disease in the sample."

#1B
dataIn = np.random.random_integers(0,1,(1000,200)) ## make 1000 rows, 200 columns, entries 0/1 
hasMinor = np.random.binomial(1,.25,900) # has minor allele in the controls 
hasMinor = np.append ( hasMinor, np.random.binomial(1,0.95,100) ) # has minor allele in the cases 
dataIn[:,0] = hasMinor # assume first SNP is causal, so the 0/1 is not randomly distributed. 

#Sanity Check
print "1B"
print "Control Group", sum ( dataIn [0:899,0] )/900. # the dot. is needed when int divides int, to cast into decimal 
print "Case Group",sum ( dataIn [900:999,0] )/100. 

#1C
def oneC(dataIn):
    maxStat = 0
    for i in range(0,200):
        pplus = sum(dataIn[900:999, i])/100.
        pminus = sum(dataIn[0:899, i]) / 900.
        pa = (pminus + pplus) / 2.
        ncp = (pplus - pminus) / math.sqrt(2*pa*(1-pa)) * math.sqrt(2*(100*900)/(100+900))
        if(ncp > maxStat):
            maxStat = ncp
    
    alpha = .05/200.
    z = stat.norm.ppf(alpha/2)
    if(maxStat > -1 * z or maxStat < z):
        print "It is significant."
    else:
        print "It is not significant."

print "1C"
oneC(dataIn)

#1D
#np.corrcoef(dataIn[:,1],dataIn[:,2]) # correlation of snps s2 and s3 (remember, in python, indexing starts at zero.)
def oneD(dataIn):
    corMatrix = np.corrcoef( dataIn.transpose() ) ## correlation of the 200 Snps in the data. 
    np.shape(corMatrix) # should be dim of 200x200
    
    (np.sum(np.abs(corMatrix)>0.1) - 200) / 2 # number of pairs that have absolute correlation over 0.1 
    
    index = np.where(np.abs(corMatrix)>0.1) # row/column index of where the abs cor is over 0.1 
    rowIndex = list(index[0])
    colIndex = list(index[1])
    
    #print len(rowIndex)
    currentIndex = 0
    #remove duplicates
    while currentIndex < len(rowIndex):
        if(rowIndex[currentIndex] > colIndex[currentIndex]):
            rowIndex.pop(currentIndex)
            colIndex.pop(currentIndex)
        else:
            currentIndex = currentIndex + 1
            
    #print len(rowIndex)
    
    tagArray = []
    
    while len(rowIndex) > 0:
        removeArray = []
        maxCorr = 0
        maxSNP = 0
        for i in range(0,200):
            numCorr = rowIndex.count(i) + colIndex.count(i)
            if(numCorr > maxCorr):
                maxCorr = numCorr
                maxSNP = i
        #print "Row",rowIndex
        #print "Col",colIndex     
        #print "Max",maxSNP
        
        tagArray.append(maxSNP)
        currentIndex = 0
        while currentIndex < len(rowIndex):
            if(rowIndex[currentIndex] == maxSNP):
                removeArray.append(colIndex[currentIndex])
                rowIndex.pop(currentIndex)
                colIndex.pop(currentIndex)
            elif(colIndex[currentIndex] == maxSNP):
                removeArray.append(rowIndex[currentIndex])
                rowIndex.pop(currentIndex)
                colIndex.pop(currentIndex)
            else:
                currentIndex = currentIndex + 1
                
        for i in removeArray:
            currentIndex = 0
            while currentIndex < len(rowIndex):
                if(rowIndex[currentIndex] == i):
                    rowIndex.pop(currentIndex)
                    colIndex.pop(currentIndex)
                elif(colIndex[currentIndex] == i):
                    rowIndex.pop(currentIndex)
                    colIndex.pop(currentIndex)
                else:
                    currentIndex = currentIndex + 1
                    
    print tagArray

print "1D"
oneD(dataIn)

#1E
def oneE():
    pplus = .95
    pminus = .25
    pa = (pplus + pminus) / 2.
    ncp = (pplus - pminus) / math.sqrt(2*pa*(1-pa)) * math.sqrt(2*(100*900)/(100+900))
    power = stat.norm.cdf(stat.norm.ppf(.05/2) + ncp) + 1 - stat.norm.cdf(-1 * stat.norm.ppf(.05/2) + ncp)
    print "Analytically:", power
    
def oneESim(dataIn):
    pplus = sum(dataIn[900:999, 0])/100.
    pminus = sum(dataIn[0:899, 0]) / 900.
    pa = (pminus + pplus) / 2.
    ncp = (pplus - pminus) / math.sqrt(2*pa*(1-pa)) * math.sqrt(2*(100*900)/(100+900))
    alpha = .05
    z = stat.norm.ppf(alpha/2)
    if(ncp > -1 * z or ncp < z):
        return True
    else:
        return False
        
print "1E"
oneE()
        
num = 0
for i in range(0, 1000):
    dataIn = np.random.random_integers(0,1,(1000,200)) ## make 1000 rows, 200 columns, entries 0/1 
    hasMinor = np.random.binomial(1,.25,900) # has minor allele in the controls 
    hasMinor = np.append ( hasMinor, np.random.binomial(1,0.95,100) ) # has minor allele in the cases 
    dataIn[:,0] = hasMinor # assume first SNP is causal, so the 0/1 is not randomly distributed. 
    if(oneESim(dataIn)):
        num = num + 1

print "Simulation:", num/1000
    
#2A
dataIn = np.random.random_integers(0,1,(1000,200)) ## make 1000 rows, 200 columns, entries 0/1 
hasMinor = np.random.binomial(1,.25,900) # has minor allele in the controls 
hasMinor = np.append ( hasMinor, np.random.binomial(1,0.95,100) ) # has minor allele in the cases 
dataIn[:,0] = hasMinor # assume first SNP is causal, so the 0/1 is not randomly distributed.

for i in range(0,180):
    dataIn[i] = dataIn[0]
    
for i in range(900, 920):
    dataIn[i] = dataIn[900]
  
print "2A" 
oneC(dataIn)
oneD(dataIn)
oneE()

num = 0
for i in range(0, 1000):
    dataIn = np.random.random_integers(0,1,(1000,200)) ## make 1000 rows, 200 columns, entries 0/1 
    hasMinor = np.random.binomial(1,.25,900) # has minor allele in the controls 
    hasMinor = np.append ( hasMinor, np.random.binomial(1,0.95,100) ) # has minor allele in the cases 
    dataIn[:,0] = hasMinor # assume first SNP is causal, so the 0/1 is not randomly distributed.
    
    for i in range(0,180):
        dataIn[i] = dataIn[0]
        
    for i in range(900, 920):
        dataIn[i] = dataIn[900]
    if(oneESim(dataIn)):
        num = num + 1

print "Simulation:", num/1000

#2B
dataIn = np.random.random_integers(0,1,(1000,200)) ## make 1000 rows, 200 columns, entries 0/1 
hasMinor = np.random.binomial(1,.25,900) # has minor allele in the controls 
hasMinor = np.append ( hasMinor, np.random.binomial(1,0.95,100) ) # has minor allele in the cases 
dataIn[:,0] = hasMinor # assume first SNP is causal, so the 0/1 is not randomly distributed.

for i in range(0,450):
    dataIn[i] = dataIn[0]
    
for i in range(900, 950):
    dataIn[i] = dataIn[900]

print "2B"
oneC(dataIn)
oneD(dataIn)
oneE()

num = 0
for i in range(0, 1000):
    dataIn = np.random.random_integers(0,1,(1000,200)) ## make 1000 rows, 200 columns, entries 0/1 
    hasMinor = np.random.binomial(1,.25,900) # has minor allele in the controls 
    hasMinor = np.append ( hasMinor, np.random.binomial(1,0.95,100) ) # has minor allele in the cases 
    dataIn[:,0] = hasMinor # assume first SNP is causal, so the 0/1 is not randomly distributed.
    
    for i in range(0,450):
        dataIn[i] = dataIn[0]
        
    for i in range(900, 950):
        dataIn[i] = dataIn[900]

    if(oneESim(dataIn)):
        num = num + 1

print "Simulation:", num/1000

#2C
def breed (dataIn,i,j): 
    ''' i mates j, make person k '''
    personI = dataIn[i,:]
    personJ = dataIn[j,:]
    ''' cross over happens '''
    flip = personI == personJ # same allele will be passed on. 
    flip = np.where(flip==False) # any position with different allele, the allele will be chosen randomly 
    # print ( np.shape(flip)[1] ) # how many times you need to do the flips 
    personK = deepcopy(personI) ## WILL NOT WORK . if you do personK = personI without deepcopy 
#!! this is a bit incorrect, the true process of 'cross-over' is not simple fliping of bits, usually, "chunks" are cross-over 
#!! you can modify this process. 
    ''' much more complicated if alleles are not independent '''
    personK[flip] = np.random.binomial( 1,.5,np.shape(flip)[1] )
    return personK
    
#offspring1 = breed(dataIn,1,2)    # first offspring of pair 1,2
#np.corrcoef(dataIn[1,:],offspring1)
#np.corrcoef(dataIn[2,:],offspring1)
#
#offspring2 = breed(dataIn,1,2)    # second offspring of pair 1,2
#np.corrcoef(dataIn[1,:],offspring2)
#np.corrcoef(dataIn[2,:],offspring2)
#
#np.corrcoef(offspring1,offspring2) # correlation between sibs. 

dataIn = np.random.random_integers(0,1,(1000,200)) ## make 1000 rows, 200 columns, entries 0/1 
hasMinor = np.random.binomial(1,.25,900) # has minor allele in the controls 
hasMinor = np.append ( hasMinor, np.random.binomial(1,0.95,100) ) # has minor allele in the cases 
dataIn[:,0] = hasMinor # assume first SNP is causal, so the 0/1 is not randomly distributed.

for i in range(2,180):
    dataIn[i] = breed(dataIn,0,1)
    
for i in range(902, 920):
    dataIn[i] = breed(dataIn,900,901)
    
print "2C"
oneC(dataIn)
oneD(dataIn)
oneE()

num = 0
for i in range(0, 1000):
    dataIn = np.random.random_integers(0,1,(1000,200)) ## make 1000 rows, 200 columns, entries 0/1 
    hasMinor = np.random.binomial(1,.25,900) # has minor allele in the controls 
    hasMinor = np.append ( hasMinor, np.random.binomial(1,0.95,100) ) # has minor allele in the cases 
    dataIn[:,0] = hasMinor # assume first SNP is causal, so the 0/1 is not randomly distributed.
    
    for i in range(2,180):
        dataIn[i] = breed(dataIn,0,1)
        
    for i in range(902, 920):
        dataIn[i] = breed(dataIn,900,901)
    if(oneESim(dataIn)):
        num = num + 1

print "Simulation:", num/1000