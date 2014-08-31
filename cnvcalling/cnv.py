#!/usr/bin/env python
from __future__ import division
import os
import argparse
from collections import defaultdict
import pysam
import numpy as np
import scipy as sp
from scipy import stats
from scipy.stats import poisson
import matplotlib.pyplot as plt

def findThreeInARow(possibleCNVList):
    cnvList = []
    for (chr, start, end, baseSum) in possibleCNVList:
        primaryElem = (chr, start, end, baseSum)
        firstIndex = possibleCNVList.index((chr, start,end, baseSum))
        for (chr2, start2, end2, baseSum2) in possibleCNVList[firstIndex:]:
            if (end == start2) and (chr == chr2):
                secondaryElem = (chr2, start2, end2, baseSum2)
                secondIndex = possibleCNVList.index((chr2, start2,end2, baseSum2))
                for (chr3, start3, end3, baseSum3) in possibleCNVList[secondIndex:]:
                    if (end2 == start3) and (chr2 == chr3):
                        tertiaryElem = (chr3, start3, end3, baseSum3)
                        cnvList += [(primaryElem,secondaryElem,tertiaryElem)]
    return cnvList

def cnv(bamfile, windowSize, probability, genomeindexfile):
    bam = pysam.Samfile(bamfile, 'rb')
    chr = None
    allBaseSums = []
    listForArray = []
    zScores = []
    zScoreIndex = []
    giantList = []
    possibleCNVList = []
    upperLimit = 0
    lowerLimit = 0
    for col in bam.pileup():
        pos = col.pos
        cov = col.n
        if bam.getrname(col.tid) != chr:
            chr = bam.getrname(col.tid)
            baseSum = cov
            numw = int(pos/windowSize) #which window
        else:
            if int(pos/windowSize) == numw:
                baseSum += cov
            else:
                if baseSum >= windowSize*3:
                    allBaseSums += [(chr, numw*windowSize, (numw+1)*windowSize, baseSum)]
                    listForArray += [baseSum]
                numw = int(pos/windowSize)
                baseSum = cov
    if baseSum >= windowSize*3:
        allBaseSums += [(chr, numw*windowSize, (numw+1)*windowSize, baseSum)]
        listForArray += [baseSum]
    average = np.mean(listForArray) #new lambda for poisson distribution
    upperLimit = np.percentile(listForArray, probability*100)
    lowerLimit = np.percentile(listForArray, (1-probability)*100)
    cutOff = poisson.ppf(probability, average)
    for (chr, start, end, baseSum) in allBaseSums:
        if baseSum > cutOff:
            possibleCNVList += [(chr, start, end, baseSum)]
    """
    stdDev = np.std(listForArray)
    
    for (chr, start, end, baseSum) in allBaseSums:
        zScoreTemp = ((baseSum - average)/float(stdDev))
        giantList += [(chr, start, end, baseSum, zScoreTemp)] #all windows' base sums and z scores
        if (zScoreTemp >= determinedZScore) or (zScoreTemp <= -determinedZScore):
            possibleCNVList += [(chr, start,end,baseSum,zScoreTemp)]
    """
    cnvList = findThreeInARow(possibleCNVList)
    return allBaseSums, cnvList, upperLimit, lowerLimit

def plotCNV(bamfile, windowSize, probability, windowNumber, genomeindexfile):
    (allBaseSums, cnvList, upperLimit, lowerLimit) = cnv(bamfile, windowSize, probability, genomeindexfile)
    chrLengths = {}
    chromList = []
    newGiantList = []
    with open(genomeindexfile) as infile:
        for (rownum, line) in enumerate(infile):
            line = line.strip().split('\t')
            chrom = line[0]
            chromList += [chrom] #list of chromosomes in order
            length = line [1]
            chrLengths[chrom] = length #dictionary of chromosome lengths
    xticks = set()
    xticks.add(0)
    fig = plt.figure(figsize=(len(chromList),4))
    ax = fig.add_axes((0.08, 0.2, 0.9, 0.7))
    cum = 0
    for (chr, start, end, baseSum) in allBaseSums:
        #baseSum = np.log(baseSum)
        indexChr = chromList.index(chr)
        for chrom in chromList[:indexChr]:
            cum += int(chrLengths[chrom])
            xticks.add(cum)
        ax.scatter(np.mean([cum+start,cum+end]),baseSum,c=(184/255,184/255,180/255),linewidth=0,alpha=0.8,s=10)
        newGiantList += [(chr,start+cum,end+cum,baseSum)]
        cum = 0
    totalCum = sum(map(int, chrLengths.values()))
    xticks.add(totalCum)
    xticks = sorted(list(xticks))
    ax.set_xticks([(xticks[i]+xticks[i+1])/2.0 for i in xrange(len(xticks)-1)])
    ax.set_xticklabels(chromList, size='medium', weight='black')
    ax.set_xlabel('Chromosome', size='xx-large', weight='black')
    ax.set_ylabel('Coverage', size='xx-large', weight='black')
    ax.set_xlim(0, totalCum)
    #ax.set_yscale('log')
    ax.set_ylim((lowerLimit,upperLimit))
    
    #moving average
    begin = 0
    average = []
    #begin and finish for encompassing window of windows in which average will be taken
    whichWindow = None
    windowList = []
    giantDict = defaultdict(list)
    for (chr, start, end, baseSum) in newGiantList: #cumulative start and end
        giantDict[chr].append((start, end, baseSum))
    for chr in chromList:
        for (start, end, baseSum) in sorted(giantDict[chr]): #cumulative start and end
            if int(start/(windowSize*windowNumber)) == whichWindow:
                baseSumSum += baseSum
                #zScoreSum += zScore
            else:
                if whichWindow != None:
                    average += [(((whichWindow*(windowSize*windowNumber))+(windowSize*windowNumber)/2),baseSumSum/windowNumber)]
                    windowList = []
                   #average += [(((whichWindow*(windowSize*windowNumber))+(windowSize*windowNumber)/2),zScoreSum/windowNumber)] #position,average
                   #windowList = []
                whichWindow = int(start/(windowSize*windowNumber))
                baseSumSum = baseSum
    ax.plot([t[0] for t in average],[t[1] for t in average],color='k',linewidth=1.5,linestyle='-')
    #ax.set_yscale('log')
    root = os.path.splitext(os.path.basename(bamfile))[0]
    plt.savefig('{}.CNVplot.png'.format(root), dpi=300)
    plt.close(fig)
    
    with open('{}.CNVvalues.txt'.format(root), 'w') as writeFile:
        for elem in cnvList:
            writeFile.write('{}\n'.format(elem))

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('bamfile', help='input BAM file')
    parser.add_argument('-w', '--windowSize', type=int, default=2000, help='window size in bases')
    parser.add_argument('-p', '--probability', type=float, default=0.95, help='probably for results') 
    parser.add_argument('-n', '--windowNumber', type=int, default=25, help='number of windows for moving average')
    parser.add_argument('-i', '--genomeindexfile',help='input genome index')
    return parser

def main(argv=None):
    parser = get_parser()
    args = parser.parse_args()
    plotCNV(args.bamfile, args.windowSize, args.probability, args.windowNumber, args.genomeindexfile)
