#!/usr/bin/env python
from __future__ import division
import os
import re
import string
import argparse
from collections import Counter, defaultdict
from itertools import izip, product, combinations
import pysam
import numpy as np
import scipy as sp
from scipy.stats import fisher_exact
import pandas as pd
from Bio import SeqIO
from pyfasta import Fasta
import snp
"""
def flipBases(section): #section: one read from readBaseList
    table = string.maketrans('ATGCN', 'TACGN')    
    flipped = string.translate(section[1], table)
    newBase = string.tranlate(section[2], table)
    return (section[0],flipped,newBase,False,section[4],section[5])
"""
def calcRange(parentGroup):
    #min:readStart of position, max:readStart+readLength from that position
    sortedParentGroup = sorted(parentGroup,key=lambda x: x[4]) #sort by starting position
    lowerBound = sortedParentGroup[0][4] #first read, start position is 5th in read info
    higherBoundRead = sortedParentGroup[-1]
    higherBound = higherBoundRead[4] + len(higherBoundRead[1])
    return lowerBound,higherBound

def compareParentGroups(parentGroup1,parentGroup2):
    (lowerBound1,higherBound1) = calcRange(parentGroup1)
    (lowerBound2,higherBound2) = calcRange(parentGroup2)
    lowerBound = min(lowerBound1,lowerBound2)
    higherBound = max(higherBound1,higherBound2)
    return (lowerBound,higherBound)

def CsWithinRange(chr,lowerBound,higherBound,openFile):
    #openFile = Fasta(genomeFile)
    range = openFile[chr][lowerBound:higherBound]
    rangeStr = str(range)
    cPositionListInString = []
    for c in re.finditer('C',rangeStr):
        cPositionListInString += [c.start()]
    cPositionList = [(i+lowerBound,'C') for i in cPositionListInString]
    return cPositionList,openFile[chr][lowerBound:higherBound+1]

def GsWithinRange(chr,lowerBound,higherBound,openFile):
    #openFile = Fasta(genomeFile)
    range = openFile[chr][lowerBound:higherBound]
    rangeStr = str(range)
    gPositionListInString = []
    for g in re.finditer('G',rangeStr):
        gPositionListInString += [g.start()]
    gPositionList = [(i+lowerBound,'G') for i in gPositionListInString]
    return gPositionList,openFile[chr][lowerBound:higherBound+1]

def getBases(parentGroup,position):
    listOfBasesAtPosition = []
    for (pos,read,base,isReverse,readStart,readLength) in parentGroup:
        #account for if position before a read starts or if cPosition after a read ends
        if (position-readStart) >= 0 and (position-readStart) < readLength:
            listOfBasesAtPosition += [read[position-readStart]]
    return listOfBasesAtPosition

def methylationCalc(listOfBasesAtPosition):
    cCount = float(listOfBasesAtPosition.count('C'))
    tCount = float(listOfBasesAtPosition.count('T'))
    if tCount + cCount != 0:
        methylationFraction = cCount/(tCount+cCount)
    else:
        methylationFraction = 0
    return cCount,tCount,round(methylationFraction,3)

def fisherTest(cCount1,cCount2,tCount1,tCount2):
    oddsratio, pvalue = fisher_exact([[cCount1, cCount2], [tCount1, tCount2]])
    return pvalue

def epimutations(bamfile,genomeFile,coverage=5,majorAlleleFreq=0.9,buffer=0.1,minRead=3):
    bam = pysam.Samfile(bamfile,'rb')
    openFile = Fasta(genomeFile)
    root = os.path.splitext(os.path.basename(bamfile))[0]
    epimutations = open('{}.epimutations.txt'.format(root), 'w')
    hetList = []
    readBaseList = []
    newReadBaseList = []
    negParentGroup1 = []
    negParentGroup2 = []
    posParentGroup1 = []
    posParentGroup2 = []
    cPositionList = []
    gPositionList = []
    for col in bam.pileup():
        chr = bam.getrname(col.tid)
        pos = col.pos
        refNuc = openFile[chr][pos]
        baseList = []
        readList = []
        for pileupRead in col.pileups:
            isReverse = pileupRead.alignment.is_reverse
            read = pileupRead.alignment.seq
            readStart = pileupRead.alignment.pos
            readLength = pileupRead.alignment.rlen
            base = pileupRead.alignment.seq[pileupRead.qpos]
            baseList += [(base,isReverse)]
            readBaseList += [(pos,read,base,isReverse,readStart,readLength)]
        (aCount,tCount,cCount,gCount,total,posCov,negCov) = snp.alleleCount(baseList,refNuc)
        if total >= coverage:
            (freqList,freqDict,allele1,allele2,snpType) = snp.snpDetermine(aCount,tCount,cCount,gCount,total,refNuc,majorAlleleFreq,buffer)
            if snpType == 'het': #epimutations-only heterozygous snps
                hetList += [pos]       
                for (pos,read,base,isReverse,readStart,readLength) in readBaseList:
                    if pos == hetList[-1] and isReverse:
                        if base == allele1:
                            negParentGroup1 += [(pos,read,base,isReverse,readStart,readLength)]
                        elif base == allele2:
                            negParentGroup2 += [(pos,read,base,isReverse,readStart,readLength)]
                    elif pos == hetList[-1] and isReverse == False:
                        if base == allele1:
                            posParentGroup1 += [(pos,read,base,isReverse,readStart,readLength)]
                        elif base == allele2:
                            posParentGroup2 += [(pos,read,base,isReverse,readStart,readLength)]
                #for finding bounds of reads at each position
                openFile = Fasta(genomeFile)
                if posParentGroup1 != [] and posParentGroup2 != []:
                    (posLowerBound,posHigherBound) = compareParentGroups(posParentGroup1,posParentGroup2)
                    (cPositionList,refNucSequence) = CsWithinRange(chr,posLowerBound,posHigherBound,openFile)
                    #print refNucSequence
                    #print cPositionList,refNucSequence
                if negParentGroup1 != [] and negParentGroup2 != []:
                    (negLowerBound,negHigherBound) = compareParentGroups(negParentGroup1,negParentGroup2)
                    (gPositionList,refNucSequence) = GsWithinRange(chr,negLowerBound,negHigherBound,openFile)
                    #print refNucSequence
                posLine = ('SNP location:  '+'%s\n')%pos
                epimutations.write(posLine)
                #get bases at each position where C or G is a reference nucleotide, calculate methylation levels
                methylationLevel1 = 0
                methylationLevel2 = 0
                pvalue = 0
                if (posParentGroup1 != [] and posParentGroup2 != []) or (negParentGroup1 != [] and negParentGroup2 != []):
                    positionList = sorted(cPositionList+gPositionList)
                    for position in positionList:
                        if position[1] == 'C':
                            listOfBasesAtCPosition1 = getBases(posParentGroup1,position[0])
                            listOfBasesAtCPosition2 = getBases(posParentGroup2,position[0])
                            #print position[0],listOfBasesAtCPosition1, listOfBasesAtCPosition2
                            if len(listOfBasesAtCPosition1) >= minRead and len(listOfBasesAtCPosition2) >= minRead:
                                (cCount1,tCount1,methylationLevel1) = methylationCalc(listOfBasesAtCPosition1)
                                (cCount2,tCount2,methylationLevel2) = methylationCalc(listOfBasesAtCPosition2)
                                pvalue = fisherTest(cCount1,cCount2,tCount1,tCount2)
                                #eachLine = ('%s\t\t%-5s\t%-5s\t%s\n')%(position[0],methylationLevel1,methylationLevel2,pvalue)
                        elif position[1] == 'G':
                            listOfBasesAtGPosition1 = getBases(negParentGroup1,position[0])
                            listOfBasesAtGPosition2 = getBases(negParentGroup2,position[0])
                            if len(listOfBasesAtGPosition1) >= minRead and len(listOfBasesAtGPosition2) >= minRead:
                                (cCount1,tCount1,methylationLevel1) = methylationCalc(listOfBasesAtGPosition1)
                                (cCount2,tCount2,methylationLevel2) = methylationCalc(listOfBasesAtGPosition2)
                                pvalue = fisherTest(cCount1,cCount2,tCount1,tCount2)
                        eachLine = ('%s\t\t%-5s\t%-5s\t%s\n')%(position[0],methylationLevel1,methylationLevel2,pvalue)
                        epimutations.write(eachLine) 
                negParentGroup1 = []
                negParentGroup2 = []
                posParentGroup1 = []
                posParentGroup2 = []
                cPositionList = []
                gPositionList = []
    epimutations.close()

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('bamfile',help='input BAM file')
    parser.add_argument('-c', '--coverage', type=int, default=5, help='coverage or minimum number of reads desired')
    parser.add_argument('-m', '--majorAlleleFreq',type=float, default=0.9, help='frequency to be considered homozygous allele')
    parser.add_argument('-b', '--buffer',type=float,default=0.1, help='buffer on either side of 0.5 to be considered heterozygous allele')
    parser.add_argument('-r', '--minRead',type=int,default=3, help='minimum number of reads for epimutation detection')
    parser.add_argument('-g', '--genomeFile')
    return parser

def main(argv=None):
#if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    epimutations(args.bamfile,args.genomeFile,args.coverage,args.majorAlleleFreq,args.buffer,args.minRead)    
