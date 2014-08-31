#!/usr/bin/env python
import os
import re
import argparse
from itertools import izip
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio import SeqIO

def main(argv=None):
#if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help='input GTF file')
    parser.add_argument('fasta', help='input FASTA file')
    parser.add_argument('cgmap', help='input CGmap file')
    parser.add_argument('-d', '--depth', type=int, default=4, help='minimum read depth')
    parser.add_argument('-p', '--promoterSize', type=int, default=1000, help='promoter size')
    parser.add_argument('-c', '--context', type=str, default='CG', help='methylation context')
    args = parser.parse_args()

    with open(args.fasta) as infile:
        fasta = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'))
        for id in fasta:
            fasta[id] = str(fasta[id].seq).upper()

    ctx_map = {}
    for id in fasta:
        ctx_map[id] = ['-']*len(fasta[id])
        cg = [match.start() for match in re.finditer(r'(?=(CG))', fasta[id])]
        for pos in cg:
            ctx_map[id][pos] = 'X'
        chg = [match.start() for match in re.finditer(r'(?=(C[ACT]G))', fasta[id])]
        for pos in chg:
            ctx_map[id][pos] = 'Y'
        chh = [match.start() for match in re.finditer(r'(?=(C[ACT][ACT]))', fasta[id])]
        for pos in chh:
            ctx_map[id][pos] = 'Z'
        rcg = [match.start()-1 for match in re.finditer(r'(?<=(CG))', fasta[id])]
        for pos in rcg:
            ctx_map[id][pos] = 'x'
        rchg = [match.start()-1 for match in re.finditer(r'(?<=(C[AGT]G))', fasta[id])]
        for pos in rchg:
            ctx_map[id][pos] = 'y'
        rchh = [match.start()-1 for match in re.finditer(r'(?<=([AGT][AGT]G))', fasta[id])]
        for pos in rchh:
            ctx_map[id][pos] = 'z'

    for id in ctx_map:
        ctx_map[id] = ''.join(ctx_map[id])

    print "Fasta processed"

    refdict = ctx_map
    with open(args.cgmap) as infile:
        cgmap_dict = {}
        for id in refdict.keys():
            cgmap_dict[id] = ['-' for _ in xrange(len(refdict[id]))]
        for line in infile:
            line = line.strip().split()
            chr = line[0]
            pos = int(line[2]) - 1 # Transfer to 0-based
            context = line[3]
            level = float(line[5])
            depth = int(line[7])
            if context in ['CG', 'CHG', 'CHH'] and depth >= args.depth:
                cgmap_dict[chr][pos] = level

    print "CG map processed" 

    ctxs = {'CG': ('x', 'X'), 'CHG': ('y', 'Y'), 'CHH': ('z', 'Z')}
    inv_ctxs = {'X': 'CG', 'Y': 'CHG', 'Z': 'CHH'}
    ref = ctx_map
    cgmap = cgmap_dict

    l = defaultdict(list)
    for chr in set(ref) & set(cgmap):
        for tag, mlevel in izip(ref[chr], cgmap[chr]):
            tag = tag.upper()
            if tag in inv_ctxs and mlevel != '-':
                l[inv_ctxs[tag]].append(mlevel)
 
    geneInfoList = []
    geneName = None
    with open(args.gtf) as infile:
        for line in infile:
            line = line.strip().split()
            (chr, start, end, direction, geneInfo) = (line[0], int(line[3]), int(line[4]), line[6], line[9].split('"')[1])
            if geneInfo != geneName: #new gene
                if geneName != None:
                    geneInfoList += [(geneName, chr, minBound, maxBound, direction)]
                minBound = start
                geneName = geneInfo
                maxBound = end
            else: #same gene
                if end > maxBound:
                    maxBound = end

    gene_mlevel = defaultdict(lambda: defaultdict(float))
    for ctx in ctxs:
        for t in geneInfoList:
            m = []
            if t[1] in ref and t[1] in cgmap:
                for tag, mlevel in izip(ref[t[1]][t[2]:t[3]], cgmap[t[1]][t[2]:t[3]]):
                    if tag in ctxs[ctx] and mlevel != '-':
                        m.append(mlevel)
                if len(m) > 0:
                    gene_mlevel[t[0]][ctx] = np.mean(m)

    promoterInfo = []
    for (geneName, chr, minBound, maxBound, direction) in geneInfoList:
        promoterSize = args.promoterSize
        if direction==('+'):
            (promMinBound, promMaxBound) = (minBound-promoterSize,minBound)
        else:
            (promMinBound, promMaxBound) = (maxBound, maxBound+promoterSize)
        promoterInfo += [(geneName,chr, promMinBound, promMaxBound, direction)]

    prom_mlevel = defaultdict(lambda: defaultdict(float)) 
    for ctx in ctxs:
        for elem in promoterInfo:
            s = []
            if elem[1] in ref and elem[1] in cgmap:
                for tag, mlevel in izip(ref[elem[1]][elem[2]:elem[3]], cgmap[elem[1]][elem[2]:elem[3]]):
                    if tag in ctxs[ctx] and mlevel != '-':
                        s.append(mlevel)
                if len(s) > 0:
                    prom_mlevel[elem[0]][ctx] = np.mean(s)

    res = []
    for geneInfo in geneInfoList:
        ctx = args.context
        res += [(geneInfo,gene_mlevel[geneInfo[0]][ctx],prom_mlevel[geneInfo[0]][ctx])]
    res = sorted(res, key=lambda t: t[-1], reverse=True)
    
    root = os.path.splitext(os.path.basename(args.cgmap))[0]
    with open('{}.methLevels.{}.txt'.format(root, args.context), 'w') as writeFile:
        writeFile.write('#CG:{}\tCHG:{}\tCHH:{}\n'.format(*[np.mean(l[ctx]) for ctx in ['CG', 'CHG', 'CHH']]))
        for t in res:
            writeFile.write('\t'.join(map(str, t[0])))
            writeFile.write('\t{:.3f}'.format(t[1]))
            writeFile.write('\t{:.3f}'.format(t[2]))
            writeFile.write('\n')
