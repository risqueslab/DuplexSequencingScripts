#!/bin/env python
'''
mutContextSpectrum.py
Version 1.22
by Brendan Kohrn
July 29, 2014

Find the mutation context for a list of mutants in mutpos format, as follows:
Chromosome    ReferenceNT    Position    Depth    number_of_mutants    mutants_to_T    mutants_to_C    mutants_to_G    mutants_to_A    insertions    deletions

'''

import sys
import re
import collections
from collections import defaultdict
import argparse
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt


class Files:
    def __init__(self, inFasta, inMutPos, minC, maxC, bufferSize = 3):
        # DEBUG sys.stderr.write('>>Starting File Initialization\n')
        self.lineNum = 0
        tmpMutPos = open(inMutPos, 'r')
        self.MutPos = mutPosFile(tmpMutPos)
        tmpFasta = open(inFasta, 'r')
        self.Fasta = fastaWin(tmpFasta, bufferSize)
        #self.Fasta.advance(endChr = self.MutPos.chrom, endPos = self.MutPos.pos)
        self.minClonality = minC
        self.maxClonality = maxC

        # DEBUG sys.stderr.write('>>File Initialization Complete\n')
    
    def __iter__(self):
        return(self)
    
    def next(self):
        # DEBUG sys.stderr.write('>>Advancing...\n')
        endTest = False
        while endTest == False:
            if self.MutPos.line.fileEnd == False:
                endTest = self.MutPos.next()
            else:
                break
        if self.MutPos.line.fileEnd == False:
            fastaTest = self.Fasta.advance(endChr = self.MutPos.chrom, endPos = self.MutPos.pos)
            # DEBUG sys.stderr.write('>>Not EOF = %s\n' % fastaTest)
            if fastaTest == False:
                raise StopIteration
            self.lineNum = self.MutPos.lineNum
            #if self.lineNum % 1000 == 0:
                #print('%s lines processed' % self.lineNum)
            # DEBUG sys.stderr.write('>>Advancing complete.\n')
            return(self.MutPos.line)
        else:
            raise StopIteration
    
    def close(self):
        # DEBUG sys.stderr.write('>>Closing Files\n')
        self.MutPos.close()
        self.Fasta.close()
        # DEBUG sys.stderr.write('>>Closing Complete\n')


class mutPosFile:
    def __init__(self, inFile):
        # DEBUG sys.stderr.write('>>Initializing MutPos\n')
        self.file = inFile
        self.chrom = ""
        self.pos = 1
        self.line = mutPosLine("-",0)
        self.lineNum = 0
        #self.next()
    
    def next(self):
        self.lineNum += 1
        self.line = mutPosLine(self.file.readline(), self.lineNum)
        if self.line.fileEnd == False:
            while self.line.fileEnd == False and self.line.Ts == 0 and self.line.Cs == 0 and self.line.Gs == 0 and self.line.As == 0:
                self.lineNum += 1
                self.line = mutPosLine(self.file.readline(), self.lineNum)
            if self.line.fileEnd == False:
                self.pos = self.line.pos
                self.chrom = self.line.chrom
                # DEBUG sys.stderr.write('>>MutPos Advanced to line %s...\n' % self.lineNum)
                return(True)
            else:
                return(False)
        else:
            sys.stderr.write('>>MutPos EOF reached\n')
            return(False)
    
    def __str__(self):
        return(str(self.line))
    
    def close(self):
        self.file.close()
        return(True)


class mutPosLine:
    def __init__(self, inLine, inLineNum):
        if inLine == "-":
            self.fileEnd = False
            self.lineNum = -1
            linebins = inLine.split()
            self.chrom = ''
            self.refBase = ''
            self.pos = -1
            self.depth = -1
            self.muts = -1
            self.Ts = -1
            self.Cs = -1
            self.Gs = -1
            self.As = -1
            self.ins = -1
            self.dels = -1
            #self.Ns = -1
            #self.clonalDepth = self.depth - self.Ns
            self.clonality = 0
        elif inLine == "":
            self.fileEnd = True
            self.lineNum = -1
            linebins = inLine.split()
            self.chrom = ''
            self.refBase = ''
            self.pos = -1
            self.depth = -1
            self.muts = -1
            self.Ts = -1
            self.Cs = -1
            self.Gs = -1
            self.As = -1
            self.ins = -1
            self.dels = -1
            #self.Ns = -1
            #self.clonalDepth = self.depth - self.Ns
            self.clonality = 0
        else:
            self.lineNum = inLineNum
            linebins = inLine.split()
            self.fileEnd = False
            self.chrom = linebins[0]
            self.refBase = linebins[1].upper()
            self.pos = int(linebins[2])
            self.depth = int(linebins[3])
            self.muts = int(linebins[4])
            self.Ts = int(linebins[5])
            self.Cs = int(linebins[6])
            self.Gs = int(linebins[7])
            self.As = int(linebins[8])
            self.ins = int(linebins[9])
            self.dels = int(linebins[10])
            #self.Ns = int(linebins[11])
            #self.clonalDepth = self.depth - self.Ns
            self.clonality = max(float(self.Ts)/self.depth, float(self.Cs)/self.depth, float(self.Gs)/self.depth, float(self.As)/self.depth)
    
    def makeMuts(self):
        outMuts = []
        if self.Ts:
            outMuts.append(Mutation(self.chrom, self.pos, self.refBase, "T"))
        if self.Cs:
            outMuts.append(Mutation(self.chrom, self.pos, self.refBase, "C"))
        if self.Gs:
            outMuts.append(Mutation(self.chrom, self.pos, self.refBase, "G"))
        if self.As:
            outMuts.append(Mutation(self.chrom, self.pos, self.refBase, "A"))
        return(outMuts)
            
    
    def __str__(self):
        return('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (self.lineNum, self.chrom, self.refBase, self.pos, self.depth, self.Ts, self.Cs, self.Gs, self.As, self.ins, self.dels))


class fastaWin:
    """
    This class enables most of the oporations carried out by the 
    program; it handles the rolling window and sequence retrieval 
    sections.  
    """
    
    def __init__(self, dataSource, bufferSize = 3):
        #~ self.debugCtr = 0
        # DEBUG sys.stderr.write('>>Initializing FASTAWIN\n')
        self.sourceFile = dataSource
        
        ##initialize the window in the first chromosome
        line = self.sourceFile.readline().strip().split(">")[1].split(" ")[0]
        self.chrom = line
        
        line = self.sourceFile.readline().strip().upper()
        self.sizeMax=2*len(line) #size of the window
        self.bsize=int(bufferSize) #size of the buffer desired
        self.eof=False #has the end of the file been reached?
        self.really=False #Are you sure?
        #how long is one line in this reference genome anyway?
        self.lLength=len(line) 

        self.data = []
        #do one line of N's to make sure that any mutations near 
        #the begining can be processed
        self.data.extend(list('N'*self.lLength)) 
        self.data.extend(list(line))#load in the first line
        #the leftmost position on the reference genome.  1-indexed
        self.minPos=1-self.lLength 
        #the rightmost position on the reference genome.  1-indexed
        self.maxPos=self.lLength 
        self.pos = 1
    
    def advance(self, endChr = None, endPos = None):
        # DEBUG sys.stderr.write('>>Advancing FASTA to %s:%s\n' % (endChr, endPos))
        if endChr == None and endPos == None:
            self.pos += 1
            if self.pos > self.maxPos - self.bsize:
                self.moveWin()
        elif endChr != None and endPos == None: 
            while self.chrom != endChr and self.chrom != '':
                self.usedChrs.append(self.chrom)
                self.chrom = self.NewChrom()
        elif endPos != None and endChr == None:
            self.pos = endPos
            while self.pos > self.maxPos - self.bsize:
                self.moveWin()
        elif endPos != None and endChr != None:
            while self.chrom != endChr and self.chrom != '':
                self.chrom = self.NewChrom()
            self.pos = endPos
            while self.pos > self.maxPos - self.bsize and self.chrom != '':
                self.moveWin()
        if self.chrom == '':
            return(False)
        else:
            return(True)
    
    def moveWin(self):
        dataSource = self.sourceFile.readline().strip().upper()
        #move the window over one line
        if self.eof == False:
            if dataSource=="":
                self.eof=True
            else:
                inSeq=dataSource
                while len(inSeq) != self.lLength:
                    inSeq+="N"
                self.maxPos+=self.lLength
                self.minPos+=self.lLength
                wPos=self.minPos%self.sizeMax-1
                for base in list(inSeq):
                    try:
                        self.data[wPos]=base
                        wPos+=1
                    except:
                        print('%s\n%s\t%s\t%s\t%s' % (self.data, wPos, base, inSeq, self.lLength))
                        raise
                
            return(True)
        else:
            if self.really == False:
                inSeq="N"*self.lLength
                self.maxPos+=self.lLength
                self.minPos+=self.lLength
                wPos=self.minPos%self.sizeMax-1
                for base in list(inSeq):
                    self.data[wPos]=base
                    wPos+=1 
                self.really=True
            else:
                return(False)
        
    def getSeq(self, pos):
        seq=""
        while int(pos) - 1 + 3 >= self.maxPos:
            fTest=self.moveWin()
            if fTest==False:
                return("")
        for i in range(int(pos) - 1, int(pos) -1 + 3):
            b=(i+self.lLength)%self.sizeMax-1
            if b == -1:
                b = self.sizeMax - 1
            seq+=str(self.data[b])
        #~ print(pos, seq)
        #~ self.debugCtr += 1
        #~ if self.debugCtr == 5:
            #~ exit()
        return(seq)

    def NewChrom(self):
        # DEBUG sys.stderr.write('>>Switching Chromosomes\n')
        # DEBUG sys.stderr.write('>>>>Old Chromosome: %s\n' % self.chrom)
        mvTest=self.sourceFile.readline().strip()
        while ">" not in mvTest and mvTest != "":
            mvTest=self.sourceFile.readline().strip()
        if mvTest != "":
            self.chrom=mvTest.split(">")[1].split(" ")[0]
            self.data = []
            newData = self.sourceFile.readline().strip().upper()
            self.lLength = len(newData)
            self.data.extend(list('N'*self.lLength))
            self.data.extend(list(newData))
            self.minPos=1-self.lLength
            self.maxPos=self.lLength
            self.sizeMax = 2 * self.lLength
        else:
            return('')
        # DEBUG sys.stderr.write('>>>>New Chromosome: %s\n' % self.chrom)
        return(self.chrom)
    
    def close(self):
        self.sourceFile.close()
        return(True)


#~ class mutType:
    #~ def __init__(self, inRef, inMut):
        #~ self.refNT = inRef
        #~ self.mutNT = inMut
    #~ 
    #~ def __str__(self):
        #~ return("%s>%s" % (self.refNT, self.mutNT))


class Mutation:
    def __init__(self, inChrom, inPos, inRef, inMut):
        # DEBUG sys.stderr.write('>>Initializing mutation with MutPos Line %s\n' % inMutPosLine.lineNum)
        self.chrom = inChrom
        self.pos = inPos
        self.Type = "%s>%s" % (inRef, inMut)
        self.context = ''
    
    def setContext(self, fasta):
        # DEBUG sys.stderr.write('>>Getting Context...\n')
        self.context = fasta.getSeq(self.pos)
    
    def __str__(self):
        strToReturn = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (self.chrom, self.pos, self.Type, self.context)
        return(strToReturn)

def revComp(inSeq):
    revCompTable = {"A": "T", "G": "C", "C": "G", "T": "A", "N": "N"} 
    revSeq = inSeq.upper()[::-1]
    revCompSeq = ""
    for base in revSeq:
        revCompSeq += revCompTable[base]
    return(revCompSeq)

def main():
    parser=ArgumentParser()
    parser.add_argument("--ref", action="store", dest="ref", 
            help="The reference genome in FASTA format", required = True)
    parser.add_argument("--mutpos", action="store", dest="mutpos", 
            help="The mutpos file", required = True)
    parser.add_argument("--label", action="store", dest="label", required = True)
    parser.add_argument("--debug", action="store_true", dest="debug", help="Print out all mutations processed with +strand sequence context.")
    parser.add_argument('-c', '--minClonality', 
                    action = 'store', 
                    type = float,
                    dest = 'minClonal',
                    default = 0,
                    help = 'The minimum clonality at which to examine a mutation. [0]'
                    )
    parser.add_argument('-C', '--maxClonality', 
                    action = 'store', 
                    type = float,
                    dest = 'maxClonal', 
                    default = 1,
                    help = 'The maximum clonality at which to examine a mutation. [1]'
                    )
    o=parser.parse_args()
    
    allFiles = Files(o.ref, o.mutpos, 0, 1, 1)
    allMuts = []
    mutSpect = {
        "C>A":{"ACA":0., "ACC":0., "ACG":0., "ACT":0., "CCA":0., "CCC":0., "CCG":0., "CCT":0., "GCA":0., "GCC":0., "GCG":0., "GCT":0., "TCA":0., "TCC":0., "TCG":0., "TCT":0.}, 
        "C>G":{"ACA":0., "ACC":0., "ACG":0., "ACT":0., "CCA":0., "CCC":0., "CCG":0., "CCT":0., "GCA":0., "GCC":0., "GCG":0., "GCT":0., "TCA":0., "TCC":0., "TCG":0., "TCT":0.}, 
        "C>T":{"ACA":0., "ACC":0., "ACG":0., "ACT":0., "CCA":0., "CCC":0., "CCG":0., "CCT":0., "GCA":0., "GCC":0., "GCG":0., "GCT":0., "TCA":0., "TCC":0., "TCG":0., "TCT":0.},  
        "T>A":{"ATA":0., "ATC":0., "ATG":0., "ATT":0., "CTA":0., "CTC":0., "CTG":0., "CTT":0., "GTA":0., "GTC":0., "GTG":0., "GTT":0., "TTA":0., "TTC":0., "TTG":0., "TTT":0.}, 
        "T>C":{"ATA":0., "ATC":0., "ATG":0., "ATT":0., "CTA":0., "CTC":0., "CTG":0., "CTT":0., "GTA":0., "GTC":0., "GTG":0., "GTT":0., "TTA":0., "TTC":0., "TTG":0., "TTT":0.}, 
        "T>G":{"ATA":0., "ATC":0., "ATG":0., "ATT":0., "CTA":0., "CTC":0., "CTG":0., "CTT":0., "GTA":0., "GTC":0., "GTG":0., "GTT":0., "TTA":0., "TTC":0., "TTG":0., "TTT":0.}
        }
    mutSpectKeys = {"all":["C>A","C>G", "C>T", "T>A", "T>C", "T>G"], "C":["ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT"], "T":["ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC", "TTG", "TTT"]}
    
    typeTrans = {"A>C":"T>G","A>G":"T>C","A>T":"T>A","G>A":"C>T","G>C":"C>G","G>T":"C>A"}
    lineNum=0#DEBUG
    print("Processing mutations...")
    for line in allFiles:
        lineNum += 1
        if lineNum % 1000 == 0:
            print('%s mutation sites processed' % lineNum)
        if o.minClonal <= line.clonality <= o.maxClonal:
            newMuts = line.makeMuts()
            for mut in newMuts:
                mut.setContext(allFiles.Fasta)
                if o.debug:
                    print(mut.chrom, mut.pos, mut.Type, mut.context)
                if mut.Type in typeTrans.keys():
                    mutSpect[typeTrans[mut.Type]][revComp(mut.context)] += 1
                else:
                    try:
                        mutSpect[mut.Type][mut.context] += 1
                    except Exception:
                        print(mut.Type, mut.context)
                        print(mut.chrom, mut.pos)
                        print("fasta dump")
                        print("Position: ", allFiles.Fasta.chrom, allFiles.Fasta.pos)
                        print(allFiles.Fasta.minPos, allFiles.Fasta.maxPos)
                        print("Line Data: ", allFiles.Fasta.lLength)
                        print(allFiles.Fasta.data)
                        raise
    
    allCounts = [0]
    allLabels = ['']
    metaLabels = ['']
    print("Preparing data table...")
    for mutType in mutSpectKeys["all"]:
        for mutKey in mutSpectKeys[mutType.split(">")[0]]:
            allCounts.append(mutSpect[mutType][mutKey])
            allLabels.append(mutKey)
            metaLabels.append(mutType)
        allCounts.append(0)
        allLabels.append('')
        metaLabels.append('')
    totalCounts = sum(allCounts)
    dataFile = open("%s.mcs.dat.txt" % o.label, "w")
    dataFile.write("Type\tContext\tCount\tProportion")
    for ind in xrange(len(allCounts)):
        dataFile.write("\n%s\t%s\t%s\t" % (metaLabels[ind],allLabels[ind],allCounts[ind]))
        allCounts[ind] /= totalCounts
        dataFile.write("%s" % allCounts[ind])
    dataFile.close()

    colTrans = {"C>A":'c', "C>G":'0.2', "C>T":'r', "T>A":'0.75', "T>C":'g', "T>G":'m', "":'w'}
    print("Building Figure...")
    plt.figure(figsize=(11,3))
    ind=np.arange(len(allCounts))
    width = 1
    rects = plt.bar(ind, allCounts, width, color=[colTrans[x] for x in metaLabels])
    plt.ylabel("Proportion of Mutations")
    plt.title("Mutation Spectrum: %s" % o.label)
    plt.xticks(ind+width/2., allLabels, rotation='vertical', fontsize=7)
    plt.xlim([0, ind.size])
    plt.yticks(fontsize=7)
    legendCreator = []
    for mutType in mutSpectKeys["all"]:
        legendCreator.append(plt.Rectangle((0,0), 1, 1, fc=colTrans[mutType]))
    plt.figlegend(legendCreator, mutSpectKeys["all"], loc=5, fontsize=7)
    print("Saving...")
    plt.savefig("%s.mcs.png" % o.label)
    

    allFiles.close()

if __name__ == "__main__":
    main()