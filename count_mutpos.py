#!/usr/bin/python

from argparse import ArgumentParser
import numpy


def read_file(inputFile):
	
	file = open(inputFile, "r")
	fileList = []
	lineCount = 0
	for line in file:
		if lineCount > 0:
			fileList.append(int(line.split()[3]))
		else:
			lineCount += 1
	file.close()
	return fileList
	
	
def calc_avg_depth(depthList):
	return sum(depthList)/ float(len(depthList))
	
def calc_med_depth(depthList):
    return numpy.median(numpy.array(depthList))
    
def sum_nts_sequenced(depthList):	
	return sum(depthList)


def main():

		parser = ArgumentParser()
		parser.add_argument('-mu', '--in_mutpos', action='store', dest = 'mu', help = 'A full mutpos file.', default = None)
		parser.add_argument('-exmu', '--in_exon_mutpos', action ='store', dest = 'exons', help = 'A mutpos file of just exons.', default = None)
		parser.add_argument('-inmu', '--in_intron_mutpos', action ='store', dest = 'introns', help = 'A mutpos file of just introns.', default = None)
		
		o = parser.parse_args()
		
		# ------------ Read mutpos files and calculate depths and # of nucleotides-------------

		if o.mu != None and o.exons != None  and o.introns != None:
			
			mutposList = read_file(o.mu)
			intronList = read_file(o.introns)
			exonList = read_file(o.exons) 
			
			avgDepth = calc_avg_depth(mutposList)
			avgExonDepth = calc_avg_depth(exonList)
			avgIntronDepth = calc_avg_depth(intronList)
			
			medianDepth = calc_med_depth(mutposList)
			medianExonDepth = calc_med_depth(exonList)
			medianIntronDepth = calc_med_depth(intronList)
			
			totalNts = sum_nts_sequenced(mutposList)
			totalExonNts = sum_nts_sequenced(exonList)
			totalIntronNts = sum_nts_sequenced(intronList)
			
			print("The average depth is: %d") % avgDepth
			print("The average exonic depth is: %d") % avgExonDepth
			print("The average intronic depth is: %d") % avgIntronDepth
			print("")
			print("The median depth is: %d") % medianDepth
			print("The median exonic depth is: %d") % medianExonDepth
			print("The median intronic depth is: %d") % medianIntronDepth
			print("")
			print("The total # of DCS nucleotides sequenced is: %d") % totalNts
			print("The total # of EXONIC DCS nucleotides sequenced is: %d") % totalExonNts
			print("The total # of INTRONIC DCS nucleotides sequenced is: %d") % totalIntronNts
			
if __name__ == "__main__":
		main()