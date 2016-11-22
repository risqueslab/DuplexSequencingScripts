#!/usr/bin/env python

# Mutation_Spectrum_Plots.py
#
# Last Modified by Dana Nachmanson 9/16/16
#
# This script takes input "Count-Muts" files which are created
# by running Count_Muts.py (v. 1.41) on a pileup file and plots
# the mutational spectrum information found in those files. 
#
# It does expect that your filename starts with a sample name separated by a period. 
# For example: "Home/Documents/DCS_Analysis/Sample1.countmuts.txt" or "Sample2.countmuts.noSNPs.unique.txt"
#
# Usage: python Mutation_Spectrum_Plots.py -dcs 1.dcs.countmuts 2.dcs.countmuts 3.dcs.countmuts -sscs 1.sscs.countmuts 2.sscs.countMuts [Optional: -dcstxt -sscstxt]
#
# Can input any number of DCS or SSCS input countmuts files and this program will plot mutational frequencies of those samples combined. This program can 
# also longitudinally combine all of the countmuts files that are inputted into a larger output file if so desired. 

from argparse import ArgumentParser
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import collections
import seaborn


def count_mutations(countMutsFile, isDCS, outFile):
	"""
	Takes an input count_muts file a boolean indicating whether the count_muts data was made from 
	DCS data and whether to make a concatenated text file of count_muts files and returns a tuple of:
	(sample name, dictionary of mutations and their frequencies).
	"""
	
	# Check if file name has path attached to beginning and remove it if so.
	countMutsFileName = countMutsFile
	if '/' in countMutsFileName:
		countMutsFileName = countMutsFile.split('/')[-1]

	sampleName = str(countMutsFileName.split('.')[0])
	
	inFile = open(countMutsFile, 'r')
	if outFile:
		outFile.write("Sample: " + sampleName)
		
	mutFreqDict = {}
	sumGCSeq = 0
	sumTASeq = 0

	for line in inFile:
		
		if outFile:
			outFile.write(line)
		
		mutsInfo = line.split()
		if len(mutsInfo) > 0:
			refBase = str(mutsInfo[0])

			# Sum the number of G's/C's sequenced and A's/T's sequenced
			if len(mutsInfo) == 3:
				if ('A' in mutsInfo[0]) or ('T' in mutsInfo[0]):
					sumTASeq += int(mutsInfo[2])
				elif ('G' in mutsInfo[0]) or ('C' in mutsInfo[0]):
					sumGCSeq += int(mutsInfo[2])

			# For SSCS data, store the mut frequency, for DCS data, store the mut count.
			elif refBase == 'A' or refBase == 'T' or refBase == 'C' or refBase == 'G':
				mutBase = str(mutsInfo[2].split(':')[0]) 
				if isDCS:
					mutFreqDict[refBase + ">" + mutBase] = float(mutsInfo[3])
				else:
					mutFreqDict[refBase + ">" + mutBase] = float(mutsInfo[4])			
	if outFile:
		outFile.write('\n')
	
	inFile.close()
	
	# If the countmuts file is DCS then calculate new frequencies
	if isDCS:
		mutFreqDict = make_dcs_frequencies(mutFreqDict, sumGCSeq, sumTASeq)

	# Sort dictionary
	oMutFreqDict = collections.OrderedDict(sorted(mutFreqDict.items()))
	
	return (sampleName, oMutFreqDict)
	
	
def make_dcs_frequencies(mutFreqDict, sumGCSeq, sumTASeq):
	"""
	Take an input of dictionary with mutation counts and the sum of G's/C's sequenced
	and A's/T's sequenced. Return dictionary with mutational frequencies calculated
	for dcs data.
	"""
	
	dcsDict = {}
	
	dcsDict['C>A\nG>T'] = float(mutFreqDict['G>T'] + mutFreqDict['C>A'])/sumGCSeq
	dcsDict['C>G\nG>C'] = float(mutFreqDict['G>C'] + mutFreqDict['C>G'])/sumGCSeq
	dcsDict['C>T\nG>A'] = float(mutFreqDict['G>A'] + mutFreqDict['C>T'])/sumGCSeq
	dcsDict['T>C\nA>G'] = float(mutFreqDict['A>G'] + mutFreqDict['T>C'])/sumTASeq
	dcsDict['T>G\nA>C'] = float(mutFreqDict['A>C'] + mutFreqDict['T>G'])/sumTASeq
	dcsDict['T>A\nA>T'] = float(mutFreqDict['A>T'] + mutFreqDict['T>A'])/sumTASeq
	
	return dcsDict
	

def plot_muts(mutDicts, isDCS):
	"""
	This function takes an input dictionary of mutation frequencies and whether
	the frequencies come from a DCS source. It then creates a bar plot with the 
	type of mutation on the x-axis and the frequency of the mutation on the y-axis.
	"""
	
	# Record Mutations
	muts = mutDicts[0][1].keys()
	n = len(mutDicts)

	fig = plt.figure()
	ax = fig.add_subplot(111)

	# Space between clusters of bars
	space = 0.2
	# Width of bars determined by # of samples
	width = (1 - space) / len(mutDicts)

	# Add each sample to the bar graph
	for i, sample in enumerate(mutDicts):
		freqs = sample[1].values()
		# Position of the bars as a function of number of samples, space and width of bars.
		pos = [j - (1 - space) / 2. + i * width for j in range(1,len(muts)+ 1)]
		# Color of bars can be changed here by changing "cm.winter" to another colormap scheme
		ax.bar(pos, freqs, width, label=sample[0], color=cm.winter(float(i) / n))

	# Create a title:
	titlePrefix = "SSCS"
	if isDCS:
		titlePrefix = "DCS"
	fig.suptitle(titlePrefix + ' Mutation Spectrum', fontsize=14, fontweight='bold')

	# Label y axis and x axis
	ax.set_ylabel("Mutation Frequency", fontweight='bold')
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	
	indeces = range(1, len(muts)+1)
	ax.set_xticks(indeces)
	ax.set_xticklabels(muts)
	ax.set_xlim([0,len(muts) + 1])
	ax.set_ylabel("Mutation Frequency", fontweight='bold')
	ax.set_xlabel("Mutation Type", fontweight='bold')
	if not isDCS:
		plt.setp(plt.xticks()[1], rotation=90)

	# Create and format legend
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles[::-1], labels[::-1],loc='best')
	plt.subplots_adjust(bottom=0.15)
	
	plt.savefig(titlePrefix + ".CountMutsPlots.pdf")


def plot_SSCS_ratio(mutDicts):
	"""
	This function takes an input dictionary of SSCS mutation frequencies. It then
	calculates and plots the ratio of G>T/C>A frequencies in all the samples contained
	in the dictionary.
	"""
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	sampleNames = []
	sscsRatios = []
	# Record sample names and calculate G>T/C>A frequency ratio
	for sample in mutDicts:
		sampleNames.append(str(sample[0]))
		sscsRatios.append(float(sample[1]["G>T"])/ float(sample[1]["C>A"]))

	# Space between bars
	space = 0.2
	# Width of bars determined by # of samples
	width = (1 - space) / len(sampleNames)
	yPos = range(1, len(sampleNames)+1)
	plt.bar(yPos, sscsRatios, width=width, align='center', color='b')
	# Create a title:
	plt.suptitle('SSCS G>T/C>A Ratio (Oxidative Damage Ratio)', fontsize=14, fontweight='bold')

	# Label y axis and x axis
	plt.ylabel("G>T/C>A Ratio", fontsize=12, fontweight='bold')
	plt.xlabel("Sample", fontsize=12, fontweight='bold')

	plt.xticks(yPos, sampleNames)
	ax.set_xlim([0,len(sampleNames) + 1 ])
	
	plt.savefig("SSCSRatio.Plot.pdf")

def main():

		parser = ArgumentParser()
		parser.add_argument('-dcs', '--in_dcs_list', nargs='+', action='store', dest = 'inDCS', help = 'A list of input DCS files.', default = None)
		parser.add_argument('-sscs', '--in_sscs_list', nargs='+', action ='store', dest = 'inSSCS', help = 'A list of input SSCS files.', default = None)
		parser.add_argument('-dcstxt', '--dcs_txt', action='store_true', dest='makeDCSTxt', help='Output a text file with concatenated dcs countmuts files')
		parser.add_argument('-sscstxt', '--sscs_txt', action='store_true', dest='makeSSCSTxt', help='Outputs  text file with concatenated scs countmuts files')
		o = parser.parse_args()

		# ------------ Read files and record mutations -------------

		# Take all DCS Files and create a list of dictionaries with the files, if files exist and there is at least 1.
		if o.inDCS != None and len(o.inDCS) > 0 :
			dcsDicts = []
			dcsOutFile = None
			if o.makeDCSTxt:
				dcsOutFile = open("DCS_CountMuts.txt", 'w')
			for dcsFile in o.inDCS:
				dcsDicts.append(count_mutations(dcsFile, True, dcsOutFile))
			if o.makeDCSTxt:
				dcsOutFile.close()
				

		# Take all SSCS Files and create a list of dictionaries with the files, if files exist and there is at least 1.
		if o.inSSCS != None and len(o.inSSCS) > 0 :
			sscsDicts = []
			sscsOutFile = None
			if o.makeSSCSTxt:
				sscsOutFile = open("SSCS_CountMuts.txt", 'w')
			for sscsFile in o.inSSCS:
				sscsDicts.append(count_mutations(sscsFile, False, sscsOutFile))
			if o.makeSSCSTxt:
				sscsOutFile.close()

		# ---------------------- Make plots ------------------------
		
		# Plot DCS mutation frequencies
		if dcsDicts:
			plot_muts(dcsDicts, True)

		# Plot SSCS mutation frequencies 
		if sscsDicts:
			plot_muts(sscsDicts, False)
			plot_SSCS_ratio(sscsDicts)


if __name__ == "__main__":
		main()