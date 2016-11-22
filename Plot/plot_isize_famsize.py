#!/usr/bin/python

from argparse import ArgumentParser
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import collections
import seaborn

# Create a scatter plot of Family Size versus Insert Size
#
# Uses sscs.read1.fastq file and mem.sscs.sort.sam file
# Usage: >> python plotISizeByFamilySize.py -fq sample1.sscs.read1.fastq -sam sample1_mem.sscs.sam
# Last modified by Dana Nachmanson 10/10/16

def record_family_sizes(fqFile):
	"""
	This function takes in a fastQ file of read 1 SSCS reads and stores the 
	family size that comprises its tag. Returns a dictionary with a tag as a 
	key and family size as value.
	"""
	
	fq = open(fqFile, "r")
	famSizeDict = {}
	count = 1 
	tag = ""
	
	for line in fq:
		if count == 1:
			tag = str(line.split("@")[1].split("#")[0])
		elif count == 3:
			fam = int(line.split("+")[1])
			famSizeDict[tag] = [fam]
		elif count == 4:
			count = 0
		count += 1
	
	fq.close()
	
	return famSizeDict
	

def add_insert_size(famSizeDict, samFile):
	"""
	This function takes in a SAM file of aligned DCS reads (from data in fastQ)
	It also takes in a dictionary of tags and their average family sizes.
	It returns the same dictionary with an iSize field added to each tag found in SAM. 
	"""
	
	sam = open(samFile, "r")
	tag = ""
	
	for line in sam:
		lineList = line.split()
		# Check if the line is a header line
		if len(lineList) > 13:
			iSize = abs(int(lineList[8]))
			# Doesn't matter if unmapped (AKA iSize == 0)
			tag = str(lineList[0].split("#")[0])
			famSizeDict[tag].append(iSize)
			
	sam.close()
	return famSizeDict

def plot_iSize_fSize(iSizeFamSizeDict):
	
	xFam = []
	yFam = []
	for tag in iSizeFamSizeDict:
		if len(iSizeFamSizeDict[tag]) > 1:
			famSize = iSizeFamSizeDict[tag][0]
			iSize = iSizeFamSizeDict[tag][1]
			if iSize != 0 and iSize < 3000:
				xFam.append(famSize)
				yFam.append(iSize)
	
	plt.scatter(xFam, yFam, alpha=0.5)
	plt.ylabel("iSize", fontsize=12, fontweight='bold')
	plt.xlabel("Family Size", fontsize=12, fontweight='bold')
	plt.show() 
		

def main():

		parser = ArgumentParser()
		parser.add_argument('-fq', '--in_dcs_fq', action='store', dest = 'fq', help = 'A fastq file of SSCS reads.', default = None)
		parser.add_argument('-sam', '--in_dcs_sam', action ='store', dest = 'sam', help = 'A sam file of aligned SSCS reads.', default = None)
		o = parser.parse_args()
		
		# ------------ Read fq file and record insert sizes and tags-------------

		# Take all DCS Files and create a list of dictionaries with the files, if files exist and there is at least 1.
		if o.fq != None and o.sam != None :
			famSizeDict = record_family_sizes(o.fq)
			iSizeFamSizeDict = add_insert_size(famSizeDict, o.sam)
			plot_iSize_fSize(iSizeFamSizeDict)
			

if __name__ == "__main__":
		main()