#!/usr/bin/python

from sys import argv

# Usage:		 
# python filter_pileup.py genomic_positions_file.bed pileupfile outputfile Y/N
#                                                           
# Filters a pileup file to include only genomic positions as indicated in bed file.
# Can also indicate to only filter out those genomic positions in bed file. By indicating
# Y instead of N at the end of the terminal command. 
# 
# Last Modified by Dana Nachmanson 9/12/16 to allow for filtering OUT regions as well. 

script, region_file, pileup_file, filtered_file, remove_reg = argv

filter = open(region_file,'r+')
pileup = open(pileup_file, 'r')
filtered_file = open(filtered_file, 'w')

if remove_reg.lower() == 'y':
    filter_out = True
else:
    filter_out = False

if not filter_out:
    print("Creating a filtered pileup using the genomic coordinates contained in %r.") %region_file
else:
    print("Creating a filtered pileup using the genomic coordinates NOT contained in %r.") %region_file


# Read in genomic regions to filter
regions = []

for line in filter:
    regions.append(line)

# Read through pileup file and find lines to keep/discard based on region filter
for line in pileup:

    test_reg = line.split()
    in_region = False

    for i in regions:
        if (test_reg[0] == i.split("chr")[1].split()[0]) and (int(test_reg[1]) >= int(i.split()[1])) and (int(test_reg[1]) <= int(i.split()[2])):
            # If you'd like to keep those regions in the bed file.
            in_region = True
            if not filter_out:
                filtered_file.write(line)

    if filter_out and not in_region:
        filtered_file.write(line)

filter.close()
pileup.close()
filtered_file.close()