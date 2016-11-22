# Written By: Alexander M. West
# Version: Feb. 11, 2016
# Usage: python <this>.py <transtable> <bed file> <reference genome> <mutpos file> <output name prefix> <minimum_clonality> <maximum_clonality>

import sys

transtable = {}
gene_list = []
gene_data = {}
chromosomes_of_interest = set()
reference_genome = {}
mutation_list = []
synonymous_tally = [0, 0]  #synonymous in 0, non-synonymous in 1
aaAbrev = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H',
           'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
           'Tyr': 'Y', 'Val': 'V', 'Asx': 'B', 'Glx': 'Z', 'Ter': 'Ter'}

output_file = open(sys.argv[5] + ".consequences.txt", "w")
nonsynonymous_out_file = open(sys.argv[5] + ".nonsynonymous.txt", "w")
SvN_output_file = open(sys.argv[5] + ".SvN.txt", "w")

# input: a Kennedy lab transtable, example line: "Met ATG,ATA, Start"
# iterates through lines of transtable and places them in a dictionary
def import_transtable(file):
    transtable_file = open(file)
    while True:
        line = transtable_file.readline()
        if not line: break
        if ">" in line:
            continue
        else:
            AA, Codons, Type = line.rstrip().split()
            for i in (Codons.split(",")):
                transtable[i] = AA
    transtable_file.close()


# input: a Kennedy lab bed file, example line: "MT2    3306    4261    ND1    0    +"
# the bed file should take into account the lack of legacy "N"
# most current Kennedy files are off by -1, probably because they were trying to account for zero index of reference.
# iterates through lines of bed file and places them in a dictionary
# makes a dictionary of relevant chromosomes to filter what is read from the genome reference
def import_gene_bed(file):
    bed_file = open(file)
    while True:
        line = bed_file.readline()
        if not line:
            break
        else:
            chromosome, start, end, gene, number, direction = line.rstrip().split()
            chromosomes_of_interest.add(chromosome)
            gene_list.append(gene)
            gene_data[gene] = [chromosome, start, end, gene, direction]
    bed_file.close()


# input: a fasta format reference file - should not contain legacy "N"
# iterates through the reference and builds a char list of selected chromosomes in the genome
# zero indexed
def import_reference_genome(file):
    reference_file = open(file)
    while True:
        line = reference_file.readline()
        if not line: break
        if ">" in line:
            chromosome = line.rstrip().split()[0].replace(">", "")
            print chromosome
        else:
            if chromosome in chromosomes_of_interest:
                sequence_fragment = line.rstrip()
                for i in range(0, len(sequence_fragment)):
                    reference_genome[chromosome].append(sequence_fragment[i])
    reference_file.close()


# input: a mutpos file
# iterates through the file to build a list of mutations in (chromosome, position, nucleotide) format.
def import_mutation_list(file, min_threshold, max_threshold):
    mutpos_file = open(file)
    while True:
        line = mutpos_file.readline()
        if not line: break
        else:
            chromosome, reference, position, reads, mutations, T_count, C_count, G_count, A_count, Ins, Del, N = line.rstrip().split()
            if float(max_threshold) >= float(T_count) / float(reads) > float(min_threshold):
                mutation_list.append((chromosome, position, "T"))
            if float(max_threshold) >= float(C_count) / float(reads) > float(min_threshold):
                mutation_list.append((chromosome, position, "C"))
            if float(max_threshold) >= float(G_count) / float(reads) > float(min_threshold):
                mutation_list.append((chromosome, position, "G"))
            if float(max_threshold) >= float(A_count) / float(reads) > float(min_threshold):
                mutation_list.append((chromosome, position, "A"))
    mutpos_file.close()


# input: the mutant position and nucleotide
# output: prints to file a string with either NC, or the gene/reference/mutant.
def report_consequence(chromosome, position, mutant):
    found_in_gene = False
    for i in gene_list:
        if chromosome == gene_data[i][0]:
            start = int(gene_data[i][1])
            end = int(gene_data[i][2])
            direction = gene_data[i][4]
            if int(position) >= start and int(position) <= end:
                found_in_gene = True
                if direction == "+":
                    distance = int(position) - start
                    reading_frame = distance % 3
                else:
                    distance = end - int(position)
                    reading_frame = abs(distance % 3 - 2)
                # all following reference_genome[chromosome] lookups are - 1 because reference is 0 indexed.
                if reading_frame == 0:
                    reference_seq = reference_genome[chromosome][position - 1] + reference_genome[chromosome][position - 1 + 1] + reference_genome[chromosome][position - 1 + 2]
                    mutant_seq = mutant + reference_genome[chromosome][position - 1 + 1] + reference_genome[chromosome][position - 1 + 2]
                elif reading_frame == 1:
                    reference_seq = reference_genome[chromosome][position - 1 - 1] + reference_genome[chromosome][position - 1] + reference_genome[chromosome][position - 1 + 1]
                    mutant_seq = reference_genome[chromosome][position - 1 - 1] + mutant + reference_genome[chromosome][position - 1 + 1]
                else:
                    reference_seq = reference_genome[chromosome][position - 1 - 2] + reference_genome[chromosome][position - 1 - 1] + reference_genome[chromosome][position - 1]
                    mutant_seq = reference_genome[chromosome][position - 1 - 2] + reference_genome[chromosome][position - 1 - 1] + mutant
                if "N" in reference_seq or "N" in mutant_seq:
                    output_file.write(chromosome + "\t" + direction + "\t" + reference_genome[chromosome][position - 1] + str(position) + mutant + "\t" + gene_data[i][3] + "\t" + reference_seq + "\t" + mutant_seq + "\t" + ".\n")
                else:
                    if direction == "+":
                        reference_AA = transtable[reference_seq]
                        mutant_AA = transtable[mutant_seq]
                    else:
                        #print reference_seq
                        reference_AA = transtable[reverse_sequence(reference_seq)]
                        mutant_AA = transtable[reverse_sequence(mutant_seq)]
                    output_file.write(chromosome + "\t" + direction + "\t" + reference_genome[chromosome][position - 1] + str(position) + mutant + "\t" + gene_data[i][3] + "\t" + reference_AA + "\t" + mutant_AA + "\t"),
                    if reference_AA == mutant_AA:
                        output_file.write(".\n")
                        synonymous_tally[0] += 1
                    else:
                        #positions 0-2 = amino acid 1, 3-6 = amino acid 2, and so forth.  int functionally rounds down.
                        output_file.write(aaAbrev[reference_AA] + str(int(float(distance)/3.0) + 1) + aaAbrev[mutant_AA] + "\n")
                        nonsynonymous_out_file.write(gene_data[i][3] + "\t" + aaAbrev[reference_AA] + str(int(float(distance)/3.0) + 1) + aaAbrev[mutant_AA] + "\n")
                        synonymous_tally[1] += 1
    if found_in_gene == False:
        output_file.write(chromosome + "\t*\t" + reference_genome[chromosome][position - 1] + str(position) + mutant + "\t.\tNC\tNC\t.\n")

#input: name of a gene imported to the gene_list via import_gene_bed
#output: a SvN summary of mutations in each amino acid of the input gene
#checks every amino acid position for a gene, and evaluates SvN ratio of all possible variants.
def SvN_tally(gene):
    chromosome, start, end, gene, direction = gene_data[gene]
    synonymous_count = 0
    nonsynonymous_count = 0
    distance = 0
    while int(start) + distance < int(end):
        if direction == "+":
            reading_frame = distance % 3
        else:
            reading_frame = abs(distance % 3 - 2)
        position = int(start) + distance
        # all following reference_genome[chromosome] lookups are - 1 because reference is 0 indexed.
        if reading_frame == 0:
            reference_seq = reference_genome[chromosome][position - 1] + reference_genome[chromosome][position - 1 + 1] + reference_genome[chromosome][position - 1 + 2]
        elif reading_frame == 1:
            reference_seq = reference_genome[chromosome][position - 1 - 1] + reference_genome[chromosome][position - 1] + reference_genome[chromosome][position - 1 + 1]
        else:
            reference_seq = reference_genome[chromosome][position - 1 - 2] + reference_genome[chromosome][position - 1 - 1] + reference_genome[chromosome][position - 1]
        if 'N' in reference_seq: continue
        if direction == "+":
            reference_AA = transtable[reference_seq]
        else:
            reference_AA = transtable[reverse_sequence(reference_seq)]
        #print "Reference: " + reference_AA + "\n"
        nucleotides = ['T', 'C', 'G', 'A']
        nucleotides.remove(reference_genome[chromosome][position - 1])
        bases = {}
        for i in (0, 1, 2):
            bases[i] = reference_seq[i]
        for i in nucleotides:
            bases[reading_frame] = i
            mutant_seq = bases[0] + bases[1] + bases[2]
            if direction == "+":
                mutant_AA = transtable[mutant_seq]
            else:
                mutant_AA = transtable[reverse_sequence(mutant_seq)]
            if mutant_AA == reference_AA:
                synonymous_count += 1
                #print "Variant: " + mutant_AA + "\t" + str(synonymous_count) + "\t" + str(nonsynonymous_count) + "\n"
            else:
                nonsynonymous_count += 1
                #print "Non-Variant: " + mutant_AA + "\t" + str(synonymous_count) + "\t" + str(nonsynonymous_count) + "\n"
        distance += 1
    SvN_output_file.write(gene + "\t" + str(synonymous_count) + "\t" + str(nonsynonymous_count) + "\t" + str(float(synonymous_count) / nonsynonymous_count) + "\n")

# input: a 3 nucleotide codon
# return: the inverted 3 nucleotide codon
def reverse_sequence(sequence):
    return_sequence = ""
    for i in (2, 1, 0):
        if sequence[i] == 'G': return_sequence += 'C'
        if sequence[i] == 'C': return_sequence += 'G'
        if sequence[i] == 'T': return_sequence += 'A'
        if sequence[i] == 'A': return_sequence += 'T'
    #print return_sequence
    return return_sequence


import_transtable(sys.argv[1])
import_gene_bed(sys.argv[2])
for i in chromosomes_of_interest:
    reference_genome[i] = []
import_reference_genome(sys.argv[3])
import_mutation_list(sys.argv[4], sys.argv[6], sys.argv[7])
#print transtable
#print gene_list
#print gene_data
#print chromosomes_of_interest
#for i in chromosomes_of_interest:
#print reference_genome
#print mutation_list
for i in mutation_list:
    report_consequence(i[0], int(i[1]), i[2])
for i in gene_list:
    SvN_tally(i)
nonsynonymous_out_file.write("\nSynonymous mutations: " + str(synonymous_tally[0]) + "\nNon-synonymous mutations: " + str(synonymous_tally[1]) + "\n")

output_file.close()
nonsynonymous_out_file.close()
SvN_output_file.close()

print "Synonymous mutations: " + str(synonymous_tally[0])
print "Non-synonymous mutations: " + str(synonymous_tally[1])

# TestCases for 41zSample
# report_consequence(73, 'G')
# report_consequence(3446, 'C')
# report_consequence(8552, 'A')
# report_consequence(14363, 'C')
