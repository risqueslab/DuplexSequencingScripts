#USAGE: python mutpos_analyzer.py <mutpos file> <consequences file> <output name> <Reference file>
#OUTPUT: analysis, consequences, haplogrep, and mitomaster files suitable for excel or website input
#Written By: Alexander M. West
#Version: November 4, 2015

import sys

#open i/o files
reference_file = open(sys.argv[4])
mutpos_file = open(sys.argv[1])
consequences_file = open(sys.argv[2])
output_file = open(sys.argv[3] + ".analysis.txt", "w")
haplogrep_out = open(sys.argv[3] + ".haplogrep.txt", "w")
mitomaster_out = open(sys.argv[3] + ".mitomaster.txt", "w")

#put headers in output file
output_file.write ("position\tbase\treference\tdloop\ttimes sequenced\ttotal mutations\tT\tC\tG\tA\tINS\tDEL\tN\tmutations frequency\tmutation type\ttriplet frequncy\ttriplet\thaplogrep\tsynonymous\tdirection\trefNt\tCompact NT\tgene name\tAA#\trefAA\tTmutAA\tCmutAA\tGmutAA\tAmutAA\taaCompact\n")

d_loop_end = 577
d_loop_start = 16023
tally_types = ['all', 'clonal', 'subclonal', 'rare', 'rare_d_loop', 'rare_not_loop']
base_tally = {'all': {'T':0, 'C':0, 'G':0, 'A':0, 'Total':0}, 'clonal': {'T':0, 'C':0, 'G':0, 'A':0, 'Total':0}, 'subclonal': {'T':0, 'C':0, 'G':0, 'A':0, 'Total':0}, 'rare': {'T':0, 'C':0, 'G':0, 'A':0, 'Total':0}, 'rare_d_loop': {'T':0, 'C':0, 'G':0, 'A':0, 'Total':0}, 'rare_not_loop': {'T':0, 'C':0, 'G':0, 'A':0, 'Total':0}}
mutation_tally = {'all': {'TA':0, 'TC':0, 'TG':0, 'AT':0, 'AC':0, 'AG':0, 'CT':0, 'CA':0, 'CG':0, 'GA':0, 'GC':0, 'GT':0, 'Total':0}, 'clonal': {'TA':0, 'TC':0, 'TG':0, 'AT':0, 'AC':0, 'AG':0, 'CT':0, 'CA':0, 'CG':0, 'GA':0, 'GC':0, 'GT':0, 'Total':0}, 'subclonal': {'TA':0, 'TC':0, 'TG':0, 'AT':0, 'AC':0, 'AG':0, 'CT':0, 'CA':0, 'CG':0, 'GA':0, 'GC':0, 'GT':0, 'Total':0}, 'rare': {'TA':0, 'TC':0, 'TG':0, 'AT':0, 'AC':0, 'AG':0, 'CT':0, 'CA':0, 'CG':0, 'GA':0, 'GC':0, 'GT':0, 'Total':0}, 'rare_d_loop': {'TA':0, 'TC':0, 'TG':0, 'AT':0, 'AC':0, 'AG':0, 'CT':0, 'CA':0, 'CG':0, 'GA':0, 'GC':0, 'GT':0, 'Total':0}, 'rare_not_loop': {'TA':0, 'TC':0, 'TG':0, 'AT':0, 'AC':0, 'AG':0, 'CT':0, 'CA':0, 'CG':0, 'GA':0, 'GC':0, 'GT':0, 'Total':0}}
consequences_reference = {}
#gene list and synonymous tally third values should be loaded from mito.bed for a generalized program.
#currently values are hard coded from the reference files Human_mito.bed.txt
gene_list = ['ND1', 'ND2', 'COX1', 'COX2', 'ATP8', 'ATP6', 'COX3', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'CYTB']
synonymous_tally = {'ND1': [0, 0, 955], 'ND2': [0, 0, 1041], 'COX1': [0, 0, 1541], 'COX2': [0, 0, 683], 'ATP8': [0, 0, 206], 'ATP6': [0, 0, 680], 'COX3': [0, 0, 783], 'ND3': [0, 0, 345], 'ND4L': [0, 0, 296], 'ND4': [0, 0, 1377], 'ND5': [0, 0, 1811], 'ND6': [0, 0, 524], 'CYTB': [0, 0, 1138]}
reference_line = 1

#input: a string with number and non-number characters
#returns: a string with only the numbers
def extract_number(string):
    output_string = ""
    for i in range(len(string)):
        if string[i].isdigit():
            output_string += string[i]
    return output_string

#import MutationConsequences.py's output file
while True:
    line = consequences_file.readline()
    if not line: break
    con_chromosome, con_direction, con_summary, con_gene, con_ref, con_mutant, con_aaCompact = line.rstrip().split()
    con_position = extract_number(con_summary)
    consequences_reference[con_position] = [con_chromosome, con_direction, con_summary, con_gene, con_ref, con_mutant, con_aaCompact]

#output genome reference information
for line in mutpos_file:
    chromosome, base, position, sequenced, total_mutations, T, C, G, A, insertion, deletion, N = line.rstrip().split()
    #print reference lines until arriving at the start position of mutpos file
    while reference_line < int(position):
        line_fields = reference_file.readline().rstrip().split()
        if int(line_fields[0]) < 3107:
            output_file.write(line_fields[0] + "\t"),
        else:
            output_file.write(str(int(line_fields[0]) + 1) + "\t"),
        output_file.write(line_fields[1] + "\t"),
        output_file.write(" ".join(line_fields[2:]) + "\t"),
        reference_line += 1
        if int(position) < d_loop_end or int(position) > d_loop_start:
            output_file.write("D-LOOP\n"),
        else:
            output_file.write("-\n"),
    line_fields = reference_file.readline().rstrip().split()
    if int(line_fields[0]) < 3107:
        output_file.write(line_fields[0] + "\t"),
    else:
        output_file.write(str(int(line_fields[0]) + 1) + "\t"),
    output_file.write(line_fields[1] + "\t"),
    output_file.write(" ".join(line_fields[2:]) + "\t"),
    reference_line += 1
    if int(position) < d_loop_end or int(position) >= d_loop_start:
        output_file.write("D-LOOP\t"),
    else:
        output_file.write("-\t"),

    #ignore repetitive positions
    if int(position) < 301 or (int(position) > 317 and int(position) < 3564) or int(position) > 3570:
        #logic using no "N", output using "N" placeholder at 3107
        if int(position) >= 3107:
            out_position = str(int(position) + 1)
        else:
            out_position = position
        #identify most prevalent mutation type and quantity of that mutation
        #identify triplet quantity - breaks with quadruplets
        mutation_quantity = 0
        mutation_nucleotide = '-'
        mutation_type = '-'
        triplet_quantity = 0
        triplet = base
        if int(total_mutations) > 0:
            if T > mutation_quantity:
                mutation_quantity = T
                mutation_nucleotide = 'T'
            if int(T) > 0:
                triplet += 'T'
            if C > mutation_quantity:
                triplet_quantity = mutation_quantity
                mutation_quantity = C
                mutation_nucleotide = 'C'
            if int(C) > 0:
                triplet += 'C'
            if G > mutation_quantity:
                triplet_quantity = mutation_quantity
                mutation_quantity = G
                mutation_nucleotide = 'G'
            if int(G) > 0:
                triplet += 'G'
            if A > mutation_quantity:
                triplet_quantity = mutation_quantity
                mutation_quantity = A
                mutation_nucleotide = 'A'
            if int(A) > 0:
                triplet += 'A'
            mutation_type = base + mutation_nucleotide

            #assign mutation category
            mutation_category = ''
            mutation_frequency = float(mutation_quantity) / float(sequenced)
            if mutation_frequency < .01:
                if int(position) < d_loop_end or int(position) > d_loop_start:
                    mutation_category = 'rare_d_loop'
                else:
                    mutation_category = 'rare_not_loop'
            elif mutation_frequency < .99 : mutation_category = 'subclonal'
            else: mutation_category = 'clonal'
            
            #tally bases and mutations
            for i in ('all', mutation_category):
                base_tally[i][base] += int(sequenced) - int(total_mutations)
                base_tally[i]['T'] += int(T)
                base_tally[i]['C'] += int(C)
                base_tally[i]['G'] += int(G)
                base_tally[i]['A'] += int(A)
                base_tally[i]['Total'] += int(sequenced)
                mutation_tally[i][mutation_type] += 1
                mutation_tally[i]['Total'] += 1
        else:
            base_tally['all'][base] += int(sequenced)
            base_tally['all']['Total'] += int(sequenced)
        if len(triplet) < 3:
            triplet = '-'

        #write to output file row information
        #output_file.write("AA\t" + str(int(position) % 3) + "\t"),
        output_file.write(sequenced + "\t" + total_mutations + "\t" + T + "\t" + C + "\t" + G + "\t" + A + "\t" + insertion + "\t" + deletion + "\t" + N + "\t"),
        s = float(sequenced)
        output_file.write(str(float(mutation_quantity) / s) + "\t" + mutation_type + "\t" + str(float(triplet_quantity) / s)+ "\t" + triplet + "\t"),
        #outputs to haplogrep and mitomaster need to conform to old standard that keeps legacy 3107
        if int(total_mutations) > 0:
            output_file.write(out_position + mutation_nucleotide + "\t")
            mitomaster_out.write (sys.argv[3] + "\t" + out_position + "\t" + base + "\t" + mutation_nucleotide + "\t" + base + out_position + mutation_nucleotide + "\n")
            if mutation_category == 'clonal' or mutation_category == 'subclonal':
                haplogrep_out.write(out_position + mutation_nucleotide + "\t")
        else:
            output_file.write(out_position + "\t" + "\t")
        if position in consequences_reference:
            if consequences_reference[position][4] != 'NC':
                if consequences_reference[position][5] == 'Ref' or consequences_reference[position][5] == 'NP' or consequences_reference[position][5] == consequences_reference[position][4]:
                    if mutation_category == 'rare_not_loop':
                        synonymous_tally[consequences_reference[position][3]][0] += 1
                    output_file.write("synonymous\t")
                else:
                    if mutation_category == 'rare_not_loop':
                        synonymous_tally[consequences_reference[position][3]][1] += 1
                    output_file.write(consequences_reference[position][4] + " -> "+ consequences_reference[position][5] + "\t")
            else:
                output_file.write("not coding\t")
            output_file.write(consequences_reference[position][6] + "\n")
        else:
            output_file.write("\n")
    else:
        output_file.write("\n")

#write to output summary statsitics
for i in base_tally['rare']:
    base_tally['rare'][i] = base_tally['rare_d_loop'][i] + base_tally['rare_not_loop'][i]
for i in mutation_tally['rare']:
    mutation_tally['rare'][i] = mutation_tally['rare_d_loop'][i] + mutation_tally['rare_not_loop'][i]

for i in tally_types:
    output_file.write("\n\n\n*** Summary of " + i + " mutations: ***\n")
    output_file.write("Nucleotide:\tQuantity:\n")
    output_file.write("T:\t" + str(base_tally[i]['T']) + "\n")
    output_file.write("C:\t" + str(base_tally[i]['C']) + "\n")
    output_file.write("G:\t" + str(base_tally[i]['G']) + "\n")
    output_file.write("A:\t" + str(base_tally[i]['A']) + "\n")
    output_file.write("Total:\t" + str(base_tally[i]['A'] + base_tally[i]['G'] + base_tally[i]['C'] + base_tally[i]['T']) + "\n")

#refactor this section
    output_file.write("\nMutation Names:\tQuantity:\tFrequency:\tProportion:\n")
    try:
        output_file.write("TC:\t" + str(mutation_tally[i]['TC']) + "\t" + str(mutation_tally[i]['TC'] / float(base_tally['all']['T']))+ "\t" + str(mutation_tally[i]['TC'] / float(mutation_tally[i]['Total'])) + "\n")
    except:
        output_file.write("TC:\t" + "0\t0\n")
    try:
        output_file.write("AG:\t" + str(mutation_tally[i]['AG']) + "\t" + str(mutation_tally[i]['AG'] / float(base_tally['all']['A']))+ "\t" + str(mutation_tally[i]['AG'] / float(mutation_tally[i]['Total'])) + "\n")
    except:
        output_file.write("AG:\t" + "0\t0\n")
    try:
        output_file.write("GA:\t" + str(mutation_tally[i]['GA']) + "\t" + str(mutation_tally[i]['GA'] / float(base_tally['all']['G']))+ "\t" + str(mutation_tally[i]['GA'] / float(mutation_tally[i]['Total'])) + "\n")
    except:
        output_file.write("GA:\t" + "0\t0\n")
    try:
        output_file.write("CT:\t" + str(mutation_tally[i]['CT']) + "\t" + str(mutation_tally[i]['CT'] / float(base_tally['all']['C']))+ "\t" + str(mutation_tally[i]['CT'] / float(mutation_tally[i]['Total'])) + "\n")
    except:
        output_file.write("CT:\t" + "0\t0\n")
    try:
        output_file.write("GC:\t" + str(mutation_tally[i]['GC']) + "\t" + str(mutation_tally[i]['GC'] / float(base_tally['all']['G']))+ "\t" + str(mutation_tally[i]['GC'] / float(mutation_tally[i]['Total'])) + "\n")
    except:
        output_file.write("GC:\t" + "0\t0\n")
    try:
        output_file.write("CG:\t" + str(mutation_tally[i]['CG']) + "\t" + str(mutation_tally[i]['CG'] / float(base_tally['all']['C']))+ "\t" + str(mutation_tally[i]['CG'] / float(mutation_tally[i]['Total'])) + "\n")
    except:
        output_file.write("CG:\t" + "0\t0\n")
    try:
        output_file.write("AC:\t" + str(mutation_tally[i]['AC']) + "\t" + str(mutation_tally[i]['AC'] / float(base_tally['all']['A']))+ "\t" + str(mutation_tally[i]['AC'] / float(mutation_tally[i]['Total'])) + "\n")
    except:
        output_file.write("AC:\t" + "0\t0\n")
    try:
        output_file.write("TG:\t" + str(mutation_tally[i]['TG']) + "\t" + str(mutation_tally[i]['TG'] / float(base_tally['all']['T']))+ "\t" + str(mutation_tally[i]['TG'] / float(mutation_tally[i]['Total'])) + "\n")
    except:
        output_file.write("TG:\t" + "0\t0\n")
    try:
        output_file.write("AT:\t" + str(mutation_tally[i]['AT']) + "\t" + str(mutation_tally[i]['AT'] / float(base_tally['all']['A']))+ "\t" + str(mutation_tally[i]['AT'] / float(mutation_tally[i]['Total'])) + "\n")
    except:
        output_file.write("AT:\t" + "0\t0\n")
    try:
        output_file.write("TA:\t" + str(mutation_tally[i]['TA']) + "\t" + str(mutation_tally[i]['TA'] / float(base_tally['all']['T']))+ "\t" + str(mutation_tally[i]['TA'] / float(mutation_tally[i]['Total'])) + "\n")
    except:
        output_file.write("TA:\t" + "0\t0\n")
    try:
        output_file.write("GT:\t" + str(mutation_tally[i]['GT']) + "\t" + str(mutation_tally[i]['GT'] / float(base_tally['all']['G']))+ "\t" + str(mutation_tally[i]['GT'] / float(mutation_tally[i]['Total'])) + "\n")
    except:
        output_file.write("GT:\t" + "0\t0\n")
    try:
        output_file.write("CA:\t" + str(mutation_tally[i]['CA']) + "\t" + str(mutation_tally[i]['CA'] / float(base_tally['all']['C']))+ "\t" + str(mutation_tally[i]['CA'] / float(mutation_tally[i]['Total'])) + "\n")
    except:
        output_file.write("CA:\t" + "0\t0\n")
    output_file.write("Total:\t" + str(mutation_tally[i]['Total']) + "\n")

    output_file.write("\nMutation Combinations:\tQuantity:\tFrequency:\n")
    try:
        output_file.write("TC + AG:\t" + str(mutation_tally[i]['TC'] + mutation_tally[i]['AG']) + "\t" + str( (mutation_tally[i]['TC'] + mutation_tally[i]['AG'])/float(base_tally['all']['T'] + base_tally['all']['A']) )+ "\n")
    except ZeroDivisionError:
        output_file.write("TC + AG:\t" + "0\t0\n")
    try:
        output_file.write("GA + CT:\t" + str(mutation_tally[i]['GA'] + mutation_tally[i]['CT']) + "\t" + str( (mutation_tally[i]['GA'] + mutation_tally[i]['CT'])/float(base_tally['all']['G'] + base_tally['all']['C']) )+ "\n")
    except ZeroDivisionError:
        output_file.write("GA + CT:\t" + "0\t0\n")
    try:
        output_file.write("GC + CG:\t" + str(mutation_tally[i]['GC'] + mutation_tally[i]['CG']) + "\t" + str( (mutation_tally[i]['GC'] + mutation_tally[i]['CG'])/float(base_tally['all']['G'] + base_tally['all']['C']) )+ "\n")
    except ZeroDivisionError:
        output_file.write("GC + CG:\t" + "0\t0\n")
    try:
        output_file.write("AC + TG:\t" + str(mutation_tally[i]['AC'] + mutation_tally[i]['TG']) + "\t" + str( (mutation_tally[i]['AC'] + mutation_tally[i]['TG'])/float(base_tally['all']['A'] + base_tally['all']['T']) )+ "\n")
    except ZeroDivisionError:
        output_file.write("AC + TG:\t" + "0\t0\n")
    try:
        output_file.write("AT + TA:\t" + str(mutation_tally[i]['AT'] + mutation_tally[i]['TA']) + "\t" + str( (mutation_tally[i]['AT'] + mutation_tally[i]['TA'])/float(base_tally['all']['A'] + base_tally['all']['T']) )+ "\n")
    except ZeroDivisionError:
        output_file.write("AT + TA:\t" + "0\t0\n")
    try:
        output_file.write("GT + CA:\t" + str(mutation_tally[i]['GT'] + mutation_tally[i]['CA']) + "\t" + str( (mutation_tally[i]['GT'] + mutation_tally[i]['CA'])/float(base_tally['all']['G'] + base_tally['all']['C']) )+ "\n")
    except ZeroDivisionError:
        output_file.write("GT + CA:\t" + "0\t0\n")

output_file.write("\n\n\n*** Summary of rare non d-loop mutation consequences by gene: ***\n")
output_file.write("Gene:\tSynonymous:\tNon-Synonymous:\tSynonymous:\tNon-Synonymous:\n")
for i in gene_list:
    output_file.write(i + "\t" + str(synonymous_tally[i][0]) + "\t" + str(synonymous_tally[i][1]) + "\t" + str(synonymous_tally[i][0]/float(synonymous_tally[i][2])) + "\t" + str(synonymous_tally[i][1]/float(synonymous_tally[i][2]))+ "\n")

output_file.write("\n\n\n*** Summary of transition/transverse mutations by rarity: ***\n")
output_file.write("Rarity:\tTC + AG:\tGA + CT:\tGC + CG:\tAC + TG:\tAT + TA:\tGT + CA:\n")
for i in ('clonal', 'subclonal', 'rare'):
    output_file.write(i + "\t" + str(mutation_tally[i]['TC'] + mutation_tally[i]['AG']) + "\t" + str(mutation_tally[i]['GA'] + mutation_tally[i]['CT']) + "\t" + str(mutation_tally[i]['GC'] + mutation_tally[i]['CG']) + "\t" + str(mutation_tally[i]['AC'] + mutation_tally[i]['TG']) + "\t" + str(mutation_tally[i]['AT'] + mutation_tally[i]['TA']) + "\t" + str(mutation_tally[i]['GT'] + mutation_tally[i]['CA']) + "\n")

 
#close i/o files
reference_file.close()
mutpos_file.close()
output_file.close()
haplogrep_out.close()
mitomaster_out.close()
                  