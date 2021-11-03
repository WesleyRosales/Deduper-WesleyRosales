#!/usr/bin/env python

import argparse
import re

# PCR duplicates can be identified by having the same chromosome (col 3), position (col 4), strand (col 2), and Unique Molecular Index (UMI). This program should check for matches in these 4 categories as well as consider soft clipping in the algorithm.

#1. Use Samtools to sort by position

#2. Argparse for input .sam file, UMI file, and output filename.
def get_args():
    my_parser = argparse.ArgumentParser(prog='deduplicate',
                description='Elimintates PCR duplicate reads in sorted sam files, retaining only the first duplicate read. Outputs a deduplicated sam file in the directory the program is called.')

    my_parser.add_argument('-f',
                        '--file',
                        action='store',
                        help='enter sam file',
                        required=True)
    
    my_parser.add_argument('-p',
                        '--paired',
                        action='store_true',
                        help='Add if sam file contains paired reads',
                        required=False)

    my_parser.add_argument('-u',
                        '--umi',
                        action='store',
                        help='enter file containing list of UMIs separated by newlines',
                        required=False)

    my_parser.add_argument('-o',
                        '--output',
                        action='store',
                        help='enter output file name',
                        required=True
    )
    return my_parser.parse_args()

args = get_args()

#If user indicates that the file is paired-end, quit the program.
if args.paired:
    print("No paired-end functionality. Quitting.")

#Function to isolate the UMI in the NAME column of a sam file
def findUMI(preUMI: str) -> str:
    import re
    umi = re.search(":[TCGA]+$", preUMI).group()[1:]
    return umi

#Function to determine strand of read using bitwise FLAG
def determineStrand(flag: int) -> str:
    if ((flag & 16) != 16):
        return "f"
    else:
        return "r"

#Function to alter start_pos to reflect read start position
def readStart(forward_strand: bool, CIGAR: str, position: int) -> int:

    #adjust forward strand position only if a soft clip is detected in beginning of CIGAR string

    if forward_strand:
        soft_clip = re.search("^\d+S", CIGAR) #search for soft clipped integer in beginning of string
        if soft_clip != None:  #change the position of soft clipped reads in the list to its real position ONLY if it will affect START position
            soft_clip = soft_clip.group()[:-1]
            adjusted_pos = position - int(soft_clip)
            return adjusted_pos
        
        else:
            return position

    #splits CIGAR string into list of chunks of integers with letters. These are the individual matches/mismatches, deletions, insertions, etc within the CIGAR string
    cigar_list = re.findall("\d+\w", CIGAR)

    #Remove softclipping in beginning of CIGAR string list if there is one
    if cigar_list[0].endswith("S"):
        del cigar_list[0]

    #add integers that are not insertions to a list to add to start position
    add_list = []
    for i in cigar_list:
        if i.endswith("I") != True:
            add_list.append(int(i[:-1]))
    
    #add the sum of the add_list to the position. return the adjusted position, this is the start position of the read itself on the reverse strand
    adjusted_pos = int(position) + sum(add_list)
    return adjusted_pos

##save the UMIs to a list
#open entire UMI file 
umi_file = open(args.umi, "r")

#read UMI file as string and split by "\n"
umi_list = umi_file.read().split("\n")
    
umi_file.close()

#open the sam input file and the sam output file
input_file = open(args.file, "r")
output_file = open(args.output, "w")

position_dict = {}
for umi in umi_list:
    position_dict[umi] = []
previous_chr = "0"
previous_ident = ""
while True:
    inp_line = input_file.readline()
    inp_list = inp_line.split("\t")
    if inp_line == "":
        break

    #3. Write all header rows to output file.
    if re.match("^@", inp_list[0]):
        output_file.writelines(inp_line)
        continue
    
    #to speed up processing of reads, reset position list after every chromosome
    if inp_list[2] != previous_chr:
        for umi in umi_list:
            position_dict[umi] = []
        previous_chr = inp_list[2]

    # isolate UMI
    umi = findUMI(inp_list[0])

    # determine strand read is on with FLAG of alignment. "f" = forward, "r" = reverse
    strnd = determineStrand(int(inp_list[1]))
    
    #create an identifier for each unique read. Any read with the same identifier will be removed
    ident = "-".join([umi,
                    strnd,
                    inp_list[2],
                    str(readStart(strnd == "f",
                                inp_list[5], 
                                int(inp_list[3])))])

    #any reads with identical position identification will be PCR duplicates because they have the same UMI, strand, chromosome, and position.
    if ident != previous_ident:
        if umi in umi_list:
            if ident not in position_dict[umi]:
                previous_ident = ident[:]
                output_file.writelines(inp_line)
                position_dict[umi].append(ident)

input_file.close()
output_file.close()