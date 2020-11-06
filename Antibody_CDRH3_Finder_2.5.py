#!/usr/bin/python

"""
Program: Antibody_CDRH3_Finder_2.3
File: Antibody_CDRH3_Finder_2.3.py
Version: 2.3
Date: 30/10/2020
Author: James Sweet-Jones
Address: Institute of Structural and Molecular Biology, Division of Biosciences, University College London
#############################################################################
Description:
===========

Takes in a number of paired Antibody light and heavy chains amino acid sequences in fasta format and outputs summary of chains

#############################################################################

Usage:
=======

Antibody_CHRH3_Finder_2.3.py [x]

where x is a fasta-formatted file where identifiers and sequences are wrapped or unwrapped
e.g.
>8E10_L|8E10 - (HUMAN) human
EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAW
YQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDF
TLTISRLEPADFAVYYCQQYGSSPSITFGQGTRLEI
KR
>8E10_H|8E10 - (HUMAN) human
QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHW
VRQAPGQRLEWMGWINAGNGNTKYSQKFQGRVTITR
DTSASTAYMELSSLRSEDTAVYYCARAMILRIGHGQ
PQGYWGEGTLVT

or

>8E10_L|8E10 - (HUMAN) human
EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPADFAVYYCQQYGSSPSITFGQGTRLEIKR
>8E10_H|8E10 - (HUMAN) human
QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHWVRQAPGQRLEWMGWINAGNGNTKYSQKFQGRVTITRDTSASTAYMELSSLRSEDTAVYYCARAMILRIGHGQPQGYWGEGTLVT


Light and heavy chains must be noted in the identifiers with "L|" or "H|" where order of light/heavy chain does not matter.
Identifiers passed "|" must be identical. All sequences in input fasta file must be paired in this format otherwise this version of the script won't work!


#############################################################################
Revision History:
================

v1.0 - took in inputs from command line
v1.1 - takes in inputs from plain text file of heavy chains, one per line
v2.0 - takes in fasta-formatted file of heavy and light chains and returns statistics about chains
v2.2 - takes in fasta-formatted file of paired light and heavy chains and writes instances of normal heavy and light chains to output file.
        statistics about chains are printed to display
v2.3 - takes in wrapped or unwrapped fasta-formatted file of paired light and heavy chains and writes pairs normal heavy and light chains to output file.
v2.4 - improved stringency on identifying paired sequences
v2.5 - takes in fasta format file of paired heavy and light chains where order does not matter
"""
#############################################################################
#Import libraries

import sys
import re
import os

#############################################################################

def correct_input_format(input):
    """
    Input: a fasta file that may be either wrapped sequence text unwrapped
    Return: file an intermediate file that will be unwrapped
    30/10/20 Original by JSJ
    """

    with open(sys.argv[1],"r") as f:
        sequences = f.read()
        sequences = re.split("^>", sequences, flags=re.MULTILINE) # Only splits string at the start of a line.
        del sequences[0]
        f.close()

    with open("correctformat.temp.txt", "w+") as correctformat:
        for fasta in sequences:
            header, sequence = fasta.split("\n", 1)
            header = ">" + header + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
            sequence = sequence.replace("\n","") + "\n"
            correctformat.write(header + sequence)
        correctformat.close()

##############################################################################
def Heavy_Chain_Identifier(x):
    """

    Input: x --- An antibody heavy chain amino acid sequence
    Return: true if chain is normal, false if not

    23/10/2020 Original by JSJ
    """

    k = len(x)
    Number_of_cysteines = 0
    CAR = 0
    First_Cys_motif_start_position = 0
    Second_Cys_motif_start_position = 0
    WG = 0
    WG_position = 0
    CDRH3_loop = ""
    Max_CDRH3_insertions = 8

    #Define number and location of cysteine residues in sequence and if the distance between residues is within expected range
    for i in range(k):
        if x[i] == "C" and Number_of_cysteines == 0:
            Number_of_cysteines += 1
            First_Cys_motif_start_position = i+1
            continue
        elif x[i] == "C" and Number_of_cysteines == 1:
            Second_Cys_motif_start_position = i+1
            Number_of_cysteines += 1
        elif x[i] == "C" and Number_of_cysteines >1:
            Number_of_cysteines += 1
    Cys_distance = Second_Cys_motif_start_position - First_Cys_motif_start_position

    #Define number and location of tryptophan residues in sequence and locate the last tryptophan

    W_positions = []
    for i in range(Second_Cys_motif_start_position-1, k):
        if x[i] == "W":
            W_positions.append(i+1)

    #Define CDRH3 loop by printing all residues between second cysteine residue +2 and final tryptophan residue

    if len(W_positions) > 0:
        WG_position = W_positions[-1]
        for i in range(Second_Cys_motif_start_position+2, WG_position-1):
            CDRH3_loop = CDRH3_loop + x[i]

        #Calculates number of inserted bases by subtracting length of CDRH3 by length of CDRH3 without insertions (8)
    len_CDRH3 = len(CDRH3_loop)
    CDRH3_insertion_length = len(CDRH3_loop) - Max_CDRH3_insertions
    CDRH3_insertion = ""
    CDRH3_insertion_modulus = len(CDRH3_insertion)
    if CDRH3_insertion_length != 0:
        for i in range(WG_position-3-CDRH3_insertion_length, WG_position-3):
            CDRH3_insertion = CDRH3_insertion + x[i]

    #Print out results

    if Number_of_cysteines != 2:
        return False
    elif CDRH3_insertion_modulus > 9:
        return False
    elif len_CDRH3 == 0:
        return False
    elif Number_of_cysteines == 2 and CDRH3_insertion_length <= 0:
        return True
    elif len(W_positions) == 0:
        return False
    elif Number_of_cysteines == 2 and 70 <= Cys_distance <= 80:
        return True
    else:
        return False

#############################################################################
def Light_Chain_Identifier(x):
    """

    Input: x --- An antibody light chain amino acid sequence in one letter format
    Return: true if chain is normal, false if not

    28/10/2020 Original by JSJ
    """
    k = len(x)
    Number_of_cysteines = 0
    First_Cys_motif_start_position = 0
    Second_Cys_motif_start_position = 0
    #Define number and location of cysteine residues in sequence and if the distance between residues is within expected range
    for i in range(k):
        if x[i] == "C" and Number_of_cysteines == 0:
            Number_of_cysteines += 1
            First_Cys_motif_start_position = i+1
            continue
        elif x[i] == "C" and Number_of_cysteines == 1:
            Second_Cys_motif_start_position = i+1
            Number_of_cysteines += 1
        elif x[i] == "C" and Number_of_cysteines >1:
            Number_of_cysteines += 1
        else:
            continue
    Cys_distance = Second_Cys_motif_start_position - First_Cys_motif_start_position
    if Number_of_cysteines == 2 and Cys_distance:
        return True
    else:
        return False

#*********************************************************
#*** Main program  ***
#*********************************************************

#if using wrapped sequences, must use converter script first!!!
#Run Normal_chain_identifier on input sequences and run tally on normal and irregular sequences
#Return results to display



Normal_antibodies = 0
Irregular_antibodies = 0
with open(sys.argv[1], 'r') as input:
    correct_input_format(input)
    f = open("correctformat.temp.txt","r")
    output = open("Initial_screening_output.txt", 'w+')
    filtered = open("Initial_screening_filtered_out.txt", 'w+')
    for line in f: # loop through temp file to locate light chains
        if line[0] == ">" and re.search(r'L\|', line):
            light_chain_identifier             = line
            light_chain_identifier_split       = light_chain_identifier.split("|")
            light_chain_identifier_split_check = str(light_chain_identifier_split[1])
            light_chain                        = f.readline() #line +1
            light_chain_removed_dels           = re.sub('X','',light_chain) #Remove unidentified/deleted amino acids (X) from sequence before writing it otherwise modelling software will reject
            heavy_chain_identifier             = f.readline() #line +2
            heavy_chain_identifier_split       = heavy_chain_identifier.split("|")
            heavy_chain_identifier_split_check = str(heavy_chain_identifier_split[1])
            heavy_chain                        = f.readline() #line +3
            heavy_chain_removed_dels           = re.sub('X','',heavy_chain) #Remove unidentified/deleted amino acids (X) from sequence before writing it
            if re.search(r'H\|', heavy_chain_identifier) and light_chain_identifier_split_check == heavy_chain_identifier_split_check:
            #If paired heavy and light chain are both "normal" then we consider them as one normal antibody
                if Light_Chain_Identifier(light_chain_removed_dels) == True and Heavy_Chain_Identifier(heavy_chain_removed_dels) == True:
                    output.write(light_chain_identifier) #Write normal antibody sequences to output in fasta format
                    output.write(light_chain_removed_dels)
                    output.write(heavy_chain_identifier)
                    output.write(heavy_chain_removed_dels)
                    light_chain_identifier = None
                    light_chain            = None
                    heavy_chain_identifier = None
                    heavy_chain            = None
                    Normal_antibodies      += 1
                else:
                    filtered.write(light_chain_identifier) #Write irregular antibody sequences to filtered output file in fasta format
                    filtered.write(light_chain)
                    filtered.write(heavy_chain_identifier)
                    filtered.write(heavy_chain)
                    Irregular_antibodies += 1
            elif light_chain_identifier_split_check != heavy_chain_identifier_split_check:
                print("WARNING: " + light_chain_identifier_split_check.strip() + " is not a paired sequence")


        elif line[0] == ">" and re.search(r'H\|', line):
            heavy_chain_identifier             = line
            heavy_chain_identifier_split       = heavy_chain_identifier.split("|")
            heavy_chain_identifier_split_check = str(heavy_chain_identifier_split[1])
            heavy_chain                        = f.readline() #line +1
            heavy_chain_removed_dels           = re.sub('X','',heavy_chain) #Remove unidentified/deleted amino acids (X) from sequence before writing it otherwise modelling software will reject
            light_chain_identifier             = f.readline() #line +2
            light_chain_identifier_split       = light_chain_identifier.split("|")
            light_chain_identifier_split_check  = str(light_chain_identifier_split[1])
            light_chain                        = f.readline() #line +3
            light_chain_removed_dels           = re.sub('X','',light_chain) #Remove unidentified/deleted amino acids (X) from sequence before writing it
            if re.search(r'L\|', light_chain_identifier) and light_chain_identifier_split_check == heavy_chain_identifier_split_check:
            #If paired heavy and light chain are both "normal" then we consider them as one normal antibody
                if Light_Chain_Identifier(light_chain_removed_dels) == True and Heavy_Chain_Identifier(heavy_chain_removed_dels) == True:
                    output.write(light_chain_identifier) #Write normal antibody sequences to output in fasta format
                    output.write(light_chain_removed_dels) #written light chain first to imput into modelling script
                    output.write(heavy_chain_identifier)
                    output.write(heavy_chain_removed_dels)
                    light_chain_identifier = None
                    light_chain            = None
                    heavy_chain_identifier = None
                    heavy_chain            = None
                    Normal_antibodies      += 1
                else:
                    filtered.write(light_chain_identifier) #Write irregular antibody sequences to filtered output file in fasta format
                    filtered.write(light_chain)
                    filtered.write(heavy_chain_identifier)
                    filtered.write(heavy_chain)
                    Irregular_antibodies += 1

            elif light_chain_identifier_split_check != heavy_chain_identifier_split_check:
                print("WARNING: " + heavy_chain_identifier_split_check.strip() + " is not a paired sequence")



f.close()
os.remove("correctformat.temp.txt") #delete temp file where lines are unwrapped
output.close()
print("You have entered", Normal_antibodies, "normal antibodies,  ", Irregular_antibodies ," irregular antibodies")
