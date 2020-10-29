#!/usr/bin/python

"""
Program: HeavyChainCDRH3Finder
File: HeavyChainCDRH3Finder.py
Version: 2.0
Date: 28/10/2020
Function: Locate CDRH3 loop within antibody heavy chains and output summary statistics about those Normal_chains
Author: James Sweet-Jones
Address: Institute of Structural and Molecular Biology, Division of Biosciences, University College London
#############################################################################
Description:
===========

Takes in a number of paired Antibody light and heavy chains amino acid sequences in fasta format and outputs summary of chains

#############################################################################

Usage:
=======

Antibody_CHRH3_Finder_2.0.py [x]

where x is a fasta-formatted file where identifiers and sequences are each on a single line
e.g.
>8E10_L|8E10 - (HUMAN) human
EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPADFAVYYCQQYGSSPSITFGQGTRLEIKR
>8E10_H|8E10 - (HUMAN) human
QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHWVRQAPGQRLEWMGWINAGNGNTKYSQKFQGRVTITRDTSASTAYMELSSLRSEDTAVYYCARAMILRIGHGQPQGYWGEGTLVT

Light and heavy chains must be noted in the identifiers with "L|" or "H|" where light chains preceed their paired heavy chain
All sequences in input fasta file must be paired in this format otherwise this version of the script won't work!


#############################################################################
Revision History:
================

v1.0 - took in inputs from command line
v1.1 - takes in inputs from plain text file of heavy chains, one per line
v2.0 - takes in fasta-formatted file of heavy and light chains and returns statistics about chains
v2.2 - takes in fasta-formatted file of paired light and heavy chains and writes instances of normal heavy and light chains to output file.
        statistics about chains are printed to display
"""
#############################################################################
#Import libraries

import sys
import re

#############################################################################


def Heavy_Chain_Identifier(x):
    """

    Input: x --- An antibody heavy chain amino acid sequence in one letter format
    Return: Print statement of statistics regarding CDRH3 loop of x

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
        print(x, "Your sequence is not typical, there are ", Number_of_cysteines ," cysteine residues where there should be 2")
        return False
    elif CDRH3_insertion_modulus > 9:
        print(x, "Your sequence has a particularly long insertion in the CDRH3 loop, which may compromise its binding ability")
        return False
    elif len_CDRH3 == 0:
        print(x, "Your sequence has no or an irregular CDRH3 loop, which may compromise its binding ability")
        return False
    elif Number_of_cysteines == 2 and CDRH3_insertion_length <= 0:
        print(x, "Your sequence is ", k, "residues long with a CDRH3 loop", CDRH3_loop)
        return True
    elif len(W_positions) == 0:
        print(x, "The CDRH3 loop could not be found on this sequence as the tryptophan closing the loop could not be found")
    elif Number_of_cysteines == 2 and 70 <= Cys_distance <= 80:
        print(x, "Your sequence is ", k, "residues long with a CDRH3 loop", CDRH3_loop,". This loop has an insertion of", CDRH3_insertion_length, "residues, which are :", CDRH3_insertion)
        return True
    else:
        print(x)
        print("This is an irregular sequence")
        return False

#############################################################################
def Light_Chain_Identifier(x):
    """

    Input: x --- An antibody light chain amino acid sequence in one letter format
    Return: Print statement of statistics regarding CDRH3 loop of x

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
    if Number_of_cysteines == 2:
        print(x, "This is a normal antibody light chain")
        return True
    else:
        print(x)
        print(x, "This antibody light chain does not have exactly 2 cysteines")
        return False

#*********************************************************
#*** Main program  ***
#*********************************************************


#Run Normal_chain_identifier on input sequences and run tally on normal and irregular sequences
#Return results to display



Normal_antibodies = 0
Irregular_antibodies = 0
with open(sys.argv[1], 'r') as f:
    output = open("output.txt", 'w+')
    for line in f:
        if line[0] == ">" and re.search(r'L\|', line):
            light_chain_identifier = line
            light_chain = f.readline() #line +1
            heavy_chain_identifier = f.readline() #line +2
            heavy_chain = f.readline() #line +3
            if re.search(r'H\|', heavy_chain_identifier):
                if Light_Chain_Identifier(light_chain) == True and Heavy_Chain_Identifier(heavy_chain) == True:
                    output.write(light_chain_identifier)
                    output.write(light_chain)
                    output.write(heavy_chain_identifier)
                    output.write(heavy_chain)
                    light_chain_identifier = None
                    light_chain = None
                    heavy_chain_identifier = None
                    heavy_chain = None
                    Normal_antibodies += 1
                else:
                    print("One of the chains of this antibody is not normal")
                    Irregular_antibodies += 1
            else:
                print("WARNING: There are sequences that are not paired")



print("You have entered", Normal_antibodies, "normal antibodies,  ", Irregular_antibodies ," irregular antibodies")
