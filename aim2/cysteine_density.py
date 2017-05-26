#!/usr/bin/env python

# Calculating the cysteine density in a moving frame window for FASTA file containing numerous sequences

#to keep the program for general use, argparse and sys used to allow the input pile to be user defined through the command line

import argparse
import sys

MAX_C_COUNT = 5;
FRAME_SIZE = 10;
FRAME_ADVANCE = 10;
CYSTEINE = 'C'

parser = argparse.ArgumentParser(description='Calculate cysteine density for protein sequences.')

parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),default=sys.stdin)
args = parser.parse_args()

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
for record in SeqIO.parse(args.infile, "fasta"):
    seq = str(record.seq)
    my_prot = ProteinAnalysis(seq)
    #print '{}\t {}\t {}\t {}'.format(record.id, my_prot.molecular_weight(), my_prot.isoelectric_point(), seq)

    knotFound = False;
    cCount = 0;

#CREATE another loop for multiple, advance by 3, 5 and 10 OR copy the loop 3 times

    #Looping through sequence, by size FRAME_ADVANCE
    for i in range(0, len(seq), FRAME_ADVANCE):
        frame = seq[i:FRAME_SIZE+i]

        #print frame

        #Count how many CYSTEINE are in frame
        for j in range(len(frame)):
            if (CYSTEINE in frame[j]):
                cCount = cCount + 1;

        #DOES A KNOT EXIST?
        if (cCount >= MAX_C_COUNT):
            knotFound = True;

    if (knotFound):
        print "knot found: " + record.id
