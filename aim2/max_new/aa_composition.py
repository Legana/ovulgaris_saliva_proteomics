#!/usr/bin/env python

# Calculate the amino acid composition in percentage for a FASTA file containing numerous sequences

#to keep the program for general use, argparse and sys used to allow the input pile to be user defined through the command line

import argparse
import sys
import math
parser = argparse.ArgumentParser(description='Calculate amino acid composition for protein sequences.')

parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),default=sys.stdin)
parser.add_argument('-c','--comparison-protein', type=int, default=1, help='Choice of protein')
args = parser.parse_args()


comparison_protein = args.comparison_protein

#add constants
cz1 = {'K':0.0238, 'H':0.0082, 'R':0.0164, 'B':0.0450, 'T':0.0612, 'S':0.0292, 'Z':0.0558, 'P':0.0364, 'G':0.0156, 'A':0.0230,
'V':0.0274, 'M':0.0068, 'I':0.0196, 'L':0.0244, 'Y':0.0186, 'F':0.0420, 'W':0.0198, 'C':0.0082, 'N':0, 'O':0.0, 'U':0.0, 'J':0.044}
cz2 = {'K':0.0306, 'H':0.0314, 'R':0.0312, 'B':0.0708, 'T':0.0332, 'S':0.0240, 'Z':0.0640, 'P':0.0341, 'G':0.0194, 'A':0.0251,
'V':0.0341, 'M':0.0134, 'I':0.0239, 'L':0.0485, 'Y':0.0388, 'F':0.0471, 'W':0.0285, 'C':0.0121, 'N':0, 'O':0.0, 'U':0.0, 'J':0.0724}

#option to use two protein options (either cz1 - default, or cz2 (-c2))
if comparison_protein == 1:
    cz=cz1
else:
    cz=cz2

#normalise values so sum becomes 1

cz_sum = sum(cz.values())


for k in cz:
  cz[k] = cz[k]/cz_sum

#to read from a FASTA file with a loop over entries using SeqIO define the FASTA sequences and analyse them by ProteinAnalysis

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
for record in SeqIO.parse(args.infile, "fasta"):
    seq = str(record.seq)
    my_prot = ProteinAnalysis(seq)
    aa_percentage = my_prot.get_amino_acids_percent()


    dot_product = 0
    for k in aa_percentage.keys():
        if k == 'D' or k == 'N':
            dot_product += aa_percentage[k]*cz['B']
        elif k == 'E' or k == 'Q':
            dot_product += aa_percentage[k]*cz['Z']
        elif k == 'I' or k == 'L':
            dot_product += aa_percentage[k]*cz['J']
        else:
            dot_product += aa_percentage[k]*cz[k]



    print '{}\t{}\t{}\t{}'.format(record.id,math.sqrt(dot_product),my_prot.molecular_weight(),my_prot.isoelectric_point())
