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
cz1 = {'K':0.01627, 'H':0.00528, 'R':0.00941, 'B':0.03380, 'T':0.05138, 'S':0.02778, 'Z':0.03816, 'P':0.03162, 'G':0.02077, 'A':0.02581,
'V':0.02339, 'M':0.00455, 'I':0.01493, 'L':0.01859, 'Y':0.01026, 'F':0.0254, 'W':0.00969, 'C':0.00676, 'N':0.0, 'O':0.0, 'U':0.0, 'J':0.03352}
cz2 = {'K':0.02093, 'H':0.02023, 'R':0.01791, 'B':0.05319, 'T':0.02787, 'S':0.02283, 'Z':0.04377, 'P':0.02962, 'G':0.02583, 'A':0.02817,
'V':0.02912, 'M':0.00898, 'I':0.018216, 'L':0.03696, 'Y':0.02141, 'F':0.02851, 'W':0.01395, 'C':0.00998, 'N':0.0, 'O':0.0, 'U':0.0, 'J':0.05517}

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
