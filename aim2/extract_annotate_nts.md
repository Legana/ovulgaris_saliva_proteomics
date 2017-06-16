#extract partial RNA sequence names to match DNA names from Trinity.FASTA

```bash
?????
```
###tricky as the IDS are 6 frame translated proteins e.g >lcl|TRINITY_DN10130_c0_g1_i1_frame_6_orf_2 as well as longest ORF e.g. >lcl|TRINITY_DN11042_c1_g1::TRINITY_DN11042_c1_g1_i1::g.5994::m.5994 and need to transform this to: >TRINITY_DN10130_c0_g1_i1 and >TRINITY_DN11042_c1_g1_i1 in order to match the Trinity file.


#place extracted names to text file (>ids.txt)

#remove the lcl| prefix from IDs

```bash
sed 's/lcl|//g' ids.txt
```

#use the ids.txt to extract corresponding nucleotide sequences (and output to a file)

```bash
xargs samtools faidx Trinity.fasta < ids.txt > 40.fasta
```

#use the ids.txt to extract corresponding transdecodor annotations

```bash
while read line; do grep $line longest_orfs.gff3; done < ids.txt > 40transdecodor.gff3
```
#do the same for known_pg and novel_pg gff3 
