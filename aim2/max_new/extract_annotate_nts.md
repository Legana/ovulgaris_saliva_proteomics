#delete lines with ID headings
sed 's/ID//g' putt_ids_u.txt > putt_ids_noheadings.txt

#delete blank lines
sed '/^$/d' putt_ids_noheadings.txt > putt_ids_noheadings_orspaces.txt

#keep only unique protein ids
cat putt_ids_noheadings_orspaces.txt | sort -u > putt_prot_ids.txt

#extract partial RNA sequence names to match DNA names from Trinity.FASTA (remove >lcl| pattern)

```bash
cat putt_prot_ids.txt | sed -E 's/.*(TRINITY_[A-Z]+[0-9]+_c[0-9]+_g[0-9]_i[0-9]+).*/>lcl|\1/' | sed 's/>lcl|//g' > putt_nucl_ids.txt
```

or

```bash
cat putt_prot_ids.txt | sed -E 's/.*(TRINITY_[A-Z]+[0-9]+_c[0-9]+_g[0-9]+_i[0-9]+).*/\1/'
```

###tricky as the IDS are 6 frame translated proteins e.g >lcl|TRINITY_DN10130_c0_g1_i1_frame_6_orf_2 as well as longest ORF e.g. >lcl|TRINITY_DN11042_c1_g1::TRINITY_DN11042_c1_g1_i1::g.5994::m.5994 and need to transform this to: >TRINITY_DN10130_c0_g1_i1 and >TRINITY_DN11042_c1_g1_i1 in order to match the Trinity file.

#remove the tr sequences as these are not found in our Transdecodor file

```bash
cat putt_nucl_ids.txt | sed '/^tr/d' > putt_nucl_ids_no_tr.txt
```

#use the ids.txt (without the 3 swissprot tr sequences) to extract corresponding nucleotide sequences (and output to a file)

```bash
xargs samtools faidx Trinity.fasta < putt_nucl_ids_no_tr.txt > putt_nt.fasta
```

#same but for protein sequences

```bash
xargs samtools faidx maxquant_fresh_as.fasta < putt_prot_ids.txt > putt_prot.fasta
```

#use the ids.txt to extract corresponding transdecodor annotations

```bash
while read line; do grep $line longest_orfs.gff3; done < putt_nucl_ids_no_tr.txt > putt_transdecodor.gff3
```
#do the same for known_pg and novel_pg gff3

```bash
while read line; do grep $line known_pg.gff3; done < putt_nucl_ids_no_tr.txt > putt_known.gff3
```

```bash
while read line; do grep $line novel_pg.gff3; done < putt_nucl_ids_no_tr.txt > putt_novel.gff3
```
