# extracting sequences corresponding to maxquant hits

#extracting the sequences obtained by maxquant from the input fasta file  

```bash
while read protein_id; do samtools faidx ../../known_novel_sequences.fsa "$protein_id"; done < <(awk -F '\t' '{printf("%s\n",$1)}' maxquant_proteins.tsv)
```
#Extract IDs obtained from various BLAST results (saved from Rstudio write_delim(blast_uniprot,"blast_uniprot.tsv",delim = "\t")) and extract according sequences from original fasta file
#(now blasting these through geneious will not take as long)

```bash
while read protein_id; do samtools faidx maxquant_clean.fasta "$protein_id"; done < <(awk '{printf("%s\n",$1)}' blast_uniprot.tsv)
```

#output it to FASTA file so it can be imported to Geneious and blasted with BLASTP.

```bash
while read protein_id; do samtools faidx maxquant_clean.fasta "$protein_id"; done < <(awk '{printf("%s\n",$1)}' blast_uniprot.tsv) > blastp_uniprot.fasta
```
