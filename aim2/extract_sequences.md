# extracting sequences corresponding to maxquant hits

#extracting the sequences obtained by maxquant from the input fasta file  

```bash
while read protein_id; do samtools faidx ../../known_novel_sequences.fsa "$protein_id"; done < <(awk -F '\t' '{printf("%s\n",$1)}' maxquant_proteins.tsv)
```
