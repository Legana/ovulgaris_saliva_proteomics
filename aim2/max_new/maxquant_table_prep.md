
# Maxquant works


#take the column with IDs from Maxquant in it "unpacked" with rstudio

```bash
cat lfq_long.tsv | awk '{print $2}' > maxq_ids.txt
```

#use ids to extract corresponding sequences from maxquant input file

```bash
while read protein_id; do samtools faidx ovulgaris_maxquantdb.fasta "$protein_id"; done < <(awk -F '\t' '{printf("%s\n",$1)}' maxq_ids.txt) > maxquant_db_new.fasta
```

#use SED to substitute the stop codon to a blank

```bash
cat maxquant.fasta | sed s/*// > maxquant_fresh.fasta
```

#remove first 2 lines as it had the Majority protein id heading

```bash
cat maxquant_fresh.fasta | tail -n +2 > maxquant_fresh_as.fasta
```

#BLAST against online database

```bash
blastp -query maxquant_fresh_as.fasta -db swissprot -remote -outfmt '6 qaccver saccver pident length evalue mismatch gapopen qstart qend sstart send stitle' -max_target_seqs 1 > blastp_swissprot.tsv
```

```bash
blastp -query maxquant_fresh_as.fasta -db nr -remote -outfmt '6 qaccver saccver pident length evalue mismatch gapopen qstart qend sstart send stitle' -max_target_seqs 1 > blastp_nr.tsv
```

#Blast against local databases downloaded from uniprot toxins (all manually reviewed venom proteins and toxins) and arachnoserver (all peptides) on 27/06/2017 with 6659 and 1838 protein sequences, respectively.

#made into a BLAST database first

```bash
makeblastdb -in uniprot-venom_toxins_metazoa.fasta -out uniprotBLAST -dbtype prot
```

```bash
blastp -query maxquant_fresh_as.fasta -db uniprotBLAST -outfmt '6 qaccver saccver pident length evalue mismatch gapopen qstart qend sstart send staxids stitle' -max_target_seqs 1 > blastp_uniprot.tsv
```

```bash
makeblastdb -in arachnoserver.fasta -out arachnoserverBLAST -dbtype prot
```

```bash
blastp -query maxquant_fresh_as.fasta -db arachnoserverBLAST -outfmt '6 qaccver saccver pident length evalue mismatch gapopen qstart qend sstart send stitle' -max_target_seqs 1 > blastp_arachnoserver.tsv
```

#truncate sequences to start with the first M

```bash
bioawk -c fastx '{printf("%s\n%s\n",">"$name, substr($seq, index($seq, "M")))}' maxquant_fresh_as.fasta > maxquant_m_seq.fasta
```

#Run the truncated Maxquant file with SignalP

```bash
./signalp -t euk -f short maxquant_m_seq.fasta > maxquant_short_out
```

#get the name and SignalP (yes/no) columns, tab separated, from the SignalP output and omit the first 3 rows as it was headings junk

```bash
cat maxquant_short_out | awk '{printf("%s\t%s\n",$1,$10)}' | tail -n +3 > maxq_sigp.tsv
```

#run cysteine program to show max number of cysteines per 30aa

```bash
./cysteine_density.py maxquant_fresh_as.fasta > maxq_cyst.tsv
```

#run mwpi program to calculate molecular weight and isoelectric point

```bash
./mwpi.py maxquant_fresh_as.fasta > maxq_mwpi.tsv
```

#get the sequence length

```bash
bioawk -c fastx '{ print $name, length($seq) }' maxquant_fresh_as.fasta > maxq_seqlen.tsv
````
