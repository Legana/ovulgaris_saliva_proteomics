# Making blastp databases and BLASTing against online and local databases

  Online data bases:

Swissprot database

```bash
blastp -query really_novel_uniqueness_pg.fasta -db swissprot -remote -outfmt '6 qaccver saccver pident length evalue mismatch gapopen qstart qend sstart send stitle' -max_target_seqs 1 > blastp_swissprot.txt
```

nr database

```bash
blastp -query really_novel_uniqueness_pg.fasta -db nr -remote -outfmt '6 qaccver saccver pident length evalue mismatch gapopen qstart qend sstart send stitle' -max_target_seqs 1 > blastp_nr.txt
```

  Local databases:

Downloaded all (6622) manually reviewed venom proteins and toxins (canocical and isoform) from Uniprot and made it into a BLAST databases

```bash
makeblastdb -in uniprot-venom_toxins_metazoa.fasta -out uniprotBLAST -dbtype prot
```

```bash
blastp -query really_novel_uniqueness_pg.fasta -db uniprotBLAST -outfmt '6 qaccver saccver pident length evalue mismatch gapopen qstart qend sstart send staxids stitle' -max_target_seqs 1 > blastp_uniprot.txt
```

Creating a arachnoserver database

```bash
makeblastdb -in arachnoserver.fasta -out arachnoserverBLAST -dbtype prot
```

```bash
blastp -query really_novel_uniqueness_pg.fasta -db arachnoserverBLAST -outfmt '6 qaccver saccver pident length evalue mismatch gapopen qstart qend sstart send stitle' -max_target_seqs 1 > blastp_arachnoserver.txt
```

Converting copied .txt from excel .csv in windows format so must convert to unix

```bash
dos2unix -c mac novel_peptide_characteristics_copy.txt
```

Use join command to merge blastp output files with novel_peptide_characteristics_copy table using the first column that is identical in both to merge

```bash
join -1 1 -2 1 swissb.txt npcc.txt
```
