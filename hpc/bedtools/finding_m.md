#Extracting sequences from fasta file that start with an M using bioawk



#displaying the name with ">" suffix and sequence underneath to match fasta format using printf

```bash
bioawk -c fastx '{printf("%s\n%s\n", ">"$name, $seq)}'  really_novel_uniqueness_pg.fasta | head
```

#index function to provide reference for substring function

```bash
bioawk -c fastx '{printf("\n" $seq) index($seq, "M")}'  really_novel_uniqueness_pg.fasta | head
```

#inserted index function into substring function
```bash
bioawk -c fastx '{printf("\n" $seq) substr($seq, index($seq, "M"))}'  really_novel_uniqueness_pg.fasta | head
```

#combining printf, substring and index to extract the novel sequences that start with M

```bash
bioawk -c fastx '{printf("%s\n%s\n",">"$name, substr($seq, index($seq, "M")))}'  really_novel_uniqueness_pg.fasta
```

#Put output in a fasta file called novel_m_seq.fasta

```bash
bioawk -c fastx '{printf("%s\n%s\n", ">"$name, substr($seq, index($seq, "M")))}'  really_novel_uniqueness_pg.fasta > novel_m_seq.fasta
```
