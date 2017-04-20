# Collecting characteristics of novel sequences

obtaining signal peptide information from the extracted protein sequences to start with M (novel_m_seq.fasta) using Signalp-4.1

changed file format from .fasta to .fsa so Signalp-4.1 could read it. Default settings were used (eukaryote type and short format )

```bash
./signalp -t euk -f short novel_m_seq.fsa > novel_m_seq.fsa_short_out
```

#displaying the 10th column (signalp Y/N )

```bash
cat novel_m_seq.fsa_short_out  | awk '{print $10}'
```

#Added a pattern stating if column 10 equals Y, print column 10

```bash
cat novel_m_seq.fsa_short_out  | awk '$10=="Y"{print $10}'
```

#| show number of lines for which this is true
#(count. resulting in 23)

```bash
cat novel_m_seq.fsa_short_out  | awk '$10=="Y"{print $10}' | wc -l
```

#count the '>' at the beginning of a line (^) (equals sequence names (353)

```bash
grep -c '^>' novel_m_seq.fsa
```

#count C for each sequence

```bash
bioawk -c fastx ' {print (split($0,a,"C")-1) }' novel_m_seq.fsa
```

#sequence length

```bash
bioawk -c fastx '{ print $name, length($seq) }' novel_m_seq.fsa
```
