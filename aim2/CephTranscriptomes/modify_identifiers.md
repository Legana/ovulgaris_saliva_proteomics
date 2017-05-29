# Modify identifiers in FASTA files

#edit the sequence identifiers assembled with Trinity to have a more meaningful name relating to the species from which it came
#in addition change the | in the identifier to _

#therefore it changes from
#>contig00004|m.4 contig00004|g.4 type:internal len:188 gc:universal contig00004:1-561(+)
#to
#>aacu_contig00004_m.5 contig00004_g.5 type:3prime_partial len:59 gc:universal contig00004:176-3(-)

```bash
cat abdopus_aculeatus_longestorfs.fasta | sed 's/>/>aacu_/' | sed s/\|/_/
```
#change the | to _ in the second part of the identifier name as well

```bash
cat abdopus_aculeatus_longestorfs.fasta | sed 's/>/>aacu_/' | sed s/\|/_/ | sed s/\|/_/
```

#OR to do the same thing but use "g" at end which changes all | into _

```bash
cat abdopus_aculeatus_longestorfs.fasta | sed 's/>/>aacu_/' | sed s/\|/_/g
```
