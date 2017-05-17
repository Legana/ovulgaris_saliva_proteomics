# Removing stopcodons from a FASTA file

#use SED to substitute one expression to another (in this case substitute * to nothing)

```bash
cat maxquant.fasta | sed s/*// > maxquant_clean.fasta
```
