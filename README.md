# ovulgaris_saliva_proteomics
Proteomics of Octopus vulgaris saliva and salivary glands. Subheadings refer to specific aspects of the analysis with links to scripts, notes and tools used.


# Running HPC Analyses

These are contained within the folder 'hpc' and require copies of large data files.  The hpc `README` for details on how to access these files

The `proteogenomics` folder within `hpc` contains just small files and a README.  See the readme within this directory for details on running that analysis.


# Distinguishing known and novel peptide with bedtools

-  [peptide_uniqueness.rb](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/hpc/bedtools/peptide_uniqueness.rb)
-  [score_by_uniqueness.rb](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/hpc/bedtools/score_by_uniqueness.rb)
-  [peptides_from_gff.rb](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/hpc/bedtools/peptides_from_gff.rb)
-  [filter_gff.rb](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/hpc/proteogenomics/filter_gff.rb)


# Extracting protein sequences to start with a Methionine

- [finding_m.md](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/hpc/bedtools/finding_m.md)

# Sequence length between known and novel proteins

- [known_novel_characteristics.Rmd](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/known_novel_characteristics.Rmd)

# CRRF analysis of cysteine rich regions in Arachnoserver, Tox-Prot and SwissProt

 - [chasing_cysteines.Rmd](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/aim2/chasing_cysteines.Rmd)

 # Obtaining protein characteristics for the protein table

- [maxquant_table_prep.md](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/aim2/max_new/maxquant_table_prep.md)

 # Creating the protein characteristics table, extracting putative toxin candidates and performing a the Fisher's Exact Test on cysteine rich regions and secreted sequences

- [fresh_maxquant_table.Rmd](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/aim2/max_new/fresh_maxquant_table.Rmd)

# Dot product tool to match sequence amino acid composition to alpha- and beta- cephalotoxin

- [aa_composition.py](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/aim2/max_new/aa_composition.py)

# Toxin domain heatmap

- [heatmap_TD.Rmd](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/aim2/max_new/heatmap_TD.Rmd)

# SSCR proteins heatmap

- [heatmap_SSCR.Rmd](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/aim2/max_new/heatmap_SSCR.Rmd)

# Cephalopod serine protease tree

- [sp_tree5.Rmd] https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/aim2/max_new/phylogenetics/sp_tree5.Rmd

# Serine protease tree including additional taxa and sequences

- [SP_tree.Rmd](https://github.com/Legana/ovulgaris_saliva_proteomics/blob/master/aim2/max_new/phylogenetics/SP_tree.Rmd)
