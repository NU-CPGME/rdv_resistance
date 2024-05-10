# The impact of remdesivir on SARS-CoV-2 evolution <i>in vivo </i>
<b>Ted Ling Hu</b>, Lacy M. Simons, Estefany Rios Guzman, Alexandre Machado de Sant'Anna Carvalho, Maria Francesca Agnes, Arghavan Alisoltanidekhordi, Egon A. Ozer, <b>Ramon Lorenzo-Redondo</b> & <b>Judd F. Hultquist</b>

<hr>

This repository contains the scripts needed to generate the figures and analysis as reported in Ling Hu et al. 2024 (JCI, under review). The script may need to be adapted to the local environment. Due to IRB constraints we are unable to share clinical data used to generate this data. We do however include GISAID accession IDs used to generate the trees in Figure 2. 


# Highlights
<hr>
<ul>
  <li>Remdesivir was the first FDA approved antiviral for the treatment of SARS-CoV-2 and has been widely administered. </li>
  <li>Genomic surveillance conducted by Northwestern University reveals that there is preferential selection of mutations in specific nucleotides that may be implicated in enhanced viral fitness in SARS-CoV-2 when treated with remdesivir.</li>
</ul>

# Summary
<hr>
In this retrospective cohort study, we characterize SARS-CoV-2 diversity and diversification over time in two hospitalized patient cohorts, one who received remdesivir and one who did not. Using viral whole genome sequencing, phylogenetics, and statistical modelling across a total of 385 patients, we find that mutations previously implicated in remdesivir resistance arose rarely in vivo. However, remdesivir administration did result in the preferential selection of mutations across the genome associated with increased fitness. This suggests that variants with enhanced fitness may be selected for in the presence of antiviral therapy as an indirect way to counteract this selective pressure.


# Dependencies
<hr>
Python
<ul>
  <li> Pandas </li>
  <li> Numpy </li>
  <li> scipy </li>
  <li> seaborn </li>
  <li> math </li>
  <li> matplotlib </li>
</ul>
R
<ul>
  <li> dplyr </li>
  <li> treeio </li>
  <li> ggtree </li>
  <li> ggtreeExtra </li>
  <li> ggplot2 </li>
  <li> RColorBrewer </li>
</ul>

MAFFT v7.453

MEGAX v10.1.8.69

IQ-Tree v2.0.5
<ul>
  <li> ModelFinder </li>
</ul>
TreeTime v0.7.6


# Phylogenetic analyses

### Alignment

mafft --auto --thread -1 --keeplength --addfragments Sequences.fasta NC_045512.fasta > Aligned.fasta

### IQtree2 ML Phylogenies

iqtree2 -s Aligned.fasta -T AUTO --alrt 1000 #for Chicago phylogeny also -B 1000 was used

### MEGA X

### iVar

### Sequence IDs from GISAID

<a href="https://github.com/tedlinghu/rdv_resistance/blob/main/Data/Supplementary%20Table%202.xlsx">Accession IDs for sequences from Northwestern</a>
