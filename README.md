# PanMAN and PanMAP tutorial (PEQG 2026 Pangenome Workshop)

[Intro summary description]

This is a tutorial for building pangenome trees and placing reads or samples onto them using Panmap.

## Environment setup

All programs used in this tutorial can be installed using conda: [Panmap](https://github.com/amkram/panmap),
[DIPPER](https://github.com/TurakhiaLab/DIPPER), and [TWILIGHT](https://github.com/TurakhiaLab/TWILIGHT). Docker is also
supported, and instructions can be found in each repository.

```bash
conda create -n pantree python=3.11 -y
conda activate pantree

conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge
conda config --env --set channel_priority strict

conda install -c conda-forge -c bioconda \
  panmap \
  panman \
  dipper \
  twilight \
  -y
```

## Build a panman

**Make a guide tree and alignment.** PanMAN can be built from several different input types: 
[Pangraph](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001034#tab2) output, 
[GFA](https://github.com/GFA-spec/GFA-spec), or MSA + nwk. For simplicity, this tutorial will build a 1000 SARS-CoV-2 
genome PanMAN from an MSA and guide tree.

```bash
git clone https://github.com/AlanZhangUCSC/PEQG_2026_pangenome_trees.git && cd PEQG_2026_pangenome_trees

mkdir -p example_output

dipper -i r -I input_data/sars.subsampled_1000.fa \
  -O example_output/sars.subsampled_1000.guide.nwk

twilight --cpu-only -i input_data/sars.subsampled_1000.fa \
  -t example_output/sars.subsampled_1000.guide.nwk \
  -o example_output/sars.subsampled_1000.guide.aln
```

**Build a PanMAN from alignment and guide tree.** Note that PanMANs built from MSA doesn't contain information on
structural variations. Currently, PanMAN doesn't infer SVs internally and uses SV information from Pangraph or GFA.

```bash
panmanUtils -M example_output/sars.subsampled_1000.guide.aln \
  -N example_output/sars.subsampled_1000.guide.nwk \
  -o sars.subsampled_1000
```

The summary stats shows that the tree contains 1969 nodes ~26,000 combined mutations (~13 mutations / branch), which is
within the expected range for 1000 randomly sampled SARS-CoV-2 genomes. This is a small example for the purpose of the
tutorial. For the rest of the demo, I will use larger, prebuilt PanMANs.

<details>
<summary>Click to see summary stats</summary>

```console
$ panmanUtils --summary sars.subsampled_1000.panman  | head -n14
starting reading panman
Data load time: 24799370 nanoseconds 
PanMAN Summary
No of Trees: 1
No of Recombinations: 0
Tree 0 Summary
Total Nodes in Tree: 1969
Total Samples in Tree: 985
Total Substitutions: 13793
Total Insertions: 6872
Total Deletions: 5225
Total Inversions: 0
Max Tree Depth: 47
Mean Tree Depth: 29.1127
```

</details>

## Build a Panmap index

Panmap sketches [syncmer](https://peerj.com/articles/10805/) seeds of all the genomes in the PanMAN and stores them 
compactly in a phylogeny-guided format. Similar to PanMAN, branches are annotated with the changes in syncmer sketches 
between parent and child nodes. 

```bash
panmap input_data/sars_20000_twilight_dipper.panman \
  --index input_data/sars_20000_twilight_dipper.panman.idx \
  -k 19 -s 8 -l 3
```

Here I used the default seed parameters (k=19, s=8, l=3), where `k` is the k-mer size, `s` is the s-mer size for syncmer 
selection, and `l` is the number of linked syncmers (same as 
[k-min-mers](https://genome.cshlp.org/content/early/2023/08/14/gr277679123) but with syncmers instead of minimizers). We
found that linked-syncmers provide slightly higher accuracy than single syncmers.

## Run Panmap

**Panmap has two modes: single sample mode and metagenomic mode.** In single sample mode, Panmap assumes the sample
contains a single haplotype and places the sample as a single unit on the tree. In metagenomic mode, Panmap assumes the
sample contains a mix of haplotypes and places each read separately on the tree.

## Place and genotype a single-haplotype sample

[Run Panmap: single sample mode]

[Run Panmap: metagenomic mode (haplotype deconvolution)]

[Run Panmap: metagenomic mode (reads assignment)]


