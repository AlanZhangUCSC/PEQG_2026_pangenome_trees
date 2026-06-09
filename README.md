# PanMAN and PanMAP tutorial (PEQG 2026 Pangenome Workshop)

This tutorial covers building pangenome trees with PanMAN and placing reads or samples onto them using Panmap.

[1. Environment setup](#environment-setup)

[2. Build a PanMAN from raw sequences](#step-1-build-a-panman)

[3. Use Panmap to place a clonal sample](#place-a-single-haplotype-sample)

[4. Use Panmap to place a metagenomic sample](#place-a-metagenomic-sample)

## Environment setup

All programs used in this tutorial can be installed using conda: [Panmap](https://github.com/amkram/panmap), 
[PanMAN](https://github.com/TurakhiaLab/panman), [DIPPER](https://github.com/TurakhiaLab/DIPPER), 
[TWILIGHT](https://github.com/TurakhiaLab/TWILIGHT), and [pangolin](https://github.com/cov-lineages/pangolin). Docker is 
also supported for Panmap, PanMAN, DIPPER, and TWILIGHT; instructions can be found in each repository.

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
  pangolin \
  -y
```

## Step 1: Build a panman

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

**Build a PanMAN from alignment and guide tree.** Note that PanMANs built from MSA do not contain information on
structural variations. Currently, PanMAN doesn't infer SVs internally and uses SV information from Pangraph or GFA.

```bash
panmanUtils -M example_output/sars.subsampled_1000.guide.aln \
  -N example_output/sars.subsampled_1000.guide.nwk \
  -o sars.subsampled_1000
```

The summary stats show that the tree contains 1969 nodes and ~26,000 combined mutations (~13 mutations / branch), which 
is within the expected range for 1000 randomly sampled SARS-CoV-2 genomes. This is a small example for tutorial 
purposes; the rest of the demo uses larger, prebuilt PanMANs.

<details>
<summary>Expand to view summary stats</summary>

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

## Step 2: Run panmap

**Panmap has two modes: single sample mode and metagenomic mode.** In single sample mode, Panmap assumes the sample
contains a single haplotype and places the sample as a single unit on the tree. In metagenomic mode, Panmap assumes the
sample contains a mix of haplotypes and places each read separately on the tree.

### Place a single-haplotype sample

**Build a Panmap index.** Panmap sketches [syncmer](https://peerj.com/articles/10805/) seeds from all genomes in the 
PanMAN and stores them compactly in a phylogeny-guided format. This step should finish in under 5 seconds.

```bash
panmap input_data/sars_20000_twilight_dipper.panman \
  --index input_data/sars_20000_twilight_dipper.panman.idx \
  -k 19 -s 8 -l 3
```

**Place and genotype a sample.** Panmap scores k-mer similarity between the sample and all the genomes in the PanMAN, 
and picks the most similar genome as the reference genome for downstream alignment and genotyping. Efficiency is 
achieved by updating the similarity scores incrementally from seed changes along the tree rather than recomputing the 
similarity scores for each genome.

```bash
panmap -i input_data/sars_20000_twilight_dipper.panman.idx \
  input_data/sars_20000_twilight_dipper.panman \
  input_data/isolate_R1.fastq.gz input_data/isolate_R2.fastq.gz \
  --stop consensus \
  -o example_output/isolate
```

Output files:
  - `isolate.ref.fa`: The top scoring reference genome used for alignment and genotyping.
  - `isolate.bam`, `isolate.mpileup`, `isolate.vcf`: Alignment, pileup, and genotyping results.
  - `isolate.consensus.fa`: The consensus genome generated from the alignment.
  - `isolate.placement.tsv`: Top scoring node(s) for the sample in different metrics.

<details>
<summary>Expand to view the placement tsv </summary>

```console
$ column -t example_output/isolate.placement.tsv 
metric                score      nodes
log_raw               74.335591  node_7618,node_7619,node_7620,node_7621,node_7622,node_7623,node_7624,node_7625,node_7626,node_7627,node_7628,node_7629,node_7630,node_7631
log_cosine            0.788133   node_7618,node_7619,node_7620,node_7621,node_7622,node_7623,node_7624,node_7625,node_7626,node_7627,node_7628,node_7629,node_7630,node_7631
containment           0.196785   node_7617,node_7618,node_7619,node_7620,node_7621,node_7622,node_7623,node_7624,node_7625,node_7626,node_7627,node_7628,node_7629,node_7630,node_7631,node_18266,node_18267
weighted_containment  1.067023   node_7617,node_7618,node_7619,node_7620,node_7621,node_7622,node_7623,node_7624,node_7625,node_7626,node_7627,node_7628,node_7629,node_7630,node_7631,node_18266,node_18267
log_containment       0.454889   node_7618,node_7619,node_7620,node_7621,node_7622,node_7623,node_7624,node_7625,node_7626,node_7627,node_7628,node_7629,node_7630,node_7631
```

</details>

### Place a metagenomic sample

**Build a Panmap metagenomic index.** Metagenomic mode uses a different index file, since it requires additional
information (seed directions, positions, etc.) for per-read placement.

```bash
panmap input_data/sars_20000_twilight_dipper.panman \
  --index-mgsr input_data/sars_20000_twilight_dipper.panman.meta.idx 
```

**Haplotype deconvolution.** Panmap can estimate the haplotype composition of a mixed-strain sample. The example here 
uses a SARS-CoV-2 wastewater amplicon sample collected at Point Loma on 1/9/2022 by [Karthikeyan et al.](https://www.nature.com/articles/s41586-022-05049-6).The sample has already been 
preprocessed (trimming, amplicon stack information, etc.). On how to preprocess wastewater samples, please refer to 
[Panmap documentation](https://amkram.github.io/panmap/metagenomic.html).

```bash
panmap input_data/sars_20000_twilight_dipper.panman \
  input_data/SRR19707934.trimmed.fastq \
  --meta \
  --index input_data/sars_20000_twilight_dipper.panman.meta.idx \
  --amplicon-depth input_data/SRR19707934.amplicon_stacks.tsv \
  --mask-reads-relative-frequency 0.01 \
  --em-delta-threshold 0.00001 \
  --output example_output/SRR19707934 
```

This outputs a .mgsr.abundance.out file containing the haplotype abundance for each sample. To get the lineage
proportions, we can post-process the file using a custom script and Pangolin.

```bash
bash scripts/get_lineages.sh \
  --abundance example_output/SRR19707934.mgsr.abundance.out \
  --panman input_data/sars_20000_twilight_dipper.panman \
  --output example_output/SRR19707934.lineages.tsv
```

The lineage abundance file shows that the sample is predominantly composed of Omicron lineages, with BA.1.1 as the major 
component and smaller fractions of BA.1.15, BA.1, and BA.1.18. This is consistent with the clinical submissions around 
the same time.

<details>
<summary>Expand to view the haplotype and lineage abundance file </summary>

Haplotype abundance:

```console
$ column -t example_output/SRR19707934.mgsr.abundance.out  | head
node_7758,node_7759                                  0.22605
node_7747                                            0.15243
node_7928                                            0.14815
node_7874,node_7732,node_7937,node_7938,node_7939    0.12265
node_7760,node_7761,node_7941,node_7942,node_7943    0.12207
USA/MI-CDC-QDX33539980/2022|OM776151.1|2022-01-31    0.04937
USA/VA-CDC-ASC210751051/2022|OM911803.1|2022-02-22   0.03502
node_8247,node_8248,node_8249,node_8250,node_8251    0.02634
node_8020                                            0.02451
USA/VI-CDC-2-5383947/2021|OM338536.1|2021-12-28      0.01887
```

Lineage abundance:

```console
$ column -t example_output/SRR19707934.lineages.tsv  
BA.1.1     0.50055
BA.1.15    0.18317
BA.1       0.13372
BA.1.18    0.050850000000000006
BA.1.20    0.04937
BA.1.17.2  0.01887
BA.1.1.4   0.01848
B.1        0.01499
BC.1       0.0136
BA.1.1.15  0.01114
BA.1.13    0.00524
```

</details>

<br>

**Reads assignment/placement.** Panmap can also place reads onto the tree using a competitive mapping like approach,
with specific applications in eDNA assignment. Here we use a prebuilt vertebrate mitochondrial genome PanMAN containing
~15k genomes, across ~8k species. The read sample contains 1,000,000 subsampled reads from permafrost and lake
sediment samples collected across the Arctic by [Wang et al. 2021](https://www.nature.com/articles/s41586-021-04016-x).
The original study reported strong megafauna signals,  particularly from mammoths.

```bash
# Build a Panmap index
panmap input_data/v_mtdna.panman \
  --index-mgsr input_data/v_mtdna.panman.meta.idx \
  -k 15 -s 8 -l 1

# Run with filter-and-assign setting
panmap input_data/v_mtdna.panman \
  input_data/subsampled.fastq.gz \
  --meta \
  --index input_data/v_mtdna.panman.meta.idx \
  --filter-and-assign \
  --discard 0.6 --dust 5 \
  --taxonomic-metadata input_data/v_mtdna.meta.tsv \
  --output example_output/subsampled
```

As expected, a distinct cluster of reads is assigned to the Elephantidae node, consistent with the original study.

<details>
<summary>Expand to view read assignment results</summary>

```console
$ sort -k3,3 -gr example_output/subsampled.mgsr.assignedReads.out  | column -t
node_14404,NC_007596.2,DQ188829.2             Elephantidae    21  1,2,3,4,6,15,18,25,26,27,28,29,30,31,32,33,34,35,36,38,39
node_14402,JF912199.1,node_14403,NC_015529.1  Elephantidae    18  1,2,3,4,15,18,25,26,27,28,29,30,31,32,33,34,35,36
node_14400                                    Elephantidae    6   1,2,3,4,6,15
node_14399                                    Elephantidae    6   1,2,3,4,6,15
node_14401,NC_005129.2,DQ316068.1             Elephantidae    5   1,3,4,6,38
node_14409,KY499555.1,NC_035230.1             Elephantidae    4   2,3,4,18
node_14407                                    Elephantidae    4   2,3,4,18
node_14406,AJ224821.1,NC_000934.1             Elephantidae    4   1,2,3,4
node_14405                                    Elephantidae    4   1,2,3,4
node_14408,NC_020759.1,JN673264.1             Elephantidae    3   2,3,4
node_14398                                    Elephantidae    3   1,2,6
node_14410,KY364233.1,NC_035800.1             Elephantidae    2   1,2
node_9733                                     Kyphosidae      1   21
node_94,MN977920.1,NC_050664.1                Scincidae       1   19
node_7761,NC_081957.1,OR066216.1              Gobiidae        1   11
node_7521,MG321595.1,NC_037248.1              Stomiidae       1   20
node_7142,NC_068841.1,ON005612.1              Myctophidae     1   14
node_6963                                     Gonorynchidae   1   16
node_6633,NC_082182.1,OM960959.1              Pegasidae       1   7
node_5832,KP013104.1,NC_028291.1              Leuciscidae     1   5
node_4765,OQ603602.1,NC_077583.1              Nemacheilidae   1   9
node_4735                                     Nemacheilidae   1   10
node_4314,NC_031600.1,AP011347.1              Botiidae        1   8
node_4313                                     Botiidae        1   8
node_4312                                     Botiidae        1   8
node_4311                                     Botiidae        1   8
node_356,NC_012443.1,EF222190.1               Chamaeleonidae  1   13
node_3252,NC_039135.1,AP018342.1              Notacanthidae   1   12
node_14456,NC_028563.1,KT818536.1             Chlamyphoridae  1   0
node_12279,NC_020154.1,AY954504.1             Herpelidae      1   37
node_12060                                    Plethodontidae  1   17
node_12031,MN259079.1,NC_044873.1             Pipidae         1   24
node_11280                                    Dicroglossidae  1   23
node_11279                                    Dicroglossidae  1   23
node_11144                                    Megophryidae    1   22
```


</details>

<br>
