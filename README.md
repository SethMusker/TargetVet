# TL;DR
TargetVet is an R-based toolkit for evolutionary studies (especially phylogenomics) using target capture (aka Hyb-seq) which allows you to use avaliable WGS raw reads, a reference genome, or assembled contigs from target capture sequence data (such as those made by **HybPiper** - https://github.com/mossmatters/HybPiper), to:
 1. Identify putative **paralogs** and **missing genes** using:

    a. WGS raw reads from one or more related species: `VetTargets_WGS.R`,

    b. Genome assemblies (draft or reference): `VetTargets_genome.R`, or

    c. Target Capture data assembled using HybPiper: `VetHybPiper.sh`.
 2. Nicely visualise the **genomic context** of your targets relative to your study group.
 3. Identify genes with **huge intron(s)** (which are often present e.g. in many mammals [[1]] and should not be assumed to be in linkage equilibrium): Only `VetTargets_genome.R`.
 4. Extract supercontigs to use as targets for phylogenetics in closely related species (or population genomics): `TargetSupercontigs.R`.

# Dependencies
The R scripts run on the command line via `Rscript` (see examples below) and are therefore platform-independent. The functions `VetHybPiper.sh`,`map_REF_to_targets.sh` and `map_WGS_to_targets.sh` use the `bash` shell programming language, the standard for linux.

Dependencies:
```
#-- R (https://cran.r-project.org/)
#-- R packages:
    optparse # for parsing opts
    tidyr, dplyr # main workhorses in all scripts
    forcats # factor releveling
    ggplot2 # does most of the plotting
    ggrepel # useful for labeling plots
    progress # for nice progress bars 

# for VetHybPiper.sh and DetectParalogs.R
    nplr and cumSeg # for automated paralog detection via breakpoint regression
    gplots # for paralogy and missingness heatmaps 
    dendextend # for plotting dendrograms and/or phylogenies
    ape # for phylogenetic trees

# for TargetSupercontigs.R
    Biostrings # for handling sequences
    Rsamtools # for indexing and extracting sequences from reference genome(s)

#  To install these, run:
    install.packages(c("tidyr","dplyr","forcats","ggplot2","ggrepel","optparse","progress","nplr","cumSeg","gplots","ape","dendextend"))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
    BiocManager::install("Biostrings")
    BiocManager::install("Rsamtools")

#-- Optional:
ncbi-BLAST+ (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) needs to be in your $PATH if you plan to run VetHybPiper.sh or map_REF_to_targets.sh
BBMap suite (https://sourceforge.net/projects/bbmap/) # for VetHybPiper.sh, optional but recommended - see below

#-- NOTE
All of the wrapper scripts (*.sh) require linux OS.
```

# Identify paralogs from HybPiper: **`VetHybPiper.sh`**
This script is designed as an easy-to-use wrapper around several TargetVet scripts with the aim of automated identification of paralogs in HybPiper-assembled target capture data. It also produces detailed summary tables (so you can apply your own filtering and paralog-detection methods), as well as nice figures that help you visualise patterns of paralogy and missingness in your assemblies, across all genes and samples.

The script runs using Bash, the standard shell-scripting language on Linux systems. You can see the available options by running `bash VetHybPiper.sh -h`. The 6 minimum required arguments are all things (mostly files) that you'll already have in your HybPiper directory, so you won't have to do anything new before running this. 

You will need to have the BLAST+ bin in your path.

Apart from the above dependencies, there is the optional added dependency of `dedupe.sh` from the BBMap suite. I have found that HybPiper contigs often contain two (or more) **near-identical copies of the same sequence for many genes, leading to anomalously high paralogy scores**. This is possibly because HybPiper is hard-coded to use the 'careful' mode of the SPAdes assembler and therefore assembles each haplotype of the same region separately, rather than collapsing them into a 'consensus' contig. Using `dedupe.sh` by setting `-d TRUE` with the default minimum percent identity threshold of 97 (which can be changed, e.g. `-m 98`) considerably reduces noise in the data sets I've tested and improves paralog detection accuracy.

By default `blastn` is used, but if your targets are quite divergent from your samples you may wish to use `tblastx` for more sensitive querying (set `-B tblastx`). However, note that because of the evaluation of all possible reading frames, this option produces many more redundant hits (i.e. within-contig matches to essentially the same part of the target) than `blastn`. I have implemented an effective 'thinning' algorithm in `VetTargets_genome.R` for removing these redundancies, but it can add a bit of time to the overall analysis (around 1 to 2 minutes per sample).

## What does `VetHybPiper.sh` do?
1. For each sample, fetch all the assembled contigs from the HybPiper output directory and collate them into a single all-target fasta.
2. For each sample, nucleotide BLAST the contigs to the target fasta (the one you used as the reference for HybPiper).
3. For each sample, run `VetTargets_genome.R` to calculate paralogy indices and more (see below). By default the script's plotting functionality is turned off as it can take a long time and is, in this context, not likely to be interesting. It can however be switched on using `-P TRUE`.
4. Across the sample set, run `DetectParalogs.R` on the previous step's output. If the reference fasta contains multiple copies of each target (in which case you need to specify `-M TRUE`), this will be done separately for each target source, with results output in separate folders named after each.

To get detailed usage info and help, run `bash VetHybPiper.sh -h`.

## How does `DetectParalogs.R` detect paralogs?
In short, it looks for shared patterns of paralogy across samples. 

Automated identification of paralogs is attempted via breakpoint regression on the raw paralogy percent estimates, with the independent variable being the order in which targets appear in the mean-paralogy-percent-sorted sequence. Hypothetically, the breakpoint should occur where the paralogy rate suddenly increases, or becomes consistently high. This is obviously quite crude, which is why `DetectParalogs.R` also makes lots of informative plots (e.g. Figure 1) to help you visualise the data and make informed decisions about which genes to call paralogs and exclude (or phase into separate loci if you can't afford or don't want to drop them). 

In cases where you have targets designed based on outgroup taxa and several outgroup taxa in the dataset, the patterns of paralogy they show may be quite different from those in your ingroup taxa. As such, you can provide a list of the sample names of the ingroup taxa, in which case another, separate breakpoint analysis will be conducted only on them. This can help a lot to reduce noise in the data and identify paralogs more accurately.

![g84089](https://user-images.githubusercontent.com/16952350/133785807-5f62faa7-a47d-428a-b937-e7a1d755afad.png)
**Figure 1**: A heatmap of paralogy estimates showing how `DetectParalogs.R` groups genes based on their shared paralogy patterns, and how paralogy estimates vary with phylogeny (at left is the species tree estimated by Johnson et al. (2016)).

# Designing a project: **`VetTargets_genome.R`** and **`VetTargets_WGS.R`**


You're planning a target capture phylogenomics project and you've designed or chosen your target set. What next? Well, if you've got genomic shotgun sequence data and/or a reference genome for one or more closely related species, then rather than moving forward with your sequencing project and hoping for decent recovery without too many paralogs or other nasty surprises, it's possible to use that WGS data to make life for your future self a little easier, while also improving the efficiency of your project.

TargetVet provides methods for using WGS data and/or a reference genome to vet (i.e. filter) your target set before bait design. 

In the event that you're already with target capture data in hand, you can still use TargetVet to help filter your assembled genes (e.g. to identify paralogs and guide their separate assembly). See the previous section.

## Usage
### `VetTargets_genome.R`
```
# to get a description of the available arguments (several are required)
Rscript VetTargets_genome.R -h

# example run
Rscript VetTargets_genome.R --blast_file blastn_yourtargets_to_yourgenome.withHeader.txt --min_fragment_length 0 --min_pident 70 --max_intron_length 1000 --max_intron_percent 50 --output_prefix targets2genome --min_display_intron 300
```

### Getting the blast result for input
You can either do it yourself (see below) or use `map_REF_to_targets.sh`. Run `bash map_REF_to_targets.sh -h` to get usage info.

#### DIY blast run
**NB:** The blast file must contain the following (tab-separated) fields, and have them as its first line:
```
qseqid	sseqid	pident	length	mismatch	gapopen	qlen	qstart	qend	slen	sstart	send	evalue	bitscore
```
To get usable blast results you need to tell blast to return these fields. Here's a blastn example:
```
makeblastdb -in yourgenome.fasta -dbtype nucl -out yourgenome_db

blastn -query yourtargets.fasta -db yourgenome_db -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore" -evalue 1e-6 -out blastn_yourtargets_to_yourgenome.txt  
```
Then just paste in the header line above.

### Running the script will generate the following files:

1. `<output_prefix>.IntronFlags.txt`: Per-target summary of introns/intergenic regions with logical flags for exceeding specified max intron length and percent of supercontig made up of introns, as well as the actual values.
   
2. `<output_prefix>.IntronStats.txt`: Position and length of each intron on each target.
   
3. `<output_prefix>.CoverageStats_AcrossChromosomes.txt`: Per-target summary of coverage across 'chromosomes' in the reference genome. Fields are `qseqid mean_paralogy_rate paralog_percent full_percent unique_percent missing_percent n_chroms`. This is informative about **paralogs**. You will be most interested in `mean_paralogy_rate`, which is the average number of segments from the genome covering each position in the target. Single-copy loci will have a value around 1, duplicates around 2, triplicates around 3, etc. This value considers the full target length, so the minumim value is 0. To check for **missing** targets, look at `missing_percent`, which is the proportion of the target without any matches (i.e. with coverage 0). Remember that if you set `--min_fragment_length` too high, `missing_percent` will be inflated. Keeping it at zero should be fine.
   
4. `<output_prefix>.CoverageStats_PerChromosome.txt`: this has the same fields as the across-chromosomes summary, but the values are reported for every target-chromosome matched pair. If the target only matches to one chromosome, it will have a single row with values identical to those in the across-chromosomes summary.

# **`TargetSupercontigs.R`**
Extract supercontigs from (draft) genome assemblies for more effective target capture sequencing in **closely related species**.

## Rationale
Phylogenomic projects in non-model organisms typically use exons as their targets, often even if their study group consists of very recently diverged, and therefore phylogenetically difficult, taxa. This results in significantly reduced coverage and poor assembly of intronic or flanking regions, which are consistently shown to improve phylogenetic inference due to their faster mutation rates.

At the same time, whole-genome shotgun sequencing is now relatively cheap, making draft genome assembly possible for us poor phylogeneticists working on non-model organisms.

`TargetSupercontigs.R` aims to locate target genes in the draft genome and extract the full exon + intron ('supercontig'). These can then be used as targets instead of (or in combination with) the original targets(if they're from a different species), provided all the species you're working with are relatively closely related. Note that at deep time scales introns will become too divergent to be useful as targets; exactly how deep the time scale limit is is very hard to say, but the probes can apparently tolerate 80% or even just 70% similarity between the probe and true sequences.

## Usage
```
# to get usage and help
Rscript TargetSupercontigs.R -h
# basic run with default parameters and blacklist from VetTargets_WGS.R
Rscript TargetSupercontigs.R -b blastn_yourtargets_to_yourgenome.withHeader.txt -g yourgenome.fasta -o test1 -x test1.removed 
```


# References
<a id="1">(1)</a> Gatesy, J. & Springer, M. S. (2014). Phylogenetic analysis at deep timescales: Unreliable gene trees, bypassed hidden support, and the coalescence/concatalescence conundrum, _Molecular Phylogenetics and Evolution_, 80, 231-266.

## Notes to self
##### TODO:
Make script to group exons from the same gene into separate target loci when introns are super big. (?)
