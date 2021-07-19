
## TL;DR
For evolutionary studies (especially phylogenomics) using target capture (aka Hyb-seq): Use avaliable WGS raw reads or a reference genome (draft genomes should work) to:
 1. identify putative **paralogs** and **missing genes** (both `VetTargets_WGS.R` and `VetTargets_genome.R`).
 2. identify genes with **huge intron(s)** (which are often present e.g. in many mammals [[(1)]] and should definitely NOT be assumed to be in linkage equilibrium). (Only `VetTargets_genome.R`).
 3. Nicely visualise the **genomic context** of your targets relative to your study group.
 4. Extract supercontigs to use as targets for phylogenetics in closely related species (or population genomics). (`TargetSupercontigs.R`.)

Currently the scripts run on the command line via Rscript. Dependencies:
```
R (https://cran.r-project.org/)
R packages:
    tidyverse
    ggrepel
    optparse
    Biostrings
    segmented
```

## `VetTargets_WGS.R` and `VetTargets_genome.R`

So you're doing a very important phylogenomics and you've designed or chosen your target set for your favorite organism. What next? Well, you've also got some genomic shotgun sequence data and/or a reference genome for one or more closely related species. Rather than moving forward with your sequencing project and hoping for decent recovery without too many nasty paralogs or other surprises, it's possible to use that WGS data to make life for your future self a little easier while also improving the efficiency of your project.

TargetVet provides methods for using WGS data and/or a reference genome to vet (i.e. filter) your target set before bait design. 

In the event that you're already with target capture data in hand, you can still use TargetVet to help filter your assembled genes (e.g. to identify paralogs and guide their separate assembly).

## Usage
#### `VetTargets_genome.R`
This takes a blast search of your targets to your genome and outputs the following:
1. `<output_prefix>.IntronFlags.txt`: Per-target summary of introns/intergenic regions with logical flags for exceeding specified max intron length and percent of supercontig made up of introns, as well as the actual values.
2. `<output_prefix>.IntronStats.txt`: Position and length of each intron on each target.
3. `<output_prefix>.CoverageStats_AcrossChromosomes.txt`: Per-target summary of coverage across 'chromosomes' in the reference genome. Fields are `qseqid mean_paralogy_rate paralog_percent full_percent unique_percent missing_percent n_chroms`. This is informative about **paralogs**. You will be most interested in `mean_paralogy_rate`, which is the average number of segments from the genome covering each position in the target. Single-copy loci will have a value around 1, duplicates around 2, triplicates around 3, etc. This value considers the full target length, so the minumim value is 0. To check for **missing** targets, look at `missing_percent`, which is the proportion of the target without any matches (i.e. with coverage 0). Remember that if you set `--min_fragment_length` too high, `missing_percent` will be inflated. Keeping it at zero should be fine.
4. `<output_prefix>.CoverageStats_PerChromosome.txt`: this has the same fields as the across-chromosomes summary, but the values are reported for every target-chromosome matched pair. If the target only matches to one chromosome, it will have a single row with values identical to those in the across-chromosomes summary.

```
# to get a description of the available arguments (several are required)
Rscript VetTargets_genome.R -h

# example run
Rscript VetTargets_genome.R --blast_file blastn_yourtargets_to_yourgenome.withHeader.txt --min_fragment_length 0 --min_pident 70 --max_intron_length 1000 --max_intron_percent 50 --output_prefix targets2genome --min_display_intron 300
```
**NB:** The blast file must contain the following (tab-separated) fields, and have them as its header:
```
qseqid	sseqid	pident	length	mismatch	gapopen	qlen	qstart	qend	slen	sstart	send	evalue	bitscore
```
To get usable blast results you need to tell blast to return these fields. Here's a blastn example:
```
makeblastdb -in yourgenome.fasta -dbtype nucl -out yourgenome_db

blastn -query yourtargets.fasta -db yourgenome_db -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore" -evalue 1e-6 -out blastn_yourtargets_to_yourgenome.txt  
```
Then just paste in the header line above.

## `TargetSupercontigs.R` (forthcoming)
Extract supercontigs from (draft) genome assemblies for more effective target capture sequencing in **closely related species**.

### Rationale
Phylogenomic projects in non-model organisms typically use exons as their targets, often even if their study group consists of very recently diverged, and therefore phylogenetically difficult, taxa. This results in significantly reduced coverage and poor assembly of intronic or flanking regions, which are consistently shown to improve phylogenetic inference due to their faster mutation rates. 

At the same time, whole-genome shotgun sequencing is now relatively cheap, making draft genome assembly possible for us poor phylogeneticists working on non-model organisms. 

`TargetSupercontigs.R` aims to locate target exons in the draft genome and extract the full exon + intron ('supercontig'). These can then be used as targets instead of (or in combination with) the original exons, provided all the species you're working with are relatively closely related (at deep time scales introns will become too divergent to be useful as targets; exactly how deep the time scale limit is is very hard to say, but the probes can apparently tolerate 80% or even just 70% similarity between the probe and true sequences).

## References
<a id="1">(1)</a> Gatesy, J. & Springer, M. S. (2014). Phylogenetic analysis at deep timescales: Unreliable gene trees, bypassed hidden support, and the coalescence/concatalescence conundrum, _Molecular Phylogenetics and Evolution_, 80, 231-266.

##### TODO:
Make script to group exons from the same gene into separate target loci when introns are super big. (?)
##### DONE 
Make sure to flag as paralogs regions on the same chromosome matching the same target.
