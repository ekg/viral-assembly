# viral assembly techniques from ONT reads

Assembly methods often aim to produce a single canonical representation of every genomic locus in their input.
This objective can sometimes frustrate analyses that depend on information about variation represented in a given biosample.
Assemblers are also often monolithic systems, where each step in the assembly process is implemented in a single program.
This reduces configurability of the method.
It may be difficult to understand all asumptions of the method, or impossible to modify certain hyperparameters without modifying the software itself.

## approach

In this repository, I explore a modular assembly method for viral genomes sequenced using long read sequencing from Oxford Nanopore (ONT).
The workflow consists of a number of independent steps:

1. select longer reads in the data
2. all-vs-all mapping
3. induction of graph with seqwish
4. graph simplification
5. read mapping to graph, coverage calculation
6. inference of viral quasispecies

Intermediate data is either in the form of sequence graphs, alignments between sequences and/or graphs, and coverage information over graphs.

## example workflow

To clarify these steps, let's work through a step by step example of such a workflow.
Here I'm using [data from a nasal swab of a SARS-CoV-2 patient](https://openresearch.labkey.com/wiki/ZEST/Ncov/page.view?name=SARS-CoV-2%20Deep%20Sequencing).

I'll use the swab Oxford Nanopore data.
We can get the data with the sra toolkit (installed in guix via `guix package -i sra-tools`).

```bash
fastq-dump https://sra-download.ncbi.nlm.nih.gov/traces/sra52/SRR/010879/SRR11140751
# symlink for convenience
ln -s SRR11140751.fastq swab_ONT.fastq
```

### (1) we filter for read length, keeping only long reads (most are very short, due the reverse transcription / amplification genesis of the data):

```bash
cat swab_ONT.fastq \
    | paste - - - - \
    | tr ' ' '_' \
    | awk 'length($2) > 2000 { print $1; print $2; print $3; print $4; }' \
    | pigz >swab_ONT.2kbp.fastq.gz
```

### (2) all versus all mapping, using all kmers to drive the mapping, and 48 mapping threads

```bash
minimap2 -c -w 1 -X -t 48 swab_ONT.2kbp.fastq.gz swab_ONT.2kbp.fastq.gz >swab_ONT.2kbp.paf
```

#### (2.1) optionally filter the alignments for read length using fpa (installable via `cargo install fpa_lr`):

```bash
fpa drop -l 1000 <swab_ONT.2kbp.paf >swab_ONT.2kbp.drop1kbp.paf
```

This will result in a more open and fragmentary graph with short tips where reads don't align full length against each other.
Different settings of length might work better given different data and applications.

### (3) induction of graph with seqwish

```bash
seqwish -t 48 -k 16 -s swab_ONT.2kbp.fastq.gz -p swab_ONT.2kbp.paf -g swab_ONT.2kbp.gfa
```

The resulting graph can be visualized with [Bandage](https://github.com/rrwick/Bandage), yielding something like this:

![Bandage rendering of swab_ONT.2kbp.gfa](https://github.com/ekg/viral-assembly/raw/master/PRJNA607948/swab_ONT.2kbp.gfa.Bandage.png)

### (4) graph simplification

```bash
odgi build -g swab_ONT.2kbp.gfa -o - \  # read in the GFA file
    | odgi prune -i - -b 3 -o - \       # run a pruning operation
    | odgi view -i - -g >swab_ONT.2kbp.odgi-prune.b3.gfa  # write GFA
```

This uses [odgi](https://github.com/vgteam/odgi) to remove edges from the graph that are not among the best three edges for any node.
It is crude, but brings out a linear structure that is expected of the SARS-CoV-2 genome.
The resulting graph still has a lot of sequence in it, presumably derived from local errors, low-quality reads, and chimeric PCR products (RT-PCR was used to generate cDNA from the input sample).

![Bandage rendering of swab_ONT.2kbp.odgi-prune.b3.gfa](https://github.com/ekg/viral-assembly/raw/master/PRJNA607948/swab_ONT.2kbp.odgi-prune.b3.gfa.Bandage.png)

### (5) read mapping to graph, coverage calculation

I haven't yet implemented this in my testing, but I would use [GraphAligner](https://github.com/maickrau/GraphAligner) or [gyeet](https://github.com/ekg/gyeet).
[vg](https://github.com/vgteam/vg) has three mappers which might work (`map`, `mpmap`, and `giraffe`).
While `vg map` can return high quality alignments, it is sometimes slow on graphs with complicated looping structures.

Graphs with complex local loops are difficult for any mapper that treats them as an actual graph (all those mentioned except for `gyeet`) because there is a risk of alignments getting stuck in small graph motifs which are "universal" models of any DNA sequence.
These can arise during the seqwish graph induction step.
Removing them may be possible with `odgi break`, but I have not been able to get this to work reliably enough to prevent GraphAligner from getting stuck sometimes.

### (6) inference of viral quasispecies

With the coverage of reads mapped back to the variation graph induced from long reads, we should be able to make a guess as to the viral quasispecies that are present.
We could apply something like [virus-vg](https://bitbucket.org/jbaaijens/virus-vg) or [vg-flow](https://bitbucket.org/jbaaijens/vg-flow) (probably the latter, as it scales better) to the variation graph.
These are designed for illumina data, but by constructing a cleaned variation graph from the ONT reads, we can effectively reduce their error rate and allow the application of these methods.
Other, cruder methods may also be possible to apply.
