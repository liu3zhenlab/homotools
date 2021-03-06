# homotools

### NOTE: This is NOT a finished project. In particular, **homograph** needs to be improved. If you are interested in contributing to this project, please contact [Sanzhen Liu](liu3zhen@ksu.edu).
## What can homotools do?

The package is dedicated to structurally compare sequences of homologous genes from different genomes of a species. Briefly, the package has been developed to perform the following analyses.

1. The script **geneseq** can extract gene sequence, gene annotation, other featured data in GTF/GFF and BED format, transcribed sequences, protein sequences, and gene sequences with transcribed sequences highlighted. The positions of gene annotation and other features are adjusted to the coordinates relative to the newly extracted sequences.
2. The script **homocomp** implements homologous searching and visualization of the comparison. With a query gene sequence, homologous sequences can be identified in a genome. The homologous gene could be found if a related GTF is supplied. The alignment between two homologous regions are plotted.
3. The script **homosnpeff** finds and annotates variants between query and reference sequences.
4. The script **homograph** visualizes alignments of homologous sequences, determines number of haplotypes, and views the haplotype graph. This module is still under development.
5. The script **homostack** visualizes alignments among multiple homologous sequences.

## Motivation

The development of **homotools** was motivated by the tedious procedure for extracting sequences and related information of a gene from a reference genome and homologous genes from other genomes. Unix shell scripts, Perl, and R were combined for the development.

## Requirements and installation

Shell scripts, Perl, and R were combined for the development. Both Perl and R are generally installed. If needed, please refer to [Perl](https://www.perl.org/) and [R](https://www.r-project.org/) for installation guides. [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), [MURmer](http://mummer.sourceforge.net/), [BEDtools](https://bedtools.readthedocs.io/en/latest/) are also required.

To run **homograph**, [cd-hit](http://weizhong-lab.ucsd.edu/cd-hit/) and [pggb](https://github.com/pangenome/pggb) are required.

Aternatively, all required packages can be installed through [conda](https://docs.conda.io/en/latest/)

```
conda create -n homotools
conda activate homotools
conda install -c anaconda perl
conda install -c r r-base
conda install -c bioconda blast mummer bedtools
conda install -c bioconda cd-hit pggb
```

#### download homotools and run scripts

Basically, the installation of homotools is just to copy all homotools files in a directory. Scripts can be directly used afterwards.

```
git clone https://github.com/liu3zhenlab/homotools.git 
perl ./homotools/homocomp
# test run
cd data
sh homocomp_testrun.sh
```
The output directory "out" will be produced with many PDF plots and other outputs.  

## geneseq

#### Usage

```
Usage: geneseq --fas <fasta> --gene <genename> --gtf <GTF> [options]
- extract fasta sequence and gtf information for the input gene
- produce a new gtf with adjusted positions relative to gene fasta
- also can extract other features, cds, cdna, protein sequences of the genes
[Options]
  --fas|f <file>      : reference fasta file (required)
  --gene|i <str>      : gene name (required)
  --gtf|g <file>      : GTF file (required)
  --prefix|p <str>    : prefix name for outputs (geneseq)
  --ext5|e <num>      : bp extension from 5' site (0)
  --ext3|x <num>      : bp extension from 3' site (0)
  --othergff|o <file> : other GTF/GFF files containing other features (optional)
                        Note: multiple GTF/GFFs can be input by reusing --othergff
  --cds|c <file>      : a fasta file containing coding sequences of the gene (optional)
  --cdna|d <file>     : a fasta file containing transcript sequences of the gene (optional)
  --protein|a <file>  : a fasta file containing protin sequences of the gene (optional)
  --help|h            : help information
```

#### Output from geneseq

```
FASTA of gene sequence         : <prefix>.1.<gene>.fasta
Length of gene sequence        : <prefix>.2.<gene>.length
Annotation BED                 : <prefix>.3.<gene>.original.bed
Annotation GTF                 : <prefix>.3.<gene>.original.gtf
Pos-adusted annotations        : <prefix>.4.pos.adjusted.gtf.bed (directory)
gene sequence w/annotation     : <prefix>.5.cdna.highlighted (directory)

(optional outputs)
GFF or GTF of other features   : <prefix>.6.other.gffs  (directory)
CDS, cDNA, protein sequences   : <prefix>.7.cds.cdna.prot (directory)
```

## homocomp

#### Usage

```
Usage: homocomp --query <fasta> --db <blast_db> [options]
[Options]
    --query <file>   fasta file containing a sequence as the query; required
    --qryadd <file>  bed file to highlight regions in query; optional 
    [format, 7 columns separated by Tab]:
    chr start(0-based) end(1-based) label height(0.01-0.1) strand(+/-) color(R compatible)
    --db <db_name>   blast database; required
    --dbacc <str>    accession name of reference database; optional
    --ref <file>     fasta file of the reference; optional
                     if not supplied, reference sequences will be extracted from --db
    --tchr <str>     targeted chromosome or contig; optional
    --tstart <num>   bp position for the region start; optional
                     if specified, the value will be used as the left end
    --tend <num>     bp position for the region end; optional
                     if specified, the value will be used as the right end
    --tgtf <file>    the gene annotation GTF file of the target; optional
    --strand <str>   plus (or +) or minus (or -) (plus)
    --evalue <str>   maximal E-value (1e-10)
    --identity <num> minimal percentage of identity from 0 to 100 (80)
    --match <num>    minimal bp match of an alignment (100)
    --coverage <num> minimal coverage of the query (0)
    --repeatdust     remove repetitive blastn alignments if specified
    --lext <num>     extended bp from the left side (0);
    --rext <num>     extended bp from the right side (0);
    --expand <num>   expand scale relative to the query length (10)
                     e.g., 1kb query, 1kb x 10 extened from both sides will be scanned for the hit range
    --prefix <str>   the output directory and the prefix for output files (hpout)
    --threads <num>  number of cpus (1)
    --noblastn       skip blastn if specified
    --bandcol <str>  a valid R color name (bisque3)
    --version        version information
    --help:          help information.
```

#### Output from homocomp

Three alignment plots and two dotplots between two sequences are output.

```
Alignment plot without alignment filtering   : <prefix>.o1.alnplot.pdf
Alignment plot with alignment filtering      : <prefix>.o2.filt.alnplot.pdf
Alignment plot with filtered unique alignment: <prefix>.o2.filt.uniq.alnplot.pdf
Dotplot without alignment filtering          : <prefix>.o3.dotplot.pdf
Dotplot with alignment filtering             : <prefix>.o4.filt.dotplot.pdf
```

## homograph

#### Usage

```
Usage: perl ./homograph --dir <dir_containing_fasta_files> [options]
[Options]
    --dir <path>      path to a collection of fasta files with suffix of .fa, .fas, or .fasta; required
    --ref <file>      fasta file of the reference sequence; required
    --mgpara <str>    minigraph parameters (-xggs -t1 -j 0.2 -l 100)
    --cdhitpara <str> parameters for cd-hit-est (-g 1 -s 0.8 -c 0.8 -r 0)
    --pggbpara <str>  parameters for -p and -s in pggb (-p 90 -s 1000)
    --gene_gtf <file> gtf file of gene structure using coordinates on --ref; optional
    --threads <num>   number of threads (1)
    --prefix <str>    prefix of outputs (hgout)
    --main_label <float>     label to add as the title of haplotype plotting figure (haplotype graph)
    --refgap_prop <float>    ratio of each gap between segments to the canvas x length (0.02)
    --nonrefgap_prop <float> ratio of each gap between segments to the canvas x length (0.05)
    --maxseg_prop <float>    proprotion of the maxium fragment that will not be shortened segments to the canvas x length (0.5)
    --version         version information
    --help:           help information
```

#### Major output from homocomp

```
Plot of alignments of multiple sequences ordered by cd-hit clustering: <prefix>.4.pggb.in.fasta.alignments.png
Layout plot based on GFA output                                      : <prefix>.4.pggb.in.fasta.gfa.layplot.png
```

## homostack

#### Usage

```
Usage: homostack --query <fasta> --db <blast_db> [options]
[Options]
    --seq <file>     fasta file containing a sequence as the query; required
                     multiple sequences can be input by using --seq multiple times
    --annot <file>   bed file to highlight regions in query; optional 
                     [format]:chr start(0-based) end(1-based) label color(R compatible) height(0.01-0.1)
                     if specified, the number --annot needs to match the number of --seq
    --identity <num> minimal percentage of identity from 0 to 100 (80)
    --match <num>    minimal bp match of an alignment (100)
    --prefix <str>   the output directory and the prefix for output files (hsout)
    --threads <num>  number of cpus (1)
    --bandcol <str>  a valid R color name (bisque3)
    --version        version information
    --help:          help information.
```

#### Output from homocomp

```
Plot of sequential alignments of multiple sequences : <prefix>.3.alnstack.pdf
```

## Bug report

Please report any bugs or suggestion on github or by email to Sanzhen Liu ([liu3zhen@ksu.edu](mailto:liu3zhen@ksu.edu)).

## License

The package homotools is distributed under MIT licence.

### 
