# homotools

## What can homotools do?

The package is dedicated to structurally compare sequences of homologous genes from different genomes of a species. Briefly, the package has been developed to perform the following analyses. Our [Wiki](https://github.com/liu3zhenlab/homotools/wiki) provides more detailed explanation.

1. **geneseq** can extract gene sequence, gene annotation, other featured data in GTF/GFF and BED format, transcribed sequences, protein sequences, and gene sequences with transcribed sequences highlighted. The positions of gene annotation and other features are adjusted to the coordinates relative to the newly extracted sequences.
2. **homocomp** implements homologous searching and visualization of the comparison. With a query gene sequence, homologous sequences can be identified in a genome. The homologous gene could be found if a related GTF is supplied. The alignment between two homologous regions are plotted.
3. **homomine** search the best homolog in a targeted genome, identify and annotate polymorphisms, including SNP, INDEL, and structural variation (SV).
4. **homostack** visualizes alignments among multiple homologous sequences.
5. **homograph** visualizes alignments of homologous sequences, determines number of haplotypes, and views the haplotype graph. This module is still under development.

## Flowchart for modules  
<img src="https://github.com/liu3zhenlab/homotools/blob/5009d08c3b78dd223db478154e3e88cc06caa156/flowcharts/homotools_flowchart_v2.png" alt="flowchart" width=550 />

Usages of these modules can refer to [Wiki](https://github.com/liu3zhenlab/homotools/wiki).

## Motivation

The development of **homotools** was motivated by the tedious procedure for extracting sequences and related information of a gene from a reference genome and homologous genes from other genomes. Unix shell scripts, Perl, and R were combined for the development.

## Requirements and installation

Shell scripts, Perl, and R were combined for the development. Both Perl and R are generally installed. If needed, please refer to [Perl](https://www.perl.org/) and [R](https://www.r-project.org/) for installation guides. [Java](https://www.java.com/en/download/) is required to run [SNPEff](https://pcingola.github.io/SnpEff/se_introduction/). [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), [MURmer4](http://mummer.sourceforge.net/), [BEDtools](https://bedtools.readthedocs.io/en/latest/) are also required.

To run **homomine**, the multiple sequence aligner [MAFFT](https://mafft.cbrc.jp/alignment/software) is required.

To run **homograph**, [cd-hit](http://weizhong-lab.ucsd.edu/cd-hit/) and one of the following aligners for multiple sequence alignment are required.  
[Clustal Omega](http://www.clustal.org/omega/)  
[MUSCLE](https://www.drive5.com/muscle/)  
[MAFFT](https://mafft.cbrc.jp/alignment/software)  

All required packages can be installed through [conda](https://docs.conda.io/en/latest/):
```
git clone https://github.com/liu3zhenlab/homotools.git
cd homotools
conda env create -f homotools.yml
# test run
cd data
sh homocomp_testrun.sh
```

*Alternatively*, separate packages can be installed separately with Conda.
```
conda create -n homotools
conda activate homotools
conda install -c anaconda perl openjdk
conda install -c r r-base r-knitr r-rmarkdown
conda install -c bioconda blast mummer4 bedtools cd-hit clustalo muscle mafft

# clone homotools to a local directory.
git clone https://github.com/liu3zhenlab/homotools.git 
perl ./homotools/homocomp
# test run
cd data
sh homocomp_testrun.sh
```

## Bug report

Please report any bugs or suggestion on github or by email to Sanzhen Liu ([liu3zhen@ksu.edu](mailto:liu3zhen@ksu.edu)).

## License

The package homotools is distributed under MIT licence.

## versions
0.3.2: updated homograph  
0.3.1: simplified homomini  
0.3.0: updated homograph outputs  
0.2.9: added the parameter of *colClasses="character"* to utils/homomine_report.Rmd to avoid recognizing nucleotide T as a logic value  
0.2.8: changed "echo -e" to Perl code to be more relialbe in different systems  
0.2.7: fixed a bug in homomine to properly display SV tables  
0.2.6: added homomine and updated homograph  
0.2.5: changed "delta-filter -1" to "delta-filter -r" in the filter for uniquely mapped alignments  
0.2.4: added the module homomine to compare homologs in two genomes and annotate polymorphisms  
0.2.3: allowed an automatic choice for strand (plus or minus) in homocomp  
