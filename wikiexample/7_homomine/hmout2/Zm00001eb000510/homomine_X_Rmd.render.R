library(rmarkdown)
library(knitr)

render('Zm00001eb000510/homomine_report.Rmd',
  params = list(
    cwd="/homes/liu3zhen/scripts2/homotools/wikiexample/7_homomine/hmout2",
    scriptdir="/homes/liu3zhen/scripts2/homotools/utils/",
    qrygene="Zm00001eb000510",
    qrybase="B73",
    qryseq="Zm00001eb000510/1_B73_Zm00001eb000510/1.Zm00001eb000510.fasta",
    qrybed="Zm00001eb000510/Zm00001eb000510_to_A188_output/data/B73.label.bed",
    tgtbase="A188",
    tgtgene="Zm00056a000060",
    tgtseq="Zm00001eb000510/2_Zm00001eb000510_to_A188/Zm00001eb000510_A188.4.target.fas",
    tgtbed="Zm00001eb000510/Zm00001eb000510_to_A188_output/data/A188.label.bed",
    qrysvte="Zm00001eb000510/Zm00001eb000510_to_A188_output/data/SV_on_B73.bed",
    tgtsvte="Zm00001eb000510/Zm00001eb000510_to_A188_output/data/SV_on_A188.bed",
    tgthit="Zm00001eb000510/Zm00001eb000510_to_A188_output/data/Zm00001eb000510_to_A188.hit.output",
    datadir="Zm00001eb000510/Zm00001eb000510_to_A188_output/data",
    figuredir="Zm00001eb000510/Zm00001eb000510_to_A188_output/figures",
    qryinput="Zm00001eb000510/Zm00001eb000510_to_A188_output/data/Zm00001eb000510_B73.search.input",
    nucmeraln="Zm00001eb000510/Zm00001eb000510_to_A188_output/data/Zm00001eb000510_to_A188.uniq.filt.delta.txt"),
  knit_root_dir="/homes/liu3zhen/scripts2/homotools/wikiexample/7_homomine/hmout2",
  output_dir="/homes/liu3zhen/scripts2/homotools/wikiexample/7_homomine/hmout2/Zm00001eb000510",
  output_format="html_document",
  output_file="Zm00001eb000510.homomine.report.html")
