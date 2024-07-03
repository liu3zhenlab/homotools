library(rmarkdown)
library(knitr)

render('Zm00001eb002700/homomine_report.Rmd',
  params = list(
    cwd="/homes/liu3zhen/scripts2/homotools/wikiexample/7_homomine/hmout",
    scriptdir="/homes/liu3zhen/scripts2/homotools/utils/",
    qrygene="Zm00001eb002700",
    qrybase="B73",
    qryseq="Zm00001eb002700/1_B73_Zm00001eb002700/1.Zm00001eb002700.fasta",
    qrybed="Zm00001eb002700/Zm00001eb002700_to_A188_output/data/B73.label.bed",
    tgtbase="A188",
    tgtgene="",
    tgtseq="Zm00001eb002700/2_Zm00001eb002700_to_A188/Zm00001eb002700_A188.4.target.fas",
    tgtbed="empty",
    qrysvte="Zm00001eb002700/Zm00001eb002700_to_A188_output/data/SV_on_B73.bed",
    tgtsvte="Zm00001eb002700/Zm00001eb002700_to_A188_output/data/SV_on_A188.bed",
    tgthit="Zm00001eb002700/Zm00001eb002700_to_A188_output/data/Zm00001eb002700_to_A188.hit.output",
    datadir="Zm00001eb002700/Zm00001eb002700_to_A188_output/data",
    figuredir="Zm00001eb002700/Zm00001eb002700_to_A188_output/figures",
    qryinput="Zm00001eb002700/Zm00001eb002700_to_A188_output/data/Zm00001eb002700_B73.search.input",
    nucmeraln="Zm00001eb002700/Zm00001eb002700_to_A188_output/data/Zm00001eb002700_to_A188.uniq.filt.delta.txt"),
  knit_root_dir="/homes/liu3zhen/scripts2/homotools/wikiexample/7_homomine/hmout",
  output_dir="/homes/liu3zhen/scripts2/homotools/wikiexample/7_homomine/hmout/Zm00001eb002700",
  output_format="html_document",
  output_file="Zm00001eb002700.homomine.report.html")
