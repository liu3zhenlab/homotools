library(rmarkdown)
library(knitr)

render('Zm00001eb214410/homomine_report.Rmd',
  params = list(
    cwd="/homes/liu3zhen/scripts2/homotools/wikiexample/7_homomine/hmout",
    scriptdir="/homes/liu3zhen/scripts2/homotools/utils/",
    qrygene="Zm00001eb214410",
    qrybase="B73-5.0",
    qryseq="Zm00001eb214410/1_B73-5.0_Zm00001eb214410/1.Zm00001eb214410.fasta",
    qrybed="Zm00001eb214410/Zm00001eb214410_to_A188v1_output/data/B73-5.0.label.bed",
    tgtbase="A188v1",
    tgtgene="Zm00056a007222",
    tgtseq="Zm00001eb214410/2_Zm00001eb214410_to_A188v1/Zm00001eb214410_A188v1.4.target.fas",
    tgtbed="Zm00001eb214410/Zm00001eb214410_to_A188v1_output/data/A188v1.label.bed",
    qrysvte="Zm00001eb214410/Zm00001eb214410_to_A188v1_output/data/SV_on_B73-5.0.bed",
    tgtsvte="Zm00001eb214410/Zm00001eb214410_to_A188v1_output/data/SV_on_A188v1.bed",
    tgthit="Zm00001eb214410/Zm00001eb214410_to_A188v1_output/data/Zm00001eb214410_to_A188v1.hit.output",
    datadir="Zm00001eb214410/Zm00001eb214410_to_A188v1_output/data",
    figuredir="Zm00001eb214410/Zm00001eb214410_to_A188v1_output/figures",
    qryinput="Zm00001eb214410/Zm00001eb214410_to_A188v1_output/data/Zm00001eb214410_B73-5.0.search.input",
    nucmeraln="Zm00001eb214410/Zm00001eb214410_to_A188v1_output/data/Zm00001eb214410_to_A188v1.uniq.filt.delta.txt"),
  knit_root_dir="/homes/liu3zhen/scripts2/homotools/wikiexample/7_homomine/hmout",
  output_dir="/homes/liu3zhen/scripts2/homotools/wikiexample/7_homomine/hmout/Zm00001eb214410",
  output_format="html_document",
  output_file="Zm00001eb214410.homomine.report.html")
