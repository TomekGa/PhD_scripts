## Files you to analyse the original output from the AmpliSAS software

1. Amplisas_output_analysis.R - semi-automatic analysis of AmpliSAS output - includes various filtering stages, alignments of allele sequences and generation of final binary matrix of genotypes
2. Functions_to_be_loaded.R - set of functions used for (almost) all other scripts
3. logit_PAF.R - plotting logit transformation of PAF values to choose a better PAF threshold
4. Duplicates_trees.R - drawing trees to compare two repeated but independent runs graphically
5. Dendograms.R - drawing dendrograms of individuals (based on the binary matrix of genotypes) or alleles (based on aligned sequences)
6. cluster_examination.R - examining similar MHC alleles to look for the patterns of improper clustering parameters for AmpliSAS
7. Individual_trees_shinyApp.R - interactive app to display the tree of all found alleles in a given run and mark alleles of a particular chosen individual. It also helps to set the proper PAF threshold. To check how it works, run the command below in R:
   
```{r}
require(shiny)
runGitHub("TomekGa/PhD_scripts",ref = "main",subdir = "AmpiSAS_resultAnalysis/Individual_trees_shinyApp.R")
```
