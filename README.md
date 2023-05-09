# Supporting code for SDS analysis and tree building

## How to reproduce the paper figures and tables for trees,  coalescence timings and selection estimates.

### Dependencies
See install_rpackages.R for a list of required packages - also see export/sessionInfo_release.txt for additional information about the environment and package versions.

### Steps to populate export directory 
Starting from trees with assigned mutations (../cache/PDD.RDS) the following produces tree, selection and timing figures together with associated tables in the export directory:

```bash
cd phyloanalysis
Rscript export_paper_analyses.R
Rscript -e "rmarkdown::render('TimingSDS8transformation.Rmd',output_file='../export/TimingSDS8transformation.html')"
```

Additional notes for going from read counts (e.g. output of vafCorrect ( https://github.com/cancerit/vafCorrect )) to filtered variants and trees are available upon request.


