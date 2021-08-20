library( wasabi )
library( here )
library( tidyverse )

h5_file <- snakemake@input[[1]]

salmon_dir <- str_replace( string = h5_file, pattern = "/quant.sf", replacement = "" )

prepare_fish_for_sleuth( salmon_dir )
