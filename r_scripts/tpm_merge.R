library(tidyverse)
library(ggplot2)
library(here)

files <- paste("output/salmon",dir("output/salmon"), sep="/")
sample.tpm.df <- read.table()


######### TPM MERGE EXAMPLE

```{r libraries}
library( tidyverse )
library( here )
```
```{r rginfo}
samples <- read_tsv( file = here( "data", "project_rgid_info.txt" ) ) %>%
  pull( RGSM ) %>%
  unique( . )
```
```{r read_salmon}
tibble_count <- tibble()
tibble_tpm <- tibble()
for ( idx in 1:length( samples ) ) {
  curr_rgsm <- samples[ idx ]

  # note that if it's compressed just add the .gz to it. read_tsv can read gzipped files also.
  # columns: Name   Length  EffectiveLength TPM NumReads
  current_salmon <- read_tsv( file = here( "results", "salmon", curr_rgsm, "quant.gene.sf" ), col_types = 'cdddd' )

  # if this was transcripts you'd use 'transcript_id' instead of 'gene_id'
  temp_tpm <- select( current_salmon, Name, TPM )
  colnames( temp_tpm ) <- c( 'gene_id', curr_rgsm )

  temp_count <- select( current_salmon, Name, NumReads )
  colnames( temp_count ) <- c( 'gene_id', curr_rgsm )

  # counts
  if ( sum( dim( tibble_count ) ) > 0 ) {
    tibble_count = merge( x = tibble_count, y = temp_count, by = "gene_id" )
  } else {
    tibble_count = temp_count
  }

  # tpm
  if ( sum( dim( tibble_tpm ) ) > 0 ) {
    tibble_tpm = merge( x = tibble_tpm, y = temp_tpm, by = "gene_id" )
  } else {
    tibble_tpm = temp_tpm
  }

  rm( curr_rgsm )
  rm( current_salmon )
  rm( temp_tpm )
  rm( temp_count )
}
```
```{r write_output}
## adding the gz automatically triggers write_tsv to compress the file with gzip
# count
write_tsv( x = tibble_count, path = here( 'results', 'merged_salmon_gene_count.tsv.gz' ) )
# tpm
write_tsv( x = tibble_tpm, path = here( 'results', 'merged_salmon_gene_tpm.tsv.gz' ) )
```
```{r session_information}
Sys.time()
getwd()

sessionInfo()
```
