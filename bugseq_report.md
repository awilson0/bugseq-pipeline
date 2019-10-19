bugseq\_report
================
Heather
2019-10-18

``` r
library(tidyverse)
```

    ## ── Attaching packages ───────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.2.1     ✔ purrr   0.3.2
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.3
    ## ✔ tidyr   1.0.0     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
# prepare ops data
#source("src/extract_performance.R")
#performance <- extract_performance()
```

``` r
resfinder_df = read_tsv("resfinder.tsv")
```

    ## Parsed with column specification:
    ## cols(
    ##   `Isolate ID` = col_character(),
    ##   Gene = col_character(),
    ##   `Predicted Phenotype` = col_character(),
    ##   `%Identity` = col_double(),
    ##   `%Overlap` = col_double(),
    ##   `HSP Length/Total Length` = col_character(),
    ##   Contig = col_character(),
    ##   Start = col_double(),
    ##   End = col_double(),
    ##   Accession = col_character()
    ## )

``` r
head(resfinder_df)
```

    ## # A tibble: 6 x 10
    ##   `Isolate ID` Gene  `Predicted Phen… `%Identity` `%Overlap`
    ##   <chr>        <chr> <chr>                  <dbl>      <dbl>
    ## 1 GCF_0019315… aac(… gentamicin              99.9      100  
    ## 2 GCF_0019315… aph(… kanamycin              100        100  
    ## 3 GCF_0019315… aph(… hygromicin             100        100  
    ## 4 GCF_0019315… blaC… ampicillin, cef…       100        100  
    ## 5 GCF_0019315… dfrA… trimethoprim            99.8      100  
    ## 6 GCF_0019315… floR  chloramphenicol         98.2       99.9
    ## # … with 5 more variables: `HSP Length/Total Length` <chr>, Contig <chr>,
    ## #   Start <dbl>, End <dbl>, Accession <chr>

``` r
detailed_summary_df = read_tsv("detailed_summary.tsv")
```

    ## Parsed with column specification:
    ## cols(
    ##   `Isolate ID` = col_character(),
    ##   `Gene/Plasmid` = col_character(),
    ##   `Predicted Phenotype` = col_character(),
    ##   `%Identity` = col_double(),
    ##   `%Overlap` = col_double(),
    ##   `HSP Length/Total Length` = col_character(),
    ##   Contig = col_character(),
    ##   Start = col_double(),
    ##   End = col_double(),
    ##   Accession = col_character(),
    ##   `Data Type` = col_character()
    ## )

``` r
head(detailed_summary_df)
```

    ## # A tibble: 6 x 11
    ##   `Isolate ID` `Gene/Plasmid` `Predicted Phen… `%Identity` `%Overlap`
    ##   <chr>        <chr>          <chr>                  <dbl>      <dbl>
    ## 1 GCF_0019315… None           <NA>                    NA           NA
    ## 2 GCF_0019315… aac(3)-IVa     gentamicin              99.9        100
    ## 3 GCF_0019315… aph(3')-Ia     kanamycin              100          100
    ## 4 GCF_0019315… aph(4)-Ia      hygromicin             100          100
    ## 5 GCF_0019315… blaCTX-M-65    ampicillin, cef…       100          100
    ## 6 GCF_0019315… dfrA14         trimethoprim            99.8        100
    ## # … with 6 more variables: `HSP Length/Total Length` <chr>, Contig <chr>,
    ## #   Start <dbl>, End <dbl>, Accession <chr>, `Data Type` <chr>
