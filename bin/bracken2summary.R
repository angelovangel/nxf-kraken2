#!/usr/bin/env Rscript

# generate a DT datatable summary for all samples
# reads all the bracken tables ( Channel.collect()? ) and generate a summary of abundance per taxon across samples

library(data.table)
library(dplyr)
library(tidyr)
library(DT)

arg <- commandArgs(trailingOnly = TRUE)

bracken2summary <- function(files) {
  lapply(files, fread) %>%
    bind_rows(.id = "id") %>% 
    dplyr::select(c(1,2,8)) %>% 
    dplyr::arrange(desc(fraction_total_reads)) %>%
    tidyr::pivot_wider(names_from = id, values_from = fraction_total_reads) %>% 
    datatable(class = 'hover row-border order-column', 
              extensions = 'Buttons',
              options = list(searchHighlight = TRUE, 
                             dom = 'Brtip', 
                             buttons = c('copy', 'csv', 'excel'),
                             pageLength = 20), 
              filter = 'top', 
              selection = 'single') %>% 
    formatPercentage(c(2,3), digits = 2) %>%
    DT::saveWidget(file = "bracken_summary.html")
}

bracken2summary(arg)