#!/usr/bin/env Rscript

# generate a DT datatable from bracken table output and save it as a standalone html

library(data.table)
library(dplyr)
library(DT)

arg <- commandArgs(trailingOnly = TRUE)

bracken2dt <- function(path, outfile) {
  dt <- fread(path) %>% 
        arrange(desc(fraction_total_reads)) %>% 
        DT::datatable()
  DT::saveWidget(dt, file = outfile)
  }


#cat(arg, "\n")
bracken2dt(arg[1], arg[2])
