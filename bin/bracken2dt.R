#!/usr/bin/env Rscript

# generate a DT datatable from bracken table output and save it as a standalone html

library(data.table)
library(dplyr)
library(DT)

arg <- commandArgs(trailingOnly = TRUE)
taxid_string <- "<a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id="
# <a href="http://rstudio.com">RStudio</a>

bracken2dt <- function(path, outfile) {
  fread(path) %>% 
    dplyr::arrange(desc(fraction_total_reads)) %>% 
    dplyr::select(c(1,2,6,7)) %>%
    dplyr::mutate(taxonomy_id = paste(taxid_string, taxonomy_id, "'>", taxonomy_id, "</a", sep = "") ) %>% # make html links to tax db
  
    DT::datatable(class = 'hover row-border order-column', 
                  escape = FALSE, # needed for the html links
                  colnames = c("Species", "Taxonomy_ID", "Estimated reads" , "Fraction of total reads"),
                  extensions = 'Buttons',
                  options = list(searchHighlight = TRUE, 
                                 dom = 'Brtip', 
                                 buttons = c('copy', 'csv', 'excel'),
                                 pageLength = 12), 
                  filter = 'top', 
                  selection = 'single', 
                  caption = paste("Bracken abundance estimation for", 
                                  tools::file_path_sans_ext(basename(outfile) )
                                  ) 
                  ) %>%
    formatPercentage(4, digits = 2) %>%
    DT::saveWidget(file = outfile)
  }


#cat(arg, "\n")
bracken2dt(arg[1], arg[2])
