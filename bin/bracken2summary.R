#!/usr/bin/env Rscript

# generate a DT datatable and heatmap summary for all samples
# reads all the bracken tables - Channel.collect()
# and generate a summary of abundance per taxon across samples
# save datatable and heatmap widgets

library(data.table)
library(dplyr)
library(tidyr)
library(DT)
library(d3heatmap)


arg <- commandArgs(trailingOnly = TRUE)

bracken2summary <- function(files) {
  names(files) <- basename(tools::file_path_sans_ext(files))
  # this is the tidy data
  df <- lapply(files, fread) %>%
      bind_rows(.id = "id") %>% 
      dplyr::select(c(1,2,8)) %>% 
      dplyr::arrange(desc(fraction_total_reads))
  df %>% fwrite("bracken_summary_long.csv")
  
  #this is a wide form, for datatable
  dfw <- df %>%
      tidyr::pivot_wider(names_from = id, 
                       values_from = fraction_total_reads)
                       #values_fill = list(fraction_total_reads = 0))  # fill 0 if no value there
  dfw %>% fwrite("bracken_summary_wide.csv")
  
  # top 10 species for each sample, saved as wide
  dfw_top10 <- df %>%
    group_by(id) %>%
    top_n(10, fraction_total_reads) %>% # take top10 of each sample
    tidyr::pivot_wider(names_from = id, 
                       values_from = fraction_total_reads) %>%
    head(50) # take max 50 rows of the table
                      # values_fill = list(fraction_total_reads = 0))
  #----
  # make and save datatable widget 
  # only if samples <= 12
  if(length(files) <= 12) {
    dfw %>%
      datatable(class = 'hover row-border order-column', 
                extensions = 'Buttons',
                options = list(searchHighlight = TRUE, 
                             dom = 'Brtip', 
                             buttons = c('copy', 'csv', 'excel'),
                             pageLength = 20), 
                filter = 'top', 
                selection = 'single') %>% 
      formatPercentage(c(2:ncol(dfw)), digits = 2) %>% #format all columns as percent
      DT::saveWidget(file = "bracken_summary_table.html")
    }
  
  #----
  # # make and save heatplot for top10 taxa for samples >=2 and < 24
   if(length(files) >= 2 & length(files) <= 24) {
     data.frame(row.names = dfw_top10$name, dfw_top10[,-1]) %>% 
       d3heatmap(colors = "YlOrRd", 
                 Rowv = FALSE, 
                 show_grid = 1,
                 dendrogram = 'none',
                 xaxis_font_size = "10px", 
                 yaxis_font_size = "10px")  %>%
       DT::saveWidget(file = "bracken_summary_heatmap.html")
     }
  
}
# execute the function
bracken2summary(arg)
