#!/usr/bin/env Rscript

# generate a DT datatable summary for all samples
# reads all the bracken tables - Channel.collect()
# and generate a summary of abundance per taxon across samples

library(data.table)
library(dplyr)
library(tidyr)
library(DT)
#library(apexcharter) # too complicated, not used

arg <- commandArgs(trailingOnly = TRUE)

bracken2summary <- function(files) {
  names(files) <- basename(tools::file_path_sans_ext(files))
  # this is the data used for apexcharter
  df <- lapply(files, fread) %>%
      bind_rows(.id = "id") %>% 
      dplyr::select(c(1,2,8)) %>% 
      dplyr::arrange(desc(fraction_total_reads))
  df %>% fwrite("bracken_summary_long.csv")
  
  #this is a wide form, for datatable
  dfw <- df %>%
      tidyr::pivot_wider(names_from = id, 
                       values_from = fraction_total_reads, 
                       values_fill = list(fraction_total_reads = 0))  # fill 0 if no value there
  dfw %>% fwrite("bracken_summary.csv")
  
  # make and save datatable widget 
  # only if samples <= 6 
  if(length(files) <= 6) {
  dfw %>%
    datatable(class = 'hover row-border order-column', 
              extensions = 'Buttons',
              options = list(searchHighlight = TRUE, 
                             dom = 'Brtip', 
                             buttons = c('copy', 'csv', 'excel'),
                             pageLength = 20), 
              filter = 'top', 
              selection = 'single') %>% 
    formatPercentage(c(2,3), digits = 2) %>%
    DT::saveWidget(file = "bracken_summary_table.html") 
    }
  #----
  # make and save apex chart (top10 species only)
  # only if samples <= 6 
  # df %>% 
  #   dplyr::top_n(10, fraction_total_reads) %>%
  #   apex(type = "bar", aes(id, fraction_total_reads, fill = name)) %>% 
  #   ax_chart(stacked = TRUE, stackType = "100%") %>% 
  #   ax_legend(horizontalAlign = "left", 
  #             position = "bottom",
  #             itemMargin = list(horizontal =5, vertical = 0)) %>%
  #   ax_colors(scales::brewer_pal(palette = "Set2")(8)) %>%
  #   ax_dataLabels(enabled = TRUE) %>%
  #   htmlwidgets::saveWidget(file = "bracken_summary_plot.html")
  
}
# execute the function
bracken2summary(arg)
