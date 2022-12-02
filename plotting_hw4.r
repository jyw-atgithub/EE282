#!#!/usr/bin/Rscript

library(ggplot2)
library(data.table)
library(scales)
# for length distribution
lenfiles <- dir('/Users/Oscar/Desktop/plotting', pattern = '*length.txt', full.names = FALSE)
for (f in lenfiles) {
  tab <- read.table(f, header=FALSE)
  len <- ggplot(data=tab)

  png(filename = paste(substring(f,0,14), ".png"),width = 800, height = 600, units = "px",)

  print(
    len + geom_histogram(mapping=aes(x=tab[,2]) , binwidth = 0.04) + scale_x_log10()+
    labs(title= paste(f), x = "log10 length (bp)", y = "count" )
  )

  dev.off()
}

# for GC content distribution
gcfiles <- dir('/Users/Oscar/Desktop/plotting', pattern = '*gc.txt', full.names = FALSE)

for (f in gcfiles) {
  tab <- read.table(f)
  gc <- ggplot(data= tab)
  png(filename = paste(substring(f,0,10), ".png"),width = 800, height = 600, units = "px")
  print(
    gc + geom_histogram(mapping=aes(x=tab[,2]), binwidth = 0.005) +
    labs(title= paste(f), x = "GC content (%)", y = "count" )
  )
  dev.off()
}