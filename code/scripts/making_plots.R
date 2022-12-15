library(ggplot2)
library(tidyverse)
library(data.table)
library(scales)
setwd('/Users/Oscar/Desktop/plotting/SNP_project')
tab1 <-read.table('dmel-all-r6.46.digSNPs.tajimasd.vk',header = TRUE)
tab2 <-read.table('dmel-all-r6.46.five_prime_UTRSNPs.tajimasd.vk',header = TRUE)
tab3 <-read.table('dmel-all-r6.46.liSNPs.tajimasd.vk',header = TRUE)
tab4 <-read.table('dmel-all-r6.46.nonsynSNPs.tajimasd.vk',header = TRUE)
tab5 <-read.table('dmel-all-r6.46.pigSNPs.tajimasd.vk',header = TRUE)
tab6 <-read.table('dmel-all-r6.46.siSNPs.tajimasd.vk',header = TRUE)
tab7 <-read.table('dmel-all-r6.46.synSNPs.tajimasd.vk',header = TRUE)
tab8 <-read.table('dmel-all-r6.46.three_prime_UTRSNPs.tajimasd.vk',header = TRUE)
tab9 <-read.table('merged.annotated.tajimasd.vk',header = TRUE)

#tabvec <- c('tab1', 'tab2', 'tab3', 'tab4', 'tab5', 'tab6', 'tab7', 'tab8')
#lp <-lapply(tabvec, summary($TajimaD)) #doesn'twork as expected
#for (i in tabvec) {
#  print(i)
#  summary(paste(i$TajimaD))
#}
summary(tab1$TajimaD)
summary(tab2$TajimaD)
summary(tab3$TajimaD)
summary(tab4$TajimaD)
summary(tab5$TajimaD)
summary(tab6$TajimaD)
summary(tab7$TajimaD)
summary(tab8$TajimaD)
summary(tab9$TajimaD)



TDfiles <- dir('/Users/Oscar/Desktop/plotting/SNP_project', pattern = '*.tajimasd.vk', full.names = FALSE)
TDfiles

for (i in seq_along(TDfiles)){
  ifile <- TDfiles[i]
  print(ifile)
  
  TD <- read.table(ifile, header=TRUE)
  TD <- subset(TD, CHROM %in% c("2L", "2R", "3L", "3R", "4", "X"))
  print(summary(TD$TajimaD))
  print(sd(TD$TajimaD))
  
  pi_plot <- ggplot(data=TD, mapping=aes(x= BIN_START, y= TajimaD))
  png(filename = paste(substring(ifile, 16) , ".png"))
  
  print(
    pi_plot + geom_point(aes(color=CHROM)) + 
      geom_hline(yintercept = mean(TD$TajimaD))+
      labs(title=paste("TajimaD, window=10000", substring(ifile, 16,)))
  )
  dev.off()
}







p <- ggplot(data=tab4, mapping=aes(x= BIN_START, y=TajimaD))
p + geom_point(aes(color=CHROM)) + geom_hline(yintercept = mean(tab4$TajimaD, na.rm=TRUE))
p + geom_point(aes(color=CHROM)) + facet_grid(facets = CHROM ~ .)
p + geom_boxplot(mapping=aes(x= CHROM, y=TajimaD, color=CHROM))

two.way <- aov(TajimaD ~ CHROM, data = tab9)
summary(two.way)

summary(pi1$PI)
pi1 <-read.table('dmel-all-r6.46.nonsynSNPs.pi',header = TRUE)
q1 <-ggplot(data=pi1, mapping=aes(x= BIN_START, y= PI))
q1 + geom_point(aes(color=CHROM)) + geom_hline(yintercept = mean( pi1$PI , na.rm=TRUE))
q1 + geom_point(aes(color=CHROM)) + facet_grid(facets = CHROM ~ .)
q1 + geom_boxplot(mapping=aes(x= CHROM, y=TajimaD, color=CHROM))

summary(pi2$PI)
pi2 <-read.table('dmel-all-r6.46.synSNPs.pi',header = TRUE)
q2 <-ggplot(data=pi2, mapping=aes(x= BIN_START, y= PI))
q2 + geom_point(aes(color=CHROM)) + geom_hline(yintercept = mean( pi2$PI , na.rm=TRUE))
q2 + geom_point(aes(color=CHROM)) + facet_grid(facets = CHROM ~ .)
q2 + geom_boxplot(mapping=aes(x= CHROM, y=TajimaD, color=CHROM))


pi2 <-read.table('parr_out.annotated.pi',header = TRUE)

summary(pi2$PI)
sd(pi2$PI)

q2 <-ggplot(data= pi2 %in% c("2L", "2R"), mapping=aes(x= BIN_START, y= PI))
q2 + geom_point(aes(color=CHROM)) + geom_hline(yintercept = mean( pi2$PI , na.rm=TRUE))
q2 + geom_point(aes(color=CHROM)) + facet_grid(facets = CHROM ~ .)

q2 + geom_boxplot(mapping=aes(x= CHROM, y= PI, color=CHROM))





frq1 <- read.table('dmel-all-r6.46.digSNPs.freq2', header=TRUE)
sfs1 <- ggplot(data=frq1)
sfs1 + geom_histogram(mapping = aes(x= N_CHR, binwidth=1)) + 
  labs(title="simple SFS of dig SNPs",x ="allele frequency", y = "count")+ 
  scale_x_continuous(breaks=c(1:30))

frq2 <- read.table('dmel-all-r6.46.five_prime_UTRSNPs.freq2', header=TRUE)
sfs2 <- ggplot(data=frq2)
sfs2 + geom_histogram(mapping = aes(x= N_CHR), binwidth=1) + 
  labs(title="simple SFS of 5'_UTR SNPs",x ="allele frequency", y = "count")+ 
  scale_x_continuous(breaks=c(1:30))

frq3 <- read.table('dmel-all-r6.46.liSNPs.freq2', header=TRUE)
sfs3 <- ggplot(data=frq3)
sfs3 + geom_histogram(mapping = aes(x= N_CHR), binwidth=1) + 
  labs(title="simple SFS of li SNPs",x ="allele frequency", y = "count")+ 
  scale_x_continuous(breaks=c(1:30))

frq4 <- read.table('dmel-all-r6.46.merged.annotated.freq2', header=TRUE)
sfs4 <- ggplot(data=frq4)
sfs4 + geom_histogram(mapping = aes(x= N_CHR), binwidth=1) + 
  labs(title="simple SFS of ALL SNPs",x ="allele frequency", y = "count")+ 
  scale_x_continuous(breaks=c(1:30))

frq5 <- read.table('dmel-all-r6.46.nonsynSNPs.freq2', header=TRUE)
sfs5 <- ggplot(data=frq5)
sfs5 + geom_histogram(mapping = aes(x= N_CHR), binwidth=1) + 
  labs(title="simple SFS of NONSYNonymous SNPs",x ="allele frequency", y = "count")+ 
  scale_x_continuous(breaks=c(1:30))

frq6 <- read.table('dmel-all-r6.46.pigSNPs.freq2', header=TRUE)
sfs6 <- ggplot(data=frq6)
sfs6 + geom_histogram(mapping = aes(x= N_CHR), binwidth=1) + 
  labs(title="simple SFS of pig SNPs",x ="allele frequency", y = "count")+ 
  scale_x_continuous(breaks=c(1:30))

frq7 <- read.table('dmel-all-r6.46.siSNPs.freq2', header=TRUE)
sfs7 <- ggplot(data=frq7)
sfs7 + geom_histogram(mapping = aes(x= N_CHR), binwidth=1) + 
  labs(title="simple SFS of si SNPs",x ="allele frequency", y = "count")+ 
  scale_x_continuous(breaks=c(1:30))

frq8 <- read.table('dmel-all-r6.46.synSNPs.freq2', header=TRUE)
sfs8 <- ggplot(data=frq8)
sfs8 + geom_histogram(mapping = aes(x= N_CHR), binwidth=1) + 
  labs(title="simple SFS of SYNonymous SNPs",x ="allele frequency", y = "count")+ 
  scale_x_continuous(breaks=c(1:30))

jpeg(filename = "3_prime_UTR.jpg",
    width = 1200, height = 1200, units = "px", pointsize = 14,
    bg = "white", res = NA, family = "", restoreConsole = TRUE)
frq9 <- read.table('dmel-all-r6.46.three_prime_UTRSNPs.freq2', header=TRUE)
sfs9 <- ggplot(data=frq9)
sfs9 + geom_histogram(mapping = aes(x= N_CHR), binwidth=1) + 
  labs(title="simple SFS of 3'UTR SNPs",x ="allele frequency", y = "count")+ 
  scale_x_continuous(breaks=c(1:30))
dev.off()

frq1$TYPE <- "dig"
frq2$TYPE <- "five_prime_UTR"
frq3$TYPE <- "li"
#data.frame(append(df1, c(x1='z'), after=1))
frqall <- rbind(frq1, frq2, frq3)
sfsall <- ggplot(data=frqall)
sfsall + geom_histogram(mapping = aes(x= N_CHR, y = (..count..)/sum(..count..), fill= TYPE), binwidth=1, position = "dodge") + 
  labs(title="simple SFS of all",x ="allele frequency", y = "count")+ 
  scale_x_continuous(breaks=c(1:30))


freq2files <- dir('/Users/Oscar/Desktop/plotting/SNP_project', pattern = '*.freq2', full.names = FALSE)
freq2files

for (i in seq_along(freq2files)){
  ifile <- freq2files[i]
  frq <- read.table(ifile, header=TRUE)
  #frq <- subset(frq, CHROM %in% c("2L", "2R", "3L", "3R", "4", "X"))
  #frq <- subset(x = frq, N_CHR %% 2 == 0)
  sfsv <- ggplot(data=frq)
  png(filename = paste(substring(ifile, 16,) , ".png"))
  print(sfsv + geom_histogram(mapping = aes(x= N_CHR/2), binwidth=1) + 
    labs(title=paste("simple SFS of", substring(ifile, 16,)), x ="allele frequency", y = "count")+ 
    scale_x_continuous(breaks=c(1:30)))
  dev.off()
}


freqlist <- lapply(X = freq2files, FUN = function(X) {read.table(file = X, header = TRUE)})
plotlist <- lapply(
  X = freqlist,
  FUN = function(X)
    ggplot(data = X, mapping = aes( x= N_CHR)) +
    geom_histogram(binwidth = 1) +
    labs(title="simple SFS of ****i",x ="allele frequency", y = "count") + 
    scale_x_continuous(breaks=c(1:30))
)




pi2 <-read.table('parr_out.annotated.pi',header = TRUE)
pi2<- subset(pi2, CHROM %in% c("2L", "2R", "3L", "3R", "4", "X"))
summary(pi2$PI)
sd(pi2$PI)
pi2 <-read.table("dmel-all-r6.46.r6.annotated.pi",header = TRUE)


q2 <-ggplot(data= pi2, mapping=aes(x= BIN_START, y= PI))


q2 + geom_point(aes(color=CHROM)) + geom_hline(yintercept = mean( pi2$PI , na.rm=TRUE))
q2 + geom_point(aes(color=CHROM)) + facet_grid(facets = CHROM ~ .)

q2 + geom_boxplot(mapping=aes(x= CHROM, y= PI, color=CHROM))

pifiles <- dir('/Users/Oscar/Desktop/plotting/SNP_project', pattern = '*.pi', full.names = FALSE)
pifiles

for (i in seq_along(pifiles)){
  ifile <- pifiles[i]
  print(ifile)

  pi <- read.table(ifile, header=TRUE)
  pi <- subset(pi, CHROM %in% c("2L", "2R", "3L", "3R", "4", "X"))
  print(summary(pi$PI))
  print(sd(pi$PI))
  
  pi_plot <- ggplot(data=pi, mapping=aes(x= BIN_START, y= PI))
  png(filename = paste(substring(ifile, 16) , ".png"))
  
  print(
    pi_plot + geom_point(aes(color=CHROM)) + 
      geom_hline(yintercept = mean(pi$PI))+
          labs(title=paste("pi, window=10000", substring(ifile, 16,)))
               )
  dev.off()
}


c = 0
length(fl)
l <- c(0:length(fl))

for (i in l) {
  print(i)
  print(paste("tab",i, sep=" "))
  a <- paste("tab",i, sep=" ")
  print(fl[i])
  #read.table(fl[i], header = TRUE)
}

load_data <- function(path) { 
  files <- dir(path, pattern = '*.vk', full.names = TRUE)
  tables <- lapply(files, read.table)
}

t <- load_data('/Users/Oscar/Downloads/')
