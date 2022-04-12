library(ggplot2)
library(tidyr)
library(dplyr)
setwd("~/Desktop/Covid/wastewater/R")


########################################################################################
#  Comparison of average variant frequency for resequenced samples
########################################################################################

### input file
dataTable <- read.csv("tallymut_line_snv2_final.tsv", sep = "\t")


### Batches
oldBatch_1 <- "20201223_HWKGTDRXX"    # these two batches contain the original sequencing data
oldBatch_2 <- "20210108_JDL4L"
newBatch   <- "20210129_JG6LR"        # this batch contains the resequencing data


### Dates of resequenced Lausanne samples
vd_dates <- c("2020-12-08","2020-12-09","2020-12-10","2020-12-13","2020-12-15","2020-12-17","2020-12-19"
              ,"2020-12-25", "2020-12-27", "2020-12-31", "2021-01-02", "2021-01-04")

### Dates of resequenced Zurich samples
zh_dates <- c("2020-12-17","2020-12-19","2020-12-20","2020-12-21", "2020-12-23", "2020-12-25",
              "2020-12-29", "2020-12-31", "2021-01-02", "2021-01-04")


### select the data rows that correspond to these samples
scatterDat <- select(dataTable, date, plantname, batch, pos, frac)       # select relevant columns
scatterDat$sample <- paste(scatterDat$date,scatterDat$plantname)
scatterDat1 <- scatterDat[is.element(scatterDat$date, vd_dates) & scatterDat$plantname=="Lausanne (VD)",]
scatterDat2 <- scatterDat[is.element(scatterDat$date, zh_dates) & scatterDat$plantname=="Zürich (ZH)",]
scatterDat <- rbind(scatterDat1, scatterDat2)


### rename the batches and merge the two original batches
scatterDat$batch <- gsub('20201223_HWKGTDRXX', 'batch1', scatterDat$batch)
scatterDat$batch <- gsub('20210108_JDL4L', 'batch1', scatterDat$batch)
scatterDat$batch <- gsub('20210129_JG6LR', 'batch2', scatterDat$batch)

### select rows that correspond to the batches
scatterDat <- scatterDat[scatterDat$batch == 'batch2' | scatterDat$batch == 'batch1',]
scatterDat <- na.omit(scatterDat)

### create wide table
scatterDatSpread <- spread(scatterDat, key=batch, value=frac) #, convert = FALSE, drop = TRUE)
scatterDatSpread <- na.omit(scatterDatSpread)


### compute mean fraction per sample (per batch)
scatterDatSpread$pos <- as.factor(scatterDatSpread$pos)
scatterDatSpread <- select(scatterDatSpread, sample, batch1, batch2)
scatterDatSpread <- aggregate(.~ sample, data=scatterDatSpread, FUN=mean)


### create scatter plot
p <- ggplot(scatterDatSpread, aes(x=batch1, y=batch2)) +
    geom_point(shape = 21, alpha = 0.5, size = 2, fill="steelblue") +
    theme_classic() +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12)) +
    #theme(legend.position="none") +
    xlab("B.1.1.7 fraction at first sequencing") +
    ylab("B.1.1.7 fraction at second sequencing") +
    geom_smooth(formula = y ~ x, method='lm')+
    #scale_color_brewer(palette = "Spectral")+
    stat_cor(label.y=0.085, method="pearson", na.rm = TRUE, size=5)

p

