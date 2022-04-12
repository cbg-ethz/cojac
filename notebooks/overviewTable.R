library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(ggnewscale)
library(tidyverse)
library(tidyr)
library(reshape2)
library(ggpubr)
library(ggh4x)
library(reshape2)
library(grid)

###########################################################################################################
selctData.function <- function(dt, usedPlants, usedMutTypes, omitSample){
  
  # filter for plants of interest
  dt <- filter(dt, is.element(dt$plantname, usedPlants))
  dt$plantname <- factor(dt$plantname)
  
  # filter out the first ski resort sample from 2020-12-20 because it is not actually from the resort but a larger area
  dt <- filter(dt, dt$plantname != 'Ski resort' | as.character(dt$date) != '2020-12-20')
  dt$date <- factor(dt$date)
  
  return(dt)
}

#####################################################################

addVariantColData.function <- function(dt, alVar, beVar, gaVar){
  # split data table by variant and add variant column
  al <- dt[is.element(dt$al, usedMutTypes)==TRUE,]
  be <- dt[is.element(dt$be, usedMutTypes)==TRUE,]
  ga <- dt[is.element(dt$ga, usedMutTypes)==TRUE,]
  
  al$variant <- alVar
  be$variant <- beVar
  ga$variant <- gaVar
  
  #uk$variant <- factor(uk$variant)
  #za$variant <- factor(za$variant)
  #br$variant <- factor(br$variant)
  
  return(rbind(al, be, ga))
}

#####################################################################

addVariantColAmpli.function <- function(dt, alVar, beVar, gaVar){
  # split data table by variant and add variant column
  al <- dt[!is.na(dt$al),]
  be <- dt[!is.na(dt$be),]
  ga <- dt[!is.na(dt$ga),]
  
  al$variant <- alVar
  be$variant <- beVar
  ga$variant <- gaVar
  
  return(rbind(al, be, ga))
}


addAmpliColData.function <- function(dt, alVar, beVar, gaVar){
  
  al <- filter(dt, dt$variant==alVar)
  be <- filter(dt, dt$variant==beVar)
  ga <- filter(dt, dt$variant==gaVar)
  
  al <- select(al, sample, pos, date, plantname, cov, var, variant)
  be <- select(be, sample, pos, date, plantname, cov, var, variant)
  ga <- select(ga, sample, pos, date, plantname, cov, var, variant)
  
  # alpha amplicon information
  a72_AL <- c(21765, 21766, 21767, 21768, 21769, 21770, 21991, 21992, 21993)
  a78_AL <- c(23604, 23709)
  a92_AL <- c(27972, 28048 ,28111)
  a93_AL <- c(28111, 28280, 28281, 28282)
  
  al$ampli <- al$pos
  al$ampli <- ifelse(al$pos<21765, "1-71", al$ampli)
  al$ampli <- ifelse(is.element(al$ampli, a72_AL), "72", al$ampli)
  al$ampli <- ifelse(al$pos>21993 & al$pos<23604, "73-77", al$ampli)
  al$ampli <- ifelse(is.element(al$ampli, a78_AL), "78", al$ampli)
  al$ampli <- ifelse(al$pos>23709 & al$pos<27972, "79-91", al$ampli)
  al$ampli <- ifelse(is.element(al$ampli, a92_AL), "92", al$ampli)
  al$ampli <- ifelse(is.element(al$ampli, a93_AL), "93", al$ampli)
  al$ampli <- ifelse(al$pos>28283, "94-98", al$ampli)
  
  # positions appears both in amplicon 92 an 93
  repeatPos <- filter(al, al$pos==28111)
  repeatPos$ampli <- "93"
  
  al <- rbind(al, repeatPos)
  
  # gamma amplicon information
  a71_GA <- c(21614, 21621 ,21638) 
  a76_GA <- c(23012, 23063)
  a95_GA <- c(28877, 28878)
  
  ga$ampli <- ga$pos
  ga$ampli <- ifelse(ga$pos<21614, "1-70", ga$ampli)
  ga$ampli <- ifelse(is.element(ga$ampli, a71_GA), "71", ga$ampli)
  ga$ampli <- ifelse(ga$pos>21638 & ga$pos<23012, "72-75", ga$ampli)
  ga$ampli <- ifelse(is.element(ga$ampli, a76_GA), "76", ga$ampli)
  ga$ampli <- ifelse(ga$pos>23063 & ga$pos<28877, "77-94", ga$ampli)
  ga$ampli <- ifelse(is.element(ga$ampli, a95_GA), "95", ga$ampli)
  ga$ampli <- ifelse(ga$pos>28878, "96-98", ga$ampli)
  
  
  
  # beta amplicon information
  a76_BE <- c(23012, 23063)
  
  be$ampli <- be$pos
  be$ampli <- ifelse(be$pos<23012, "1-75", be$ampli)
  be$ampli <- ifelse(is.element(be$ampli, a76_BE), "76", be$ampli)
  be$ampli <- ifelse(be$pos>23063, "77-98", be$ampli)
  
  return(rbind(al, be, ga))
}


#####################################################################

aggregateMultiBaseEvents.function <- function(dt, multiBaseEventsList){
  
  # update pos column by multi-base-events
  for(intv in multiBaseEventsList){
    multiBaseEvent <- paste(head(unlist(intv),1), "-", tail(unlist(intv),1), sep='')
    dt$pos <- ifelse(is.element(dt$pos, intv), multiBaseEvent, dt$pos)
    
  }  
  
  # aggregate rows that describe the same multi-base event, get average frac
  dataTable_agg <- aggregate(. ~  sample + date + pos + plantname + ampli + variant, data = dt, sum)
  
  
  
  #print(names(dataTable_agg))
  
  dataTable_agg$frac <- ifelse(dataTable_agg$cov == 0, NA, dataTable_agg$var/dataTable_agg$cov)
  
  
  #return(dataTable_agg[, !names(dataTable_agg) %in% c("event")] )  # return all but the event column
  return(dataTable_agg)  # return all but the event column
}
###########################################################################################################


setwd("~/Desktop/Covid/wastewater/R")

### read input files
dataTable <- read.csv("tallymut_line_snv2_final.tsv", sep = "\t")
ampliTable <- read.csv("cooc-output-lines.csv", sep = ",", na.strings=c(""))  # amplicon table

### get rid of Umlaut for Zurich
levels(dataTable$plantname) <- c(levels(dataTable$plantname), "Zurich (ZH)")
dataTable$plantname[as.character(dataTable$plantname)=='Zürich (ZH)'] <- 'Zurich (ZH)'

### select
dataTable  <- select(dataTable, batch, sample, pos, date, plantname, cov, var, frac, ga, al, be)
dataTable$date <- as.Date(dataTable$date, format="%Y-%m-%d")


### harmonize column names
connector <- unique(select(dataTable, sample, date, plantname))
ampliTable <- merge(x = ampliTable, y = connector, by = "sample", all.x = TRUE)
ampliTable$pos <- paste("Amplicon", ampliTable$amplicon, sep=" ")
names(ampliTable)[names(ampliTable) == 'amplicon'] <- 'ampli'
names(ampliTable)[names(ampliTable) == 'count'] <- 'cov'
names(ampliTable)[names(ampliTable) == 'mut_all'] <- 'var'


####################################################################
##        data selection                                           #
####################################################################

##### which plants to use
usedPlants <- c('Ski resort', 'Lausanne (VD)', 'Zurich (ZH)')

##### which mutation types to use (mut, extra, shared)
usedMutTypes <- c('mut', 'extra')

##### omitted samples
omitSample <- c('Ski resort','2020-12-20')

##### mutations pre-existing in population
preexMuts <- c(26801)


##### the variants of interest
alVar = "B.1.1.7"
beVar = "B.1.351"
gaVar = "P.1"

#multi-base events
multiBaseBE <- list(seq(11288, 11296))
multiBaseGA <- list(seq(11288, 11296))
multiBaseAL <- list(seq(11288, 11296), seq(21765, 21770), seq(21991,21993), seq(28280, 28282))
####################################################################

### hand pick amplicons
ampliTable <- filter(ampliTable, ampliTable$ampli!="73")
ampliTable <- filter(ampliTable, !is.na(ampliTable$date))

# remove wrong B.1.1.7 mutation
dataTable <- filter(dataTable, dataTable$pos != '26801')

##### add mutations that are relevant for co-occurrence
dataTable$BR[dataTable$pos=='28877'] <- "extra"
dataTable$BR[dataTable$pos=='28878'] <- "extra"



# select data
dataTable <- selctData.function(dataTable, usedPlants, usedMutTypes, omitSample)
ampliTable <- selctData.function(ampliTable, usedPlants, usedMutTypes, omitSample)

# add variant columns
dataTable <- addVariantColData.function(dataTable, alVar, beVar, gaVar)
ampliTable <- addVariantColAmpli.function(ampliTable, alVar, beVar, gaVar)

sampleList <- unique(dataTable$sample)

#write(levels(sampleList), ncolumns=1, file = "testSampleList.tsv")

####################################################################
##        data selection                                           #
####################################################################


# add amplicon column to dataTable (amplicon table already has this info)
dataTable <- addAmpliColData.function(dataTable, alVar, beVar, gaVar)

dataTable$ga[dataTable$pos=='28878']

# aggregate multi-base events
al <- aggregateMultiBaseEvents.function(filter(dataTable, dataTable$variant==alVar), multiBaseAL)
be <- aggregateMultiBaseEvents.function(filter(dataTable, dataTable$variant==beVar), multiBaseBE)
ga <- aggregateMultiBaseEvents.function(filter(dataTable, dataTable$variant==gaVar), multiBaseGA)
dataTable <- rbind(al, be, ga)


ampliTable  <- select(ampliTable, sample, date, pos, plantname, ampli, variant, cov, var, frac)

names(dataTable)
names(ampliTable)

#filter(ampliTable, is.na(dataTable$date))

# merge tables

plot_data <- rbind(dataTable, ampliTable)
#plot_data <- select(plot_data, sample, date, pos, plantname, ampli, frac)




plot_data$variant <- factor(plot_data$variant, levels=c('B.1.351', 'P.1', 'B.1.1.7'))
plot_data$plantname <- factor(plot_data$plantname, levels=c('Ski resort', 'Lausanne (VD)', 'Zurich (ZH)'))
plot_data$date <- as.Date(plot_data$date, format="%Y-%m-%d")
#plot_data$date <- as.Date(plot_data$date, format="%m-%d")
plot_data$month <- format(plot_data$date, "%b")
plot_data$monthPrint <- as.character(plot_data$month)
plot_data$monthPrint <- ifelse(as.character(plot_data$plantname)=='Lausanne (VD)' & is.element(plot_data$monthPrint, c('Sep','Oct','Nov')), 'Sep-Nov', plot_data$monthPrint)
plot_data$monthPrint <- ifelse(as.character(plot_data$plantname)=='Zurich (ZH)' & is.element(plot_data$monthPrint, c('Jul', 'Aug','Sep', 'Oct', 'Nov')), 'Jul-Nov', plot_data$monthPrint)
plot_data$monthPrint <- factor(plot_data$monthPrint, levels=c('Jul-Nov', 'Sep-Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'))

unique(plot_data$date)
plot_data$date <- factor(plot_data$date)                       # factor dates
unique(plot_data$date)  

table(plot_data$pos) 

names(plot_data)
names(ampliTable)

#l_month_merge <- 


table(plot_data$variant)
table(plot_data$monthPrint)
table(plot_data$month)
names(plot_data)












mutOrder <- rev(c("733", "913", "1059", "2749", "3267", "3828", "5230", "5388", "5648", "5986", "6954", "10323", "11288-11296", "12778", "13860", "14408", "14676", "15279", "16176", "17259", "21614", "21621", "21638", "Amplicon 71", "21765-21770", "21801", "21974", "21991-21993", "Amplicon 72", "Amplicon 73", "22812", "22813", "23012", "23063", "Amplicon 76", "23271", "23403", "23525", "23604", "23664", "23709", "Amplicon 78", "24506", "24642", "24914", "25563", "25904", "26456", "26801", "27972", "28048", "28111", "Amplicon 92", "28167", "28280-28282", "Amplicon 93", "28512", "28877","28878" , "28887", "28977", "Amplicon 95"))

## color amplicon labels
a <- c("red", "black", "black", "black", "black", "black", "black", "black", "black", "green",
       "black", "black", "black", "black", "black", "black", "black", "black", "black", "green", 
       "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", 
       "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", 
       "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", 
       "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", 
       "black", "black", "black", "black", "black", "black", "black", "black", "black", "black")
table(a)



plot_data$pos <- factor(plot_data$pos, levels = mutOrder)

length(plot_data$pos) 


#######################################################################
#   heat map colors
#######################################################################
na_color  <- 'white'
min_color <- "#FCE9EE"
#min_color <- 'lavenderblush1'
#min_color <- "#F7EBFA"
mid_color <- "#92C5DE" 
max_color <- "#2166AC"




p <- ggplot(plot_data, aes(x = pos, y = date, fill = frac, colour = "")) +    
  geom_tile() +
  ylab("position") +
  xlab("") + 
  #facet_nested( variant+ampli ~ plantname + monthPrint, scales = "free", space='free', nest_line = TRUE) +
  facet_nested( plantname + monthPrint ~ variant+ampli, scales = "free", space='free', nest_line = TRUE) +
  
  scale_fill_gradientn(colours = c(min_color, mid_color, max_color), na.value = "white", 
                       values = c(0,0.0005,1), guide = "legend")+
  scale_colour_manual(values=NA) +
  guides(colour=guide_legend("No data", override.aes=list(fill="white", color="black")))+ 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=14)) +
  theme(axis.text.y = element_text(size=13)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  #theme(panel.spacing = unit(1, "lines"))
  theme(panel.spacing=unit(.3, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1))

p
pdf(file = "overviewTable.pdf", width = 16, height = 22)




#p
g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("#CFFCED", "gainsboro","lemonchiffon", "gainsboro","#CFFCD6","gainsboro","lemonchiffon", "gainsboro", "lemonchiffon", "gainsboro", "lemonchiffon",  "#E9FCCF", "gainsboro", "lemonchiffon", "gainsboro", "lemonchiffon", "gainsboro", "lemonchiffon", "lemonchiffon", "gainsboro", "#ffd99f","#CFD6FC", "#ffcb9f", "#CFE2FC", "#ffd99f",  "#ffe39f", "#fff19f", "#ffcb9f", "#CFEDFC", "#ffd99f",  "#ffe39f", "#fff19f")
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

dev.off()

