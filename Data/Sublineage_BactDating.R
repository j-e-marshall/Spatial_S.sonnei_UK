library(ape)
library(BactDating)

#read in data
tree <- ape::read.tree("XDR_English_isolates_only.treefile")

#prepare tree
tree1= ladderize(multi2di(phy = tree, random = F))
genome_length = 4971117
alignment_length = 924
tree1$edge.length = tree1$edge.length * alignment_length

#metadata <- read_excel("metadata.xlsx")

#convert to data format
#datesfile<- metadata %>% dplyr::select(Accession, `Date of isolation`) %>% mutate(DateUse = as.Date(`Date of isolation`, format="%Y/%m/%d"))

#convert to decimal form
#datesfile<- datesfile %>% mutate(Date =decimal_date(datesfile$DateUse))

#reduce data set
#datesfiledec <- datesfile %>% dplyr::select(Accession, Date)

#prepare date data
datevec<- datesfiledec$Date
names(datevec)<- datesfiledec$Accession

## Time the tree
res1 = bactdate(tree1, datevec, showProgress = T,
                initMu = 4.6e-08*genome_length, updateMu = T, nbIts = 5E8)


## Save results
saveRDS(object = res, file = 'bact_object.rds')
write.tree(res1$tree, file = 'timedtree.rds.tree')
