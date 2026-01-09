########################################################################################################
## Library
########################################################################################################
library(ape)
library(BactDating)

########################################################################################################
## shigella sonnei tree
########################################################################################################
tree = read.tree("Ssonnei_tree.nwk")
tree1= ladderize(multi2di(phy = tree, random = F))
datevec<- datesdat$Date
names(datevec)<- datesdat$Accession

genome_length = 4971117
alignment_length = 32294 
tree1$edge.length = tree1$edge.length * alignment_length

## Time the tree 
res = bactdate(tree1, datevec, showProgress = T, 
               initMu = 4.6e-08*genome_length, updateMu = T, nbIts = 1E8)

# plot(ladderize(res$tree,F), show.tip.label = F, tip.color = 'red')
# axisPhylo(backward = F)

## Save results
saveRDS(object = res, file = 'ssonneigeospatialbacdate030524.rds')
write.tree(res$tree, file = 'HPCssonneigeospatial030524.rds.tree')

## Look at ESS values, and print them
idx = which(!is.na(colnames(res$record)))
idx_names = colnames(res$record)[idx]
ESS_val = data.frame('Var' = idx_names, 'ESS' = rep(NA, length(idx_names)))
for(i in 1:length(idx)){ESS_val[i,2] = coda::effectiveSize(res$record[,idx[i]])}
print(ESS_val)

########################################################################################################