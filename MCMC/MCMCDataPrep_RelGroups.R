
library(tidyverse)
## Set wd 
#setwd()

## Load data
datahr <- data_file 
datahr <- datahr %>% filter(Group != "non-pMSM" & year > 2014) 
datahr<- datahr %>% mutate(GroupUse = if_else(Group =="pMSM", 1, 0)) 
table(datahr$GroupUse)
datahr$year<- as.numeric(datahr$year)
datahr$GroupUse[is.na(datahr$GroupUse)] <- 0

## Add label bucket 
datahr$rel_type = NA
datahr$rel_type = rep('0', length(datahr$GroupUse))

datahr$rel_type[which(is.na(match(datahr$GroupUse, 1)) == F)] = 'Non-Ref'
datahr$rel_type[which(is.na(match(datahr$GroupUse, 0)) == F)] = 'Ref'
table(datahr$rel_type)
datahr$count = 1

################################################################################
## Prepare data: each group
################################################################################
nb_clades = 2 
nb_years = max(datahr$year)-min(datahr$year)+1
nb_groups = length(levels(as.factor(datahr$count)))
groups = levels(as.factor(datahr$count))
first_year = min(datahr$year)

ref_clade = 1 ## Ref group

################################################################################

################################################################################
## Size arrays
################################################################################
clade_number_array = array(0, dim = c(nb_clades, nb_years, nb_groups))
clade_number_array_freq_ref = array(0, dim = c(nb_clades-1, nb_years, nb_groups))
clade_number_array_freq_non_zeros = array(0, dim = c(nb_clades, nb_years, nb_groups))
clade_number_array_freq_non_zeros_ref = array(0, dim = c(nb_clades-1, nb_years, nb_groups))

non_zero_group_year = matrix(0, nb_years, nb_groups)
non_zero_group_year_genotype = array(0, dim = c(nb_clades-1, nb_years, nb_groups))
clade_number_array_freq_first_non_zeros = rep(0, nb_groups)
total_number = matrix(0, nrow = nb_groups, ncol = nb_years)

################################################################################

################################################################################
## Fill with Shigella data
################################################################################
for(i in 1:nb_groups){
  for(j in 1:nb_years){
    clade_number_array[,j,i] = table(factor(datahr$rel_type[which(datahr$year == first_year+(j-1) & datahr$count == groups[i])], levels = c('Ref', 'Non-Ref')))
  }
}

## Compute freq with respect to a ref clade
for(kkk in 1:nb_groups){
  tmp = 1
  for(i in 1:nb_clades){
    if(i != ref_clade){
      clade_number_array_freq_ref[tmp,,kkk] = clade_number_array[i,,kkk]/(clade_number_array[ref_clade,,kkk]+clade_number_array[i,,kkk]);
      tmp = tmp+1
    }
  }
}
clade_number_array_freq_ref[which(is.na(clade_number_array_freq_ref) == T)] = 0

## Counts, per group, per year, for the ref clade
clade_number_array_ref = clade_number_array[ref_clade,,]
clade_number_array_ref1 <- t(t(clade_number_array_ref))
clade_number_array_ref = array(clade_number_array_ref, dim=c(nb_clades-1, nb_years, nb_groups))

## Counts, per group, per year, per clade (without ref clade)
clade_number_array_paired = clade_number_array[-ref_clade,,]
clade_number_array_paired = array(clade_number_array_paired, dim=c(nb_clades-1, nb_years, nb_groups))

## Number of sequences per group
for(kkk in 1:nb_groups){
  total_number[kkk,] =  apply(clade_number_array[,,kkk], MARGIN = 2, sum)
}

## Which group-year-clades are not 0s (used to not compute the likelihood on those points later on)
for(kkk in 1:nb_groups){
  for(i in 1:nb_years){
    non_zero_group_year_genotype[which(clade_number_array_paired[,i,kkk]+clade_number_array_ref[,i,kkk] > 0),i,kkk] = 1
  }
}

dim(non_zero_group_year_genotype)

## Which group-year are not 0s
for(kkk in 1:nb_groups){
  for(i in 1:nb_years){
    if(sum(c(clade_number_array_paired[,i,kkk], clade_number_array_ref[,i,kkk])) > 0){
      non_zero_group_year[i, kkk] = 1
    }
  }
}



## Year of f0 introduction, per group
f0_introduction = rep(2015, nb_groups) - first_year + 1


data.MCMCreltrav = list(nb_genotypes = nb_clades,
                        nb_years = nb_years,
                        nb_groups = nb_groups,
                        data_genotype_non_ref = clade_number_array_paired,
                        data_genotype_ref = clade_number_array_ref1,
                        data_total_number = t(total_number),
                        non_zero_group_year = non_zero_group_year,
                        non_zero_group_year_genotype = non_zero_group_year_genotype,
                        number_zeros_group_year = length(which(non_zero_group_year == 0)),
                        number_zeros_group_year_genotype = length(which(non_zero_group_year_genotype == 0)),
                        
                        f0_introduction = f0_introduction)


################################################################################
## Save data 
#saveRDS(data.MCMCreltrav, 'data.MCMCreltrav.rds')

################################################################################
################################################################################


################################################################################
################################################################################

data <- data_file 
data <- data %>% filter(Group != "pMSM" & year > 2014) 
data<- data %>% mutate(GroupUse = if_else(Group =="non-pMSM", 1, 0)) 

data$year<- as.numeric(data$year)
data$GroupUse[is.na(data$GroupUse)] <- 0
## Add label bucket 
data$rel_type = NA
data$rel_type = rep('0', length(data$GroupUse))

data$rel_type[which(is.na(match(data$GroupUse, 1)) == F)] = 'Non-Ref'
data$rel_type[which(is.na(match(data$GroupUse, 0)) == F)] = 'Ref'

data$count = 1

################################################################################
## Prepare data: each group
################################################################################
nb_clades = 2 
nb_years = max(data$year)-min(data$year)+1
nb_groups = length(levels(as.factor(data$count)))
groups = levels(as.factor(data$count))
first_year = min(data$year)

ref_clade = 1 ## Ref group

################################################################################
## Size arrays
################################################################################
clade_number_array = array(0, dim = c(nb_clades, nb_years, nb_groups))
clade_number_array_freq_ref = array(0, dim = c(nb_clades-1, nb_years, nb_groups))
clade_number_array_freq_non_zeros = array(0, dim = c(nb_clades, nb_years, nb_groups))
clade_number_array_freq_non_zeros_ref = array(0, dim = c(nb_clades-1, nb_years, nb_groups))

non_zero_group_year = matrix(0, nb_years, nb_groups)
non_zero_group_year_genotype = array(0, dim = c(nb_clades-1, nb_years, nb_groups))
clade_number_array_freq_first_non_zeros = rep(0, nb_groups)
total_number = matrix(0, nrow = nb_groups, ncol = nb_years)

################################################################################

################################################################################
## Fill with Shigella data
################################################################################
for(i in 1:nb_groups){
  for(j in 1:nb_years){
    clade_number_array[,j,i] = table(factor(data$rel_type[which(data$year == first_year+(j-1) & data$count == groups[i])], levels = c('Ref', 'Non-Ref')))
  }
}

## Compute freq with respect to a ref clade
for(kkk in 1:nb_groups){
  tmp = 1
  for(i in 1:nb_clades){
    if(i != ref_clade){
      clade_number_array_freq_ref[tmp,,kkk] = clade_number_array[i,,kkk]/(clade_number_array[ref_clade,,kkk]+clade_number_array[i,,kkk]);
      tmp = tmp+1
    }
  }
}
clade_number_array_freq_ref[which(is.na(clade_number_array_freq_ref) == T)] = 0
## Counts, per group, per year, for the ref clade

clade_number_array_ref = clade_number_array[ref_clade,,]
clade_number_array_ref1 <- t(t(clade_number_array_ref))
clade_number_array_ref = array(clade_number_array_ref, dim=c(nb_clades-1, nb_years, nb_groups))

## Counts, per group, per year, per clade (without ref clade)
clade_number_array_paired = clade_number_array[-ref_clade,,]
clade_number_array_paired = array(clade_number_array_paired, dim=c(nb_clades-1, nb_years, nb_groups))

## Number of sequences per group
for(kkk in 1:nb_groups){
  total_number[kkk,] =  apply(clade_number_array[,,kkk], MARGIN = 2, sum)
}

## Which group-year-clades are not 0s (used to not compute the likelihood on those points later on)
for(kkk in 1:nb_groups){
  for(i in 1:nb_years){
    non_zero_group_year_genotype[which(clade_number_array_paired[,i,kkk]+clade_number_array_ref[,i,kkk] > 0),i,kkk] = 1
  }
}

dim(non_zero_group_year_genotype)

## Which group-year are not 0s
for(kkk in 1:nb_groups){
  for(i in 1:nb_years){
    if(sum(c(clade_number_array_paired[,i,kkk], clade_number_array_ref[,i,kkk])) > 0){
      non_zero_group_year[i, kkk] = 1
    }
  }
}


## Year of f0 introduction, per group
f0_introduction = rep(2015, nb_groups) - first_year + 1


data.MCMCrelnon = list(nb_genotypes = nb_clades,
                       nb_years = nb_years,
                       nb_groups = nb_groups,
                       data_genotype_non_ref = clade_number_array_paired,
                       data_genotype_ref = clade_number_array_ref1,
                       data_total_number = t(total_number),
                       non_zero_group_year = non_zero_group_year,
                       non_zero_group_year_genotype = non_zero_group_year_genotype,
                       number_zeros_group_year = length(which(non_zero_group_year == 0)),
                       number_zeros_group_year_genotype = length(which(non_zero_group_year_genotype == 0)),
                       f0_introduction = f0_introduction)


################################################################################
## Save data 
#saveRDS(data.MCMCrelnon, 'data.MCMCrelnon.rds')

################################################################################
################################################################################
