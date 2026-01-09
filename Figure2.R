library(ggplot2)
library(sf)
library(tidyverse)
library(viridis)
library(MASS)
library(RColorBrewer)
library(geodist)
library(lineup2)


#~~~~~~~~~~~~~~~~~~SPATIAL DISTANCE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lat_vec <- hptgeocode3$Lat
long_vec <- hptgeocode3$Long
sra_vec <- hptgeocode3$`Accession`
x <- cbind(long_vec,lat_vec)
y <- x

#pariwise spatial distance matrix in meters
pairwise_geodist <- data.frame(geodist(x,y,paired = FALSE,measure = "haversine") /1000)
colnames(pairwise_geodist) <- sra_vec
rownames(pairwise_geodist) <- sra_vec
pairwise_geodist<- round(pairwise_geodist, 4)
diag(pairwise_geodist)<-NA


################prepare matrices for analysis 

hptdat <- hptgeocode3 

#pMSM
pmsmall <- hptgeocode3 %>% filter(Group == "pMSM")
idsmsm <- unique(pmsmall$Accession)
treerowcolfilter1msmbact <- bacdist.mat[idsmsm, idsmsm]
treerowcolfilter1msmspat <- pairwise_geodist[idsmsm, idsmsm]

#non-pMSM
nonpmsmall <- hptdat %>% filter(Group == "non-pMSM")
idsnonmsm <- unique(nonpmsmall$Accession)
treerowcolfilter1msmbactnon <- bacdist.mat[idsnonmsm, idsnonmsm]
treerowcolfilter1msmspatnon <- pairwise_geodist[idsnonmsm, idsnonmsm]

#domestic groups 
nonandmsmpmsmall <- hptdat %>% filter(Group != "Travel")
idsnonmsmandmsm <- unique(nonandmsmpmsmall$Accession)
treerowcolfilter1msmbactnonandmsm <- bacdist.mat[idsnonmsmandmsm, idsnonmsmandmsm]

##time matrices (simple indicator of all pairs)
time_matmsm = matrix(1, nrow=nrow(treerowcolfilter1msmspat), ncol =ncol(treerowcolfilter1msmspat))
diag(time_matmsm) <- NA
time_matnonmsm = matrix(1, nrow= nrow(treerowcolfilter1msmspatnon), ncol =  ncol(treerowcolfilter1msmspatnon))
diag(time_matnonmsm) <- NA
time_matmsmandnon = matrix(1, nrow=nrow(treerowcolfilter1msmbactnonandmsm), ncol =ncol(treerowcolfilter1msmbactnonandmsm))
diag(time_matmsmandnon) <- NA
colnames(time_matmsmandnon)<- rownames(time_matmsmandnon)<-colnames(treerowcolfilter1msmbactnonandmsm)


#date receipt matricies
dates_msm <- outer(as.Date(pmsmall$Receipt.Date), as.Date(pmsmall$Receipt.Date), 
                   FUN = function(x, y) abs(as.numeric(x - y)))
diag(dates_msm) <- NA
colnames(dates_msm)<-rownames(dates_msm)<- pmsmall$Accession
dates_nonmsm <- outer(as.Date(nonpmsmall$Receipt.Date), as.Date(nonpmsmall$Receipt.Date), 
                      FUN = function(x, y) abs(as.numeric(x - y)))
diag(dates_nonmsm) <- NA
colnames(dates_nonmsm)<-rownames(dates_nonmsm)<- nonpmsmall$Accession


#aligning matrices 
#pMSM
msmalignedrowsall<- align_matrix_rows(treerowcolfilter1msmspat, treerowcolfilter1msmbact)
spatrowstest <- as.matrix(msmalignedrowsall[[1]])
evolrowstest <- as.matrix(msmalignedrowsall[[2]])
msmalignedall<- align_matrix_cols(spatrowstest, evolrowstest)
spataligned <- as.matrix(msmalignedall[[1]])
evolalignedfin<- as.matrix(msmalignedall[[2]])
msmaligneddates<- align_matrix_rows(spataligned, dates_msm)
spatrows_2<- as.matrix(msmaligneddates[[1]])
datesrows <- as.matrix(msmaligneddates[[2]])
msmaligneddatescol<- align_matrix_cols(spatrows_2, datesrows)
spatalignedfin <- as.matrix(msmaligneddatescol[[1]])
datesalignedfin<- as.matrix(msmaligneddatescol[[2]])


#non-pMSM 
msmalignedrowsallnon<- align_matrix_rows(treerowcolfilter1msmspatnon, treerowcolfilter1msmbactnon)
spatrowstestnon <- as.matrix(msmalignedrowsallnon[[1]])
evolrowstestnon <- as.matrix(msmalignedrowsallnon[[2]])
msmalignedallnon<- align_matrix_cols(spatrowstestnon,evolrowstestnon)
spatalignednon <- as.matrix(msmalignedallnon[[1]])
evolalignedfinnon<- as.matrix(msmalignedallnon[[2]])
nonmsmaligneddates<- align_matrix_rows(spatalignednon, dates_nonmsm)
nonspatrows_2<- as.matrix(nonmsmaligneddates[[1]])
nondatesrows <- as.matrix(nonmsmaligneddates[[2]])
nonmsmaligneddatescol<- align_matrix_cols(nonspatrows_2, nondatesrows)
nonspatalignedfin <- as.matrix(nonmsmaligneddatescol[[1]])
nondatesalignedfin<- as.matrix(nonmsmaligneddatescol[[2]])


#######sampling seq for bootstrapping

nseqmsm = length(colnames(treerowcolfilter1msmbact))
nseqnon = length(colnames(treerowcolfilter1msmbactnon))

#Functions############

#probability of <2 years evolutionary time by space
prob_bootstrapspatcum <- function(x, geo_mat, evol_mat, time_mat, date_mat){
  #same location mat
  geo_mat.tmp = geo_mat[x,x]
  time_mat.tmp = time_mat[x,x]
  date_mat.tmp = date_mat[x,x]
  evol_mat.tmp2 = evol_mat[x,x]
  tmp = time_mat.tmp * (evol_mat.tmp2< 2)*(date_mat.tmp< (365*2)) #collected no more than two years apart
  tmp[which(tmp == 0)] = NA
  tmp2 = (time_mat.tmp)*(date_mat.tmp< (365*2)) #collected no more than two years apart
  tmp2[which(tmp2 == 0)] = NA
  a = cumsum(hist(tmp*geo_mat.tmp, breaks = c(0,Pmax,1E10), plot = F)$counts)
  c = cumsum(hist(tmp2*geo_mat.tmp, breaks = c(0,Pmax,1E20), plot = F)$counts)
  prob = (a/c)
  return(prob[-length(prob)])
}


##########################Probability of <2 years evolutionary time by space, cumulative
set.seed(111)

nboot=10
nsim=10

#set windows
Pmax = c(seq(0, 20, 1))
win = 1
Pmin = Pmax-win
Pmin[which(Pmin<0)]<-0
pmid1<-Pmax

boot.outallspatpmsm = matrix(NA, nrow = length(pmid1), ncol = nboot*nsim)
boot.outallspatnonpmsm= matrix(NA, nrow = length(pmid1), ncol = nboot*nsim)

#pMSM
for(j in 1:nsim){
  print(paste0('nsim : ', j, '/', nsim))
  geo_mat = spataligned
  evol_mat = evolalignedfin
  time_mat = time_matmsm
  date_mat = datesalignedfin
  for (i in  (1:(nboot))){
    tmp = sample(nseqmsm, nseqmsm, replace = T)
    prob=prob_bootstrapspatcum(tmp, geo_mat, evol_mat,time_mat, date_mat)
    boot.outallspatpmsm[,(j-1)*nboot + i] = prob
  }}


ci1_m_tallspatpmsm = apply(boot.outallspatpmsm, 1, quantile, probs = c(0.025, 0.975, 0.5), na.rm = T)
spatpmsmprob<- data.frame(t(ci1_m_tallspatpmsm))
spatpmsmprob$pmid<- pmid1
colnames(spatpmsmprob)[3]<- "Prob"
spatpmsmprob$Group <- "pMSM"

#non-pMSM
for(j in 1:nsim){
  print(paste0('nsim : ', j, '/', nsim))
  geo_mat = spatalignednon
  evol_mat = evolalignedfinnon
  time_mat = time_matnonmsm
  date_mat = nondatesalignedfin
  for (i in  (1:(nboot))){
    tmp = sample(nseqnon, nseqnon, replace = T)
    prob=prob_bootstrapspatcum(tmp, geo_mat, evol_mat,time_mat, date_mat)
    boot.outallspatnonpmsm[,(j-1)*nboot + i] = prob
  }}


boot.ci.m_tallspattnonpmsm = apply(boot.outallspatnonpmsm, 1, quantile, probs = c(0.025, 0.975, 0.5), na.rm = T)
spatnonpmsmprob<- data.frame(t(boot.ci.m_tallspattnonpmsm))
spatnonpmsmprob$pmid<- pmid1
colnames(spatnonpmsmprob)[3]<- "Prob"
spatnonpmsmprob$pmid<- pmid1
spatnonpmsmprob$Group <- "non-pMSM"
totspat<- rbind(spatpmsmprob, spatnonpmsmprob)

#generate plot
pprob<-ggplot(aes(y=Prob, x=pmid,colour = Group), data = totspat) + geom_line()+ scale_color_manual(values = c("#b21819", "#088978"))
pprob<- pprob +geom_ribbon(aes(ymin=totspat$X2.5., ymax=totspat$X97.5., fill = Group), linetype=2, alpha=0.1)+ scale_fill_manual(values = c("#b21819", "#088978"))
probbyspat<- pprob + theme_classic() + scale_x_continuous(name = "Pairwise Spatial Distance (km)")+labs(y="Prob of <2 evol", limits = c(0,NA),,expand = c(0, 0))+scale_y_continuous(limits = c(0,NA),expand = c(0, 0)) 

#estimate number of effective transmission chains
totspat$Chains<- 1/totspat$Prob
totspat$ChainsX2.5<- 1/totspat$X97.5.
totspat$ChainsX97.5<- 1/totspat$X2.5.


##########################################################################
#prepare population counts data from worldpop
library(sp)
library(raster)
library(terra)

worldpopdata = read.csv("ppp_GBR_2017_1km_Aggregated.csv", as.is=T)
worldpoppoints = st_as_sf(worldpopdata, coords = c("X", "Y"), crs = 4326)
worldpopraster <- rasterFromXYZ(cbind(worldpopdata[,c("X","Y","Z")]),
                                crs=4326)
rs  <- rast(worldpopraster)


############################################The average population size per km in each radius 
set.seed(111)
nboot=10
nsim=10

#set windows
Pmax = c(seq(0, 20, 1))
win = 1
Pmin = Pmax-win
Pmin[which(Pmin<0)]<-0
pmid1<-Pmax

meanpopsize<- c()

#convert to km
pmax1000<- Pmax*1000

nonasdat2usered <-hptgeocode3

#pMSM
#convert long, lat coordinates into spatial points
nonasdat2pmsm <-nonasdat2usered%>% filter(Group == "pMSM" )
nonasdat2pmsm.sf <- st_as_sf(nonasdat2pmsm,coords=c("Long","Lat"),crs=4326)
p_mrcinanalysispmsm <- st_transform(
  nonasdat2pmsm.sf,
  CRS("+proj=utm +zone=30 +datum=WGS84")
)

matemptypopmsm<-matrix(NaN,nrow=nseqmsm, ncol=nseqmsm)
colnames(matemptypopmsm) <- rownames(spataligned)
rownames(matemptypopmsm) <- rownames(spataligned)

allpopmsm <- data.frame(matemptypopmsm)
meanpopsizemsm<-c()

for(i in 1:length(pmax1000)){
  buffer <- st_buffer(p_mrcinanalysispmsm, dist = pmax1000[i]) #create buffer around each case equal to spatial max specified
  print(i)
  buffer_st <- st_transform(buffer, CRS("+init=epsg:4326")) 
  poly<-vect(buffer_st$geometry)
  PopSizevval <- terra::extract(rs,poly, fun = sum, na.rm=TRUE) #calculate population size within buffer
  allpop<- data.frame(PopSizevval$Z)
  for (j in 1:nseqmsm){ 
    allpopmsm[j,] = allpop[j,1]
  } 
  diag(allpopmsm) = NA
  #ensure matrices are aligned
  tmsmalignedrowsall<- align_matrix_rows(evolalignedfin,allpopmsm)
  evolalignedfin2 <- as.matrix(tmsmalignedrowsall[[1]])
  tmsmrowstest <- as.matrix(tmsmalignedrowsall[[2]])
  tmsmalignedall<- align_matrix_cols(evolalignedfin2,tmsmrowstest)
  tmsmrowstest2 <- as.matrix(tmsmalignedall[[2]])
  spatalignedalluse= (spataligned <= Pmax[i]) #within the spatial max specified
  spatalignedalluse[which(spatalignedalluse ==0)] = NA
  rangevalsmean<- mean(tmsmrowstest2*spatalignedalluse, na.rm=T)
  meanpopsizemsm<- append(meanpopsizemsm,rangevalsmean)
  
}

msmpopsize<- data.frame(meanpopsizemsm)
msmpopsize$Group<- "pMSM"
colnames(msmpopsize)[1]<- "meanpopsize"


#non-pMSM
#convert long, lat coordinates into spatial points
set.seed(111)
nonasdat2nonpmsm <-nonasdat2usered%>% filter(Group == "non-pMSM" )
nonasdat2nonpmsm.sf <- st_as_sf(nonasdat2nonpmsm,coords=c("Long","Lat"),crs=4326)
p_mrcinanalysisnonpmsm <- st_transform(
  nonasdat2nonpmsm.sf,
  CRS("+proj=utm +zone=30 +datum=WGS84")
)
matemptypopnonmsm<-matrix(NaN,nrow=nseqnon, ncol=nseqnon)
colnames(matemptypopnonmsm) <- rownames(spatalignednon)
rownames(matemptypopnonmsm) <- rownames(spatalignednon)

allpopnonmsm <- data.frame(matemptypopnonmsm)
meanpopsizenonmsm<-c()

for(i in 1:length(pmax1000)){
  buffer <- st_buffer(p_mrcinanalysisnonpmsm, dist = pmax1000[i])  #create buffer around each case equal to spatial max specified
  print(i)
  buffer_st <- st_transform(buffer, CRS("+init=epsg:4326")) 
  poly<-vect(buffer_st$geometry)
  PopSizevval <- terra::extract(rs,poly, fun = sum, na.rm=TRUE)  #calculate population size within buffer
  allpop<- data.frame(PopSizevval$Z)
  for (j in 1:nseqnon){ 
    allpopnonmsm[j,] = allpop[j,1]
  } 
  diag(allpopnonmsm) = NA
  #ensure matrices are aligned
  tnonmsmalignedrowsall<- align_matrix_rows(evolalignedfinnon,allpopnonmsm)
  nonevolalignedfin2 <- as.matrix(tnonmsmalignedrowsall[[1]])
  tnonmsmrows <- as.matrix(tnonmsmalignedrowsall[[2]])
  tnonmsmalignedall<- align_matrix_cols(nonevolalignedfin2,tnonmsmrows)
  tnonmsmall<- as.matrix(tnonmsmalignedall[[2]])
  spatalignedalluse= (spatalignednon <= Pmax[i])  #within the spatial max specified
  spatalignedalluse[which(spatalignedalluse ==0)] = NA
  rangevalsmean<- mean(tnonmsmall*spatalignedalluse, na.rm=T)
  meanpopsizenonmsm<- append(meanpopsizenonmsm,rangevalsmean)
  
}
nonmsmpopsize<- data.frame(meanpopsizenonmsm)
nonmsmpopsize$Group<- "non-pMSM"
colnames(nonmsmpopsize)[1]<- "meanpopsize"
nonmsmpopsize$pmid<- pmid1
nonmsmpopsize$pmax<- Pmax
msmpopsize$pmid<- pmid1
msmpopsize$pmax<- Pmax
popsizes<-rbind(nonmsmpopsize,msmpopsize)

#merge mean population size data with effective transmission chains estimates
totspatpop<- left_join(totspat, popsizes,  by = c("pmid", "Group"))

#plotting
ppop<-ggplot(aes(y=Chains, x=meanpopsize,colour = Group), data = totspatpop) + geom_line()+ scale_color_manual(values = c("#b21819", "#088978"))
ppop<- ppop +geom_ribbon(aes(ymin=totspatpop$ChainsX2.5, ymax=totspatpop$ChainsX97.5, fill = Group), linetype=2, alpha=0.1)+ scale_fill_manual(values = c("#b21819", "#088978"))
popsizesplot<- ppop + theme_classic() + scale_x_log10(name = "Mean Population Size")+labs(y="Effective Transmission Chains")+ geom_smooth(method = "lm", fill = NA)+scale_y_continuous(expand = c(0, 0))

#by area 
totspatpop$area<- pi * (totspatpop$pmax)^2 
parea<-ggplot(aes(y=Chains, x=area,colour = Group), data = totspatpop) + geom_line()+ scale_color_manual(values = c("#b21819", "#088978"))
parea<- parea +geom_ribbon(aes(ymin=totspatpop$ChainsX2.5, ymax=totspatpop$ChainsX97.5, fill = Group), linetype=2, alpha=0.1)+ scale_fill_manual(values = c("#b21819", "#088978"))
popsizesareaplot<- parea + theme_classic() + scale_x_log10(name = "Area")+labs(y="Effective Transmission Chains")+ geom_smooth(method = "lm", fill = NA)+ scale_y_continuous(expand = c(0, 0))

#linear regression of number of effective transmission chains over population size
totspatpop<- totspatpop %>% mutate(non = if_else(Group == "non-pMSM",0,1))
linearreg<- lm(Chains ~ log(meanpopsize) + non, totspatpop)
confint(linearreg)
summary(linearreg)

#########################Odds ratio of same type of transmission

#generate matrices indicating pMSM pairs or non-pMSM pairs
allpmsmdef <- hptgeocode3 %>% filter(Group != "Travel")%>% dplyr::select(Accession, Group)
allpmsmdef$Group <- as.character(allpmsmdef$Group)

matemptyalldef<-data.frame(matrix(NaN,nrow=length(allpmsmdef$Accession),ncol=length(allpmsmdef$Accession)))
colnames(matemptyalldef) <- allpmsmdef$Accession
rownames(matemptyalldef) <- allpmsmdef$Accession

t1alldef <-matemptyalldef
t21alldef <- matemptyalldef

for (j in 1:length(allpmsmdef$Accession)){ 
  t1alldef[j,] =  allpmsmdef[j,2]
  t21alldef[,j] =  allpmsmdef[j,2]
} 


#concatenate columns from two df's which will give the pairs of groups for all isolates  
t31alldef <- data.frame(Map(paste, setNames(t1alldef, paste0(names(t1alldef))), t21alldef, 
                            MoreArgs = list(sep = "_")))

colnames(t31alldef) <- colnames(matemptyalldef)
rownames(t31alldef) <- rownames(matemptyalldef)

a1alldef <- t31alldef

c11alldef <- matrix(NA,ncol=length(allpmsmdef$Accession),nrow=length(allpmsmdef$Accession))
c12alldef <- data.frame(c11alldef)

for(j in 1:ncol(a1alldef)) {
  col <- colnames(a1alldef)[j]
  sep <- NA  
  sep <- separate_wider_delim(a1alldef, col, delim = "_", names = c("iso1", "iso2"))
  sep1 <- sep %>% mutate(msmpair = if_else(iso1 == "pMSM"& iso2=="pMSM", 1,
                                           if_else(iso1== "non-pMSM"&iso2=="non-pMSM",0,
                                                   NA)))
  sep2 <- sep1 %>% dplyr::select(msmpair)
  val <- sep2[,1]
  c12alldef[,j] = val
}
diag(c12alldef)<-NA

colnames(c12alldef) <- colnames(t31alldef)
rownames(c12alldef) <- rownames(t31alldef)

pairwise_matrix_same<-c12alldef

#generate matrices different type pairs pairs
pairwise_matrix_diff <- outer(allpmsmdef$Group, allpmsmdef$Group, FUN = "==") * 1
pairwise_matrix_diff<- 1-pairwise_matrix_diff
diag(pairwise_matrix_diff)<-NA
colnames(pairwise_matrix_diff) <- allpmsmdef$Accession
rownames(pairwise_matrix_diff) <- allpmsmdef$Accession

#ensure matricies are aligned
alignedrowsmsmandnon<- align_matrix_rows(pairwise_matrix_same, treerowcolfilter1msmbactnonandmsm)
geoalignedrowsmsmandnon <- as.matrix(alignedrowsmsmandnon[[1]])
evolalignedrowsmsmandnon <- as.matrix(alignedrowsmsmandnon[[2]])
alignedmsmandnon <- align_matrix_cols(geoalignedrowsmsmandnon,evolalignedrowsmsmandnon)
geoalignedmsmandnon<- as.matrix(alignedmsmandnon[[1]])
evolalignedmsmandnon <- as.matrix(alignedmsmandnon[[2]])
alignedrowsnonandmsm<- align_matrix_rows(geoalignedmsmandnon, pairwise_matrix_diff)
geoalignedrowsmsmandnon<- as.matrix(alignedrowsnonandmsm[[1]])
geodiffalignedrowsmsmandnon <- as.matrix(alignedrowsnonandmsm[[2]])
alignedmsmandnon1 <- align_matrix_cols(geoalignedrowsmsmandnon,geodiffalignedrowsmsmandnon)
geoalignedmsmandnon1 <- as.matrix(alignedmsmandnon1[[1]])
geodiffalignedmsmandnon <- as.matrix(alignedmsmandnon1[[2]])

nseqnonandmsm<-length(colnames(pairwise_matrix_same))

#function for OR of same group type in non-overlapping windows
oddsratio_bootstrapbywindows <- function(x, geo_matsame,geomatdiff, evol_mat){
  geo_mat.tmpsame = geo_matsame[x,x]
  geo_mat.tmpdiff = geomatdiff[x,x]
  evol_mat.tmp2 = evol_mat[x,x]
  tmp = geo_mat.tmpsame 
  tmp[which(tmp == 0)] = NA
  tmpdiff = geo_mat.tmpdiff 
  tmpdiff[which(tmpdiff == 0)] = NA
  
  #same type in window
  a1 = cumsum(hist(tmp*evol_mat.tmp2, breaks = c(0,Pmax,1E10), plot = F)$counts)
  a2 = cumsum(hist(tmp*evol_mat.tmp2, breaks = c(0,Pmin,1E10), plot = F)$counts)
  a = a1 - a2
  
  #different type in window
  b1 = cumsum(hist(tmpdiff*evol_mat.tmp2, breaks = c(0,Pmax,1E10), plot = F)$counts)
  b2 = cumsum(hist(tmpdiff*evol_mat.tmp2, breaks = c(0,Pmin,1E10), plot = F)$counts)
  b = b1 - b2
  
  #same type outside window
  c1 =cumsum(hist(tmp*evol_mat.tmp2, breaks = c(0,GroupMax,1E10), plot = F)$counts)
  c =c1-a
  
  # different type outside window
  d1 =cumsum(hist(tmpdiff*evol_mat.tmp2, breaks = c(0,GroupMax,1E10), plot = F)$counts)
  d = d1-b
  
  or = (a/b)/(c/d) 
  
  return(or[-length(or)])
}


#set windows
set.seed(111)

nboot = 10
nsim=10
Pmin = c(seq(0, 17.5,2.5), 20)
Pmax<- c(NA,NA,NA,NA,NA,NA,NA,NA)
windows = 2.5-.000001
Pmax[1:8] = Pmin[1:8] + windows

#last window of evaluation is >= 20 + 
GroupMax<- max(evolalignedmsmandnon, na.rm=TRUE)
Pmax[9] = Pmin[9] + GroupMax 
Pmin[which(Pmin<0)]<-0
pmid1<-(Pmin+Pmax)/2
pmid1[9]<-21

boot.outmsmandnonodd = matrix(NA, length(pmid1), nsim*nboot)
boot.outmsmandnonnonnumodd = matrix(NA, length(pmid1), nsim*nboot)

#pMSM

for(j in 1:nsim){
  print(paste0('nsim : ', j, '/', nsim))
  geo_matmsmandnon = geoalignedmsmandnon1
  geo_matmsmandnondiff = geodiffalignedmsmandnon
  evol_matmsmandnon = evolalignedmsmandnon
  for (i in (1:(nboot))){
    tmp= sample(nseqnonandmsm, nseqnonandmsm, replace = T) # Resample isolates
    or = oddsratio_bootstrapbywindows(tmp, geo_matmsmandnon, geo_matmsmandnondiff, evol_matmsmandnon) #Compute OR
    boot.outmsmandnonodd[,(j-1)*nboot + i] = or  #save result
  }
}
boot.ci1_msmandnonevol = apply(boot.outmsmandnonodd, 1, quantile, probs = c(0.025, 0.975), na.rm = T)
boot.ci.m_msmandnonevol = apply(boot.outmsmandnonodd, 1, quantile, probs = c(0.5), na.rm = T)

#non-pMSM
nonnas1<- 1- geoalignedmsmandnon1

for(j in 1:nsim){
  print(paste0('nsim : ', j, '/', nsim))
  geo_matmsmandnon = nonnas1
  geo_matmsmandnondiff = geodiffalignedmsmandnon
  evol_matmsmandnon = evolalignedmsmandnon
  for (i in (1:(nboot))){
    tmp= sample(nseqnonandmsm, nseqnonandmsm, replace = T) ## Resample isolates
    or = oddsratio_bootstrapbywindows(tmp, geo_matmsmandnon, geo_matmsmandnondiff, evol_matmsmandnon) #Compute OR
    boot.outmsmandnonnonnumodd[,(j-1)*nboot + i] = or #save result
  }
}

boot.ci1_msmandnonevolnon = apply(boot.outmsmandnonnonnumodd, 1, quantile, probs = c(0.025, 0.975), na.rm = T)
boot.ci.m_msmandnonevolnon = apply(boot.outmsmandnonnonnumodd, 1, quantile, probs = c(0.5), na.rm = T)

#combine data for plotting
msmo <- data.frame(t(boot.ci1_msmandnonevol))
nonmsmo <- data.frame(t(boot.ci1_msmandnonevolnon))
msmo$time <- pmid1
msmo$Transmission <-"pMSM-pMSM"
nonmsmo$time <-  pmid1
nonmsmo$Transmission <-"non-pMSM-non-pMSM"
cimevolpairs <- rbind(msmo,nonmsmo)

msmm <- data.frame(boot.ci.m_msmandnonevol)
colnames(msmm)[1]<- "Odds"
msmm$time <-pmid1
msmm$Transmission <-"pMSM-pMSM"
nonmsmm <- data.frame(boot.ci.m_msmandnonevolnon)
nonmsmm$time <-pmid1
colnames(nonmsmm)[1]<- "Odds"
nonmsmm$Transmission <-"non-pMSM-non-pMSM"
mmpairs <- rbind(msmm,nonmsmm)

p1pairs<-ggplot(data=mmpairs, aes(x=time, y=Odds, colour=Transmission))+geom_hline(yintercept = 1,linetype='dashed', color='grey')+ geom_point() 
p1pairs<-p1pairs+geom_errorbar(aes(ymin=cimevolpairs$X2.5., ymax=cimevolpairs$X97.5.,fill = Transmission), width = 0.5)
odds<- p1pairs +theme_classic()+ labs(x= "Evolutionary Time (Years)", y = "Odds Ratio", title = "Odds Ratio of Within Case-Group Transmission Pairs by Evolutionary Time") + scale_color_manual(values = c("#b21819", "#088978")) + scale_y_log10(limits = c(0.1,30))

####################################################################Mean spatial distance by evolutionary time 
nsim=10
nboot=10

#set windows
Pmax = seq(0,100,0.25)
windows = 0.25
Pmin = Pmax - windows
Pmin[which(Pmin<0)]<-0
pmid1<- Pmax

boot.outmsmspat = matrix(NA, length(pmid1), nsim*nboot)
boot.outmsmspatnon = matrix(NA, length(pmid1), nsim*nboot)

#pMSM
for(k in 1:nsim){
  print(paste0('K : ', k, '/', nsim))
  for(j in 1:length(Pmax)){
    evol_mat2 = evolalignedfin
    evol_mat2.tmp = (evol_mat2<=Pmax[j]) # Select only pairs that verify the evol distance
    evol_mat2.tmp[which(evol_mat2.tmp==0)]<-NA
    for(i in 1:nboot){
      tmp = sample(nseqmsm, nseqmsm, replace = T) # Resample isolates
      mean.out = mean(spataligned[tmp,tmp]*evol_mat2.tmp[tmp,tmp], na.rm = T) # Compute mean geo distance
      boot.outmsmspat[j,(k-1)*nboot + i] = mean.out # Store result
    }
  }
}
boot.ci1_msmspat = apply(boot.outmsmspat, 1, quantile, probs = c(0.025, 0.975), na.rm = T)
boot.ci.m_msmspat = apply(boot.outmsmspat, 1, quantile, probs = c(0.5), na.rm = T)

#non-pMSM
for(k in 1:nsim){
  print(paste0('K : ', k, '/', nsim))
  for(j in 1:length(Pmax)){
    evol_mat2 = evolalignedfinnon
    evol_mat2.tmp = (evol_mat2<=Pmax[j]) # Select only pairs that verify the evol distance
    evol_mat2.tmp[which(evol_mat2.tmp==0)]<-NA
    for(i in 1:nboot){
      tmp = sample(nseqnon, nseqnon, replace = T) # Resample isolates
      mean.out = mean(spatalignednon[tmp,tmp]*evol_mat2.tmp[tmp,tmp], na.rm = T) # Compute mean geo distance
      boot.outmsmspatnon[j,(k-1)*nboot + i] = mean.out# Store result
    }
  }
}

boot.ci1_nonspat = apply(boot.outmsmspatnon, 1, quantile, probs = c(0.025, 0.975), na.rm = T)
boot.ci.m_nonspat = apply(boot.outmsmspatnon, 1, quantile, probs = c(0.5), na.rm = T)


#prepare data for plotting
msmspat <- data.frame(t(boot.ci1_msmspat))
nonspat <- data.frame(t(boot.ci1_nonspat))
msmspat$time <-  pmid1
msmspat$Group <-"pMSM"
nonspat$time <- pmid1
nonspat$Group <-"non-pMSM"
cimevolspat <- rbind(msmspat,nonspat)

msmspatci <- data.frame(boot.ci.m_msmspat)
colnames(msmspatci)[1]<- "MeanSpat"
msmspatci$time <-  pmid1
msmspatci$Group <-"pMSM"
nonspatci <- data.frame(boot.ci.m_nonspat)
nonspatci$time <-  pmid1
colnames(nonspatci)[1]<- "MeanSpat"
nonspatci$Group <-"non-pMSM"
mmspat <- rbind(msmspatci,nonspatci)

mmspatunder20<- mmspat %>% filter(time<=20)
cimevolspatunder20<- cimevolspat %>% filter(time<=20)

#plot
p20<-ggplot(data=mmspatunder20, aes(x=time, y=MeanSpat, colour=Group)) + geom_line()+ scale_color_manual(values = c("#b21819", "#088978"))
p20<-p20+geom_ribbon(aes(ymin=cimevolspatunder20$X2.5., ymax=cimevolspatunder20$X97.5.,fill = Group), linetype=2, alpha=0.1)+ scale_fill_manual(values = c("#b21819", "#088978"))+scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
spatevol<- p20 +theme_classic()+ labs(x= "Evolutionary Time (Years)", y = "Mean Cumulative Spatial Distance (km)")

p100<-ggplot(data=mmspat, aes(x=time, y=MeanSpat, colour=Group)) + geom_line()+ scale_color_manual(values = c("#b21819", "#088978"))
p100<-p100+geom_ribbon(aes(ymin=cimevolspat$X2.5., ymax=cimevolspat$X97.5.,fill = Group), linetype=2, alpha=0.1)+ scale_fill_manual(values = c("#b21819", "#088978"))
spatevol<- p100 +theme_classic()+ labs(x= "Evolutionary Time (Years)", y = "Mean Cumulative Spatial Distance (km)")
