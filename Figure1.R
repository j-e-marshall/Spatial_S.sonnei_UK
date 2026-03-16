library(binom)
library(ggtree)
library(tidyverse)
library(ape)
library(TreeTools)
devtools::install_github("xavierdidelot/BactDating")
library(BactDating)
library(coda)
library(posterior)
library(bayesplot)

#read in base data file
hptgeocode3<-  data

#~~~~~~~~~~~~~~~~~~EVOLUTIONARY TIME ESTIMATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#convert to decimals for bactdating 
#alldates <- hptgeocode3 %>% dplyr::select(Accession, Receipt.Date)
#alldates$Date <- decimal_date(alldates$Receipt.Date)
#alldatesdec <- alldates %>% dplyr::select(Accession, Date)

#read in tree
bigtree<- read.tree("phylogenetictree.tree")

#adjust length of branches 
genome_length = 4971117
alignment_length = 32294 
bigtree$edge.length = bigtree$edge.length * alignment_length

#save copy of isolate data set
groups <- hptgeocode3

#loading date data used in bactdating
datesdat <- alldatesdec

#creating date vector
datevec<- datesdat$Date
names(datevec)<- datesdat$Accession

#root to tip linear regression function -- modified to adjust axis of figure (zoomed in on period from 1980-2020)
roottotipfilt = function(tree,date,rate=NA,permTest=10000,showFig=T,colored=T,showPredInt='gamma',showText=T,showTree=T)
{
  if (!is.rooted(tree)) warning('Warning: roottotip was called on an unrooted input tree. Consider using initRoot first.\n')
  if (sum(tree$edge.length)<5) warning('Warning: input tree has small branch lengths. Make sure branch lengths are in number of substitutions (NOT per site).\n')
  #Rerranging of dates, if needed
  if (!is.null(names(date))) date=findDates(tree,date)
  
  if (var(date,na.rm=T)==0 && is.na(rate)) {warning('Warning: All dates are identical.\n');return(list(rate=NA,ori=NA,pvalue=NA))}
  n=length(date)
  ys=leafDates(tree)
  if (is.na(rate)) {
    res=lm(ys~date)
  }
  else {
    res=lm(I(ys-rate*date)~1)
    res$coefficients=c(res$coefficients,rate)
  }
  ori=-coef(res)[1]/coef(res)[2]
  rate=coef(res)[2]
  r2=summary(res)$r.squared
  correl=cor(date,ys,use='complete.obs')
  #pvalue=summary(res)$coefficients[,4][2]
  #print(c(r2,correl^2))#Equal
  
  pvalue=0
  for (i in 1:permTest) {
    date2=sample(date,n,replace=F)
    correl2=cor(date2,ys,use='complete.obs')
    if (correl2>=correl) pvalue=pvalue+1/permTest
  }
  
  if (rate<0) {warning('The linear regression suggests a negative rate.')}
  if (showFig==F) return(list(rate=rate,ori=ori,pvalue=pvalue))
  old.par=par(no.readonly = T)
  par(xpd=NA,oma = c(0, 0, 2, 0))
  if (colored) {
    normed=(date-min(date,na.rm=T))/(max(date,na.rm=T)-min(date,na.rm=T))
    cols=grDevices::rgb(ifelse(is.na(normed),0,normed),ifelse(is.na(normed),0.5,0),1-ifelse(is.na(normed),1,normed),0.5)
  } else cols='black'
  if (showTree) {
    par(mfrow=c(1,2))
    plot(tree,show.tip.label = F)
    if (colored) tiplabels(col=cols,pch=19)
    axisPhylo(1,backward = F)
  }
  plot(date,ys,col=cols,xlab=ifelse(showText,'Sampling date',''),ylab=ifelse(showText,'Root-to-tip distance',''),xaxs='i',yaxs='i',pch=19,ylim=c(0,max(ys)),xlim=c(ifelse(rate>0,1980,1980),max(date,na.rm = T)))
  #text(date,ys,labels=1:length(date))
  par(xpd=F)
  abline(res,lwd=2)
  if (rate<0) {par(old.par);return(list(rate=rate,ori=ori,pvalue=pvalue))}
  xs=seq(ori,max(date,na.rm = T),0.1)
  plim=0.05
  if (showPredInt=='poisson') {
    lines(xs,qpois(  plim/2,(xs-ori)*rate),lty='dashed')
    lines(xs,qpois(1-plim/2,(xs-ori)*rate),lty='dashed')
  }
  if (showPredInt=='gamma') {
    lines(xs,qgamma(  plim/2,shape=(xs-ori)*rate,scale=1),lty='dashed')
    lines(xs,qgamma(1-plim/2,shape=(xs-ori)*rate,scale=1),lty='dashed')
  }
  if (showText) {
    if (pvalue==0) mtext(sprintf('Rate=%.2e,MRCA=%.2f,R2=%.2f,p<%.2e',rate,ori,r2,1/permTest), outer = TRUE, cex = 1.5)
    else           mtext(sprintf('Rate=%.2e,MRCA=%.2f,R2=%.2f,p=%.2e',rate,ori,r2,pvalue), outer = TRUE, cex = 1.5)
  }
  par(old.par)
  return(list(rate=rate,ori=ori,pvalue=pvalue))
}

#' Initial tree rooting based on best root-to-tip correlation
#' @param phy An unrooted phylogenetic tree
#' @param date Dates of sampling
#' @param mtry Average number of rooting attempts per branch
#' @param useRec Whether or not to use results from previous recombination analysis
#' @return Rooted tree
#' @export
initRoot = function(phy,date,mtry=10,useRec=F) {
  #Rerranging of dates, if needed
  if (!is.null(names(date))) date=findDates(phy,date)
  
  n=length(date)
  bestcorrel=-Inf
  denom=mean(phy$edge.length)
  for (w in c(1:Ntip(phy),Ntip(phy)+(2:Nnode(phy)))) {
    if (w<=Ntip(phy)) tree=root(phy,outgroup=w,resolve.root = T) else tree=root(phy,node=w,resolve.root = T)
    wi=which(tree$edge[,1]==Ntip(tree)+1)
    toshare=sum(tree$edge.length[wi])
    attempts=max(1,ceiling(mtry*toshare/denom))
    for (a in 1:attempts) {
      tree$edge.length[wi]=toshare*c(a,attempts+1-a)/(attempts+1)
      #ys=leafDates(tree)
      ys=allDates(tree)[1:n]#This is faster
      correl=suppressWarnings(cor(date,ys,use='complete.obs'))
      if (is.na(correl)) correl=-Inf
      if (correl>bestcorrel) {bestcorrel=correl;best=c(w,a,attempts)}
    }
  }
  if (correl==-Inf) {#This happens for example if all dates are identical
    best=c(1,1,1)
  }
  
  w=best[1];a=best[2];attempts=best[3]
  if (useRec==F) {
    #Rooting without recombination
    if (w<=Ntip(phy)) tree=root(phy,outgroup=w,resolve.root = T) else tree=root(phy,node=w,resolve.root = T)
    wi=which(tree$edge[,1]==Ntip(tree)+1)
    tree$edge.length[wi]=sum(tree$edge.length[wi])*c(a,attempts+1-a)/(attempts+1)
    
  } else {
    #Rooting with recombination - need to be careful to keep correct unrec values on correct branch
    phy$node.label=sprintf('n%d',1:Nnode(phy))
    edgenames=cbind(c(phy$tip.label,phy$node.label)[phy$edge[,1]],c(phy$tip.label,phy$node.label)[phy$edge[,2]])
    unrec=phy$unrec
    unrecbest=unrec[which(phy$edge[,2]==w)]
    if (w<=Ntip(phy)) tree=root(phy,outgroup=w,resolve.root=T,edgelabel=F) else tree=root(phy,node=w,resolve.root=T,edgelabel=F)
    wi=which(tree$edge[,1]==Ntip(tree)+1)
    tree$edge.length[wi]=sum(tree$edge.length[wi])*c(a,attempts+1-a)/(attempts+1)
    tree$unrec=rep(NA,nrow(tree$edge))
    tree$unrec[wi]=unrecbest
    for (i in 1:nrow(tree$edge)) {
      nams=c(tree$tip.label,tree$node.label)[tree$edge[i,]]
      for (j in 1:nrow(edgenames)) if (setequal(nams,edgenames[j,])) {tree$unrec[i]=unrec[j];break}
    }
  }
  return(tree)
}

#' Compute dates of leaves for a given tree
#' @param phy Tree
#' @return Dates of leaves
#' @export
leafDates = function (phy) {
  rootdate=phy$root.time
  if (is.null(rootdate)) rootdate=0
  nsam=length(phy$tip.label)
  dates=rep(rootdate,nsam)
  for (i in 1:nsam) {
    w=i
    while (1) {
      r=which(phy$edge[,2]==w)
      if (length(r)==0) break
      dates[i]=dates[i]+phy$edge.length[r]
      w=phy$edge[r,1]
    }
  }
  return(dates)
}

#' Compute dates of leaves and internal nodes for a given tree
#' @param phy Tree
#' @return Dates of leaves and internal nodes
#' @export
allDates = function (phy) {
  rootdate=phy$root.time
  if (is.null(rootdate)) rootdate=0
  return(rootdate+unname(dist.nodes(phy)[Ntip(phy)+1,]))
  #o=rev(postorder(phy))#preorder
  #n=Ntip(phy)+Nnode(phy)
  #dates=rep(NA,n)
  #dates[Ntip(phy)+1]=rootdate
  #for (i in o) {
  #  dates[phy$edge[i,2]]=dates[phy$edge[i,1]]+phy$edge.length[i]
  #}
  #return(dates)
}

#' Compute dates of internal nodes for a given tree
#' @param phy Tree
#' @return Dates of internal nodes
#' @export
nodeDates = function (phy) {
  rootdate=phy$root.time
  if (is.null(rootdate)) rootdate=0
  #return(rootdate+dist.nodes(phy)[Ntip(phy)+1,])#This is not faster
  nsam=length(phy$tip.label)
  dates=rep(rootdate,nsam-1)
  for (i in 2:(nsam-1)) {
    w=i+nsam
    while (1) {
      r=which(phy$edge[,2]==w)
      if (length(r)==0) break
      dates[i]=dates[i]+phy$edge.length[r]
      w=phy$edge[r,1]
    }
  }
  return(dates)
}

findDates = function(tree,dates)
{
  date2=rep(NA,Ntip(tree))
  for (i in 1:Ntip(tree)) {
    wi=which(names(dates)==tree$tip.label[i])
    if (length(wi)>1) wi=wi[1]
    if (length(wi)==1) date2[i]=dates[wi]
  }
  return(date2)
}

#root to tip on full tree
#pdf("roottotipfull_plot.pdf", width=10, height=5)
roottotip(bigtree,date=datevec)
dev.off()

#root to tip post 1980
#pdf("roottotipfilt_plot.pdf", width=10, height=5)
roottotipfilt(bigtree,date=datevec, showTree = F)
dev.off()


#load bacdating object (run on cluster)
dat = readRDS(file = 'timephylogenetictree.rds')

#trace plots
plot(dat,'trace')

#check effective sample size of parameters
mcmc=as.mcmc.resBactDating(dat, burnin = 0.2)
summary(mcmc)
effectiveSize(mcmc)

#median and 95% HPD mutation rate
3.1744/4971117
3.0873/4971117
3.2732/4971117

#posterior trace plots of paramters
mcmcdraws <- as_draws(mcmc)
mcmc_trace(mcmcdraws)

#pairwise evolutionary time distances
bacdist.mat<-cophenetic.phylo(dat$tree)
diag(bacdist.mat) <- NA

#labeling cases that are in analysis 
groups$inanalysis<-1

#joining list of cases to tree tip list
comboacc <- left_join(listoftips, groups, by = "Accession")

#generating labels for tree
groups2 <-groups %>% dplyr::select("Accession", "Group","inanalysis")

#filtering for only isolates within the tree
groups2combo <- left_join(listoftips,groups2, by = "Accession")

#the tip to drop are those not in the analysis but are in the tree 
groups2combononadrop<- groups2combo%>% filter(is.na(inanalysis))

#dropping cases in the tree that are not in the analysis 
notinanalysisids <- unique(groups2combononadrop$Accession)

#drop tips we don't have data for
inanalysistreeuse <- drop.tip(dat$tree, notinanalysisids)

#plot tree figure
groups2combofilt<- groups2combo %>% filter(!is.na(inanalysis))
groups2combouse<- groups2combofilt %>% dplyr::select(Group)
rownames(groups2combouse)<- groups2combofilt$Accession

#make tree figure
p <- ggtree(inanalysistreeuse,size =.05)+ 
  theme_tree2() 
heatmap.colours <- c( "#b21819","#088978", "#f85c06")
p1<- gheatmap(p, groups2combouse, offset = 10, color=NULL, 
              colnames_position="top", 
              colnames_angle=90, colnames_offset_y = 1, 
              hjust=0, font.size=2, width=0.1)+ scale_fill_manual(values=heatmap.colours)

#Figure 1A
totbucketsinc2004years <- as.data.frame(hptgeocode3)  %>% group_by(Group, year) %>%
  summarize(cases = n())
totbucketsunder2008years<-totbucketsinc2004years %>%filter(year<2008)
totbucketsinc2004years$year <- factor(totbucketsinc2004years$year)

plot1 <- ggplot(data=totbucketsinc2004years, aes(x=year, y=cases, fill = Group)) + theme_classic() + geom_bar(position="stack", stat="identity") + scale_y_continuous("Isolates (n)", expand = c(0, 0)) + scale_fill_manual(values = c("#b21819", "#088978","#f85c06")) +labs(x="Year")+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot2 <- ggplot(data=totbucketsunder2008years, aes(x=year, y=cases, fill = Group)) + theme_classic() + geom_bar(position="stack", stat="identity") + scale_y_continuous("Isolates (n)", expand = c(0, 0)) + scale_fill_manual(values = c("#088978","#f85c06","#b21819")) +labs(x="Year")+ theme(panel.border = element_blank(), panel.grid.major = element_blank())


#Supp figs.
hptgeocode3 <- hptgeocode3 %>% mutate(pMSMfinal = if_else(Group == "pMSM",1,0), highrisktravel = if_else(Group == "Travel",1,0),domesticnonpmsm = if_else(Group == "non-pMSM",1,0))
hptgeocode3_sum <- hptgeocode3 %>%
  group_by(year) %>%
  summarise( pMSM = sum(pMSMfinal),Travel = sum(highrisktravel),nonpMSM = sum(domesticnonpmsm)) %>%
  pivot_longer(-year, names_to = "Group", values_to = "Tot") 

hptgeocode3_sum_proplong <- hptgeocode3_sum %>%
  group_by(year) %>%
  mutate(prop = Tot / sum(Tot))

propcases<- ggplot(hptgeocode3_sum_proplong, aes(x = factor(year), y = prop, fill = Group)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Year",
    y = "Proportion of Cases",
    fill = "Group"
  ) + theme_minimal() + scale_fill_manual(values = c( "#b21819","#088978","#f85c06"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#relative proportion of Cip, CRO, AZM, XDR
hptgeocode3 <- hptgeocode3 %>% mutate(XDR = if_else(AZM ==1 & CIP==1 & CRO==1,1,0))
hptgeocode3_relpropres <- hptgeocode3 %>% group_by(year, Group) %>% summarise(total = n(), azm = sum(AZM), cip = sum(CIP), cro = sum(CRO), xdr = sum(XDR)) 

hptgeocode3_relpropres <- hptgeocode3_relpropres %>% mutate(propazm = azm/total, propcip = cip/total, propcro = cro/total, propxdr = xdr/total)

hptgeocode3_relpropres_long <- hptgeocode3_relpropres %>%
  select(year, Group, propazm, propcip, propcro, propxdr) %>%
  pivot_longer(
    cols = starts_with("prop"),
    names_to = "Res",
    values_to = "Proportion"
  )

plotprop<- ggplot(hptgeocode3_relpropres_long,
                  aes(x = year, y = Proportion, color = Res)) +
  geom_point() +
  geom_line()+
  facet_wrap(~Group, ncol = 1) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Year",
    y = "Resistant (%)",
    color = "AMR Determinant"
  ) +
  theme_classic()+scale_color_brewer(palette = "RdYlBu")


library(tidyverse)
library(sf)
library(tmap)
library(raster)
library(ggspatial)
library(binom)

#~~~~~~~~~~~~~~~~~~MAP GENERATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#remove duplicate postcodes
tabl1datt2<- hptgeocode3[!duplicated(hptgeocode3$postcode), ]

#summarize case counts by postcode
tabl1datsum<- hptgeocode3%>% group_by(postcode) %>% summarize(sumcases = n(), sumpMSM = sum(pMSMfinal))

#merge geographic data to summary data 
tabl1datsummerge<- left_join(tabl1datsum, tabl1datt2, by = "postcode")

#transform to sf object
tabsf <- st_as_sf(tabl1datsummerge, coords = c("Long", "Lat"),crs=4326)

#jitter locations of specimen collection by 2.5 km 
set.seed(111)
tabsftransform<- st_transform(tabsf, CRS("+proj=utm +zone=30 + datum=WGS84"))
dtrjitter<- st_jitter(tabsftransform,2500)

#save copy of df
dtrjitterprop<-dtrjitter

#change total case summary col name
colnames(dtrjitter)[2]<- "Cases (N)"

#generate categorical cases variable
dtrjitter<- dtrjitter %>% mutate(BubblesCat = if_else(`Cases (N)` <= 5, 5,
                                                      if_else(`Cases (N)` >5 &`Cases (N)` <= 50, 50,
                                                              if_else(`Cases (N)` >50 & `Cases (N)`<= 100, 100, 
                                                                      if_else(`Cases (N)` >100,200, NA)))))
dtrjitter$BubblesCat = factor(dtrjitter$BubblesCat)

#read in national geographies shape file 
countiesultagen <- st_read("Countries_December_2021_UK_BUC_2022/CTRY_DEC_2021_UK_BUC.shp")

#select only geographies in England 
countiesultagen <- countiesultagen %>% filter(CTRY21NM == "England")
ukshp <- st_transform(countiesultagen, crs = 4326)

#plot national map
mapbig<- ggplot() +
  geom_sf(data = ukshp,fill = NA)+
  geom_sf(data = dtrjitter, aes(size = BubblesCat, alpha = 0.8)) + theme_classic()

#filtering cases to those identified in London 
dtrjitterlont <-dtrjitter %>% mutate(londoncheck = ifelse(grepl("London",HPT),1,0))
lonpoints <- dtrjitterlont %>% filter(londoncheck == 1)

#creating HPT-specific map to selection only London HPTs
counties <- st_read("Counties_and_Unitary_Authorities_December_2021_UK_BFE_2022/CTYUA_DEC_2021_UK_BFE.shp")
hptcon <- read.csv("Counties_and_Unitary_Authorities_December_2021_EN_BUC.csv")

hptconfilt <- hptcon %>% dplyr::select(CTYUA21NM, HPT)
mjoin <- left_join(counties, hptconfilt, by='CTYUA21NM')

#remove geographies that haven't been geocoded to HPTs (i.e. redundant geographies)
mjoin2 <- mjoin %>% filter(!is.na(HPT))

#merge with geocoded case data 
mjoin2red<- mjoin2 %>% dplyr::select(LONG, LAT, HPT, geometry)

#unify geographies by HPT
HPTs <- mjoin2red %>% 
  group_by(HPT) %>% 
  summarise()

#select only London HPTs
lon<-HPTs[9:11,]
lon$loncheck<-1

#merge geographies
lonmerg <- lon %>% 
  group_by(loncheck) %>% 
  summarise()

#map of cases in london hpts
mapsmall<- ggplot() +
  geom_sf(data = lonmerg,fill = NA)+
  geom_sf(data = lonpoints, aes(size = BubblesCat, alpha = 0.8)) + theme_classic()

#proportion pMSM of isolates in London/Manchester & elsewhere
dtrjitterlontman <-dtrjitter %>% mutate(londonmancheck = ifelse(grepl("London|Manch",HPT),1,0))
dtrjitterlontman<- dtrjitterlontman%>% group_by(londonmancheck) %>% summarize(pmsm = sum(sumpMSM), ncases = sum(`Cases (N)`)) %>% mutate(propmsm = binom.confint(x = pmsm,n = ncases, conf.level = 0.95, method = c("wilson"))) 
dtrjitterlontman$londonmancheck<- factor(dtrjitterlontman$londonmancheck)
dtrjitterlontman <- dtrjitterlontman %>% mutate(Location = if_else(londonmancheck == 1, "London & Manchester", "Other"))
propmsmfig<- ggplot(aes(y=dtrjitterlontman$propmsm$mean, x=Location), data = dtrjitterlontman)  +geom_point() + theme_classic()+geom_errorbar(aes(ymin=dtrjitterlontman$propmsm$lower, ymax=dtrjitterlontman$propmsm$upper, width = 0.2 ))+ scale_x_discrete(name = "HPT of Isolate Collection")+labs(y="Proportion pMSM")+scale_y_continuous(limits=c(0, .5))


#Alluvial plot of isolates
totbucketssum <- data.frame(hptgeocode3) %>% dplyr::select(Group,CIP, AZM, CRO) 
totbucketssum$CRO[totbucketssum$CRO ==0] <- "S"
totbucketssum$CRO[totbucketssum$CRO ==1] <- "R"
totbucketssum$AZM[totbucketssum$AZM ==0] <- "S"
totbucketssum$AZM[totbucketssum$AZM ==1] <- "R"
totbucketssum$CIP[totbucketssum$CIP ==0] <- "S"
totbucketssum$CIP[totbucketssum$CIP ==1] <- "R"
totbucketssum$Group<- factor(totbucketssum$Group, levels = c("pMSM", "non-pMSM", "Travel"))

col_vector = c("#088978","#b21819","#f85c06")
colnames(totbucketssum)[2:4]<- c("Azithromycin","Ceftriaxone","Ciprofloxacin")
p<- alluvial_wide( dplyr::select(totbucketssum,Group,Ciprofloxacin, Azithromycin, Ceftriaxone), fill_by = 'first_variable', stratum_labels  = F, col_vector_flow = col_vector) + theme_classic() 

alluv<- p+ geom_text(stat = "stratum", aes(label = after_stat(stratum)),
                     angle=c(90), size = 2, color = "white")+scale_y_continuous(name ="Frequency (n)")+scale_x_discrete(name ="Category")


#Supp figs
#relative proportion of cases over time by demographic group
hptgeocode3_sum <- hptgeocode3 %>%
  group_by(year) %>%
  summarise( pMSM = sum(pMSMfinal),Travel = sum(highrisktravel),nonpMSM = sum(domesticnonpmsm)) %>%
  pivot_longer(-year, names_to = "Group", values_to = "Tot") 

hptgeocode3_sum_proplong <- hptgeocode3_sum %>%
  group_by(year) %>%
  mutate(prop = Tot / sum(Tot))

propcases<- ggplot(hptgeocode3_sum_proplong, aes(x = factor(year), y = prop, fill = Group)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Year",
    y = "Cases (%)",
    fill = "Group"
  ) + theme_minimal() + scale_fill_manual(values = c( "#b21819","#088978","#f85c06"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#relative proportion of CIP, CRO, AZM, XDR
hptgeocode3 <- hptgeocode3 %>% mutate(XDR = if_else(AZM ==1 & CIP==1 & CRO==1,1,0))
hptgeocode3_relpropres <- hptgeocode3 %>% group_by(year, Group) %>% summarise(total = n(), azm = sum(AZM), cip = sum(CIP), cro = sum(CRO), xdr = sum(XDR)) 

hptgeocode3_relpropres <- hptgeocode3_relpropres %>% mutate(propazm = azm/total, propcip = cip/total, propcro = cro/total, propxdr = xdr/total)

hptgeocode3_relpropres_long <- hptgeocode3_relpropres %>%
  select(year, Group, propazm, propcip, propcro, propxdr) %>%
  pivot_longer(
    cols = starts_with("prop"),
    names_to = "Res",
    values_to = "Prop"
  )

plotprop<- ggplot(hptgeocode3_relpropres_long,
                  aes(x = year, y = Prop, color = Res)) +
  geom_point() +
  geom_line()+
  facet_wrap(~Group, ncol = 1) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Year",
    y = "Resistant (%)",
    color = "AMR Determinants"
  ) +
  theme_classic()+scale_color_brewer(palette = "RdYlBu")


