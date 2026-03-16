## MCMC fitness
## By group, by genotype 
install.packages("questionr")
library(questionr)
library(ape)
library(rstan)
library(gridExtra)
library(grid)
library(gtable)
library(tidyverse)
library(ggplot2)
install.packages("doParallel")
library(doParallel)  
rstan_options(auto_write = TRUE)
library(devtools)

############# Model
model.MCMC <- stan_model(file = "Model_fitness_1f0_1_switch_reals_independent1221.stan")

#############CIP##################################################################

############# data for MCMC ######################################################

data.MCMC_CIP = readRDS(file = 'Data_model_independent_CIP.rds')

############# parameters #############################################

## introduction
data.MCMC_CIP$yearF0 =array(rep(1, data.MCMC_CIP$nb_groups))
data.MCMC_CIP$yearIntroduction =array(rep(12, data.MCMC_CIP$nb_groups))
nb_geno = data.MCMC_CIP$nb_genotypes

#shared f0
f0_init = function(nb_geno){
  res = c(0, 0)
  res = rnorm(nb_geno, mean = 0.08, sd = 0.01)
  sum = sum(res)
  for(i in 1:2){
    res[i] = res[i]/sum
  }
  return(res)
}


##################################################################################
## Run MCMC 
##################################################################################
set.seed(111)

no_cores = 3 
registerDoParallel(cores=no_cores)  
cl = makeCluster(no_cores) 


#CIP
name_file = 'Output_per_allgroups_CIP'

set.seed(111)

for(i in 1:3) {
  print(paste0('Running chain n = ', i))
  fit_delay <- sampling(model.MCMC, data = data.MCMC_CIP, 
                        show_messages = TRUE, 
                        chains = 1,cores = 1,iter= 10000, chain_id = i,
                        control = list(adapt_delta = 0.97, max_treedepth = 10), init = list(list(f0 = f0_init( data.MCMC_CIP$nb_genotypes),fitness_genotypes_pre_switch = array(rep(rnorm(1,0,0.01),3)),fitness_genotypes_post_switch = array(rep(rnorm(1,0,0.01),3)))))
  fit = list(fit=fit_delay,
             data= data.MCMC_CIP)
  Chains=rstan::extract(fit$fit)
  saveRDS(Chains, file = paste0(name_file, '_chains_', i, '.rds'))
  saveRDS(data.MCMC_CIP, file = paste0(name_file, '_data_', i, '.rds'))
  
  m = monitor(fit$fit, print = F)
  fit$monitor = m
  saveRDS(fit, file = paste0(name_file, '_fit_', i, '.rds'))
}

print('Reading 1')
fit1 = readRDS(paste0(name_file, '_fit_', 1, '.rds'))
print('Reading 2')
fit2 = readRDS(paste0(name_file, '_fit_', 2, '.rds'))
print('Reading 3')
fit3 = readRDS(paste0(name_file, '_fit_', 3, '.rds'))

print('Fit')

fit = NULL
fit$fit = sflist2stanfit(list(fit1$fit, fit2$fit, fit3$fit))
fit$data = fit1$data

print('Chains')
Chains = rstan::extract(fit$fit)

print('Writing fit')
saveRDS(fit, paste0(name_file, '_fit_all.rds'))

print('Writing chains')
saveRDS(Chains, paste0(name_file, '_chains_all.rds'))


################################################################################
#look at WAIC
library(loo)
log_lik_1 <- extract_log_lik(fit$fit, merge_chains = F)
r_eff <- relative_eff(exp(log_lik_1), cores = 2)
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
waic_1 <- waic(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1$estimates)
print(waic_1$estimates)

#look at trace plots and Rhat
library(bayesplot)
library(coda)
library(ggplot2)
tracevars <- As.mcmc.list(fit$fit, pars = c("fitness_genotypes_post_switch[1]", "fitness_genotypes_post_switch[2]","fitness_genotypes_post_switch[3]","fitness_genotypes_vector_pre_switch[1]", "fitness_genotypes_vector_pre_switch[2]" , "fitness_genotypes_vector_pre_switch[3]", "f0[1]", "f0[2]"))
library(posterior)
traceplot<- mcmc_trace(tracevars) + theme(
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12))+ theme_classic()

traceplot
summary(as_draws(tracevars))$rhat
effectiveSize(tracevars)

################################################################################

## Load fit
## Chains
nonpmsm<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_post_switch[1]"))
pmsm<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_post_switch[2]"))
trav<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_post_switch[3]"))

nonpmsmpre<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_pre_switch[1]"))
pmsmpre<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_pre_switch[2]"))
travpre<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_pre_switch[3]"))

datindep<- data.frame("Group", "Mean", "q2.5","q972.5","q50","Switch")
colnames(datindep)[1:6]<- c("Group", "Mean", "q2.5","q972.5","q50","Switch")

datindep[1,1:6] <-c("Non-pMSM",mean(nonpmsm), c(quantile(nonpmsm, c(0.025,0.975, 0.5))),"Post")
datindep[2,1:6] <-c("pMSM",mean(pmsm), c(quantile(pmsm, c(0.025,0.975, 0.5))),"Post")
datindep[3,1:6] <-c("Travel",mean(trav), c(quantile(trav, c(0.025,0.975, 0.5))),"Post")

datindep[4,1:6] <-c("Non-pMSM",mean(nonpmsmpre), quantile(nonpmsmpre, c(0.025,0.975, 0.5)),"Pre")
datindep[5,1:6] <-c("pMSM",mean(pmsmpre), quantile(pmsmpre, c(0.025,0.975, 0.5)),"Pre")
datindep[6,1:6] <-c("Travel",mean(travpre), quantile(travpre, c(0.025,0.975, 0.5)),"Pre")

datindep$Switch <- factor(datindep$Switch, levels = c("Pre", "Post"))
datindep$Mean <-as.numeric(datindep$Mean)
datindep$q2.5 <-as.numeric(datindep$q2.5)
datindep$q972.5 <-as.numeric(datindep$q972.5)

plot1 <- ggplot(datindep, aes(x = Group, y = Mean, color = Switch)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(
    aes(ymin = q2.5, ymax = q972.5),
    width = 0.2,
    position = position_dodge(width = 0.5),linetype = "dashed"
  ) +
  labs(x = "Group", y = "CIP Relative Fitness") +
  theme_classic()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) 

datindeppost_pre<- datindep[4:6,]
colnames(datindeppost_pre)<- c("Group","Mean_pre", "q2.5_pre", "q972.5_pre", "q50_pre", "Switch")
datindeppost_post<- datindep[1:3,]
colnames(datindeppost_post)<- c("Group","Mean_post", "q2.5_post", "q972.5_post", "q50_post", "Switch")
datindeplong<- left_join(datindeppost_post,datindeppost_pre, by = "Group")
datindeplong$deltamean = datindeplong$Mean_post-datindeplong$Mean_pre
datindeplong$deltauci = datindeplong$q972.5_post-datindeplong$q2.5_pre
datindeplong$deltalci = datindeplong$q2.5_post-datindeplong$q972.5_pre
datindeplong_filt<- datindeplong %>% dplyr::select(Group, Mean_pre, q2.5_pre, q972.5_pre,Mean_post, q2.5_post, q972.5_post)

datindeplong_filt$changeinmean <- exp(datindeplong_filt$Mean_pre) - exp(datindeplong_filt$Mean_post)
datindeplong_filt$changeupperci <- exp(datindeplong_filt$q972.5_pre) - exp(datindeplong_filt$q2.5_post)
datindeplong_filt$changelowerci <- exp(datindeplong_filt$q2.5_pre) - exp(datindeplong_filt$q972.5_post)

install.packages("kable")
library(knitr)

datindeplong_filt <- datindeplong_filt %>%
  mutate(
    estimate_CI_pre_exp = sprintf("%.2f (%.2f, %.2f)", exp(Mean_pre), exp(q2.5_pre), exp(q972.5_pre)),
    estimate_CI_post_exp = sprintf("%.2f (%.2f, %.2f)", exp(Mean_post), exp(q2.5_post), exp(q972.5_post))
  )

kable(
  datindeplong_filt %>% dplyr::select(Group, estimate_CI_pre_exp, estimate_CI_post_exp),
  digits = 2,
  col.names = c("Group", "Pre-2015 (95% CI)","Post-2015 (95% CI)"))


################################################################################
## Evaluate Different F0 and Switch Years
##################################################################################
library(rstan)
library(RColorBrewer)
library(binom)
library(loo)

## Load fit
set.seed(111)
years<-seq(1,17,by=1)
dfyears<- data.frame(years)
catchtt<- data.frame(matrix(NA, nrow = 17, ncol = 18))
catchtt$X1<- years
catchttwaic<-catchttpwaic<- as.matrix(catchtt)

data.MCMC_CIP = readRDS(file = 'Data_model_independent_CIP.rds')
name_file = 'Output_per_cipswitch'

for(l in 1:17) {
  print(paste0('F0 = ', l))
  data.MCMC_CIP$yearF0 = array(rep(years[l], data.MCMC_CIP$nb_groups)) 
  for(y in 1:17) {
    if(l < y){
      print(paste0('First switch = ',l, y))
      data.MCMC_CIP$yearIntroduction = array(rep(years[y], data.MCMC_CIP$nb_groups)) 
      for(i in 1:3) {
        print(paste0('Running chain n = ', i))
        fit_delay <- sampling(model.MCMC, data = data.MCMC_CIP, 
                              show_messages = TRUE, 
                              chains = 1,cores = 1,iter= 10000, chain_id = i,
                              control = list(adapt_delta = 0.97, max_treedepth = 15), init = list(list(f0 = f0_init( data.MCMC_CIP$nb_genotypes),fitness_genotypes_pre_switch = array(rep(rnorm(1,0,0.01),3)),fitness_genotypes_post_switch = array(rep(rnorm(1,0,0.01),3)))))
        fit = list(fit=fit_delay,
                   data= data.MCMC_CIP)
        Chains=rstan::extract(fit$fit)
        saveRDS(Chains, file = paste0(name_file, '_chains_', i, '.rds'))
        saveRDS(data.MCMC_CIP, file = paste0(name_file, '_data_', i, '.rds'))
        
        m = monitor(fit$fit, print = F)
        fit$monitor = m
        saveRDS(fit, file = paste0(name_file, '_fit_', i, '.rds'))
      }
      print(y)
      fit1 = readRDS(paste0(name_file, '_fit_', 1, '.rds'))
      print('Reading 2')
      fit2 = readRDS(paste0(name_file, '_fit_', 2, '.rds'))
      print('Reading 3')
      fit3 = readRDS(paste0(name_file, '_fit_', 3, '.rds'))
      print('Fit')
      fit = NULL
      fit$fit = sflist2stanfit(list(fit1$fit, fit2$fit, fit3$fit))
      fit$data = fit1$data
      print('Chains')
      Chains = rstan::extract(fit$fit)
      log_lik_1 <- extract_log_lik(fit$fit, merge_chains = F)
      r_eff <- relative_eff(exp(log_lik_1), cores = 2)
      loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
      waic_1 <- waic(log_lik_1, r_eff = r_eff, cores = 2)
      pwaic= waic_1$estimates[2]
      waic= waic_1$estimates[3]
      catchttwaic[l,1+y]=waic
      catchttpwaic[l,1+y]=pwaic
    }}}

longwaic <- data.frame(catchttwaic) %>% tibble::rownames_to_column(var = "F0") %>% pivot_longer(cols = starts_with(c("X")), names_to = "Switch", values_to = "WAIC")
longwaic<- longwaic %>% filter(Switch != "X1" & Switch != "X2")
longwaic$switchclean <- as.numeric(gsub('X','',longwaic$Switch))-1
longwaic$startyear <- as.numeric(longwaic$F0)+2003
longwaic$switchyear <- as.numeric(longwaic$switchclean)+2003

longwaic$switches <- paste0(longwaic$startyear, sep = " , ",longwaic$switchyear)
longwaic<- longwaic %>% filter(!is.na(WAIC) & !is.na(switchclean))
waicplotall<-ggplot(aes(x = reorder(switches,WAIC), y = WAIC), data = longwaic) +geom_line(linetype= 2, color = "grey")+ geom_point(size=2)+theme_classic()

longwaic2<- longwaic %>%slice_min(WAIC, n = 25)
longwaic2[1,]
waicplotallcipfilt<-ggplot(aes(x = reorder(switches,WAIC), y = WAIC), data = longwaic2) +geom_line(linetype= 2, color = "grey")+ geom_point(size=2)+theme_classic()+ theme(axis.text.x=element_text(angle=90, hjust=1))

############################################################################################
## Plot fits for model, per group
############################################################################################

## Functions
mean.and.ci <-function(v){
  return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))
}

cexlab = 1
cexaxis = 1
cexdots= 0.9
cexmain = 0.6

group = c("Non-pMSM", "pMSM", "Travel")

nb_groups = fit$data$nb_groups
nb_genotypes = fit$data$nb_genotypes
nb_years = fit$data$nb_years

numbers_obs_all_clades = array(NA, dim = c(nb_genotypes, nb_years, nb_groups)) 
for(i in 1:nb_groups){
  numbers_obs_all_clades[1:(nb_genotypes-1),,i] = fit$data$data_genotype_non_ref[,,i]
  numbers_obs_all_clades[nb_genotypes,,i] = fit$data$data_genotype_ref[,i]
}

min_date = 2004
nb_chains = length(Chains$lp__)

pdf(width = 7, height = 15, file = "Figure_fits_pergroup_sharedf0_switch_cip.pdf", onefile = T)
par(mfrow = c(8,2), oma = c(1,1,1,0), mai = c(0.7,0.7,0.3,0.3))
colors = c("#b21819", "#088978","#f85c06")

for(c in 1:fit$data$nb_groups){
  titles = c('(Proportion CIP)',
             '(1-Proportion CIP)')
  count = 1
  d = fit$data$data_genotype_non_ref[count,,c]
  total_m = rep(0, length(d))
  total_cimin = rep(0, length(d))
  total_cimax = rep(0, length(d))
  for(i in 1:(nb_genotypes)){
    t = fit$data$data_total_number[,c]
    d_m = rep(0, length(d))
    d_ci = matrix(0, ncol = length(d), nrow = 2)
    if(i < nb_genotypes){
      d = fit$data$data_genotype_non_ref[count,,c]
      for(j in 1:nb_years){
        if(t[j]>0) {
          tmp = binom.confint(x = d[j], n = t[j], method = c("wilson"), type="central")
          d_m[j] = tmp$mean
          d_ci[1,j] = tmp$lower
          d_ci[2,j] = tmp$upper
        }
      }
      count = count +1
    }
    if(i == nb_genotypes){
      d = fit$data$data_genotype_ref[,c]
      for(j in 1:nb_years){
        if(t[j]>0) {
          tmp = binom.confint(x = d[j], n = t[j], method = c("wilson"), type="central")
          d_m[j] = tmp$mean
          d_ci[1,j] = tmp$lower
          d_ci[2,j] = tmp$upper
        }
      }
    }
    zeros = fit$data$non_zero_group_year[,c]
    
    d_m[which(zeros==0)] = NA
    d_ci[1,which(zeros==0)] = NA
    d_ci[2,which(zeros==0)] = NA
    
    f = Chains$pred_absolute_freq[,c,i,1:nb_years]
    
    f[which(is.infinite(f)==T)] = NA
    f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
    
    pch_times = fit$data$f0_introduction[c]
    
    ylims = c(0,1)
    if(pch_times>nb_years){pch_times=nb_years-1}
    if(pch_times<0){pch_times=1}  
    if(i == 1 & c < nb_groups)  { print(' yay')
      plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
           xlab="", ylab = 'Proportion', ylim = ylims, xlim = c(0,18),cex = cexdots,
           main = paste0(group[c], ' ',titles[i]), cex.main = cexmain,
           col = adjustcolor('grey30', alpha.f = 0.6), 
           yaxt = 'n', xaxt = 'n', cex.lab = cexlab, xaxs = "i",yaxs="i")}
    if(i == 1 & c == nb_groups)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                                       xlab="Time", ylab = 'Proportion', ylim = ylims,xlim = c(0,18),cex = cexdots,
                                       main = paste0(group[c], ' ',titles[i]), cex.main = cexmain,
                                       col = adjustcolor('grey30', alpha.f = 0.6), 
                                       yaxt = 'n', xaxt = 'n', cex.lab = cexlab, xaxs = "i",yaxs="i")}
    if(i > 1 & c < nb_groups)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                                     xlab="", ylab = '', ylim = ylims,xlim = c(0,18),cex = cexdots,
                                     main = paste0(group[c], ' ',titles[i]), cex.main = cexmain,
                                     col = adjustcolor('grey30', alpha.f = 0.6), 
                                     yaxt = 'n', xaxt = 'n', cex.lab = cexlab, xaxs = "i",yaxs="i")}
    if(i > 1 & c == nb_groups)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                                      xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(0,18), cex = cexdots,
                                      main = paste0(group[c], ' ',titles[i]), cex.main = cexmain,
                                      col = adjustcolor('grey30', alpha.f = 0.6), 
                                      yaxt = 'n', xaxt = 'n', cex.lab = cexlab, xaxs = "i",yaxs="i")}
    
    arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
           col = adjustcolor('grey30', alpha.f = 0.8))
    
    axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
    # axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
    axis(1, at = c(0,6,11,16, 21)+1,labels = c(0,6,11,16, 21)+min_date, cex.axis = cexaxis, lwd = 0.9,  title(xlab = "Year", ylab = "Proportion"))
    
    lines(1:nb_years,f_mean_ci[1,], lwd = 1.25,
          cex=1.4,xlab="",ylab="", 
          col = adjustcolor(colors[c], alpha.f = 0.7))
    polygon(x = c(1:nb_years, rev(1:nb_years)),
            y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colors[c], alpha.f = 0.4), border = F)
    
    # abline(v = fit$data$yearF0[c], lty = 3)
    abline(v = fit$data$yearIntroduction[c], lty = 3)
    
    total_m = total_m + f_mean_ci[1,]
    total_cimin = total_cimin + f_mean_ci[2,]
    total_cimax = total_cimax + f_mean_ci[3,]
  }
}
dev.off()



rd = c(1) ## order plot group
cols = RColorBrewer::brewer.pal('Set1', n = 1)
threshold = 0
cexdots= 0.8
pdf(width = 5, height = 9, file = "Figure_observed_predicted_loose_priors_independent_cipswitch.pdf", onefile = T)
# windows(width = 3/2.54, height = 9/2.54)

par(mfrow = c(3,1), oma = c(1,1,1,0), mai = c(0.7,0.7,0.3,0.3))
## Numbers

numbers_obs = fit$data$data_genotype_non_ref
numbers_pred_m = numbers_pred_cimax = numbers_pred_cimin = array(0, dim = c(nb_genotypes-1, nb_years, nb_groups))
for(i in 1:nb_groups){
  for(j in 1:(nb_genotypes-1)){
    numbers_pred_m[j,,i] = apply(Chains$pred_number_non_ref[,i,j,], MARGIN = 2, function(x)mean.and.ci(x)[1])
    numbers_pred_cimin[j,,i] = apply(Chains$pred_number_non_ref[,i,j,], MARGIN = 2, function(x)mean.and.ci(x)[2])
    numbers_pred_cimax[j,,i] = apply(Chains$pred_number_non_ref[,i,j,], MARGIN = 2, function(x)mean.and.ci(x)[3])
  }
  numbers_obs[,which(colSums(numbers_obs_all_clades[,,i]) < threshold),i] = NA
}
plot(numbers_obs[,,rd[1]], numbers_pred_m[,,rd[1]], col = adjustcolor(cols[rd[1]], alpha.f = 0.9), pch = 20, xlab = 'Observed', ylab = 'Predicted', yaxt = 'n', bty = 'n', xaxt = 'n',
     main = "Numbers",
     cex.axis = cexaxis, cex.lab = cexlab, cex = cexdots,
     xlim = c(0.5,max(Chains$pred_number_non_ref)), ylim = c(0.5,max(fit$data$data_genotype_non_ref, na.rm=T)), 
     log= 'xy')
abline(b=1, a=0, col = 'black', lty = 2)
axis(1, las = 1, cex.axis = cexaxis, at = c(0, 1, 10, 50, 500), labels = c(0,1, 10, 50, 500))
axis(2, las = 2, cex.axis = cexaxis, at = c(0,1, 10, 50, 500), labels = c(0,1, 10, 50, 500))
for(i in 2:nb_groups){
  points(numbers_obs[,,rd[i]], numbers_pred_m[,,rd[i]], col = adjustcolor(cols[rd[i]], alpha.f = 0.9), pch = 20, cex = cexdots)
}
cor.test(numbers_obs, numbers_pred_m) 
## Freqs absolute
freqs_obs = array(NA, dim = c(nb_genotypes, nb_years, nb_groups))
for(i in 1:nb_groups){
  for(j in 1:(nb_genotypes-1)){
    freqs_obs[j,,i] = fit$data$data_genotype_non_ref[j,,i]/fit$data$data_total_number[,i]
  }
  freqs_obs[nb_genotypes,,i] = fit$data$data_genotype_ref[,i]/fit$data$data_total_number[,i]
  
  freqs_obs[,which(colSums(numbers_obs_all_clades[,,i]) < threshold),i] = NA
}
freqs_obs[which(is.nan(freqs_obs))] = NA
freqs_pred_m = array(0, dim = c(nb_genotypes, nb_years, nb_groups))
freqs_pred_cimax = array(0, dim = c(nb_genotypes, nb_years, nb_groups))
freqs_pred_cimin = array(0, dim = c(nb_genotypes, nb_years, nb_groups))
for(i in 1:nb_groups){
  count = 1
  t = fit$data$data_total_number[,i]
  for(j in 1:(nb_genotypes-1)){
    tmp = t(apply(Chains$pred_number_non_ref[,i,count,], MARGIN = 1, function(x)x/t))
    freqs_pred_m[j,,i] = apply(tmp, MARGIN = 2, function(x)mean.and.ci(x)[1])
    freqs_pred_cimin[j,,i] = apply(tmp, MARGIN = 2, function(x)mean.and.ci(x)[2])
    freqs_pred_cimax[j,,i] = apply(tmp, MARGIN = 2, function(x)mean.and.ci(x)[3])
    count = count+1
  }
  f = t(apply(Chains$pred_number_ref[,i,], MARGIN = 1, function(x)x/t))
  freqs_pred_m[nb_genotypes,,i] = apply(f, MARGIN = 2, function(x)mean.and.ci(x)[1])
  freqs_pred_cimin[nb_genotypes,,i] = apply(f, MARGIN = 2, function(x)mean.and.ci(x)[2])
  freqs_pred_cimax[nb_genotypes,,i] = apply(f, MARGIN = 2, function(x)mean.and.ci(x)[3])
}
plot(freqs_obs[,,rd[1]], freqs_pred_m[,,rd[1]], col = adjustcolor(cols[rd[1]], alpha.f = 0.9), 
     pch = 20, xlab = 'Observed', ylab = 'Predicted', 
     yaxt = 'n', xaxt = 'n', bty='n',
     main = "Frequency",
     cex.axis = cexaxis, cex.lab = cexlab, xlim = c(0,1), ylim = c(0,1))
abline(b=1, a=0, col = 'black', lty = 2)
axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
axis(1, las = 1, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
for(i in 2:nb_groups){
  points(freqs_obs[,,rd[i]], freqs_pred_m[,,rd[i]], col = adjustcolor(cols[rd[i]], alpha.f = 0.9), pch = 20, cex = cexdots)
}

cor.test(freqs_obs, freqs_pred_m)
View(freqs_obs)

###################################################### CRO ###################
############# data for MCMC CRO######################################################

data.MCMC_CRO = readRDS(file = 'Data_model_independent_CRO.rds')

############# parameters #############################################

## introduction
data.MCMC_CRO$yearF0 =array(rep(9, data.MCMC_CRO$nb_groups))
data.MCMC_CRO$yearIntroduction =array(rep(15, data.MCMC_CRO$nb_groups))
nb_geno = data.MCMC_CRO$nb_genotypes

#shared f0
f0_init = function(nb_geno){
  res = c(0, 0)
  res = rnorm(nb_geno, mean = 0.08, sd = 0.01)
  sum = sum(res)
  for(i in 1:2){
    res[i] = res[i]/sum
  }
  return(res)
}


##################################################################################
## Run MCMC 
##################################################################################
set.seed(111)

no_cores = 3 
registerDoParallel(cores=no_cores)  
cl = makeCluster(no_cores) 


#CRO
name_file = 'Output_per_allgroups_CRO'

set.seed(111)

for(i in 1:3) {
  print(paste0('Running chain n = ', i))
  fit_delay <- sampling(model.MCMC, data = data.MCMC_CRO, 
                        show_messages = TRUE, 
                        chains = 1,cores = 1,iter= 10000, chain_id = i,
                        control = list(adapt_delta = 0.97, max_treedepth = 10), init = list(list(f0 = f0_init( data.MCMC_CRO$nb_genotypes),fitness_genotypes_pre_switch = array(rep(rnorm(1,0,0.01),3)),fitness_genotypes_post_switch = array(rep(rnorm(1,0,0.01),3)))))
  fit = list(fit=fit_delay,
             data= data.MCMC_CRO)
  Chains=rstan::extract(fit$fit)
  saveRDS(Chains, file = paste0(name_file, '_chains_', i, '.rds'))
  saveRDS(data.MCMC_CRO, file = paste0(name_file, '_data_', i, '.rds'))
  
  m = monitor(fit$fit, print = F)
  fit$monitor = m
  saveRDS(fit, file = paste0(name_file, '_fit_', i, '.rds'))
}

print('Reading 1')
fit1 = readRDS(paste0(name_file, '_fit_', 1, '.rds'))
print('Reading 2')
fit2 = readRDS(paste0(name_file, '_fit_', 2, '.rds'))
print('Reading 3')
fit3 = readRDS(paste0(name_file, '_fit_', 3, '.rds'))

print('Fit')

fit = NULL
fit$fit = sflist2stanfit(list(fit1$fit, fit2$fit, fit3$fit))
fit$data = fit1$data

print('Chains')
Chains = rstan::extract(fit$fit)

print('Writing fit')
saveRDS(fit, paste0(name_file, '_fit_all.rds'))

print('Writing chains')
saveRDS(Chains, paste0(name_file, '_chains_all.rds'))


################################################################################
#look at WAIC
log_lik_1 <- extract_log_lik(fit$fit, merge_chains = F)
r_eff <- relative_eff(exp(log_lik_1), cores = 2)
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
waic_1 <- waic(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1$estimates)
print(waic_1$estimates)

#look at trace plots and Rhat
tracevars <- As.mcmc.list(fit$fit, pars = c("fitness_genotypes_post_switch[1]", "fitness_genotypes_post_switch[2]","fitness_genotypes_post_switch[3]","fitness_genotypes_vector_pre_switch[1]", "fitness_genotypes_vector_pre_switch[2]" , "fitness_genotypes_vector_pre_switch[3]", "f0[1]", "f0[2]"))
traceplot<- mcmc_trace(tracevars) + theme_classic() + theme(
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12))

traceplot
summary(as_draws(tracevars))$rhat
effectiveSize(tracevars)

################################################################################

## Load fit
## Chains
nonpmsm<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_post_switch[1]"))
pmsm<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_post_switch[2]"))
trav<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_post_switch[3]"))

nonpmsmpre<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_pre_switch[1]"))
pmsmpre<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_pre_switch[2]"))
travpre<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_pre_switch[3]"))

datindep<- data.frame("Group", "Mean", "q2.5","q972.5","q50","Switch")
colnames(datindep)[1:6]<- c("Group", "Mean", "q2.5","q972.5","q50","Switch")

datindep[1,1:6] <-c("Non-pMSM",mean(nonpmsm), c(quantile(nonpmsm, c(0.025,0.975, 0.5))),"Post")
datindep[2,1:6] <-c("pMSM",mean(pmsm), c(quantile(pmsm, c(0.025,0.975, 0.5))),"Post")
datindep[3,1:6] <-c("Travel",mean(trav), c(quantile(trav, c(0.025,0.975, 0.5))),"Post")


datindep[4,1:6] <-c("Non-pMSM",mean(nonpmsmpre), quantile(nonpmsmpre, c(0.025,0.975, 0.5)),"Pre")
datindep[5,1:6] <-c("pMSM",mean(pmsmpre), quantile(pmsmpre, c(0.025,0.975, 0.5)),"Pre")
datindep[6,1:6] <-c("Travel",mean(travpre), quantile(travpre, c(0.025,0.975, 0.5)),"Pre")

datindep$Switch <- factor(datindep$Switch, levels = c("Pre", "Post"))
datindep$Mean <-as.numeric(datindep$Mean)
datindep$q2.5 <-as.numeric(datindep$q2.5)
datindep$q972.5 <-as.numeric(datindep$q972.5)

plot1 <- ggplot(datindep, aes(x = Group, y = Mean, color = Switch)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(
    aes(ymin = q2.5, ymax = q972.5),
    width = 0.2,
    position = position_dodge(width = 0.5),linetype = "dashed"
  ) +
  labs(x = "Group", y = "CRO Relative Fitness") +
  theme_classic()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) 

ggsave("CIP_switch_2015est_031626.pdf",plot1, height = 4, width = 5)

datindeppost_pre<- datindep[4:6,]
colnames(datindeppost_pre)<- c("Mean_pre", "q2.5_pre", "q972.5_pre", "q50_pre")
datindeppost_post<- datindep[1:3,]
colnames(datindeppost_post)<- c("Mean_post", "q2.5_post", "q972.5_post", "q50_post")

datindeplong<- left_join(datindeppost_post,datindeppost_pre, by = "Group")
datindeplong$deltamean = datindeplong$Mean_post-datindeplong$Mean_pre
datindeplong$deltauci = datindeplong$q972.5_post-datindeplong$q2.5_pre
datindeplong$deltalci = datindeplong$q2.5_post-datindeplong$q972.5_pre
datindeplong_filt<- datindeplong %>% dplyr::select(Group, Mean_pre, q2.5_pre, q972.5_pre,Mean_post, q2.5_post, q972.5_post)

datindeplong_filt$changeinmean <- exp(datindeplong_filt$Mean_pre) - exp(datindeplong_filt$Mean_post)
datindeplong_filt$changeupperci <- exp(datindeplong_filt$q972.5_pre) - exp(datindeplong_filt$q2.5_post)
datindeplong_filt$changelowerci <- exp(datindeplong_filt$q2.5_pre) - exp(datindeplong_filt$q972.5_post)

datindeplong_filt <- datindeplong_filt %>%
  mutate(
    estimate_CI_pre_exp = sprintf("%.2f (%.2f, %.2f)", exp(Mean_pre), exp(q2.5_pre), exp(q972.5_pre)),
    estimate_CI_post_exp = sprintf("%.2f (%.2f, %.2f)", exp(Mean_post), exp(q2.5_post), exp(q972.5_post))
  )

kable(
  datindeplong_filt %>% dplyr::select(Group, estimate_CI_pre_exp, estimate_CI_post_exp),
  digits = 2,
  col.names = c("Group", "Pre-2018 (95% CI)","Post-2018 (95% CI)"), digits = 2
)

################################################################################
## Evaluate Different F0 and Switch Years
##################################################################################
library(rstan)
library(RColorBrewer)
library(binom)
library(loo)


## Load fit
set.seed(111)
years<-seq(1,17,by=1)
dfyears<- data.frame(years)
catchtt<- data.frame(matrix(NA, nrow = 17, ncol = 18))
catchtt$X1<- years
catchttwaic<-catchttpwaic<- as.matrix(catchtt)

data.MCMC_CRO = readRDS(file = 'Data_model_independent_CRO.rds')
name_file = 'Output_per_croswitch'

for(l in 1:17) {
  print(paste0('F0 = ', l))
  data.MCMC_CRO$yearF0 = array(rep(years[l], data.MCMC_CRO$nb_groups)) 
  for(y in 1:17) {
    if(l < y){
      print(paste0('First switch = ',l, y))
      data.MCMC_CRO$yearIntroduction = array(rep(years[y], data.MCMC_CRO$nb_groups)) 
      for(i in 1:3) {
        print(paste0('Running chain n = ', i))
        fit_delay <- sampling(model.MCMC, data = data.MCMC_CRO, 
                              show_messages = TRUE, 
                              chains = 1,cores = 1,iter= 10000, chain_id = i,
                              control = list(adapt_delta = 0.97, max_treedepth = 15), init = list(list(f0 = f0_init( data.MCMC_CRO$nb_genotypes),fitness_genotypes_pre_switch = array(rep(rnorm(1,0,0.01),3)),fitness_genotypes_post_switch = array(rep(rnorm(1,0,0.01),3)))))
        fit = list(fit=fit_delay,
                   data= data.MCMC_CRO)
        Chains=rstan::extract(fit$fit)
        saveRDS(Chains, file = paste0(name_file, '_chains_', i, '.rds'))
        saveRDS(data.MCMC_CRO, file = paste0(name_file, '_data_', i, '.rds'))
        
        m = monitor(fit$fit, print = F)
        fit$monitor = m
        saveRDS(fit, file = paste0(name_file, '_fit_', i, '.rds'))
      }
      print(y)
      fit1 = readRDS(paste0(name_file, '_fit_', 1, '.rds'))
      print('Reading 2')
      fit2 = readRDS(paste0(name_file, '_fit_', 2, '.rds'))
      print('Reading 3')
      fit3 = readRDS(paste0(name_file, '_fit_', 3, '.rds'))
      print('Fit')
      fit = NULL
      fit$fit = sflist2stanfit(list(fit1$fit, fit2$fit, fit3$fit))
      fit$data = fit1$data
      print('Chains')
      Chains = rstan::extract(fit$fit)
      log_lik_1 <- extract_log_lik(fit$fit, merge_chains = F)
      r_eff <- relative_eff(exp(log_lik_1), cores = 2)
      loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
      waic_1 <- waic(log_lik_1, r_eff = r_eff, cores = 2)
      pwaic= waic_1$estimates[2]
      waic= waic_1$estimates[3]
      catchttwaic[l,1+y]=waic
      catchttpwaic[l,1+y]=pwaic
    }}}

longwaic <- data.frame(catchttwaic) %>% tibble::rownames_to_column(var = "F0") %>% pivot_longer(cols = starts_with(c("X")), names_to = "Switch", values_to = "WAIC")

longwaic<- longwaic %>% filter(Switch != "X1" & Switch != "X2")
longwaic$switchclean <- as.numeric(gsub('X','',longwaic$Switch))-1
longwaic$startyear <- as.numeric(longwaic$F0)+2003
longwaic$switchyear <- as.numeric(longwaic$switchclean)+2003

longwaic$switches <- paste0(longwaic$startyear, sep = " , ",longwaic$switchyear)
longwaic<- longwaic %>% filter(!is.na(WAIC) & !is.na(switchclean))
waicplotall<-ggplot(aes(x = reorder(switches,WAIC), y = WAIC), data = longwaic) +geom_line(linetype= 2, color = "grey")+ geom_point(size=2)+theme_classic()
longwaic2<- longwaic %>%slice_min(WAIC, n = 25)
longwaic2[1,]
waicplotallfilt<-ggplot(aes(x = reorder(switches,WAIC), y = WAIC), data = longwaic2) +geom_line(linetype= 2, color = "grey")+ geom_point(size=2)+theme_classic()+ theme(axis.text.x=element_text(angle=90, hjust=1))

############################################################################################
## Plot fits for model, per group
############################################################################################

## Functions
mean.and.ci <-function(v){
  return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))
}

cexlab = 1
cexaxis = 1
cexdots= 0.9
cexmain = 0.6

group = c("Non-pMSM", "pMSM", "Travel")

nb_groups = fit$data$nb_groups
nb_genotypes = fit$data$nb_genotypes
nb_years = fit$data$nb_years

numbers_obs_all_clades = array(NA, dim = c(nb_genotypes, nb_years, nb_groups)) 
for(i in 1:nb_groups){
  numbers_obs_all_clades[1:(nb_genotypes-1),,i] = fit$data$data_genotype_non_ref[,,i]
  numbers_obs_all_clades[nb_genotypes,,i] = fit$data$data_genotype_ref[,i]
}

min_date = 2012
nb_chains = length(Chains$lp__)

pdf(width = 7, height = 15, file = "Figure_fits_pergroup_sharedf0_switch_cro.pdf", onefile = T)
par(mfrow = c(8,2), oma = c(1,1,1,0), mai = c(0.7,0.7,0.3,0.3))
colors = c("#b21819", "#088978","#f85c06")

for(c in 1:fit$data$nb_groups){
  titles = c('(Proportion CRO)',
             '(1-Proportion CRO)')
  count = 1
  d = fit$data$data_genotype_non_ref[count,,c]
  total_m = rep(0, length(d))
  total_cimin = rep(0, length(d))
  total_cimax = rep(0, length(d))
  for(i in 1:(nb_genotypes)){
    t = fit$data$data_total_number[,c]
    d_m = rep(0, length(d))
    d_ci = matrix(0, ncol = length(d), nrow = 2)
    if(i < nb_genotypes){
      d = fit$data$data_genotype_non_ref[count,,c]
      for(j in 1:nb_years){
        if(t[j]>0) {
          tmp = binom.confint(x = d[j], n = t[j], method = c("wilson"), type="central")
          d_m[j] = tmp$mean
          d_ci[1,j] = tmp$lower
          d_ci[2,j] = tmp$upper
        }
      }
      count = count +1
    }
    if(i == nb_genotypes){
      d = fit$data$data_genotype_ref[,c]
      for(j in 1:nb_years){
        if(t[j]>0) {
          tmp = binom.confint(x = d[j], n = t[j], method = c("wilson"), type="central")
          d_m[j] = tmp$mean
          d_ci[1,j] = tmp$lower
          d_ci[2,j] = tmp$upper
        }
      }
    }
    zeros = fit$data$non_zero_group_year[,c]
    
    d_m[which(zeros==0)] = NA
    d_ci[1,which(zeros==0)] = NA
    d_ci[2,which(zeros==0)] = NA
    
    f = Chains$pred_absolute_freq[,c,i,1:nb_years]
    
    f[which(is.infinite(f)==T)] = NA
    f_mean_ci = apply(f, MARGIN = 2, function(x)mean.and.ci(x))
    
    pch_times = fit$data$f0_introduction[c]
    
    ylims = c(0,1)
    if(pch_times>nb_years){pch_times=nb_years-1}
    if(pch_times<0){pch_times=1}  
    if(i == 1 & c < nb_groups)  { print(' yay')
      plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
           xlab="", ylab = 'Proportion', ylim = ylims, xlim = c(9,18),cex = cexdots,
           main = paste0(group[c], ' ',titles[i]), cex.main = cexmain,
           col = adjustcolor('grey30', alpha.f = 0.6), 
           yaxt = 'n', xaxt = 'n', cex.lab = cexlab, xaxs = "i",yaxs="i")}
    if(i == 1 & c == nb_groups)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                                       xlab="Time", ylab = 'Proportion', ylim = ylims,xlim = c(9,18),cex = cexdots,
                                       main = paste0(group[c], ' ',titles[i]), cex.main = cexmain,
                                       col = adjustcolor('grey30', alpha.f = 0.6), 
                                       yaxt = 'n', xaxt = 'n', cex.lab = cexlab, xaxs = "i",yaxs="i")}
    if(i > 1 & c < nb_groups)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                                     xlab="", ylab = '', ylim = ylims,xlim = c(9,18),cex = cexdots,
                                     main = paste0(group[c], ' ',titles[i]), cex.main = cexmain,
                                     col = adjustcolor('grey30', alpha.f = 0.6), 
                                     yaxt = 'n', xaxt = 'n', cex.lab = cexlab, xaxs = "i",yaxs="i")}
    if(i > 1 & c == nb_groups)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                                      xlab="Time (years)", ylab = '', ylim = ylims,xlim = c(9,18), cex = cexdots,
                                      main = paste0(group[c], ' ',titles[i]), cex.main = cexmain,
                                      col = adjustcolor('grey30', alpha.f = 0.6), 
                                      yaxt = 'n', xaxt = 'n', cex.lab = cexlab, xaxs = "i",yaxs="i")}
    
    arrows(1:nb_years,d_ci[1,], 1:(nb_years),d_ci[2,], length=0, angle=0, code=3, lwd = 0.8,
           col = adjustcolor('grey30', alpha.f = 0.8))
    
    axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
    # axis(2, las = 2, cex.axis = cexaxis, lwd = 0.9)
    axis(1, at = c(0,6,11,16, 21)+1,labels = c(0,6,11,16, 21)+min_date, cex.axis = cexaxis, lwd = 0.9,  title(xlab = "Year", ylab = "Proportion"))
    
    lines(1:nb_years,f_mean_ci[1,], lwd = 1.25,
          cex=1.4,xlab="",ylab="", 
          col = adjustcolor(colors[c], alpha.f = 0.7))
    polygon(x = c(1:nb_years, rev(1:nb_years)),
            y = c(f_mean_ci[2,], rev(f_mean_ci[3,])), col = adjustcolor(colors[c], alpha.f = 0.4), border = F)
    
    # abline(v = fit$data$yearF0[c], lty = 3)
    abline(v = fit$data$yearIntroduction[c], lty = 3)
    
    total_m = total_m + f_mean_ci[1,]
    total_cimin = total_cimin + f_mean_ci[2,]
    total_cimax = total_cimax + f_mean_ci[3,]
  }
}
dev.off()



rd = c(1) ## order plot group
cols = RColorBrewer::brewer.pal('Set1', n = 1)
threshold = 0
cexdots= 0.8
pdf(width = 5, height = 9, file = "Figure_observed_predicted_loose_priors_independent_croswitch.pdf", onefile = T)
# windows(width = 3/2.54, height = 9/2.54)

par(mfrow = c(3,1), oma = c(1,1,1,0), mai = c(0.7,0.7,0.3,0.3))
## Numbers

numbers_obs = fit$data$data_genotype_non_ref
numbers_pred_m = numbers_pred_cimax = numbers_pred_cimin = array(0, dim = c(nb_genotypes-1, nb_years, nb_groups))
for(i in 1:nb_groups){
  for(j in 1:(nb_genotypes-1)){
    numbers_pred_m[j,,i] = apply(Chains$pred_number_non_ref[,i,j,], MARGIN = 2, function(x)mean.and.ci(x)[1])
    numbers_pred_cimin[j,,i] = apply(Chains$pred_number_non_ref[,i,j,], MARGIN = 2, function(x)mean.and.ci(x)[2])
    numbers_pred_cimax[j,,i] = apply(Chains$pred_number_non_ref[,i,j,], MARGIN = 2, function(x)mean.and.ci(x)[3])
  }
  numbers_obs[,which(colSums(numbers_obs_all_clades[,,i]) < threshold),i] = NA
}
plot(numbers_obs[,,rd[1]], numbers_pred_m[,,rd[1]], col = adjustcolor(cols[rd[1]], alpha.f = 0.9), pch = 20, xlab = 'Observed', ylab = 'Predicted', yaxt = 'n', bty = 'n', xaxt = 'n',
     main = "Numbers",
     cex.axis = cexaxis, cex.lab = cexlab, cex = cexdots,
     xlim = c(0.5,max(Chains$pred_number_non_ref)), ylim = c(0.5,max(fit$data$data_genotype_non_ref, na.rm=T)), 
     log= 'xy')
abline(b=1, a=0, col = 'black', lty = 2)
axis(1, las = 1, cex.axis = cexaxis, at = c(0, 1, 10, 50, 500), labels = c(0,1, 10, 50, 500))
axis(2, las = 2, cex.axis = cexaxis, at = c(0,1, 10, 50, 500), labels = c(0,1, 10, 50, 500))
for(i in 2:nb_groups){
  points(numbers_obs[,,rd[i]], numbers_pred_m[,,rd[i]], col = adjustcolor(cols[rd[i]], alpha.f = 0.9), pch = 20, cex = cexdots)
}
cor.test(numbers_obs, numbers_pred_m) 
## Freqs absolute
freqs_obs = array(NA, dim = c(nb_genotypes, nb_years, nb_groups))
for(i in 1:nb_groups){
  for(j in 1:(nb_genotypes-1)){
    freqs_obs[j,,i] = fit$data$data_genotype_non_ref[j,,i]/fit$data$data_total_number[,i]
  }
  freqs_obs[nb_genotypes,,i] = fit$data$data_genotype_ref[,i]/fit$data$data_total_number[,i]
  
  freqs_obs[,which(colSums(numbers_obs_all_clades[,,i]) < threshold),i] = NA
}
freqs_obs[which(is.nan(freqs_obs))] = NA
freqs_pred_m = array(0, dim = c(nb_genotypes, nb_years, nb_groups))
freqs_pred_cimax = array(0, dim = c(nb_genotypes, nb_years, nb_groups))
freqs_pred_cimin = array(0, dim = c(nb_genotypes, nb_years, nb_groups))
for(i in 1:nb_groups){
  count = 1
  t = fit$data$data_total_number[,i]
  for(j in 1:(nb_genotypes-1)){
    tmp = t(apply(Chains$pred_number_non_ref[,i,count,], MARGIN = 1, function(x)x/t))
    freqs_pred_m[j,,i] = apply(tmp, MARGIN = 2, function(x)mean.and.ci(x)[1])
    freqs_pred_cimin[j,,i] = apply(tmp, MARGIN = 2, function(x)mean.and.ci(x)[2])
    freqs_pred_cimax[j,,i] = apply(tmp, MARGIN = 2, function(x)mean.and.ci(x)[3])
    count = count+1
  }
  f = t(apply(Chains$pred_number_ref[,i,], MARGIN = 1, function(x)x/t))
  freqs_pred_m[nb_genotypes,,i] = apply(f, MARGIN = 2, function(x)mean.and.ci(x)[1])
  freqs_pred_cimin[nb_genotypes,,i] = apply(f, MARGIN = 2, function(x)mean.and.ci(x)[2])
  freqs_pred_cimax[nb_genotypes,,i] = apply(f, MARGIN = 2, function(x)mean.and.ci(x)[3])
}
plot(freqs_obs[,,rd[1]], freqs_pred_m[,,rd[1]], col = adjustcolor(cols[rd[1]], alpha.f = 0.9), 
     pch = 20, xlab = 'Observed', ylab = 'Predicted', 
     yaxt = 'n', xaxt = 'n', bty='n',
     main = "Frequency",
     cex.axis = cexaxis, cex.lab = cexlab, xlim = c(0,1), ylim = c(0,1))
abline(b=1, a=0, col = 'black', lty = 2)
axis(2, las = 2, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
axis(1, las = 1, cex.axis = cexaxis, at = seq(0,1,0.25), labels = c("0", "", "0.5", "", "1"), lwd = 0.9)
for(i in 2:nb_groups){
  points(freqs_obs[,,rd[i]], freqs_pred_m[,,rd[i]], col = adjustcolor(cols[rd[i]], alpha.f = 0.9), pch = 20, cex = cexdots)
}

cor.test(freqs_obs, freqs_pred_m)
View(freqs_obs)
