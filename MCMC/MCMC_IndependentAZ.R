## MCMC fitness
## By group, by genotype 
install.packages("questionr")
library(questionr)
library(ape)
library(rstan)
library(gridExtra)
library(grid)
library(gtable)
install.packages("doParallel")
library(doParallel)  
rstan_options(auto_write = TRUE)
library(devtools)

############# Model

model.MCMC <- stan_model(file = "Model_fitness_no_switch_reals_independent_1221.stan")

#########load data 
data.MCMC = readRDS(file = 'Data_model_independent_AZ.rds')
#read in CRO or CIP data here instead to run analysis on other determinants
#data.MCMC = readRDS(file = 'Data_model_independent_CIP.rds')
#data.MCMC = readRDS(file = 'Data_model_independent_CRO.rds')

data.MCMC$yearIntroduction =array(rep(1, data.MCMC$nb_groups))
data.MCMC$yearF0 =array(rep(1, data.MCMC$nb_groups))

f0_init = function(nb_groups, nb_geno){
  res = matrix(0, ncol = nb_geno, nrow = nb_groups)
  for(i in 1:nb_groups){
    res[i,] = rnorm(nb_geno, mean = 0.08, sd = 0.01)
    res[i,] = res[i,]/sum(res[i,])
  }
  return(res)
}

##################################################################################
## Run MCMC 
##################################################################################
no_cores = 3 
registerDoParallel(cores=no_cores)  
cl = makeCluster(no_cores) 

set.seed(111)
name_file = 'AZ_independent'

for(i in 1:3) {
  print(paste0('Running chain n = ', i))
  fit_delay <- sampling(model.MCMC, data = data.MCMC, 
                        show_messages = TRUE, 
                        chains = 1,cores = 1,iter= 10000, chain_id = i,
                        control = list(adapt_delta = 0.97, max_treedepth = 10), init = list(list(f0 = f0_init(data.MCMC$nb_groups, data.MCMC$nb_genotypes),fitness_genotypes_post_vacc = array(rep(rnorm(1,0,0.01),3)))))
  fit = list(fit=fit_delay,
             data= data.MCMC)
  Chains=rstan::extract(fit$fit)
  saveRDS(Chains, file = paste0(name_file, '_chains_', i, '.rds'))
  saveRDS(data.MCMC, file = paste0(name_file, '_data_', i, '.rds'))
  
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
library(posterior)

tracevars <- As.mcmc.list(fit$fit, pars = c("fitness_genotypes_post_vacc[1]", "fitness_genotypes_post_vacc[2]","fitness_genotypes_post_vacc[3]", "f0[1,1]", "f0[1,2]", "f0[2,1]", "f0[2,2]", "f0[3,1]", "f0[3,2]"))


traceplot<- mcmc_trace(tracevars) + theme_classic() + theme(
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12))

traceplot
summary(as_draws(tracevars))
effectiveSize(tracevars)

################################################################################
#plotting

library(rstan)
library(RColorBrewer)
library(binom)

## Load fit
## Chains
Chains=rstan::extract(fit$fit)
non<- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_post_vacc[1]"))
pmsm <- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_post_vacc[2]"))
trav <- as.matrix(fit$fit, pars = c("fitness_genotypes_vector_post_vacc[3]"))

datindep<- data.frame("Group")
datindep[1,1] <-"pMSM"
datindep[1,2] <-mean(pmsm)
datindep[1,3:5] <-quantile(pmsm,c(0.025,0.975,0.5))

colnames(datindep)[1:5]<- c("Group", "Mean","q2.5","q972.5","q50")

datindep[2,1] <-"non-pMSM"
datindep[2,2] <-mean(non)
datindep[2,3:5] <-quantile(non, c(0.025,0.975,0.5))

datindep[3,1] <-"Travel"
datindep[3,2] <-mean(trav)
datindep[3,3:5] <-quantile(trav, c(0.025,0.975,0.5))
datindep2<- datindep
datindep2[,3:6]<- exp(datindep2[,3:6])


############################################################################################
## Plot fits for model, per group
############################################################################################

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

## Functions
mean.and.ci <-function(v){
  return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))
}

cexlab = 1
cexaxis = 1
cexdots= 0.9
cexmain = 0.6


group = c("Non-pMSM", "pMSM", "Travel")

pdf(width = 7, height = 15, file = "Figure_fits_pergroupindep_AZ.pdf", onefile = T)
par(mfrow = c(8,2), oma = c(1,1,1,0), mai = c(0.7,0.7,0.3,0.3))
colors = c("#b21819", "#088978","#f85c06")


for(c in 1:fit$data$nb_groups){
  titles = c('(Proportion AZM)',
             '(1-Proportion AZM)')
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
    
    pch_times = -1
    
    ylims = c(0,1)
    if(pch_times>nb_years){pch_times=nb_years-1}
    if(pch_times<0){pch_times=1}  
    if(i == 1 & c < 18)  { print(' yay')
      plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
           xlab="", ylab = 'Proportion', ylim = ylims, xlim = c(0,18),cex = cexdots,
           main = paste0(group[c], ' ',titles[i]), cex.main = cexmain,
           col = adjustcolor('grey30', alpha.f = 0.6), 
           yaxt = 'n', xaxt = 'n', cex.lab = cexlab, xaxs = "i",yaxs="i")}
    if(i == 1 & c == 18)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                                xlab="Time", ylab = 'Proportion', ylim = ylims,xlim = c(0,18),cex = cexdots,
                                main = paste0(group[c], ' ',titles[i]), cex.main = cexmain,
                                col = adjustcolor('grey30', alpha.f = 0.6), 
                                yaxt = 'n', xaxt = 'n', cex.lab = cexlab, xaxs = "i",yaxs="i")}
    if(i > 1 & c < 18)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
                              xlab="", ylab = '', ylim = ylims,xlim = c(0,18),cex = cexdots,
                              main = paste0(group[c], ' ',titles[i]), cex.main = cexmain,
                              col = adjustcolor('grey30', alpha.f = 0.6), 
                              yaxt = 'n', xaxt = 'n', cex.lab = cexlab, xaxs = "i",yaxs="i")}
    if(i > 1 & c == 18)  {plot(1:nb_years, d_m, type="p", pch=c(rep(17, pch_times-1), rep(17, nb_years-pch_times+1)), bty = 'n',
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
pdf(width = 5, height = 9, file = "Figure_observed_predicted_loose_priors_indepndent_AZ.pdf", onefile = T)

# windows(width = 3/2.54, height = 9/2.54)
par(mfrow = c(3,1), oma = c(1,1,1,0), mai = c(0.7,0.7,0.3,0.3))
## Numbers
fit$data$data_total_number
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
# freqs_obs[which(freqs_obs == 0)] = NA
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
