


library(Delaporte)
#library(boot)


offspring.list = as.null()


load(file = '[proc]_xxk_offspring_array_.RData')
offspring.list[[1]] = offspring.array
load(file = '[proc]_lim_offspring_array_[period_2].RData')
offspring.list[[2]] = offspring.array
load(file = '[proc]_lim_offspring_array_[period_1].RData')
offspring.list[[3]] = offspring.array
load(file = '[proc]_hku_offspring_array_.RData')#;  offspring.array = offspring.array[-c(1:46)] # sporadic = 46
offspring.list[[4]] = offspring.array
load(file = '[proc]_zyj_offspring_array_.RData')
offspring.list[[5]] = offspring.array
load(file = '[proc]_shen_offspring_array_.RData')
offspring.list[[6]] = offspring.array



Dela.R.fixed.array = c(0.26,0.38,0.17,0.17,0.00,0.06)
Dela.R.vary.array = c(0.43,0.30,0.65,0.42,0.72,2.00)
Dela.dispersion.array = c(0.24,0.11,0.09,0.16,0.23,0.05)
#
NB.R.fixed.array = rep(0,6)
NB.R.vary.array = c(0.69,0.68,0.81,0.58,0.71,0.96)
NB.dispersion.array = c(0.74,0.85,0.23,0.43,0.28,0.10)




x.array = 0:99


# dev.off()
#        pdf(file = 'fitting_outcomes.pdf', width = 11, height = 7)
par(las = 1, mfrow = c(2,3), oma = c(1.5,3,0,1))
#
for(dataset.j in 1:length(offspring.list)){#             dataset.j = 4
  offspring.array = offspring.list[[dataset.j]]
  #
  par(mar = c(3,2,2,1))
  hist(offspring.array +0.5, col = 'lightgrey', border = 'darkgrey', ann = F, axes = F, freq = F, xlim = c(0,13), ylim = c(0,0.8), breaks = c(0:99), xaxs = 'i', yaxs = 'i')
  #
  NB.lambda = NB.R.fixed.array[dataset.j]; NB.lambda = ifelse(NB.lambda == 0, 1e-99, NB.lambda)
  NB.shape.alpha = NB.dispersion.array[dataset.j]; NB.shape.alpha = ifelse(NB.shape.alpha == 0, 1e-99, NB.shape.alpha)
  NB.shape.beta = NB.R.vary.array[dataset.j] / NB.shape.alpha; NB.shape.beta = ifelse(NB.shape.beta == 0, 1e-99, NB.shape.beta)
  NB.PMF.array = ddelap(
    x = x.array, log = F, 
    alpha = NB.shape.alpha, beta = NB.shape.beta, lambda = NB.lambda
  )
  lines(x.array +0.45, NB.PMF.array, col = 'royalblue', type = 'b', lwd = 1, lty = 4, pch = 2)
  #
  Dela.lambda = Dela.R.fixed.array[dataset.j]; Dela.lambda = ifelse(Dela.lambda == 0, 1e-99, Dela.lambda)
  Dela.shape.alpha = Dela.dispersion.array[dataset.j]; Dela.shape.alpha = ifelse(Dela.shape.alpha == 0, 1e-99, Dela.shape.alpha)
  Dela.shape.beta = Dela.R.vary.array[dataset.j] / Dela.shape.alpha; Dela.shape.beta = ifelse(Dela.shape.beta == 0, 1e-99, Dela.shape.beta)
  Dela.PMF.array = ddelap(
    x = x.array, log = F, 
    alpha = Dela.shape.alpha, beta = Dela.shape.beta, lambda = Dela.lambda
  )
  lines(x.array +0.55, Dela.PMF.array, col = 'purple', type = 'b', lwd = 1, lty = 2, pch = 5)
  #
  axis(1, at = c(-1:99)+0.5, labels = c(-1:99))
  axis(2)
  mtext(side = 3, adj = 0, text = paste0('(', LETTERS[dataset.j], ') dataset #', c(1,'2a','2b',3,4,5)[dataset.j]), las = 0)
  #
  if(dataset.j == 1){
    legend('right', lty = c(4,2), pch = c(2,5), col = c('royalblue','purple'), legend = c('NB','Delaporte'), bty = 'n')
  }
}
mtext(side = 1, line = 0, '# of secondary cases', outer = T)
mtext(side = 2, line = 1, outer = T, 'relative frequency', las = 0)
#      dev.off()
#      dev.off()
#      dev.off()












