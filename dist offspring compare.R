




library(Delaporte)







x.array = 0:99
comb.R.array = c(0.5,1,2, 5)
variable.k = 0.5
#col.for.dist = c('pink', 'red', 'skyblue')


# dev.off()
#        pdf(file = 'dist_offspring_compare.pdf', width = 11, height = 3)
par(las = 1, mfrow = c(1,4), oma = c(1.5,3,0,1))
#
for (comb.R.j in 1:length(comb.R.array)) {#           comb.R.j = 1
  comb.R = comb.R.array[comb.R.j]
  #comb.R = 5
  
  prop.fixed.R = 0.3
  fixed.R = comb.R * prop.fixed.R
  variable.R = comb.R - fixed.R
  delap.para.alpha = variable.k; delap.para.alpha = ifelse(delap.para.alpha == 0, 1e-99, delap.para.alpha)
  delap.para.beta = variable.R / variable.k; delap.para.beta = ifelse(delap.para.beta == 0, 1e-99, delap.para.beta)
  delap.para.lambda = fixed.R; delap.para.lambda = ifelse(delap.para.lambda == 0, 1e-99, delap.para.lambda)
  delap.PMF.array = ddelap(
    x = x.array, log = F, 
    alpha = delap.para.alpha, beta = delap.para.beta, lambda = delap.para.lambda
  )
  #
  prop.fixed.R = 0.0
  fixed.R = comb.R * prop.fixed.R
  variable.R = comb.R - fixed.R
  delap.para.alpha = variable.k; delap.para.alpha = ifelse(delap.para.alpha == 0, 1e-99, delap.para.alpha)
  delap.para.beta = variable.R / variable.k; delap.para.beta = ifelse(delap.para.beta == 0, 1e-99, delap.para.beta)
  delap.para.lambda = fixed.R; delap.para.lambda = ifelse(delap.para.lambda == 0, 1e-99, delap.para.lambda)
  NB.PMF.array = ddelap(
    x = x.array, log = F, 
    alpha = delap.para.alpha, beta = delap.para.beta, lambda = delap.para.lambda
  )
  #
  prop.fixed.R = 1.0
  fixed.R = comb.R * prop.fixed.R
  variable.R = comb.R - fixed.R
  delap.para.alpha = variable.k; delap.para.alpha = ifelse(delap.para.alpha == 0, 1e-99, delap.para.alpha)
  delap.para.beta = variable.R / variable.k; delap.para.beta = ifelse(delap.para.beta == 0, 1e-99, delap.para.beta)
  delap.para.lambda = fixed.R; delap.para.lambda = ifelse(delap.para.lambda == 0, 1e-99, delap.para.lambda)
  Pois.PMF.array = ddelap(
    x = x.array, log = F, 
    alpha = delap.para.alpha, beta = delap.para.beta, lambda = delap.para.lambda
  )
  #
  par(mar = c(3,2,2,1))
  plot(1,1, type = 'n', xlim = c(0, 10), ylim = c(0, 0.7), xaxs = 'r', yaxs = 'r', ann = F, axe = F)
  lines(x.array, Pois.PMF.array, col = 'orange', type = 'b', lty = 3, pch = 6)
  lines(x.array, NB.PMF.array, col = 'skyblue', type = 'b', lty = 4, pch = 2)
  lines(x.array, delap.PMF.array, col = 'purple', type = 'b', lty = 2, pch = 5)
  #
  axis(1); axis(2)
  # mtext(side = 3, adj = 0, '(A) LR test profile using offspring dataset', las = 0)
  mtext(side = 3, adj = 0, text = paste0('(', LETTERS[comb.R.j], ') reproduction number = ', comb.R), las = 0)
  #
  if(comb.R.j == 1){
    legend('topright', lty = c(4,2,3), pch = c(2,5,6), col = c('skyblue','purple','orange'), legend = c('NB','Delaporte', 'Poisson'), bty = 'n')
  }
}
mtext(side = 1, line = 0, '# of secondary cases', outer = T)
mtext(side = 2, line = 1, outer = T, 'relative frequency', las = 0)
#      dev.off()
#      dev.off()
#      dev.off()



















