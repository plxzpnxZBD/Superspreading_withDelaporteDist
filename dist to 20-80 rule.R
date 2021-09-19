
library(Delaporte)







prop.fixed.R.array = c(seq(0,0.1,length.out = 51), seq(0.1,0.9,length.out = 401), seq(0.9,1.0,length.out = 51))
int.array = 0:99
comb.R.array = c(0.5,1,2, 5)
#
#variable.k.array = 10^seq(-2, 2, by = 1)
variable.k.array = c(0.01,0.1,0.3,1,10)
col.for.k = c('gold', 'pink', 'darkorange', 'red', 'darkred')


# dev.off()
#        pdf(file = '20-80_rule_.pdf', width = 11, height = 4)
par(las = 1, mfrow = c(1,4), oma = c(2,4,0,1))
#
for (comb.R.j in 1:length(comb.R.array)) {#           comb.R.j = 1
  comb.R = comb.R.array[comb.R.j]
  #
  par(mar = c(3,2,2,1))
  plot(1,1, type = 'n', xlim = c(0, 1), ylim = c(0, 0.7), xaxs = 'i', yaxs = 'i', ann = F, axe = F)
  #
  for (variable.k.j in 1:length(variable.k.array)) {#       variable.k.j = 1
    #variable.k = 0.01
    variable.k = variable.k.array[variable.k.j]
    #
    sol.index.prop.array = NULL
    for (prop.j in 1:length(prop.fixed.R.array)) {#        prop.j = 101
      prop.fixed.R = prop.fixed.R.array[prop.j]
      #
      fixed.R = comb.R * prop.fixed.R
      variable.R = comb.R - fixed.R
      #
      delap.para.alpha = variable.k
      delap.para.alpha = ifelse(delap.para.alpha == 0, 1e-99, delap.para.alpha)
      delap.para.beta = variable.R / variable.k
      delap.para.beta = ifelse(delap.para.beta == 0, 1e-99, delap.para.beta)
      delap.para.lambda = fixed.R
      delap.para.lambda = ifelse(delap.para.lambda == 0, 1e-99, delap.para.lambda)
      
      #
      pm.array = ddelap(
        x = int.array, log = F, 
        alpha = delap.para.alpha, beta = delap.para.beta, lambda = delap.para.lambda
      )
      expectation.array = ddelap(
        x = int.array, log = F, 
        alpha = delap.para.alpha, beta = delap.para.beta, lambda = delap.para.lambda
      ) *int.array
      # offspring.prop.array = c(expectation.array) / comb.R
      # offspring.cumprop.array = cumsum(expectation.array) / comb.R
      #
      minor.infectee.num = comb.R*0.2
      sol.int.index = which(cumsum(expectation.array) >= minor.infectee.num)[1]
      #
      excess.infectee.num = cumsum(expectation.array)[sol.int.index] - minor.infectee.num
      excess.infector.prop = (excess.infectee.num / expectation.array[sol.int.index]) * pm.array[sol.int.index]
      threshold.int = int.array[sol.int.index]
      sol.index.prop = pdelap(
        q = threshold.int -0, lower.tail = F,
        alpha = delap.para.alpha, beta = delap.para.beta, lambda = delap.para.lambda
      ) + excess.infector.prop
      #   print(round(sol.index.prop, 2))
      #sol.index.prop = ifelse(is.na(sol.index.prop),0,sol.index.prop)
      sol.index.prop.array = c(sol.index.prop.array, sol.index.prop)
    }
    lines(prop.fixed.R.array, sol.index.prop.array, col = col.for.k[variable.k.j], lwd = 2)
  }
  axis(1); axis(2)
  axis(1, line = 1, at = c(0,1), labels = c('NB', 'Poisson'), tick = F)
 # mtext(side = 3, adj = 0, '(A) LR test profile using offspring dataset', las = 0)
  mtext(side = 3, adj = 0, text = paste0('(', LETTERS[comb.R.j], ') reproduction number = ', comb.R), las = 0)
  #
  if(comb.R.j == 1){
    legend('topleft', lwd = 2, col = rev(col.for.k), legend = rev(paste0('dispersion = ', variable.k.array)), bty = 'n')
  }
  #
}
mtext(side = 1, line = 1, 'fraction of the fixed components', outer = T)
mtext(side = 2, line = 1, outer = T, 'proportion of the most infectious cases \nthat cause 80% of secondary cases', las = 0)
#      dev.off()
#      dev.off()
#      dev.off()










# dnbinom(x = int.array, mu = 2, size = 0.2)
# pnbinom(q = int.array, mu = 2, size = 0.2, lower.tail = F)
# dnbinom(x = int.array, mu = 2, size = 0.2)*int.array











