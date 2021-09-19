
library(Delaporte)





prop.sourcecase.array = c(seq(0.01,0.1,length.out = 11), seq(0.1,0.9,length.out = 21), seq(0.9,0.99,length.out = 11))

#prop.fixed.R.array = c(seq(0,0.1,length.out = 51), seq(0.1,0.9,length.out = 401), seq(0.9,1.0,length.out = 51))
prop.Delap.fixed.R.array = c(0.2,0.4,0.6,0.8)
prop.fixed.R.array = c(1,0.3,0)
int.array = 0:999
#comb.R.array = c(0.5,1,2,5)
#
#variable.k.array = 10^seq(-2, 2, by = 1)
variable.k.array = c(0.1,0.4,1,10)
col.for.dist = c(Poisson = 'darkorange', Delaporte = 'purple', NB = 'royalblue')





comb.R = 5
prop.fixed.R = 0.3

# dev.off()
#        pdf(file = 'dist_QQ_.pdf', width = 11, height = 11)
par(las = 1, mfrow = c(4,4), oma = c(2,3,0,1))
#
for (prop.Delap.fixed.R.j in 1:length(prop.Delap.fixed.R.array)) {#   prop.Delap.fixed.R.j = 1
  prop.Delap.fixed.R = prop.Delap.fixed.R.array[prop.Delap.fixed.R.j]
  prop.fixed.R.array = c(1,prop.Delap.fixed.R,0)
  
  for (variable.k.j in 1:length(variable.k.array)) {#       variable.k.j = 1
    variable.k = variable.k.array[variable.k.j]
    #variable.k = 0.5
    
    par(mar = c(2.5,2,1.5,1))
    plot(1,1, type = 'n', xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i', frame = F, ann = F, axe = F)
    lines(c(0,1), c(0,1), lty = 1, col = 'lightgrey')
    #
    for (prop.fixed.R.j in 1:length(prop.fixed.R.array)) {#           prop.fixed.R.j = 1
      prop.fixed.R = prop.fixed.R.array[prop.fixed.R.j]
      # prop.fixed.R = 0.3
      
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
      expectation.array = pm.array * int.array#; cumsum(expectation.array)
      #
      target.index.array = apply(X = as.matrix(prop.sourcecase.array), MARGIN = 1, FUN = function(z){
        which.min(z > cumsum(pm.array))[1]
      })
      target.weight.array = 1 - (cumsum(pm.array)[target.index.array] - prop.sourcecase.array) / pm.array[target.index.array]
      #   data.frame(target.index = target.index.array, target.weight = target.weight.array)
      #
      prop.secondcase.array = 1 - rev(apply(X = as.matrix(1:length(target.index.array)), MARGIN = 1, FUN = function(z){
        y = target.index.array[z]
        sum(expectation.array[1:y]) - expectation.array[y]*(1-target.weight.array[z])
      })) / comb.R
      # prop.secondcase.array = c(apply(X = as.matrix(1:length(target.index.array)), MARGIN = 1, FUN = function(z){
      #   y = target.index.array[z]
      #   sum(expectation.array[y:length(expectation.array)]) - expectation.array[y]*(1-target.weight.array[z])
      # })) / comb.R
      #
      lines(c(0,prop.sourcecase.array,1), c(0,prop.secondcase.array,1), lwd = 2, lty = c(3,2,4)[prop.fixed.R.j], col = col.for.dist[prop.fixed.R.j])
    }
    #
    axis(1); axis(2)
    mtext(side = 3, line = 0.22, adj = 0, text = paste0('(', LETTERS[(prop.Delap.fixed.R.j -1)*4 + variable.k.j], ') fixed frac. = ', prop.Delap.fixed.R, ';  disp. = ', variable.k.array[variable.k.j]), las = 0)
    #
    if(prop.Delap.fixed.R.j == 1 & variable.k.j == 1){
      legend('bottomright', lty = c(4,2,3), col = rev(col.for.dist), legend = c('NB','Delaporte', 'Poisson'), bty = 'n')
    }
  }
}
mtext(side = 2, line = 1, 'proportion of secondary cases', outer = T, las = 0)
mtext(side = 1, line = 0, outer = T, 'proportion of the most infectious cases', las = 0)
#      dev.off()
#      dev.off()
#      dev.off()














