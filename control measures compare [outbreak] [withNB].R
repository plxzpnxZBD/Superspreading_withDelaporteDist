


library(Delaporte)




# sum(dnbinom(x = 1:999 -1, mu = 1 *1:999, size = 0.5 *1:999) / 1:999)
# final.cluster.size = 3
# dnbinom(x = final.cluster.size -1, mu = 1 *final.cluster.size, size = 0.5 *final.cluster.size) / final.cluster.size



prop.fixed.R.array = c(0.0,0.1,0.3, 0.5)
comb.R.array = c(0.6,1,1.2,1.6,2)
variable.k = 0.2
#
num.pnt = 101
control.effort.array = seq(0,100, length.out = num.pnt) /100
outbreak.size.array = seq(1,100, length.out = 100)


# dev.off()
#        pdf(file = 'control_efficacy_[outbreak]_[withNB].pdf', width = 15, height = 9)
par(las = 1, mfrow = c(4,5), oma = c(2,3,0,1))
#
for (prop.fixed.R.j in 1:length(prop.fixed.R.array)) {#          prop.fixed.R.j = 1
  prop.fixed.R = prop.fixed.R.array[prop.fixed.R.j]
  #
  for (comb.R.j in 1:length(comb.R.array)) {#               comb.R.j = 1
    comb.R = comb.R.array[comb.R.j]
    #sse.threshold = qpois(p = 0.99, lambda = comb.R* (1 -0))# 
    
    fixed.R = comb.R *(1-control.effort.array) * prop.fixed.R
    variable.R = comb.R *(1-control.effort.array) - fixed.R #-0.1
    delap.para.alpha = variable.k; delap.para.alpha = ifelse(delap.para.alpha <= 0, 1e-99, delap.para.alpha)
    delap.para.beta = variable.R / variable.k; delap.para.beta = ifelse(delap.para.beta <= 0, 1e-99, delap.para.beta)
    delap.para.lambda = fixed.R; delap.para.lambda = ifelse(delap.para.lambda <= 0, 1e-99, delap.para.lambda)
    universal_based.delap.outbreak.prob = NULL
    for(control.effort.j in 1:length(control.effort.array)){#      control.effort.j = 1
      temp.outbreak.prob = 1 - sum(ddelap(
        x = outbreak.size.array -1, 
        alpha = c(delap.para.alpha * outbreak.size.array), 
        beta = rep(delap.para.beta[control.effort.j], length(outbreak.size.array)), 
        lambda = c(delap.para.lambda[control.effort.j] * outbreak.size.array)
      ) / outbreak.size.array)
      universal_based.delap.outbreak.prob = c(universal_based.delap.outbreak.prob, temp.outbreak.prob)
    }
    universal_based.control.efficacy = 1 - universal_based.delap.outbreak.prob / universal_based.delap.outbreak.prob[1]
    
    
    fixed.R = comb.R * prop.fixed.R
    variable.R = comb.R *(1-control.effort.array) - fixed.R
    fixed.R = ifelse(variable.R <= 0, fixed.R +variable.R, fixed.R)
    variable.R = ifelse(variable.R <= 0, 1e-99, variable.R)
    delap.para.alpha = variable.k; delap.para.alpha = ifelse(delap.para.alpha <= 0, 1e-99, delap.para.alpha)
    delap.para.beta = variable.R / variable.k; delap.para.beta = ifelse(delap.para.beta <= 0, 1e-99, delap.para.beta)
    delap.para.lambda = fixed.R; delap.para.lambda = ifelse(delap.para.lambda <= 0, 1e-99, delap.para.lambda)
    target_based.delap.outbreak.prob = NULL
    for(control.effort.j in 1:length(control.effort.array)){#      control.effort.j = 1
      temp.outbreak.prob = 1 - sum(ddelap(
        x = outbreak.size.array -1, 
        alpha = c(delap.para.alpha * outbreak.size.array), 
        beta = rep(delap.para.beta[control.effort.j], length(outbreak.size.array)), 
        lambda = c(delap.para.lambda[control.effort.j] * outbreak.size.array)
      ) / outbreak.size.array)
      target_based.delap.outbreak.prob = c(target_based.delap.outbreak.prob, temp.outbreak.prob)
    }
    target_based.control.efficacy = 1 - target_based.delap.outbreak.prob / target_based.delap.outbreak.prob[1]
    
    #
    
    par(mar = c(2,2,2,1))
    plot(1,1, type = 'n', xlim = c(0, 1), ylim = c(0, 1), xaxs = 'r', yaxs = 'i', ann = F, axe = F)
    lines(control.effort.array, target_based.control.efficacy, col = 'darkorange', lty = 1, lwd = 2)
    lines(control.effort.array, universal_based.control.efficacy, col = 'darkcyan', lty = 2, lwd = 1)
    rect(xleft = 1-prop.fixed.R, ybottom = -99, xright = 99, ytop = 99, col = '#00000066', border = 'black')
    #abline(v = 1-prop.fixed.R, lty = 2, col = 'darkgrey')
    axis(1); axis(2)
    mtext(side = 3, adj = 0, text = paste0('(', LETTERS[(prop.fixed.R.j -1) * 5 + comb.R.j], ') R = ', comb.R, '; fixed frac. = ', prop.fixed.R), las = 0)
    #
    if((prop.fixed.R.j + comb.R.j) == 2){
      legend('bottom', lty = c(1,2), col = c('darkorange', 'darkcyan'), legend = c('population-wide','high-risk-specific'), bty = 'n')
    }
  }
}
mtext(side = 1, line = 0.5, '(relative) reduction in reproduction number', outer = T)
mtext(side = 2, line = 1, outer = T, 'relative risk reduction of outbreak with size > 100', las = 0)
#      dev.off()
#      dev.off()
#      dev.off()





























