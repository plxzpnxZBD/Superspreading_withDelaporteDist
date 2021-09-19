


library(Delaporte)







control.effort.array = c(0.1,0.2,0.4, 0.6)
comb.R.array = c(0.5,1,1.5, 3)
variable.k = 0.2
#
num.pnt = 101
prop.fixed.R.array = seq(1,99, length.out = num.pnt) /100


# dev.off()
#        pdf(file = 'control_efficacy_withinDelaporte.pdf', width = 12, height = 8)
par(las = 1, mfrow = c(4,4), oma = c(2,3,0,1))
#
for (control.effort.j in 1:length(control.effort.array)) {#          control.effort.j = 1
  control.effort = control.effort.array[control.effort.j]
  #
  for (comb.R.j in 1:length(comb.R.array)) {#               comb.R.j = 1
    comb.R = comb.R.array[comb.R.j]
    sse.threshold = qpois(p = 0.99, lambda = comb.R* (1 -0))# 
    
    
    fixed.R = comb.R *(1-control.effort) * prop.fixed.R.array
    variable.R = comb.R *(1-control.effort) - fixed.R #-0.1
    delap.para.alpha = variable.k; delap.para.alpha = ifelse(delap.para.alpha <= 0, 1e-99, delap.para.alpha)
    delap.para.beta = variable.R / variable.k; delap.para.beta = ifelse(delap.para.beta <= 0, 1e-99, delap.para.beta)
    delap.para.lambda = fixed.R; delap.para.lambda = ifelse(delap.para.lambda <= 0, 1e-99, delap.para.lambda)
    universal_based.delap.SSE.prob = pdelap(q = rep(sse.threshold,num.pnt) -0.01, lower.tail = F, alpha = rep(delap.para.alpha,num.pnt), beta = delap.para.beta, lambda = delap.para.lambda)
    #
    fixed.R = comb.R * prop.fixed.R.array
    variable.R = comb.R *(1-control.effort) - fixed.R
    fixed.R = ifelse(variable.R <= 0, fixed.R +variable.R, fixed.R)
    variable.R = ifelse(variable.R <= 0, 1e-99, variable.R)
    delap.para.alpha = variable.k; delap.para.alpha = ifelse(delap.para.alpha <= 0, 1e-99, delap.para.alpha)
    delap.para.beta = variable.R / variable.k; delap.para.beta = ifelse(delap.para.beta <= 0, 1e-99, delap.para.beta)
    delap.para.lambda = fixed.R; delap.para.lambda = ifelse(delap.para.lambda <= 0, 1e-99, delap.para.lambda)
    target_based.delap.SSE.prob = pdelap(q = rep(sse.threshold,num.pnt) -0.01, lower.tail = F, alpha = rep(delap.para.alpha,num.pnt), beta = delap.para.beta, lambda = delap.para.lambda)
    #
    control.efficacy = 1 - target_based.delap.SSE.prob / universal_based.delap.SSE.prob
    
    par(mar = c(3,2,2,1))
    plot(1,1, type = 'n', xlim = c(0, 1), ylim = c(0, 1), xaxs = 'r', yaxs = 'i', ann = F, axe = F)
    lines(prop.fixed.R.array, control.efficacy, col = 'orange', lty = 2, lwd = 2)
    sel.index = intersect(which(prop.fixed.R.array <= 0.5), which(prop.fixed.R.array >= 0.1))
    lines(prop.fixed.R.array[sel.index], control.efficacy[sel.index], col = 'red', lwd = 3)
    rect(xleft = 1-control.effort, ybottom = -99, xright = 99, ytop = 99, col = '#00000066', border = 'black')
    #abline(v = 1-control.effort, lty = 2, col = 'darkgrey')
    axis(1); axis(2)
    axis(1, at = 0, labels = 'NB', tick = F, line = 1)
    axis(1, at = 1, labels = 'Poisson', tick = F, line = 1)
    mtext(side = 3, adj = 0, text = paste0('(', LETTERS[(control.effort.j -1) * 4 + comb.R.j], ') R = ', comb.R, '; reduction in R = ', control.effort*100, '%'), las = 0)
    #
  }
}
mtext(side = 1, line = 0.5, 'fraction of the fixed component', outer = T)
mtext(side = 2, line = 1, outer = T, 'relative risk reduction of superspreading events', las = 0)
#      dev.off()
#      dev.off()
#      dev.off()





























