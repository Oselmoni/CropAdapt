## custom function to rapidly get p-value from many kendall taus
kendall_pvalue <- function(T, n) {
  # Z-score under H0: tau = 0
  z <- T * sqrt(9 * n * (n - 1) / (2 * (2 * n + 5)))
  
  p <- 2 * (1 - pnorm(abs(z)))
  
  return(p)
}



## custom plot
cplot = function(x,y, col=1, pch=16, cex=1, main='', xlab='', ylab='', signif=2, adj=0) {
  
  
  plot(x,y, col=col, pch=pch, cex=cex,axes=F, xlab=NA, ylab=NA)
  
  # add axes
  axis(1, at=seq(min(x), max(x), length.out=3), labels=signif(seq(min(x), max(x), length.out=3),signif))
  axis(2, at=seq(min(y), max(y), length.out=3), labels=signif(seq(min(y), max(y), length.out=3),signif))
  
  # add labels
  title(xlab=xlab, ylab=ylab, line=2)
  
  # add title
  title(main=main, adj=adj)
  
  
}
