signif_ceiling <- function(x, digits=6){
    pow <- ceiling(log10(x))
    y <- ceiling(x*10^(digits-pow))/10^(digits-pow)
    return(y)
}

plot_cov <- function(df, bgcols=c("lightgreen", "lightblue"), linecols=c("black", "red"), showMean=T, showSD=T){
    dfs <- split(df, df[,1])
    dfs <- dfs[order(unlist(lapply(dfs, nrow)), decreasing=T)]
    breakpoints <- cumsum(unlist(lapply(dfs, nrow)))+0.5
    ymax <- signif_ceiling(max(df[,-1]), 2)

    plot(0, xlim=c(0, nrow(df)), ylim=c(0, ymax), type="n")
    xl <- c(par("usr")[1], breakpoints[-length(breakpoints)])
    yb <- rep(par("usr")[3], length(dfs))
    xr <- c(breakpoints[-1], par("usr")[2])
    yt <- rep(par("usr")[4], length(dfs))
    rect(xl, yb, xr, yt, col=bgcols, border=NA)
    x = 0
    for(df in dfs){
        for(i in 2:ncol(df)){
            if(showMean){
                lines(c(x+1, x+nrow(df)), rep(mean(df[df[,i]>0,i]),2), col=linecols[i-1], lty=2)
            }
            if(showSD){
                lines(c(x+1, x+nrow(df)), rep(mean(df[df[,i]>0,i])-2*sd(df[df[,i]>0,i]),2), col=linecols[i-1], lty=3)
                lines(c(x+1, x+nrow(df)), rep(mean(df[df[,i]>0,i])+2*sd(df[df[,i]>0,i]),2), col=linecols[i-1], lty=3)
            }
            lines(x+1:nrow(df), df[,i], col=linecols[i-1])
        }
        x = x+nrow(df)
    }
}
