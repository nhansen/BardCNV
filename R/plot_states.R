plot_states <- function(df, plotchr, xlow=0, xhigh=0, muratio=0.0, sample=c(""), minscore=0, contam=0.0) {
    attach(df);
    newdf <- df[chr==plotchr & score>=minscore,]
    attach(newdf);

    if (xlow==0) {
        xlow <- min(pos)
    }

    if (xhigh==0) {
        xhigh <- max(pos)
    }

    par(mfrow=c(3,1));
    plot(pos, depthratio, pch=20, col=rgb(255,0,0,100, maxColorValue=255), xlim=c(xlow, xhigh), main=paste(sample, "Depth of Coverage Ratio on", plotchr));
    if (muratio != 0) {
       expratio <- muratio*total*0.5*(1.0 - contam) + contam * muratio;
       points(pos, expratio, col=c("darkgreen"), type="s")
    }
    ymax <- max(total) + 1
    explow <- 1.0*minor/total*(1.0 - contam) + contam*0.5
    exphigh <- 1.0 - explow
    plot(pos, pialt, pch=20, col=rgb(100,100,100,100, maxColorValue=255), xlim=c(xlow, xhigh), main="Percentage Alternate Allele Among Reads");
    points(pos, explow, col=c("darkgreen"), type="s")
    points(pos, exphigh, col=c("darkgreen"), type="s")
    plot(pos, total, ylim=c(0,ymax), type="s", col=c("darkblue"), xlim=c(xlow, xhigh), main="Predicted States");
    points(pos, minor, type="s", col=c("orange"));
    legend(xlow+(xhigh-xlow)*6.5/10, 0.53*ymax, c("total copies", "minor allele copies"), col=c("darkblue", "orange"), lwd=c(1,1), bty="n")
}

read_states <- function(filename) {
    statesframe <- read.table(filename, sep="\t", header=TRUE);
    return(statesframe);
}
