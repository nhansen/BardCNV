plot_states <- function(df, plotchr, xlow=0, xhigh=0, muratio=0.0, sample=c(""), minscore=0, contam=0.0, nolines=FALSE, segs=0) {
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
    plot(pos, depthratio, pch=20, col=rgb(150,150,150,50, maxColorValue=255), xlim=c(xlow, xhigh), main=paste(sample, "Depth of Coverage Ratio on", plotchr));
    if (!nolines) {
        if (muratio != 0) {
           expratio <- muratio*total*0.5*(1.0 - contam) + contam * muratio;
           points(pos, expratio, col=c("darkgreen"), type="s")
        }
    }
    if (length(segs) > 1) {
        add_segs(segs, muratio=muratio, contam=contam, segcol=c("darkblue"), chr=plotchr)
    }
    ymax <- max(total) + 1
    explow <- 1.0*minor/total*(1.0 - contam) + contam*0.5
    exphigh <- 1.0 - explow
    plot(pos, pialt, pch=20, col=rgb(150,150,150,50, maxColorValue=255), xlim=c(xlow, xhigh), main="Percentage Alternate Allele Among Reads");
    if (!nolines) {
        points(pos, explow, col=c("darkgreen"), type="s")
        points(pos, exphigh, col=c("darkgreen"), type="s")
    }
    else {
        if (length(segs) > 1) {
            add_segs(segs, muratio=muratio, contam=contam, segcol=c("red"), chr=plotchr, pi=1)
        }
    }
    plot(pos, total, ylim=c(0,ymax), type="s", col=c("darkblue"), xlim=c(xlow, xhigh), main="Predicted States");
    points(pos, minor, type="s", col=c("orange"));
    legend(xlow+(xhigh-xlow)*6.5/10, 0.53*ymax, c("total copies", "minor allele copies"), col=c("darkblue", "orange"), lwd=c(1,1), bty="n")
    detach(df);
    detach(newdf);
}

read_states <- function(filename) {
    statesframe <- read.table(filename, sep="\t", header=TRUE);
    return(statesframe);
}

read_segs <- function(filename) {
    segs <- read.table(filename, sep=" ", header=FALSE);
    return(segs);
}

plot_segs <- function(segs, segcol=c("black"), segchr=c("1")) {
    chrsegs <- segs[segs[[1]]==segchr, ];
    x0 <- chrsegs[[2]];
    y0 <- chrsegs[[4]];
    x1 <- chrsegs[[3]];
    y1 <- chrsegs[[4]];
    plot(x0[[1]], y0[[1]], xlim=range(chrsegs[[2]]), ylim=range(chrsegs[[4]]), col=c("white"));
    segments(x0, y0, x1, y1, col=segcol);
}

add_segs <- function(segs, segcol=c("black"), addchr=c("1"), muratio=1.0, contam=0.01, pi=0) {
    chrsegs <- segs[segs[[1]]==addchr, ];
    minorcount <- chrsegs[[4]];
    minorcount[minorcount=="."] <- 1;
    minorcount <- as.numeric(as.character(minorcount));
    majorcount <- chrsegs[[5]];
    x0 <- chrsegs[[2]];
    x1 <- chrsegs[[3]];
    if (pi == 0) {
        y0 <- muratio*majorcount*0.5*(1.0-contam) + contam*muratio;
        y1 <- y0;
        segments(x0, y0, x1, y1, col=segcol);
        y0 <- muratio*minorcount*0.5*(1.0-contam) + contam*muratio*0.5;
        y1 <- y0;
        segments(x0, y0, x1, y1, col=rgb(0,0,100,100,maxColorValue=255), lwd=2.5);
    }
    else {
        #explowstart <- 1.0*minorcount/(majorcount+0.01)*(1.0 - contam) + contam*0.5;
        explowstart <- (1.0*minorcount*(1.0-contam) + 1.0*contam)/(majorcount*(1.0-contam) + 2.0*contam);
        explowend <- explowstart;
        exphighstart <- 1.0 - explowstart;
        exphighend <- exphighstart;
        segments(x0, explowstart, x1, explowend, col=c("black"));
        segments(x0, exphighstart, x1, exphighend, col=c("black"));
    }
}

cell_contam <- function(read_contam, tumor_ploidy) {
    contam <- read_contam*tumor_ploidy/(read_contam*tumor_ploidy + 2.0 - 2.0*read_contam);
    return(contam);
}

mu_1143_n0 <- 0.454294;
mu_1143_n5 <- 0.547722;
mu_1143_n20 <- 0.615508;
mu_1143_n40 <- 0.708369;
mu_1143_n60 <- 0.795230;
mu_1143_n80 <- 0.874237;
ct_1143_n0 <- 0.01;
ct_1143_n5 <- 0.087657;
ct_1143_n20 <- 0.313364;
ct_1143_n40 <- 0.548940;
ct_1143_n60 <- 0.732495;
ct_1143_n80 <- 0.879547;
mu_1954_n0 <- 0.363883;
mu_1954_n5 <- 0.460782;
mu_1954_n20 <- 0.547224;
mu_1954_n40 <- 0.916383;
mu_1954_n60 <- 0.919437;
mu_1954_n80 <- 0.943711;
ct_1954_n0 <- 0.01;
ct_1954_n5 <- 0.105082;
ct_1954_n20 <- 0.358049;
ct_1954_n40 <- 0.597963;
ct_1954_n60 <- 0.769930;
ct_1954_n80 <- 0.899234;
