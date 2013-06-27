read_segs <- function(filename) {
    segs <- read.table(filename, sep=" ", header=FALSE);
    return(segs);
}

plot_segs <- function(segs, segcol=c("black"), chr=c("1")) {
    chrsegs <- segs[segs[[1]]==chr, ];
    x0 <- chrsegs[[2]];
    y0 <- chrsegs[[4]];
    x1 <- chrsegs[[3]];
    y1 <- chrsegs[[4]];
    plot(x0[[1]], y0[[1]], xlim=range(chrsegs[[2]]), ylim=range(chrsegs[[4]]), col=c("white"));
    segments(x0, y0, x1, y1, col=segcol);
}

add_segs <- function(segs, segcol=c("black"), chr=c("1"), muratio=1.0, contam=0.01, pi=0) {
    chrsegs <- segs[segs[[1]]==chr, ];
    minorcount <- chrsegs[[4]];
    minorcount[minorcount=="."] <- 1;
    minorcount <- as.numeric(as.character(minorcount));
    majorcount <- chrsegs[[5]];
    x0 <- chrsegs[[2]];
    x1 <- chrsegs[[3]];
    if (pi == 0) {
        y0 <- muratio*majorcount*0.5*(1.0-contam) + contam*muratio*0.5;
        y1 <- y0;
        segments(x0, y0, x1, y1, col=segcol);
        y0 <- muratio*minorcount*0.5*(1.0-contam) + contam*muratio*0.5;
        y1 <- y0;
        segments(x0, y0, x1, y1, col=rgb(0,0,100,100,maxColorValue=255), lwd=2.5);
    }
    else {
        explowstart <- 1.0*minorcount/(majorcount+0.01)*(1.0 - contam) + contam*0.5;
        explowend <- explowstart;
        exphighstart <- 1.0 - explowstart;
        exphighend <- exphighstart;
        segments(x0, explowstart, x1, explowend, col=c("black"));
        segments(x0, exphighstart, x1, exphighend, col=c("black"));
    }
}

