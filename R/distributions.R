poisratios <- function(lambda1, lambda2, min1=10, nopoints=100) {
	x <- rpois(nopoints, lambda1);
	y <- rpois(nopoints, lambda2);
	ratio <- ifelse(x < min1, NA, y/x);
	xy <- cbind(x, y, ratio);
	return(xy);
}

altfreq <- function(coverage, pi) {
	if (runif(1) > 0.5) {
		pi <- 1.0 - pi;
	}
	altreads <- rbinom(1, coverage, pi);
	return(altreads/coverage);
}
