getSignal <-
function(chr.pileup, window) {
start <- window[1]
end <- window[2]
cov.window <- as.numeric( chr.pileup[start : end] )
    return( cov.window ) #numeric vector
}

leftBound <-
function(peak.x,pileup) {
    pile.derl1 <- c(0, diff(pileup))
    count <- peak.x
    fx2 <- pile.derl1[peak.x]
    x2 <- peak.x
    fx1 <- 0
    x1 <- 0
    while(count >= 2){
        count = count - 1
        if(count==1)return(x2)
        if(pile.derl1[count] != 0){
            fx1 <- pile.derl1[count]
            x1 <- count 
            if( (fx1 < 0) & (fx2 > 0) )return(x2)
            else{#if((fx1>=0)&(fx2>0)){
                fx2 <- fx1
                x2 <- x1
            }
        }
    }
}

rightBound <-
function(peak.x, pileup) {
    pile.der1 <- c(0, diff(pileup))
    count <- peak.x
    fx1 <- pile.der1[peak.x]
    x1 <- peak.x
    fx2 <- 0
    x2 <- 0
    support.bound <- length(pileup)
    while(count < support.bound){
        count <- count + 1
        if(count == support.bound)return(x1)
        if(pile.der1[count] != 0) {
            fx2 <- pile.der1[count]
            x2 <- count
            if((fx1 < 0) & (fx2 > 0))return(x1)
            else{
                fx1 <- fx2
               x1 <- x2
                }
            }
}
}


wavCWTPeaksCF <-
function (x, snr.min = 3, scale.range = NULL, length.min = 10,
    noise.span = NULL, noise.fun = "quantile", noise.min = NULL) {
    if (!is(x, "wavCWTTree"))
        stop("Input must be an object of class wavCWTTree")
    xatt <- attributes(x)
    endtimes <- attr(x, "endtime")
    times <- attr(x, "time")
    scale <- attr(x, "scale")
    noise <- attr(x, "noise")
    wavelet <- attr(x, "wavelet")
    series <- attr(x, "series")
    branch.hist <- attr(x, "branch.hist")
    sampling.interval <- abs(diff(times[1:2]))
    if (!is.element(wavelet, "gaussian2"))
        stop("Only CWT developed using the Mexican hat (gaussian2) filter are supported")
    if (is.null(noise.min))
        noise.min <- quantile(abs(attr(x, "noise")), prob = 0.05)
    if (is.null(scale.range))
        scale.range <- scale[range(which(branch.hist > quantile(branch.hist,
            prob = 0.8)))]
    if (is.null(noise.span))
        noise.span <- max(0.01 * diff(range(times)), 5 * sampling.interval)
    noise.levels <- unlist(lapply(endtimes, function(x, noise.fun,
        times, times.range, noise, noise.min, noise.span) {
        time.start <- x - noise.span
        if (time.start < times.range[1])
            time.start <- times.range[1]
        time.end <- x + noise.span
        if (time.end < times.range[2])
            time.end <- times.range[2]
        ix <- which(times >= time.start & times <= time.end)
        noise.local <- noise.fun(abs(noise[ix]))
        if (noise.local < noise.min)
            noise.local <- noise.min
        noise.local
    }, noise.fun = switch(noise.fun, quantile = function(x) {
        quantile(x, probs = 0.95)
    }, sd = sd, mad = function(x) {
        mad(x, center = 0)
    }), times = times, times.range = range(times), noise = noise,
        noise.min = noise.min, noise.span = noise.span))
    tmpargs <- lapply(x, function(x) unlist(lapply(x, function(x,
        imax) x[imax], imax = which.max(x$extrema))))
    peaks <- data.frame(do.call("rbind", tmpargs))
    peaks <- cbind(data.frame(branch = row.names(peaks)), peaks,
        data.frame(iendtime = attr(x, "iendtime")))
    peak.snr <- peaks[["extrema"]]/noise.levels
    peak.scale <- peaks[["scale"]]
    branch.lengths <- unlist(lapply(x, function(x, scale.range) length(which(x$scale >=
        scale.range[1] & x$scale <= scale.range[2])), scale.range = scale.range))
    good.snr <- peak.snr >= snr.min
    good.scale <- peak.scale >= scale.range[1]
    good.length <- branch.lengths >= length.min
    iendtime.min <- max(as.integer(noise.span/sampling.interval/4),
        3)
    iendtime.max <- length(times) - iendtime.min + 1
    good.end <- peaks[["iendtime"]] > iendtime.min & peaks[["iendtime"]] <
        iendtime.max
    peaks <- peaks[which(good.snr & good.scale & good.length &
        good.end), ]
    if(nrow(peaks) == 0)
      return(c())
    row.names(peaks) <- as.character(seq(nrow(peaks)))
    z <- list(x = times[peaks$iendtime], y = series[peaks$iendtime])
    attr(z, "peaks") <- peaks
    attr(z, "snr.min") <- snr.min
    attr(z, "scale.range") <- scale.range
    attr(z, "length.min") <- length.min
    attr(z, "noise.span") <- noise.span
    attr(z, "noise.fun") <- noise.fun
    attr(z, "noise.min") <- noise.min
    return( z )
}
