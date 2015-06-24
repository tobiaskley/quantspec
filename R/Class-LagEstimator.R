#' @include generics.R
#' @include Class-Weight.R
#' @include Class-LagOperator.R
#' 
NULL
################################################################################
#' Class for a lag-window type estimator.
#'
#' For a given time series Y a lag-window estimator of the Form
#' \deqn{\hat{f}(\omega) = \sum_{|k|< n-1 } K_n(k) \Gamma(Y_0,Y_k) \exp(-i \omega k)} 
#' will be calculated on initalization. The \code{LagKernelWeight} K_n is determined
#' by the slot \code{weight} and the \code{LagOperator} \eqn{\Gamma(Y_0,Y_k)} is defined 
#' by the slot lagOp.
#' 
#'
#' @name LagEstimator-class
#' @aliases LagEstimator
#' @exportClass LagEstimator
#'
#' @keywords S4-classes
#'
#' @slot Y the time series where the lag estimator was applied one
#' @slot weight a \code{\link{Weight}} object to be used as lag window
#' @slot lagOp a \code{\link{LagOperator}} object that determines which
#'             kind of bivariate structure should be calculated.
#'
################################################################################
setClass(
  Class = "LagEstimator",
  representation=representation(
    Y = "numeric",
    weight = "Weight",
    lagOp = "LagOperator"   
  ),
  contains = "QSpecQuantity"
)

setMethod(
  f = "initialize",
  signature = "LagEstimator",
  definition = function(.Object, Y,frequencies,lagOp,weight){
   
    levels.1 = getLevels(lagOp,1)
    levels.2 = getLevels(lagOp,2)
    
    .Object@Y = Y
    .Object@frequencies = frequenciesValidator(frequencies, length(Y), steps=1:6)
    .Object@levels[[1]] = levels.1
    .Object@levels[[2]] = levels.2
    .Object@weight = weight
    .Object@lagOp = lagOp
    
    levels.all = union(levels.1,levels.2)
    ln1 = length(levels.1)
    ln2 = length(levels.2)
    ln = length(levels.all)
    Q <- Y
    Q <- quantile(Y,probs = levels.all)
    
    
    res <- array(dim = c(length(Y),ln1,ln2))
    
    Fhat = getValues(lagOp,levels.1 = levels.1,levels.2 = levels.2)
    n = length(Fhat)/(ln1*ln2)
    Kernel = getValues(weight)
    if(length(Kernel)<n){
      warning("number of 'weights' is too small, filled up with zeroes")
      Kernel = c(Kernel,rep(0,length(Fhat[,1,1])-length(Kernel)))
    }
    if(length(Kernel)>length(Fhat[,1,1])){
      Kernel = Kernel[1:length(Fhat[,1,1])]
    }
    
    for (l1 in 1:ln1) {
      for (l2 in 1:ln2) {  #\tau_2,\tau_1 can be calculated from \tau_1,\tau_2 
        if((l1 > l2)&&!(levels.1[l1]==levels.2[l2])&&is.element(levels.1[l1],levels.2)&&is.element(levels.2[l2],levels.1)){ 
          res[,l1,l2] <- Conj(res[,l2,l1])
          next
        }
        A1 <- Fhat[,l1,l2]
        A1 <- A1# - levels.1[l1] * levels.2[l2]
        A1 <- A1*Kernel
        
        A2 <- Fhat[,l2,l1]
        A2 <- A2# - levels.1[l2] * levels.2[l1]  				 # same as above!
        A2 <- A2*Kernel
        
        A0 <- Fhat[1,l2,l1]
        A0 <- A0# - levels.1[l2] * levels.2[l1]					 # same as above!
        A0 <- A0*Kernel[1]
        res[,l1,l2] <- 1/(2*pi) * (fft(A1) + Conj(fft(A2)) - A0)
      }
    }
    .Object@values = res
    return(.Object)
  })

################################################################################
#' Get values from a lag-window type estimator.
#'
#' The returned array of \code{values} is of dimension \code{[J,K1,K2]},
#' where \code{J=length(frequencies)}, \code{K1=length(levels.1)} and
#' \code{K2=length(levels.2))}. At position \code{(j,k1,k2)}
#' the returned value is the one corresponding to \code{frequencies[j]},
#' \code{levels.1[k1]} and \code{levels.2[k2]} that are closest to the
#' \code{frequencies}, \code{levels.1} and \code{levels.2}
#' available in \code{object}; \code{\link{closest.pos}} is used to determine
#' what closest to means. 
#' 
#' @name getValues-LagEstimator
#' @aliases getValues,LagEstimator-method
#'
#' @keywords Access-functions
#'
#' @param object \code{LagEstimator} of which to get the values
#' @param frequencies a vector of frequencies for which to get the values
#' @param levels.1 the first vector of levels for which to get the values
#' @param levels.2 the second vector of levels for which to get the values
#'
#' @return Returns data from the array \code{values} that's a slot of
#'          \code{object}.
#'
#' @seealso
#' An example on how to use this function is analogously to the example given in
#' \code{\link{getValues-QuantilePG}}.
################################################################################

setMethod(f = "getValues",
          signature = signature(
            "LagEstimator"),
          definition = function(object,
                                frequencies=2*pi*(0:(length(object@Y)-1))/length(object@Y),
                                levels.1=getLevels(object,1),
                                levels.2=getLevels(object,2)) {
            
            # workaround: default values don't seem to work for generic functions?
            if (!hasArg(frequencies)) {
              frequencies <- 2*pi*(0:(length(object@Y)-1))/length(object@Y)
            }
            if (!hasArg(levels.1)) {
              levels.1 <- object@levels[[1]]
            }
            if (!hasArg(levels.2)) {
              levels.2 <- object@levels[[2]]
            }
            # end: workaround
            
            ##############################
            ## (Similar) Code also in Class-FreqRep!!!
            ## (Similar) Code also in Class-QuantileSD!!!
            ##############################
            
            # Transform all frequencies to [0,2pi)
            frequencies <- frequencies %% (2*pi)
            
            # Create an aux vector with all available frequencies
            oF <- object@frequencies
            f <- frequencies
            
            # returns TRUE if x c y
            subsetequal.approx <- function(x,y) {
              X <- round(x, .Machine$double.exponent-2)
              Y <- round(y, .Machine$double.exponent-2)
              return(setequal(X,intersect(X,Y)))
            }
            
            C1 <- subsetequal.approx(f[f <= pi], oF)
            C2 <- subsetequal.approx(f[f > pi], 2*pi - oF[which(oF != 0 & oF != pi)])
            
            if (!(C1 & C2)) {
              warning("Not all 'values' for 'frequencies' requested were available. 'values' for the next available Fourier frequencies are returned.")
            }
            
            # Select columns
            c.1.pos <- closest.pos(object@levels[[1]],levels.1)
            c.2.pos <- closest.pos(object@levels[[2]],levels.2)
            
            if (!subsetequal.approx(levels.1, object@levels[[1]])) {
              warning("Not all 'values' for 'levels.1' requested were available. 'values' for the next available level are returned.")
            }
            
            if (!subsetequal.approx(levels.2, object@levels[[2]])) {
              warning("Not all 'values' for 'levels.2' requested were available. 'values' for the next available level are returned.")
            }
            
            
            J <- length(frequencies)
            K1 <- length(levels.1)
            K2 <- length(levels.2)
            res <- array(dim=c(J, K1, K2))
            
            
            if (class(object@weight) == "LagKernelWeight") {
              
              # Select rows
              r1.pos <- closest.pos(oF, f[f <= pi])
              r2.pos <- closest.pos(-1*(2*pi-oF),-1*f[f > pi])
              
              if (length(r1.pos) > 0) {
                res[which(f <= pi),,] <- object@values[r1.pos,c.1.pos,c.2.pos]
              }
              if (length(r2.pos) > 0) {
                res[which(f > pi),,] <- Conj(object@values[r2.pos,c.1.pos,c.2.pos])
              }
              
            }
            
            return(res)
          }
)

  ################################################################################
  #' Create an instance of the \code{LagEstimator} class.
  #'
  #' A \code{LagEstimator} object can be created from \code{numeric}, a \code{ts},
  #' or a \code{zoo} object. Also a \code{\link{LagOperator}} and a 
  #' \code{\link{Weight}} object can be used to create different types of 
  #' estimators.
  #'
  #'
  #' @name LagEstimator-constructor
  #' @aliases lagEstimator
  #' @export
  #'
  #' @keywords Constructors
  #'
  #' @param Y a time series (\code{numeric}, \code{ts}, or \code{zoo} object) or a 
  #'          \code{\link{LagOperator}} from which to determine the \code{LagEstimator}
  #' @param frequencies A vector containing (Fourier-)frequencies at which to determine the
  #'                     smoothed periodogram.
  #' @param levels.1 the first vector of levels for which to compute the LagEstimator
  #' @param levels.2 the second vector of levels for which to compute the LagEstimator               
  #' @param weight Object of type \code{\link{Weight}} to be used for smoothing.
  #' @param type if \code{Y} is a time series, this indicates which LagOperator will be used 
  #' @return Returns an instance of \code{LagEstimator}.
  #'
  #' @examples
  #' Y <- rnorm(100)
  #' levels.1 <- c(0.1,0.5,0.9)
  #' weight <- lagKernelWeight(W = WParzen,  bw = 10, K = length(Y))
  #' lagOp <- clippedCov(Y,levels.1 = levels.1)
  #' lagEst <- lagEstimator(lagOp,weight = weight)
################################################################################
lagEstimator <- function(Y,
                         frequencies=2*pi/length(Y) * 0:(length(Y)-1),
                         levels.1 = .5,
                         levels.2 = levels.1,
                         weight = lagKernelWeight(K = length(Y),bw = 25),
                         type = c("clippedCov")
                         ){
  
  #ToDo Check if params are okay
  #     Transform to Fourier frequencies
  if (class(Y) == "numeric") {
    versConstr <- 1
    Y <- Y
  } else if (class(Y) == "ts") {
    versConstr <- 1
    Y <- Y[1:length(Y)]
  } else if (class(Y) == "zoo") {
    versConstr <- 1
    Y <- coredata(Y)
  } else if ( is(Y,"LagOperator")) {
    LagOp = Y
    versConstr <- 2
    if (!hasArg(frequencies)) {
      Z <- Y@Y
      frequencies <- 2*pi/length(Z) * 0:(length(Z)-1)
    }
    
    if (!hasArg(levels.1)) {
      levels.1 <- getLevels(Y,1)
    }
    if (!hasArg(levels.2)) {
      levels.2 <- getLevels(Y,2)
    }
    Y = Z
  } else {
    stop("object is neither 'numeric', 'ts', 'zoo', nor 'LagOperator'.")
  }
  if(versConstr == 1){
    if (type == "clippedCov"){
      lagOp = clippedCov(Y,levels.1 = levels.1, levels.2 = levels.2)  
    }
    
  }

   return(new(Class = "LagEstimator",Y = Y,frequencies = frequencies,weight = weight, lagOp = lagOp))
  
  
}
  
################################################################################
#' Plot the values of a \code{\link{LagEstimator}}.
#'
#' Creates a \code{K} x \code{K} plot displaying all levels combinations from the
#' argument \code{levels}.  
#' In each of the subplots either the real part (on and below the diagonal;
#' i. e., \eqn{\tau_1 \leq \tau_2}{tau1 <= tau2}) or the imaginary parts
#' (above the diagonal; i. e., \eqn{\tau_1 > \tau_2}{tau1 > tau2}) of
#' the lag-window estimator, for the combination of levels \eqn{\tau_1}{tau1}
#'  and \eqn{\tau_2}{tau2} denoted on the left and bottom margin of the plot are displayed.
#'
#' @name plot-LagEstimator
#' @aliases plot,LagEstimator,ANY-method
#' @export
#'
#' @importFrom abind abind
#'
#' @param x  The \code{\link{LagEstimator}} object to plot
#' @param ptw.CIs the confidence level for the confidence intervals to be
#'                 displayed; must be a number from [0,1]; if null, then no
#'                 confidence intervals will be plotted.
#'                 TO BE DONE!
#' @param ratio quotient of width over height of the subplots; use this
#'               parameter to produce landscape or portrait shaped plots.
#' @param widthlab width for the labels (left and bottom); default is
#'                  \code{lcm(1)}, cf. \code{\link[graphics]{layout}}.
#' @param xlab label that will be shown on the bottom of the plots; can be
#'              an expression (for formulas), characters or \code{NULL} to
#'              force omission (to save space).
#' @param ylab label that will be shown on the left side of the plots;
#'              can be an expression (for formulas), characters or
#'              \code{NULL} to force omission (to save space).
#' @param type.scaling a method for scaling of the subplots; currently there
#'                      are three options: \code{"individual"} will scale each of the
#'                      \code{K^2} subplots to minimum and maximum of the values
#'                      in that plot, \code{"real-imaginary"} will scale each of the
#'                      subplots displaying real parts and each of the subplots
#'                      displaying imaginary parts to the minimum and maximum of
#'                      the values display in these subportion of plots. The
#'                      option \code{"all"} will scale the subplots to the minimum and
#'                      maximum in all of the subplots.
#' @param frequencies a set of frequencies for which the values are to be
#'                    plotted.
#' @param levels a set of levels for which the values are to be plotted.
#'
#' @return Returns the plot described in the Description section.
#' 
#' See Birr et al. (2015)
#' @references
#' Birr, S., Volgushev, S., Kley, T., Dette, H. & Hallin, M. (2015).
#' Quantile Spectral Analysis for Locally Stationary Time Series.
#' \url{http://arxiv.org/abs/1404.4605}.
################################################################################

setMethod(f = "plot",
          signature = signature("LagEstimator"),
          definition = function(x,
                                ptw.CIs = 0, 
                                ratio = 3/2, widthlab = lcm(1), xlab = expression(omega/2*pi), ylab = NULL,
                                type.scaling = c("individual", "real-imaginary", "all"),
                                frequencies=x@frequencies,
                                levels=intersect(x@levels[[1]], x@levels[[2]])) {
            
            def.par <- par(no.readonly = TRUE) # save default, for resetting...
            
            # workaround: default values don't seem to work for generic functions?
            if (!hasArg(frequencies)) {
              frequencies <- x@frequencies
            }
            if (!hasArg(levels)) {
              levels <- intersect(x@levels[[1]], x@levels[[2]])
            }
            if (!hasArg(ptw.CIs)) {
              ptw.CIs <- 0
            }
            if (!hasArg(ratio)) {
              ratio <- 3/2
            }
            if (!hasArg(widthlab)) {
              widthlab <- lcm(1)
            }
            if (!hasArg(xlab)) {
              xlab <- expression(omega/2*pi)
            }
            if (!hasArg(ylab)) {
              ylab <- NULL
            }
            if (!hasArg(type.scaling)) {
              type.scaling <- c("individual", "real-imaginary", "all")
            }
            # end: workaround
            
            if (length(levels) == 0) {
              stop("There has to be at least one level to plot.")
            }
            
            tryCatch({
              
              K <- length(levels)
              values <- getValues(x, frequencies = frequencies,
                                  levels.1=levels, levels.2=levels)
            #text.headline <- x@weight@descr
              if (ptw.CIs > 0) {
                CI <- getPointwiseCIs(x, frequencies = frequencies,
                                      alpha=ptw.CIs, type=type.CIs,
                                      levels.1=levels, levels.2=levels)
                lowerCIs  <- CI$lowerCIs
                upperCIs  <- CI$upperCIs
                #text.headline <- (paste(text.headline, ", includes ",1-ptw.CIs,"-CI (ptw. of type '",type.CIs,"')",sep=""))
              }
              
              X <- frequencies/(2*pi)
              
              allVals <- array(values[,,], dim=c(length(X), K, K))
             
              if (ptw.CIs > 0) {
                allVals <- abind(allVals, lowerCIs, upperCIs, along=1)
              }
              type.scaling <- match.arg(type.scaling)[1]
              
              p <- K
              M <- matrix(1:p^2, ncol=p)
              M.heights <- rep(1,p)
              M.widths  <- rep(ratio,p)
              
              # Add places for tau labels
              M <- cbind((p^2+1):(p^2+p),M)
              M.widths <- c(widthlab,M.widths)
              M <- rbind(M,c(0,(p^2+p+1):(p^2+2*p)))
              M.heights <- c(M.heights, widthlab)
              
              i <- (p^2+2*p+1)
              # Add places for labels
              if (length(xlab)>0) {
                M.heights <- c(M.heights, widthlab)
                M <- rbind(M,c(rep(0,length(M.widths)-p),rep(i,p)))
                i <- i + 1
              }
              
              if (length(ylab)>0) {
                M <- cbind(c(rep(i,p),rep(0,length(M.heights)-p)),M)
                M.widths <- c(widthlab,M.widths)
              }
              
              nf <- layout(M, M.widths, M.heights, TRUE)
              
              par(mar=c(2,2,1,1))
              
              for (i1 in 1:K) {
                for (i2 in 1:K) {
                  if (i2 >= i1) {
                    switch(type.scaling,
                           "individual" = {
                             y.min <- min(Re(allVals[,i1,i2]))
                             y.max <- max(Re(allVals[,i1,i2]))},
                           "real-imaginary" = {
                             y.min <- min(Re(allVals))
                             y.max <- max(Re(allVals))},
                           "all" = {
                             y.min <- min(Re(allVals),Im(allVals))
                             y.max <- max(Re(allVals),Im(allVals))}
                    )
                    
                    plot(x=0,y=0, type="n", xlab="", ylab="", #xlab=xl, ylab=yl,
                         xlim=c(min(X), max(X)), ylim=c(y.min, y.max))
                    if (ptw.CIs > 0) {
                      polygon(x=c(X,rev(X)), y=c(Re(lowerCIs[,i1,i2]),rev(Re(upperCIs[,i1,i2]))),
                              col="lightgray", border=NA)
                    }
                    lines(x=X, y=Re(values[,i1,i2]),
                          ylim=c(min(Re(allVals)), max(Re(allVals))),
                          type="l", col="blue")
                  } else {
                    switch(type.scaling,
                           "individual" = {
                             y.min <- min(Im(allVals[,i1,i2]))
                             y.max <- max(Im(allVals[,i1,i2]))},
                           "real-imaginary" = {
                             y.min <- min(Im(allVals))
                             y.max <- max(Im(allVals))},
                           "all" = {
                             y.min <- min(Re(allVals),Im(allVals))
                             y.max <- max(Re(allVals),Im(allVals))}
                    )
                    plot(x=0,y=0, type="n", xlab="", ylab="", #xlab=xl, ylab=yl,
                         xlim=c(min(X), max(X)), ylim=c(y.min, y.max))
                    if (ptw.CIs > 0) {
                      polygon(x=c(X,rev(X)), y=c(Im(lowerCIs[,i1,i2]),rev(Im(upperCIs[,i1,i2]))),
                              col="lightgray", border=NA)
                    }
                    lines(x=X, y=Im(values[,i1,i2]),
                          ylim=c(min(Im(allVals)), max(Im(allVals))),
                          type="l", col="blue")
                  }
                }
              }
              
              par(mar=c(0,0,0,0))
              for (i in 1:p) {
                plot.new()
                text(0.5,0.5,substitute(paste(tau[1],"=",k),list(k=levels[i])), srt=90)
              }
              
              for (i in 1:p) {
                plot.new()
                text(0.5,0.5,substitute(paste(tau[2],"=",k),list(k=levels[i])))
              }
              if (length(xlab)>0) {
                plot.new()
                text(0.5, 0.5, xlab)
              }
              if (length(ylab)>0) {
                plot.new()
                text(0.5, 0.5, ylab, srt=90)
              }
              
            },  error = function(e) e,
            warning = function(w) w,
            finally = {
              par(def.par)  #- reset to default
            })
          }
)

  
  
