#' @include generics.R
NULL

################################################################################
#' Interface Class to access different types of operators on time series.
#'
#' \code{LagOperator} is an S4 class that provides a common interface to
#' implementations of an operator \eqn{\Gamma(Y)}{Gamma(Y)} which is calculated on 
#' all pairs of observations \eqn{(Y_0,Y_k)}{(Y0,Yk)} with lag smaller than maxLag 
#'
#' Currently one implementation is available:
#'     (1) \code{\link{ClippedCov}}.
#'
#' @name LagOperator-class
#' @aliases LagOperator
#'
#' @keywords S4-classes
#'
#' @slot values an array of dimension \code{c(maxLag,length(levels.1),length(levels.2))} containing the values of the operator.
#' @slot Y is the time series the operator shall be applied to
#' @slot maxLag maximum lag between two Observations
#' @slot levels a vector of numerics that determines the levels of the operator
#' @slot isRankBased A flag that is \code{FALSE} if the determined \code{values}
#'                     are based on the original time series and \code{TRUE} if it
#'                     is based on the ranks.
#'
################################################################################


setClass(
  Class = "LagOperator",
  representation=representation(
    Y = "numeric",
    values = "array",
    maxLag = "numeric",
    levels.1 = "numeric",
    levels.2 = "numeric",
    isRankBased = "logical"
    )
)

################################################################################
#' Get attribute \code{values} from a \code{LagOperator}.
#'
#' @name getValues-LagOperator 
#' @aliases getValues,LagOperator-method
#'
#' @keywords Access-functions
#'
#' @param object \code{LagOperator} from which to get the \code{values}.
#' @param levels.1 the first vector of levels for which to get the values
#' @param levels.2 the second vector of levels for which to get the values    
#' @return Returns the \code{values} attribute.
################################################################################
setMethod(f = "getValues",
          signature = "LagOperator",
          definition = function(object,levels.1,levels.2) {
            
            # workaround: default values don't seem to work for generic functions?
            if (!hasArg(levels.1)) {
              levels.1 <- object@levels.1
            }
            if (!hasArg(levels.2)) {
              levels.2 <- object@levels.2
            }
            # end: workaround
            
            
            # Select columns
        
            c.1.pos <- match(levels.1,object@levels.1)
            c.2.pos <- match(levels.2,object@levels.2)
            
            if(is.na(c.1.pos[1])){
              stop("no 'values' for 'levels.1' requested were found")
            }
            if(is.na(c.2.pos[1])){
              stop("no 'values' for 'levels.2' requested were found")
            }
            
            if(!(length(c.1.pos)==length(levels.1))){
              warning("not all requested 'levels.1' were found")
            }
            
            if(!(length(c.2.pos)==length(levels.2))){
              warning("not all requested 'levels.2' were found")
            }
            
            return(object@values[,c.1.pos,c.2.pos])
          }
)
################################################################################
#' Get attribute \code{levels} from a \code{LagOperator}.
#'
#' If the optional parameter \code{j} is supplied, then the \code{j}th vector of
#' levels will be returned, a list with all vectors otherwise.
#'
#' @name getLevels-LagOperator
#' @aliases getLevels,LagOperator-method
#'
#' @keywords Access-functions
#'
#' @param object \code{LagOperator} from which to get the \code{levels}.
#' @param j Index pointing to a set of levels in the list; optional.
#'
#' @return Returns levels attribute, as a vector of real numbers.
################################################################################
setMethod(f = "getLevels",
          signature = "LagOperator",
          definition = function(object,j) {
            if (missing("j")) {
              levels = list(object@levels.1,object@levels.2)
              names(levels) = c("levels.1","levels.2")
              return(levels)
            } else {
              if (!(j==1 | j==2)) {
                error("Index needs to be either 1 or 2.")
              } else {
                if(j==1){return(object@levels.1)}
                if(j==2){return(object@levels.2)}
              }
            }
          }
)
