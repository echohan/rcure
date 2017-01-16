#' Retrieves the estimated coefficients from rcure object
#' @description Retrieves the estimated coefficients from rcure object
#' @param x an object of rcure
#' @param ... further arguments to be passed to the coefrcure function
#' @references Cai, C., Zou, Y., Peng, Y., & Zhang, J. (2012). smcure: An R-Package for estimating semiparametric mixture cure models. Computer methods and programs in biomedicine, 108(3), 1255-1260
#' @export

coefrcure <-
function(x, ...)
{
	coef <- c(x$b,x$beta)
  names(coef) <- c(x$bnm,x$betanm)
  coef
}

