#' Plots predicted rcure object
#'
#' @description Plots predicted survival curve(s) from the estimated robust cure model
#' @param object an object of the predictrcure function
#' @param type type of plot. "S" means steps plot.
#' @param xlab label for the x axis
#' @param ylab label for the y axis
#' @param model specifies survival part in cure model, "ph" or "aft"
#' @param ... further arguments to be passed to the plotpredictrcure function
#' @references Cai, C., Zou, Y., Peng, Y., & Zhang, J. (2012). smcure: An R-Package for estimating semiparametric mixture cure models. Computer methods and programs in biomedicine, 108(3), 1255-1260
#' @importFrom graphics lines matplot plot
#' @export
#'
plotpredictrcure <-
function(object, type="S", xlab="Time",ylab="Predicted Survival Probability",model=c("ph","aft"), ...)
{
  pred <- object$prediction
  if(model=="ph"){
    pdsort <- pred[order(pred[,"Time"]),]
    if(length(object$newuncureprob)==1) plot(pdsort[,"Time"],pdsort[,1], type="S")
    else
      matplot(pdsort[,"Time"],pdsort[,1:(ncol(pred)-1)],col=1,type="S",lty=1:(ncol(pred)-1),xlab=xlab,ylab=ylab)
    }
  if(model=="aft"){
    nplot=ncol(pred)/2
    pdsort <- pred[order(pred[,1+nplot]),c(1,1+nplot)]
    plot(pdsort[,2],pdsort[,1],xlab=xlab,ylab=ylab,col=1,type="S",ylim=c(0,1))
    if(nplot>1){
      for(i in 2:nplot){
        pdsort<- pred[order(pred[,i+nplot]),c(i,i+nplot)]
        lines(pdsort[,2],pdsort[,1],lty=i,type="S")
        }
      }
    }
  }

