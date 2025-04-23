# server.R
library(shiny)
library(shinyjs)
library(dplyr)
#library(yogiroc)
#library(maveLLR)
library(rsconnect)
library(rmarkdown)
library(knitr)
library(tinytex)
library(httr)
library(jsonlite)
library(DT)
library(ggplot2)
library(shinycssloaders)
library(data.table)
library(future)
library(promises)




#' confidence interval for precision or recall
#'
#' @param is the number of correct calls
#' @param ns the number of total calls
#' @param p the confidence interval probability range e.g. 0.025;0.975
#' @param res the resolution at which to sample the call rates
#' 
#' @export
#'
#' @return the confidence interval (numerical vector)
prcCI <- function(is, ns, p=c(0.025,0.975),res=0.001) {
  stopifnot(length(is)==length(ns))
  do.call(rbind,mapply(function(i,n) {
    rates <- seq(0,1,res)
    dens <- dbinom(i,n,rates)
    cdf <- c(0,cumsum(dens[-1]*res)/sum(dens*res))
    setNames(sapply(p,function(.p) rates[max(which(cdf < .p))]),p)
  },is,ns,SIMPLIFY=FALSE))
}

#' Quick and dirty sampling of rate parameters of a binomial distribution
#' 
#' This function is used internally by samplePRCs() (below)
#' 
#' @param i numerator (number of successful events)
#' @param n denominator (total number of events)
#' @param N desired number of samples to generate
#' @param minQ a minimum constraint for the samples 
#'    Rates can't be smaller than this. This is useful to include prior information.
#' @param maxQ a maximum constarint for the samples
#' @return a vector of sampled rates
sampleRatesQD <- function(i,n,N=1000,minQ=0,maxQ=1) {
  #Simple rejection sampling
  rateSamples <- runif(N,min=minQ,max=maxQ)
  #rates are accepted if a uniform RV falls below their probability (dictated by binomial distribution)
  accept <- runif(N,min=0,max=dbinom(i,n,i/n)) < dbinom(i,n,rateSamples)
  #store accepted values
  out <- rateSamples[accept]
  #At this point we're likely short of the original N outputs we need.
  #How much more do we need sample to make our quota (N)?
  #Make 2x as many attempts as we expect to require
  M <- 2*ceiling(N/(sum(accept)/N))
  #generate the remaining samples as before
  rateSamples <- runif(M,min=rep(minQ,length.out=M),max=rep(maxQ,length.out=M))
  accept <- runif(M,min=0,max=dbinom(i,n,i/n)) < dbinom(i,n,rateSamples)
  out <- c(out,rateSamples[accept])
  return(head(out,N))
}

# rejection sampling method for single rates with constraints
rejSam <- function(i,n,minQ=0,maxQ=1) {
  x <- runif(1,min=minQ,max=maxQ)
  #add some shortcuts for extreme constraints
  if (minQ > i/n && dbinom(i,n,minQ) < 0.05) {
    return(minQ)
  }
  if (maxQ < i/n && dbinom(i,n,maxQ) < 0.05) {
    return(maxQ)
  }
  while (runif(1,0,dbinom(i,n,i/n)) > dbinom(i,n,x)) {
    x <- runif(1,min=minQ,max=maxQ)
  }
  x
}

#' Slower, but order-preserving sampling of rate parameters of a binomial distribution
#' 
#' This function is used internally by samplePRCs() (below)
#' 
#' @param i numerator (number of successful events)
#' @param n denominator (total number of events)
#' @param N desired number of samples to generate
#' @param minQ a minimum constraint for the samples 
#'    Rates can't be smaller than this. This is useful to include prior information.
#' @param maxQ a maximum constarint for the samples
#' @return a vector of sampled rates
sampleRates <- function(i,n,N=1000,minQ=NA, maxQ=NA) {
  if (all(is.na(minQ)) && all(is.na(maxQ))) {
    return(rbeta(N,i,n-i))
  } else {
    if (!all(is.na(minQ))) {
      sapply(minQ, function(mq) {
        rejSam(i,n,minQ=mq)
      })
    } else if (!all(is.na(maxQ))) {
      sapply(maxQ, function(mq) {
        rejSam(i,n,maxQ=mq)
      })
    }
  }
}

#' Sample a distribution of PRC paths based, based on likelihood dictated by data
#' 
#' @param data a data table from a yr2 object
#' @param N the number of samples 
#' @param monotonized whether to use monotonization
#' @return a list of tables containin N samples of precision and recall each. 
#'    Each table corresponds to one row in the 'data' input
samplePRCs <- function(data,N=1000,monotonized=TRUE,sr=sampleRates) {
  pb <- txtProgressBar(max=nrow(data)-1,style=3)
  randomPaths <- list(
    cbind(
      precision=sr(data[1,"tp"],data[1,"tp"]+data[1,"fp"]),
      recall=sr(data[1,"tp"],data[1,"tp"]+data[1,"fn"])
    )
  )
  setTxtProgressBar(pb,1)
  for (k in 2:(nrow(data)-1)) {
    if (monotonized) {
      randomPaths[[k]] <- cbind(
        precision=sr(data[k,"tp"],data[k,"tp"]+data[k,"fp"],minQ = randomPaths[[k-1]][,"precision"]),
        recall=sr(data[k,"tp"],data[k,"tp"]+data[k,"fn"],maxQ = randomPaths[[k-1]][,"recall"])
      )
    } else {
      randomPaths[[k]] <- cbind(
        precision=sr(data[k,"tp"],data[k,"tp"]+data[k,"fp"]),
        recall=sr(data[k,"tp"],data[k,"tp"]+data[k,"fn"])
      )
    }
    setTxtProgressBar(pb,k)
  }
  return(randomPaths)
}

#' Use the output of samplePRCs to infer confidence intervals for a PRC curve
#' 
#' @param randomPaths the output of samplePRCs
#' @param nbins the number of bins along the recall axis to use
#' @return a table listing the confidence interval at each recall bin
inferPRCCI <- function(randomPaths,nbins=50) {
  randomSamples <- do.call(rbind,randomPaths)
  q5 <- yogitools::runningFunction(
    randomSamples[,"recall"],randomSamples[,"precision"],
    nbins=nbins,fun=function(xs)quantile(xs,.025)
  )
  q95 <- yogitools::runningFunction(
    randomSamples[,"recall"],randomSamples[,"precision"],
    nbins=nbins,fun=function(xs)quantile(xs,.975)
  )
  out <- cbind(q5,q95[,2])
  rownames(out) <- NULL
  colnames(out) <- c("recall","0.025","0.975")
  return(out)
}

# prcCI <- function(i,n,p=c(0.025,0.975),res=0.001) {
#   rates <- seq(0,1,res)
#   dens <- dbinom(i,n,rates)
#   cdf <- c(0,cumsum(dens[-1]*res)/sum(dens*res))
#   # plot(rates,cdf,type="l")
#   sapply(p,function(.p) rates[max(which(cdf < .p))])
# }

#' Helper function to monotonize precision
#'
#' @param xs numerical input vector, representing precision ordered according to increasing t
#'
#' @return the monotonized equivalent vector
monotonize <- function(xs) {
  for (i in 2:length(xs)) {
    if (xs[[i]] < xs[[i-1]]) {
      xs[[i]] <- xs[[i-1]]
    }
  }
  xs
}

#Balancing concept by Yingzhou Wu and Fritz Roth (Wu et al, unpublished) 
balance.prec <- function(ppv.prec,prior) {
  ppv.prec*(1-prior)/(ppv.prec*(1-prior)+(1-ppv.prec)*prior)
}

configure.prec <- function(sheet,monotonized=TRUE,balanced=FALSE) {
  ppv <- sheet[,"ppv.prec"]
  if (balanced) {
    prior <- sheet[1,"tp"]/(sheet[1,"tp"]+sheet[1,"fp"])
    ppv <- balance.prec(ppv,prior)
  } 
  if (monotonized) {
    ppv <- monotonize(ppv)
  }
  return(ppv)
}


#' YogiRoc2 object constructor
#'
#' @param truth a boolean vector indicating the classes of the reference set
#' @param scores a matrix of scores, with rows for each entry in truth, and one column for each predictor
#' @param names the names of the predictors
#' @param high a boolean vector indicating for each predictor whether its scoring high-to-low (or low-to-high)
#'
#' @return a yogiroc2 object
#' @export
#'
#' @examples
#' #generate fake data
#' truth <- c(rep(TRUE,10),rep(FALSE,8))
#' scores <- cbind(
#'   pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
#'   pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #draw PRC curve
#' draw.prc(yrobj)
#' #calculate recall at 90% precision
#' recall.at.prec(yrobj,0.9)
yr2 <- function(truth, scores, names=colnames(scores), high=TRUE) {
  
  #make sure all input is of correct datatype
  stopifnot(is.logical(truth), 
            is.data.frame(scores) || is.matrix(scores), 
            is.numeric(scores[1,1]), 
            is.character(names),
            is.logical(high)
  )
  #make sure all input is of correct size
  stopifnot(length(truth) == nrow(scores),
            length(names) == ncol(scores),
            length(high) == 1 || length(high) == ncol(scores)
  )
  
  #apply flipping to any scores that are not high-to-low
  if (length(high) == 1 && !high) {
    scores <- -scores
  } else if (any(!high)) {
    scores[,which(!high)] <- -scores[,which(!high)]
  }
  
  #calculate and return the roc/prc tables for each score
  tables <- setNames(lapply(1:ncol(scores), function(coli) {
    #the sample prior is the share of true cases out of all cases
    appl <- which(!is.na(scores[,coli]))
    prior <- sum(truth[appl])/length(truth[appl])
    #build a table for the ROC/PRC curves by iterating over all possible score thresholds
    ts <- na.omit(c(-Inf,sort(scores[,coli]),Inf))
    data <- do.call(rbind,lapply(ts, function(t) {
      #which scores fall above the current threshold?
      calls <- scores[,coli] >= t
      #calculate True Positives, True Negatives, False Positives and False Negatives
      tp <- sum(calls & truth,na.rm=TRUE)
      tn <- sum(!calls & !truth,na.rm=TRUE)
      fp <- sum(calls & !truth,na.rm=TRUE)
      fn <- sum(!calls & truth,na.rm=TRUE)
      #calculate PPV/precision, TPR/sensitivity/recall, and FPR/fallout
      ppv.prec <- tp/(tp+fp)
      tpr.sens <- tp/(tp+fn)
      fpr.fall <- fp/(tn+fp)
      # ppv.prec.balanced <- balance.prec(ppv.prec,prior)
      #return the results
      c(
        thresh=t,tp=tp,tn=tn,fp=fp,fn=fn,
        ppv.prec=ppv.prec,tpr.sens=tpr.sens,fpr.fall=fpr.fall
      )
    }))
    #set precision at infinite score threshold based on penultimate value
    # data[nrow(data),c("ppv.prec","ppv.prec.balanced")] <- data[nrow(data)-1,c("ppv.prec","ppv.prec.balanced")]
    data[nrow(data),"ppv.prec"] <- data[nrow(data)-1,"ppv.prec"]
    
    return(data)
    
  }),names)
  
  return(structure(tables,class="yr2"))
}


#' print method for yogiroc2 objects
#'
#' @param yr2 the object
#'
#' @return nothing, just prints a description
#' @export
#'
#' 
print.yr2 <- function(yr2) {
  cat("YogiROC object\n")
  cat("Reference set size:",nrow(yr2[[1]]-2),"\n")
  cat("Predictors:",paste(names(yr2),collapse=", "),"\n")
}


#' Draw a ROC curve
#'
#' @param yr2 an underlying yogiroc2 object
#' @param col the colors to use for the predictors
#' @param legend the positioning of the legend (e.g. "bottomright). NA to disable legend.
#' @param ... additional graphical parameters (see \code{par})
#'
#' @return nothing, draws a plot
#' @export
#'
#' @examples
#' #generate fake data
#' truth <- c(rep(TRUE,10),rep(FALSE,8))
#' scores <- cbind(
#'   pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
#'   pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #draw PRC curve
#' draw.roc(yrobj)
draw.roc <- function(yr2,col=seq_along(yr2),lty=1,legend="bottomright",...) {
  stopifnot(inherits(yr2,"yr2"))
  if (length(lty) < length(yr2)) {
    lty <- rep(lty,length(yr2))
  }
  plot(
    100*yr2[[1]][,"fpr.fall"],100*yr2[[1]][,"tpr.sens"],
    type="l",
    xlab="False positive rate (%)\n(= 100%-specificity)", ylab="Sensitivity or True positive rate (%)",
    xlim=c(0,100),ylim=c(0,100),col=col[[1]], lty=lty[[1]], ...
  )
  if(length(yr2) > 1) {
    for (i in 2:length(yr2)) {
      lines(
        100*yr2[[i]][,"fpr.fall"],100*yr2[[i]][,"tpr.sens"],
        col=col[[i]], lty=lty[[i]], ...
      )
    }
  }
  if (!is.na(legend)) {
    legend(legend,sprintf("%s (AUROC=%.02f)",names(yr2),auroc(yr2)),col=col,lty=lty)
  }
}


#' Draw Precision-Recall Curve (PRC)
#' 
#' Balancing concept by Yingzhou Wu and Fritz Roth (Wu et al, unpublished) 
#'
#' @param yr2 the yogiroc2 object
#' @param col vector of colors to use for the predictors
#' @param monotonized whether or not to monotonized the curve
#' @param balanced whether or not to use prior-balancing
#' @param legend the position of the legend, e.g. "bottomleft". NA disables legend
#' @param ... additional graphical parameters (see \code{par})
#'
#' @return nothing. draws a plot
#' @export
#'
#' @examples
#' #generate fake data
#' truth <- c(rep(TRUE,10),rep(FALSE,8))
#' scores <- cbind(
#'   pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
#'   pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #draw PRC curve
#' draw.prc(yrobj)
#' #draw non-monotonized PRC curve
#' draw.prc(yrobj,monotonized=FALSE)
#' #draw balanced PRC curve
#' draw.prc(yrobj,balanced=TRUE)
draw.prc <- function(yr2,col=seq_along(yr2),lty=1,monotonized=TRUE,balanced=FALSE,legend="bottomleft",...) {
  stopifnot(inherits(yr2,"yr2"))
  if (length(lty) < length(yr2)) {
    lty <- rep(lty,length(yr2))
  }
  ppv <- function(i) {
    configure.prec(yr2[[i]],monotonized,balanced)
    # raw <- if (balanced) yr2[[i]][,"ppv.prec.balanced"] else yr2[[i]][,"ppv.prec"]
    # if (monotonized) monotonize(raw) else raw
  }
  plabel <- ifelse(balanced,"Balanced precision (%)","Precision (%)")
  plot(
    100*yr2[[1]][,"tpr.sens"],100*ppv(1),
    type="l",
    xlab="Recall (%)", ylab=plabel,
    xlim=c(0,100),ylim=c(0,100),col=col[[1]], lty=lty[[1]], ...
  )
  if(length(yr2) > 1) {
    for (i in 2:length(yr2)) {
      lines(
        100*yr2[[i]][,"tpr.sens"],100*ppv(i),
        col=col[[i]], lty=lty[[i]], ...
      )
    }
  } 
  if (!is.na(legend)) {
    legend(legend,sprintf("%s (AUBPRC=%.02f;R90BP=%.02f)",
                          names(yr2),auprc(yr2,monotonized,balanced),recall.at.prec(yr2,0.9,monotonized,balanced)
    ),col=col,lty=lty)
  }
}


#' Draw Precision-Recall Curve (PRC) with confidence intervals
#' 
#' Balancing concept by Yingzhou Wu and Fritz Roth (Wu et al, unpublished) 
#' Confidence interval concept by Jochen Weile
#'
#' @param yr2 the yogiroc2 object
#' @param col vector of colors to use for the predictors
#' @param monotonized whether or not to monotonize the curve
#' @param balanced whether or not to use prior-balancing
#' @param legend the position of the legend, e.g. "bottomleft". NA disables legend
#' @param ... additional graphical parameters (see \code{par})
#'
#' @return nothing. draws a plot
#' @export
#'
#' @examples
#' #generate fake data
#' N <- 100
#' M <- 80
#' truth <- c(rep(TRUE,N),rep(FALSE,M))
#' scores <- cbind(
#'   pred1=c(rnorm(N,1,0.2),rnorm(M,.9,0.1)),
#'   pred2=c(rnorm(N,1.1,0.2),rnorm(M,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #draw PRC curve
#' draw.prc.CI(yrobj)
#' #draw non-monotonized PRC curve
#' draw.prc.CI(yrobj,monotonized=FALSE)
draw.prc.CI <- function(yr2,col=seq_along(yr2),lty=1,
                        monotonized=TRUE,balanced=FALSE,legend="bottomleft",
                        sampling=c("accurate","quickDirty"),nsamples=1000L,monotonizedSampling=FALSE,
                        ...) {
  
  stopifnot(inherits(yr2,"yr2"))
  sr <- switch(match.arg(sampling,c("accurate","quickDirty")),
               quickDirty=sampleRatesQD,accurate=sampleRates
  )
  if (length(lty) < length(yr2)) {
    lty <- rep(lty,length(yr2))
  }
  # mon <- function(xs) if (monotonized) monotonize(xs) else xs
  ppv <- function(i) configure.prec(yr2[[i]],monotonized,balanced)
  plabel <- ifelse(balanced,"Balanced precision (%)","Precision (%)")
  plot(
    100*yr2[[1]][,"tpr.sens"],100*ppv(1),
    type="l",
    xlab="Recall (%)", ylab=plabel,
    xlim=c(0,100),ylim=c(0,100),col=col[[1]], lty=lty[[1]], ...
  )
  if(length(yr2) > 1) {
    for (i in 2:length(yr2)) {
      lines(
        100*yr2[[i]][,"tpr.sens"],100*ppv(i),
        col=col[[i]], lty=lty[[i]], ...
      )
    }
  } 
  for (i in 1:length(yr2)) {
    # x <- 100*yr2[[i]][,"tpr.sens"]
    prior <- yr2[[i]][1,"tp"]/(yr2[[i]][1,"tp"]+yr2[[i]][1,"fp"])
    # precCI <- prcCI(yr2[[i]][,"tp"],yr2[[i]][,"tp"]+yr2[[i]][,"fp"])
    precCI <- inferPRCCI(samplePRCs(yr2[[i]],N=nsamples,monotonized=monotonizedSampling,sr=sr))
    precCI[,-1] <- apply(precCI[,-1],2,function(column) {
      if (balanced) {
        column <- balance.prec(column,prior)
      } 
      # if (monotonized) {
      #   column <- monotonize(column)
      # }
      column
    })
    polygon(100*c(precCI[,1],rev(precCI[,1])),
            c(100*precCI[,"0.025"],rev(100*precCI[,"0.975"])),
            col=yogitools::colAlpha(col[[i]],0.1),border=NA
    )
  }
  if (!is.na(legend)) {
    legend(legend,sprintf("%s (AUBPRC=%.02f;R90BP=%.02f)",
                          names(yr2),auprc(yr2,monotonized,balanced),recall.at.prec(yr2,0.9,monotonized,balanced)
    ),col=col,lty=lty)
  }
}

auprc.signif2 <- function(yr2,monotonized=TRUE) {
  aucDistrs <- lapply(1:length(yr2),function(i) {
    paths <- samplePRCs(yr2[[i]],monotonized=monotonized,sr=sampleRates)
    paths <- lapply(1:nrow(paths[[1]]),function(j) {
      do.call(rbind,lapply(1:length(paths), function(row) {
        paths[[row]][j,]
      }))
    })
    sapply(paths,function(path) {
      path <- path[order(path[,2]),]
      calc.auc(path[,2],path[,1])
    })
  })
}

#' Assess the significance of AUPRC differences
#' 
#' The list returned by this functions contains four elements:
#' \describe{
#' \item{auprc}{is simply the empirical area under the precision recall curve 
#' for each predictor.}
#' \item{ci}{is a matrix listing the lower and upper end of the 95% confidence 
#' interval for the AUPRC of each predictor.}
#' \item{llr}{is a matrix with columns and rows corresponding to each predictor.
#' It lists the log likelihood ratio of how much more (or less) likely the row-wise
#' predictor is to have a greater AUPRC than the column-wise predictor.}
#' \item{pval}{is a matrix with columns and rows corresponding to each predictor.
#' It lists the p-value of how likely it would be to observe the AUPRC of the row-wise
#' predictor under the distribution of the column-wise predictor.}
#' }
#'
#' @param yr2 the yogiroc2 object
#' @param monotonized whether or not to monotonize the curve
#' @param res the resolution at which to sample the probability function 
#' (defaults to 0.001)
#'
#' @return a list containing 4 elements: "auprc" (the empirical area under the 
#' precision recall curve), "ci" (the 95%confidence interval around the auprc), 
#' "llr" (the log likelihood ratio matrix, see details), and "pval" (the p-value 
#' of each auprc against each other)
#' @export
#'
#' @examples
#' #generate fake data
#' N <- 100
#' M <- 80
#' truth <- c(rep(TRUE,N),rep(FALSE,M))
#' scores <- cbind(
#'   pred1=c(rnorm(N,1,0.2),rnorm(M,.9,0.1)),
#'   pred2=c(rnorm(N,1.1,0.2),rnorm(M,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' auprc.signif(yrobj)
auprc.signif <- function(yr2,monotonized=TRUE,res=0.001) {
  #probability range
  ps <- seq(res,1-res,res)
  #calculate the AUPRC for each probability
  #i.e. the quantiles corresponding to ps
  auprcs <- do.call(cbind,lapply(1:length(yr2),function(i) {
    precCI <- prcCI(yr2[[i]][,"tp"],yr2[[i]][,"tp"]+yr2[[i]][,"fp"],p=ps)
    apply(precCI,2,function(ppv) {
      if (monotonized) {
        ppv <- monotonize(ppv)
      }
      calc.auc(yr2[[i]][,"tpr.sens"],ppv)
    })
  }))
  
  #empirical AUCs of the predictors
  empAUCs <- auprc(yr2,monotonized=monotonized)
  
  #1. build a reverse-lookup table that returns p for a given auprc
  aucRange <- do.call(seq,as.list(c(round(range(auprcs),digits=2),res)))
  aucPs <- do.call(rbind,lapply(aucRange,function(a){
    apply(auprcs,2,function(ladder){
      c(0,ps)[[sum(ladder < a)+1]]
    })
  }))
  
  confInts <- auprcs[c("0.025","0.975"),]
  colnames(confInts) <- names(yr2)
  
  pvals <- do.call(rbind,lapply(1:ncol(aucPs),function(i) {
    sapply(1:ncol(aucPs),function(j) {
      if (i==j) NA else {
        1-c(0,ps)[[sum(auprcs[,j] < empAUCs[[i]])+1]]
      }
    })
  }))
  dimnames(pvals) <- list(names(yr2), names(yr2))
  
  llrs <- do.call(rbind,lapply(1:ncol(aucPs),function(i) {
    sapply(1:ncol(aucPs),function(j) {
      if (i==j) NA else {
        #2. iterate over range of auprcs and calculate p_A(x) * (1-p_B(x)) 
        #  (i.e. the probability that area A is smaller than x AND area B is greater than x)
        log10(calc.auc(aucRange,aucPs[,j]*(1-aucPs[,i]))/calc.auc(aucRange,aucPs[,i]*(1-aucPs[,j])))
      }
    })
  }))
  dimnames(llrs) <- list(names(yr2), names(yr2))
  
  # plot(NA,type="n",xlim=c(0,1),ylim=c(0,1),xlab="AUPRC",ylab="CDF")
  # for (i in 1:ncol(auprcs)) {
  #   lines(c(0,auprcs[,i],1),c(0,ps,1),col=i)
  # }
  # abline(v=empAUCs,col=1:length(yr2),lty="dashed")
  
  return(list(auprc=empAUCs,ci=confInts,llr=llrs,pval=pvals))
}


#' Assess the significance of AUPRC against random guessing
#' #'
#' @param yr2 the yogiroc2 object
#' @param monotonized whether or not to monotonize the curves
#' @param cycles the sample size of the null distribution to use
#'
#' @return The emprical p-values of the AUPRC against random guessing
#' @export
#'
#' @examples
#' #generate fake data
#' N <- 10
#' M <- 8
#' truth <- c(rep(TRUE,N),rep(FALSE,M))
#' scores <- cbind(
#'   pred1=c(rnorm(N,1,0.2),rnorm(M,.9,0.1)),
#'   pred2=c(rnorm(N,1.1,0.2),rnorm(M,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #call pvrandom function
#' auprc.pvrandom(yrobj)
auprc.pvrandom <- function(yr2,monotonized=TRUE,cycles=10000) {
  
  empAUCs <- auprc(yr2,monotonized=monotonized)
  
  #reconstruct truth table
  real <- yr2[[1]][1,"tp"]
  nreal <- yr2[[1]][1,"fp"]
  truth <- c(rep(TRUE,real),rep(FALSE,nreal))
  
  nullAucs <- replicate(cycles,{
    scores <- runif(real+nreal,0,1)
    ts <- na.omit(c(-Inf,sort(scores),Inf))
    pr <- t(sapply(sort(scores), function(t) {
      calls <- scores >= t
      tp <- sum(calls & truth,na.rm=TRUE)
      # tn <- sum(!calls & !truth,na.rm=TRUE)
      fp <- sum(calls & !truth,na.rm=TRUE)
      fn <- sum(!calls & truth,na.rm=TRUE)
      prec <- tp/(tp+fp)
      recall <- tp/(tp+fn)
      # fpr.fall <- fp/(tn+fp)
      c(prec,recall)
    }))
    if (monotonized){
      pr[,1] <- monotonize(pr[,1])
    }
    calc.auc(pr[,2],pr[,1])
  })
  
  pvals <- sapply(empAUCs, function(eauc) {
    sum(nullAucs >= eauc)/cycles
  })
  
  return(pvals)
}
# auprc.CI <- function(yr2,monotonized=TRUE) {
#   do.call(rbind,lapply(1:length(yr2),function(i) {
#     precCI <- prcCI(yr2[[i]][,"tp"],yr2[[i]][,"tp"]+yr2[[i]][,"fp"])
#     auprcs <- apply(precCI,2,function(ppv) {
#       if (monotonized) {
#         ppv <- monotonize(ppv)
#       }
#       calc.auc(yr2[[i]][,"tpr.sens"],ppv)
#     })
#   }))
# }

#' Helper function to calculate area under curve
#'
#' @param xs the x values of the graph
#' @param ys the corresponding y values of the graph
#'
#' @return the area under the curve
calc.auc <- function(xs,ys) {
  #calculate the sum of the areas of individual x-segments
  sum(sapply(1:(length(xs)-1),function(i) {
    #calculate interval width between datapoints on x-axis
    delta.x <- abs(xs[[i]]-xs[[i+1]])
    #calculate the average height of the two points on the y-axis
    y <- (ys[[i]]+ys[[i+1]])/2
    #area = x * y ; geometrically, this works out to be the same as the area of the polygon
    delta.x * y
  }))
}

#' Calculate area under precision recall curve (AUPRC)
#' 
#' Balancing concept by Yingzhou Wu and Fritz Roth (Wu et al, unpublished) 
#'
#' @param yr2 the yogiroc2 object
#' @param monotonized whether or not use a monotonized PRC curve
#' @param balanced whether or not to use prior-balancing
#'
#' @return a numerical vector with the AUPRC values for each predictor
#' @export
#'
#' @examples
#' #generate fake data
#' truth <- c(rep(TRUE,10),rep(FALSE,8))
#' scores <- cbind(
#'   pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
#'   pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #calculate AUPRC
#' auprc(yrobj)
#' #calculate non-monotonized AUPRC
#' auprc(yrobj,monotonized=FALSE)
#' #calculate balanced AUPRC
#' auprc(yrobj,balanced=TRUE)
auprc <- function(yr2, monotonized=TRUE, balanced=FALSE) {
  stopifnot(inherits(yr2,"yr2"))
  # ppv <- function(data) {
  #   raw <- if (balanced) data[,"ppv.prec.balanced"] else data[,"ppv.prec"]
  #   if (monotonized) monotonize(raw) else raw
  # }
  sapply(yr2,function(data) {
    calc.auc(data[,"tpr.sens"],configure.prec(data,monotonized,balanced))
  })
}


#' Calculate the area under the ROC curve
#'
#' @param yr2 the yogiroc2 object
#'
#' @return a numerical vector with the AUROC for each predictor
#' @export
#'
#' @examples
#' #generate fake data
#' truth <- c(rep(TRUE,10),rep(FALSE,8))
#' scores <- cbind(
#'   pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
#'   pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #calculate AUROC
#' auroc(yrobj)
auroc <- function(yr2) {
  stopifnot(inherits(yr2,"yr2"))
  sapply(yr2,function(data) {
    calc.auc(data[,"fpr.fall"],data[,"tpr.sens"])
  })
}

#' Calculate maximum recall at given minimum precision
#'
#' @param yr2 the yogiroc2 object
#' @param x the precision cutoff (default 0.9)
#' @param monotonized whether or not to use monotonized PRC
#' @param balanced whether or not to use prior-balancing
#'
#' @export
#'
#' @examples
#' #generate fake data
#' truth <- c(rep(TRUE,10),rep(FALSE,8))
#' scores <- cbind(
#'   pred1=c(rnorm(10,1,0.2),rnorm(8,.9,0.1)),
#'   pred2=c(rnorm(10,1.1,0.2),rnorm(8,.9,0.2))
#' )
#' #create yogiroc2 object
#' yrobj <- yr2(truth,scores)
#' #calculate R90P
#' recall.at.prec(yrobj)
#' #calculate non-monotonized R90P
#' recall.at.prec(yrobj,monotonized=FALSE)
#' #calculate balanced R90P
#' recall.at.prec(yrobj,balanced=TRUE)
recall.at.prec <- function(yr2,x=0.9,monotonized=TRUE,balanced=FALSE) {
  stopifnot(inherits(yr2,"yr2"))
  # ppv <- function(data) {
  #   raw <- if (balanced) data[,"ppv.prec.balanced"] else data[,"ppv.prec"]
  #   if (monotonized) monotonize(raw) else raw
  # }
  sapply(yr2,function(data) {
    ppv <- configure.prec(data,monotonized,balanced)
    if (any(ppv > x)) {
      max(data[which(ppv > x),"tpr.sens"])
    } else NA
  })
}

#' Extract thresholds
#'
#' @param yr2 the yogiroc2 object
#' @param x the precision cutoff (default 0.9)
#' @param monotonized whether or not to use monotonized PRC
#' @param balanced whether or not to use prior-balancing
#'
#' @return threshold ranges
#' @export
calculate_thresh_range <- function(yr2, x=0.9, monotonized=TRUE, balanced=FALSE) {
  stopifnot(inherits(yr2, "yr2"))
  
  # list to store ranges
  thresh_ranges <- vector("list", length(yr2))
  names(thresh_ranges) <- names(yr2)
  
  for (i in seq_along(yr2)) {
    data <- yr2[[i]]
    ppv <- configure.prec(data, monotonized=monotonized, balanced=balanced)
    
    # which thresholds meet cutoff
    hits <- which(ppv > x)
    
    if (length(hits) > 0) {
      max_thresh <- data[hits[1], "thresh"] # right after reaching precision cutoff
      min_thresh <- if (hits[1] > 1) {
        data[hits[1] - 1, "thresh"] # default - take threshold right before cutoff
      } else {
        -Inf
      }
      
      # Combine into a single string "min-max"
      thresh_ranges[[i]] <- sprintf("%.3f-%.3f", min_thresh, max_thresh)
      
    } else {
      thresh_ranges[[i]] <- NA_character_
    }
  }
  return(thresh_ranges)
}


calculate_thresh_range_flipper <- function(yr2, x = 0.9, monotonized = TRUE, balanced = FALSE, high = rep(TRUE, length(yr2))) {
  stopifnot(inherits(yr2, "yr2"))
  stopifnot(length(high) == length(yr2))  # ensure one 'high' per predictor
  
  thresh_ranges <- vector("list", length(yr2))
  names(thresh_ranges) <- names(yr2)
  
  for (i in seq_along(yr2)) {
    data <- yr2[[i]]
    ppv <- configure.prec(data, monotonized = monotonized, balanced = balanced)
    
    hits <- which(ppv > x)
    
    if (length(hits) > 0) {
      if (high[i]) {
        max_thresh <- data[hits[1], "thresh"]
        min_thresh <- if (hits[1] > 1) data[hits[1] - 1, "thresh"] else -Inf
      } else {
        min_thresh <- data[hits[1], "thresh"]
        max_thresh <- if (hits[1] > 1) data[hits[1] - 1, "thresh"] else Inf
      }
      
      thresh_ranges[[i]] <- sprintf("%.3f-%.3f", min_thresh, max_thresh)
    } else {
      thresh_ranges[[i]] <- NA_character_
    }
  }
  return(thresh_ranges)
}




# Copyright (C) 2020  Jochen Weile, Roth Lab
# Edited by Cindy Zhang, Roth Lab (2024)
#
# This file is part of maveLLR
#
# maveLLR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# maveLLR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with maveLLR  If not, see <https://www.gnu.org/licenses/>.

#' Build LLR function using kernel density estimation
#'
#' @param posScores fitness scores in the positive reference set
#' @param negScores fitness scores in the negative reference set
#' @param bw the bandwith to use for kernel density estimation (see kdensity package for more info).
#'   This can either be a numerical value or the name of the algorithm used to automatically choose one.
#' @param kernel the type of kernel to use (see kdensity package for more info)
#'
#' @return a list containing the LLR function in log10 scale and the positive and negative density functions
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.kernel(posScores,negScores,bw=0.1,kernel="gaussian")
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
buildLLR.kernel <- function(posScores, negScores, bw=0.1, kernel="gaussian", outlierSuppression=0.0001) {
  
  library(kdensity)
  
  posDens <- kdensity(posScores,bw=bw,kernel=kernel)
  negDens <- kdensity(negScores,bw=bw,kernel=kernel)
  
  #generate a dummy outlier point far away from the rest of the distribution
  # so we can use it to measure it's hypothetical density
  pseudo <- min(posScores) - 10*ifelse(is.numeric(bw),bw,0.1)
  #measure the hypothetical density of an outlier point
  minDensPos <- kdensity(c(posScores[-1],pseudo),bw=bw,kernel=kernel)(pseudo)
  #do the same for the negative distribution
  pseudo <- min(negScores) - 10*ifelse(is.numeric(bw),bw,0.1)
  minDensNeg <- kdensity(c(negScores[-1],pseudo),bw=bw,kernel=kernel)(pseudo)
  #whichever outlier density is higher will serve as our uniform prior
  minDens <- outlierSuppression * max(minDensPos,minDensNeg)
  
  llrFun <- function(score) sapply(score,function(s)
    # log10( max(posDens(s),minDens) / max(negDens(s),minDens) )
    log10(max(posDens(s),minDens)) - log10(max(negDens(s),minDens))
  )
  
  return(list(llr=llrFun,posDens=posDens,negDens=negDens))
}

#' Build LLR function using kernel density estimation - experimental version
#'
#' @param posScores fitness scores in the positive reference set
#' @param negScores fitness scores in the negative reference set
#' @param bw the bandwith to use for kernel density estimation (see kdensity package for more info).
#'   This can either be a numerical value or the name of the algorithm used to automatically choose one.
#' @param kernel the type of kernel to use (see kdensity package for more info)
#' 
#'
#' @return a list containing the LLR function in log10 scale and the positive and negative density functions
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.kernel(posScores,negScores,bw=0.1,kernel="gaussian")
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
buildLLR.kernelExperimental <- function(posScores, negScores, bw=0.1, kernel="gaussian", outlierSuppression=0.0001) {
  
  library(kdensity)
  
  posDens <- kdensity(posScores,bw=bw,kernel=kernel)
  negDens <- kdensity(negScores,bw=bw,kernel=kernel)
  refDens <- kdensity(c(posScores,negScores), bw=bw, kernel=kernel)
  
  # generate a dummy outlier point far away from the rest of the distribution
  # so we can use it to measure it's hypothetical density
  pseudo <- min(c(posScores, negScores)) - 10*ifelse(is.numeric(bw),bw,0.1)
  # measure the hypothetical density of an outlier point
  pseudoDens <- kdensity(c(posScores[-1],negScores,pseudo),bw=bw,kernel=kernel)(pseudo)
  priorWeight = pseudoDens * outlierSuppression
  
  llrFun <- function(score) sapply(score,function(s) {
    rawLLR = log10(posDens(s)) - log10(negDens(s))
    weight = refDens(s) / (refDens(s)+priorWeight)
    finalLLR = rawLLR * weight
  })
  
  return(list(llr=llrFun,posDens=posDens,negDens=negDens))
}


#' Build LLR function using Gaussian densities.
#'
#' @param posScores fitness scores in the positive reference set
#' @param negScores fitness scores in the negative reference set
#' @param spline logical; whether or not to apply spline monotonization. TRUE by default.
#'
#' @return a list containing the LLR function in log10 scale and the positive and negative density functions
#' @export
#'
#' @examples
#'
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.gauss(posScores,negScores)
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
buildLLR.gauss <- function(posScores, negScores, spline=TRUE) {
  
  mpos <- mean(posScores)
  spos <- sd(posScores)
  mneg <- mean(negScores)
  sneg <- sd(negScores)
  
  posDens <- function(x) dnorm(x,mpos,spos)
  negDens <- function(x) dnorm(x,mneg,sneg)
  
  llrFun <- function(score) log10(posDens(score)/negDens(score))
  
  if (spline) {
    minPoint <- optimize(llrFun,interval=c(mpos, qnorm(0.999,mneg,sneg)),maximum=FALSE)
    if (minPoint$minimum < mneg) {
      minPoint$minimum <- Inf
    }
    maxPoint <- optimize(llrFun,interval=c(qnorm(0.001,mpos,spos),mneg),maximum=TRUE)
    if (maxPoint$maximum > mpos) {
      maxPoint$maximum <- -Inf
    }
    
    llrSpline <- function(scores) {
      sapply(scores,function(score) {
        if (is.na(score)) {
          NA
        } else if (score > minPoint$minimum) {
          minPoint$objective
        } else if (score < maxPoint$maximum) {
          maxPoint$objective
        } else {
          llrFun(score)
        }
      })
    }
    
    return(list(llr=llrSpline,posDens=posDens,negDens=negDens))
    
  } else {
    return(list(llr=llrFun,posDens=posDens,negDens=negDens))
  }
  
}


#' Draw a plot that summarizes the densities and LLR calculated by other functions
#'
#' @param scores the numerical scores in the Mave map
#' @param llrFun the LLR function
#' @param posDens the density function of the positive reference set
#' @param negDens the density function of the negative reference set
#' @param posScores the numerical scores in the postive reference set
#' @param negScores the numerical scores in the negative reference set
#'
#' @return nothing
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.gauss(posScores,negScores)
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
drawDensityLLR <- function(scores, llrFun, posDens, negDens, posScores, negScores, prior=0.1) {
  
  llrTs <- llrThresholds(optiLLR(prior))
  
  opar <- par(mfrow=c(2,1))
  xlim <- range(scores,na.rm=TRUE,finite=TRUE)
  ymax <- max(c(
    posDens(seq(xlim[[1]],xlim[[2]],length.out = 100)),
    negDens(seq(xlim[[1]],xlim[[2]],length.out = 100))
  ))
  par(mar=c(.1,4,1,1))
  xs <- seq(xlim[[1]],xlim[[2]],length.out=200)
  ys <- llrFun(xs)
  ylim <- range(ys,na.rm=TRUE,finite=TRUE)
  plot(NA,type="n",xlim=xlim,ylim=ylim,
       axes=FALSE,xlab="",ylab="LLR"
  )
  drawThresh <- function(t,col,label) {
    if (t >= 0) {
      if (t < ylim[[2]]) {
        rect(xlim[[1]],t,xlim[[2]],ylim[[2]],col=col,border=NA)
        text(.1*xlim[[2]]+.9*xlim[[1]],t,label,pos=3,cex=0.8)
      }
    } else {
      if (t > ylim[[1]]) {
        rect(xlim[[1]],ylim[[1]],xlim[[2]],t,col=col,border=NA)
        text(.1*xlim[[2]]+.9*xlim[[1]],t,label,pos=1,cex=0.8)
      }
    }
  }
  drawThresh(llrTs[["patho.support"]],"lemonchiffon","Patho. support.")
  drawThresh(llrTs[["patho.moderate"]],"lightgoldenrod1","Patho. moderate")
  drawThresh(llrTs[["patho.strong"]],"goldenrod1","Patho. strong")
  drawThresh(llrTs[["patho.vstrong"]],"indianred1","Patho. very str.")
  drawThresh(llrTs[["benign.support"]],"lemonchiffon","Benign support.")
  drawThresh(llrTs[["benign.strong"]],"goldenrod1","Benign strong")
  
  lines(xs,ys)
  # plot(llrFun,from=xlim[[1]],to=xlim[[2]],xlim=xlim,axes=FALSE,xlab="",ylab="LLR")
  abline(h=0,col="gray",lty="dashed")
  axis(2)
  par(mar=c(5,4,.1,1))
  hist(scores,col="gray90",border=NA,freq=FALSE,main="",xlim=xlim,ylim=c(0,ymax),breaks=50)
  plot.function(posDens,from=xlim[[1]],to=xlim[[2]],add=TRUE,col="firebrick3",lwd=2)
  plot.function(negDens,from=xlim[[1]],to=xlim[[2]],add=TRUE,col="darkolivegreen3",lwd=2)
  abline(v=posScores,col="firebrick3")
  abline(v=negScores,col="darkolivegreen3")
  par(opar)
  
  return(invisible(NULL))
}

#' @export
optiLLR <- function(prior=0.1,posterior=0.9) {
  log10(posterior*(1-prior)/(prior*(1-posterior)))*4/3
}

#' @export
llrThresholds <- function(LLRpvst=optiLLR(0.1),X=2) {
  c(
    patho.vstrong=LLRpvst,
    patho.strong=LLRpvst/X,
    patho.moderate=LLRpvst/(X^2),
    patho.support=LLRpvst/(X^3),
    benign.support=-LLRpvst/(X^3),
    benign.strong=-LLRpvst/(X^1)
  )
}

#' Like drawDensityLLR, but with fixed y range for better comparison across different LLRs
#'
#' @param scores the numerical scores in the Mave map
#' @param llrFun the LLR function
#' @param posDens the density function of the positive reference set
#' @param negDens the density function of the negative reference set
#' @param posScores the numerical scores in the postive reference set
#' @param negScores the numerical scores in the negative reference set
#'
#' @return nothing
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.gauss(posScores,negScores)
#' drawDensityLLR_fixedRange(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
drawDensityLLR_fixedRange <- function(scores, llrFun, posDens, negDens, posScores, negScores, prior=0.1) {
  
  llrTs <- llrThresholds(optiLLR(prior))
  
  opar <- par(mfrow=c(2,1))
  xlim <- range(scores,na.rm=TRUE,finite=TRUE)
  ymax <- max(c(
    posDens(seq(xlim[[1]],xlim[[2]],length.out = 100)),
    negDens(seq(xlim[[1]],xlim[[2]],length.out = 100))
  ))
  par(mar=c(.1,4,1,1))
  xs <- seq(xlim[[1]],xlim[[2]],length.out=200)
  ys <- llrFun(xs)
  #ylim <- range(ys,na.rm=TRUE,finite=TRUE)
  ylim <- c(-3,4) # fixed range
  plot(NA,type="n",xlim=xlim,ylim=ylim,
       axes=FALSE,xlab="",ylab="LLR"
  )
  drawThresh <- function(t,col,label) {
    if (t >= 0) {
      if (t < ylim[[2]]) {
        rect(xlim[[1]],t,xlim[[2]],ylim[[2]],col=col,border=NA)
        text(.1*xlim[[2]]+.9*xlim[[1]],t,label,pos=3,cex=0.8)
      }
    } else {
      if (t > ylim[[1]]) {
        rect(xlim[[1]],ylim[[1]],xlim[[2]],t,col=col,border=NA)
        text(.1*xlim[[2]]+.9*xlim[[1]],t,label,pos=1,cex=0.8)
      }
    }
  }
  drawThresh(llrTs[["patho.support"]],"#FFC0CB","Patho. support.")
  drawThresh(llrTs[["patho.moderate"]],"#FF9999","Patho. moderate")
  drawThresh(llrTs[["patho.strong"]],"#FF6666","Patho. strong")
  drawThresh(llrTs[["patho.vstrong"]],"red","Patho. very str.")
  drawThresh(llrTs[["benign.support"]],"lightblue","Benign support.")
  drawThresh(llrTs[["benign.strong"]],"dodgerblue","Benign strong")
  
  lines(xs,ys)
  # plot(llrFun,from=xlim[[1]],to=xlim[[2]],xlim=xlim,axes=FALSE,xlab="",ylab="LLR")
  abline(h=0,col="gray",lty="dashed")
  axis(2)
  par(mar=c(5,4,.1,1))
  hist(scores,col="gray90",border=NA,freq=FALSE,main="",xlim=xlim,ylim=c(0,ymax),breaks=50)
  plot.function(posDens,from=xlim[[1]],to=xlim[[2]],add=TRUE,col="firebrick3",lwd=2)
  plot.function(negDens,from=xlim[[1]],to=xlim[[2]],add=TRUE,col="deepskyblue2",lwd=2)
  abline(v=posScores,col="firebrick3")
  abline(v=negScores,col="deepskyblue2")
  par(opar)
  
  return(invisible(NULL))
}

#' Find x-values at which the given LLR function crosses specified thresholds.
#'
#' @param llrFun the LLR Function
#' @param thresholds thresholds for each ACMG category
#' @param xlim range of x values i.e. scores
#' @param nPoints an integer specifying number of points to sample across xlim
#'
#' @return a table with threshold labels, values, and x when crossing
#' @export
#'
#' @examples
#' # Using a previously defined LLR function
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llrList <- buildLLR.kernel(posScores, negScores, bw=0.1, kernel="gaussian")
#' llrFun <- llrList$llr
#'
#' # Define thresholds of interest
#' llrTs <- llrThresholds(optiLLR(0.1))
#'
#' # Find the crossings in the range of observed scores
#' x_range <- range(c(posScores, negScores))
#' findLLRcrossings(llrFun, llrTs, xlim = x_range)
findLLRcrossings <- function(llrFun, thresholds, xlim, nPoints = 1000) {
  xs <- seq(xlim[1], xlim[2], length.out = nPoints)
  ys <- llrFun(xs)
  
  # Get indices
  patho_idx <- which(names(thresholds) == "patho.support")
  benign_support_idx <- which(names(thresholds) == "benign.support")
  benign_strong_val <- thresholds[which(names(thresholds) == "benign.strong")]
  
  # Insert 'none' as category between patho.support and benign.support
  benign_support_val <- thresholds[benign_support_idx]
  thresholds <- append(
    thresholds,
    setNames(benign_support_val, "none"),
    after = patho_idx
  )
  
  # Assign the lower value of none as the upper value of benign.support, but without changing display thresholds
  display_thresholds <- thresholds
  thresholds[which(names(thresholds) == "benign.support")] <- benign_strong_val
  display_thresholds["benign.support"] <- thresholds["none"]
  
  thresholdNames <- names(thresholds)
  
  # Interval range formatting for thresholds
  format_range <- function(name, i) {
    if (name == "none") return("[-0.32, 0.32)")
    val <- display_thresholds[i]
    if (is.na(val)) return("")
    
    if (grepl("^patho\\.", name)) {
      upper <- if (i > 1 && !is.na(display_thresholds[i - 1])) display_thresholds[i - 1] else Inf
      return(sprintf("[%.2f, %s)", val, ifelse(is.infinite(upper), "∞", sprintf("%.2f", upper))))
    } else if (grepl("^benign\\.", name)) {
      lower <- if (i < length(display_thresholds) && !is.na(display_thresholds[i + 1])) display_thresholds[i + 1] else -Inf
      bracket <- if (is.infinite(lower)) "(" else "["
      return(sprintf("%s%s, %.2f)", bracket, ifelse(is.infinite(lower), "-∞", sprintf("%.2f", lower)), val))
    }
    return("")
  }
  
  formatted_thresholds <- setNames(sapply(seq_along(thresholds), function(i) {
    format_range(names(thresholds)[i], i)
  }), names(thresholds))
  
  # Determine crossings and direction of crossing
  crossings <- list()
  for (i in seq_along(thresholds)) {
    thr <- thresholds[i]
    
    fvals <- ys - thr
    sign_changes <- which(fvals[-length(fvals)] * fvals[-1] < 0)
    
    for (j in sign_changes) {
      interval <- c(xs[j], xs[j + 1])
      if (!is.finite(fvals[j]) || !is.finite(fvals[j + 1])) next
      
      root <- tryCatch({
        uniroot(function(z) llrFun(z) - thr, interval = interval)$root
      }, error = function(e) NA)
      
      if (!is.na(root)) {
        slope <- ys[j + 1] - ys[j]
        cur_idx <- i
        
        if (slope > 0 && cur_idx + 1 <= length(thresholds)) {
          to <- thresholdNames[cur_idx]
          from <- thresholdNames[cur_idx + 1]
        } else if (slope < 0 && cur_idx + 1 <= length(thresholds)) {
          from <- thresholdNames[cur_idx]
          to <- thresholdNames[cur_idx + 1]
        } else {
          next
        }
        
        crossings[[length(crossings) + 1]] <- list(
          from = from,
          to = to,
          score = round(root, 2)
        )
      }
    }
  }
  
  # Build crossing data frame
  crossing_df <- do.call(rbind, lapply(crossings, as.data.frame))
  crossing_df$score <- as.numeric(as.character(crossing_df$score))
  crossing_df$key <- paste0(pmin(crossing_df$from, crossing_df$to), "|", pmax(crossing_df$from, crossing_df$to))
  
  # Build crossing blocks by key (sorted by score)
  crossing_blocks_by_key <- list()
  if (nrow(crossing_df) > 0) {
    crossing_df <- crossing_df[order(crossing_df$score), ]
    for (i in seq_len(nrow(crossing_df))) {
      from <- crossing_df$from[i]
      to <- crossing_df$to[i]
      score <- crossing_df$score[i]
      key <- crossing_df$key[i]
      
      block <- list()
      block[[1]] <- data.frame(
        LLR_category = "---------------------",
        LLR_threshold = "---------------------",
        score_crossing = score
      )
      block[[2]] <- data.frame(
        LLR_category = to,
        LLR_threshold = formatted_thresholds[[to]],
        score_crossing = ""
      )
      
      if (!key %in% names(crossing_blocks_by_key)) crossing_blocks_by_key[[key]] <- list()
      crossing_blocks_by_key[[key]][[length(crossing_blocks_by_key[[key]]) + 1]] <- list(
        from = from,
        to = to,
        score = score,
        rows = block
      )
    }
  }
  
  # Output starts with the first category in the threshold sequence
  output <- list()
  # Determine first_from: check if a crossing exists for the first threshold pair
  first_pair_key <- paste0(pmin(thresholdNames[1], thresholdNames[2]), "|", pmax(thresholdNames[1], thresholdNames[2]))
  
  if (first_pair_key %in% names(crossing_blocks_by_key)) {
    # Use the first actual 'from' involved in that crossing
    first_from <- crossing_blocks_by_key[[first_pair_key]][[1]]$from
  } else {
    # Fallback to first threshold name
    first_from <- thresholdNames[1]
  }
  
  output[[length(output) + 1]] <- data.frame(
    LLR_category = first_from,
    LLR_threshold = formatted_thresholds[[first_from]],
    score_crossing = ""
  )
  
  seen <- first_from
  
  # Walk through threshold pairs in original order
  for (i in seq_len(length(thresholds) - 1)) {
    from <- thresholdNames[i]
    to <- thresholdNames[i + 1]
    key <- paste0(pmin(from, to), "|", pmax(from, to))
    
    if (key %in% names(crossing_blocks_by_key)) {
      blocks <- crossing_blocks_by_key[[key]]
      for (block in blocks) {
        output[[length(output) + 1]] <- block$rows[[1]]  # crossing line
        output[[length(output) + 1]] <- block$rows[[2]]  # to row
        seen <- c(seen, block$to)
      }
    } else {
      # Insert NA crossing row and 'to' category
      output[[length(output) + 1]] <- data.frame(
        LLR_category = "---------------------",
        LLR_threshold = "---------------------",
        score_crossing = NA
      )
      output[[length(output) + 1]] <- data.frame(
        LLR_category = to,
        LLR_threshold = formatted_thresholds[[to]],
        score_crossing = ""
      )
      seen <- c(seen, to)
    }
  }
  
  # Add any remaining categories not seen
  for (name in names(formatted_thresholds)) {
    if (!(name %in% seen)) {
      output[[length(output) + 1]] <- data.frame(
        LLR_category = name,
        LLR_threshold = formatted_thresholds[[name]],
        score_crossing = ""
      )
    }
  }
  
  
  
  
  
  # Point at where LLR = 0
  fvals <- ys
  zero_crossings <- which(fvals[-length(fvals)] * fvals[-1] < 0)
  for (j in zero_crossings) {
    interval <- c(xs[j], xs[j + 1])
    root <- tryCatch({
      uniroot(function(z) llrFun(z), interval = interval)$root
    }, error = function(e) NA)
    if (!is.na(root)) {
      output[[length(output) + 1]] <- data.frame(
        LLR_category = "Crossing 0",
        LLR_threshold = 0,
        score_crossing = round(root, 2)
      )
    }
  }
  
  do.call(rbind, output)
}





        


#
# mthfr <- read.csv("~/projects/mthfr/folate_response_model5.csv")
# variants <- mthfr[mthfr$type=="substitution","hgvs"]
# scores <- mthfr[mthfr$type=="substitution","m25.score"]
#
# idb <- 337
# prsnrs <- read.csv("~/projects/mthfr/mthfr_prsNrs.csv")
#
# posref <- prsnrs[prsnrs$reference=="positive" & prsnrs$start < idb,"hgvs"]
# negref <- prsnrs[prsnrs$reference=="negative" & prsnrs$start < idb,"hgvs"]
# posScores <- na.omit(scores[variants %in% posref])
# negScores <- na.omit(scores[variants %in% negref])
# llr <- buildLLR.kernel(posScores, negScores,bw=0.1,kernel="gaussian")
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
# llr <- buildLLR.gauss(posScores, negScores)
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#
# posref <- prsnrs[prsnrs$reference=="positive" & prsnrs$start > idb,"hgvs"]
# negref <- prsnrs[prsnrs$reference=="negative" & prsnrs$start > idb,"hgvs"]
# posScores <- na.omit(scores[variants %in% posref])
# negScores <- na.omit(scores[variants %in% negref])
# llr <- buildLLR.kernel(posScores, negScores)
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
# llr <- buildLLR.gauss(posScores, negScores)
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)




plan(multisession)

server <- function(input, output, session) {
  plot_data <- reactiveVal(NULL)
  selected_variants <- reactiveVal(NULL)  # Store user-selected variants
  threshold_data <- reactiveVal(NULL) # for storing PRC thresholds - REFACTOR: include as part of prc data?
  
  preprocessed_df <- read.table("preprocessed.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
  prcdata <- reactiveVal(preprocessed_df)
  
  rv <- reactiveValues(full_df = NULL, loadingFull = TRUE)
  
  # Read full.csv async
  future({
    df <- fread("full.csv")
    as.data.frame(df)
  }) %...>% (function(big_data) {
    rv$full_df <- big_data
    rv$loadingFull <- FALSE
  }) %...!% (function(err) {
    rv$loadingFull <- FALSE
    showModal(modalDialog(
      title = "Error loading large data",
      paste("Could not load full.csv:", conditionMessage(err))
    ))
  })
  
  # Giant if/else block to handle separate logic for Main App and Advanced
  observe({
    if (input$navbar == "Basic") {
      # logic for Main App
      # Update gene names based on what exists
      observe({
        df <- prcdata()
        if (!is.null(df)) {
          gene_names <- unique(df$base__gene)
          updateSelectizeInput(session, "Main_gene", choices = gene_names, selected = character(0), server = TRUE)
        }
      })
      
      # Update scores - show only if a predictor has at least 1 P/LP and 1 B/LB
      observe({
        req(input$Main_gene)
        df <- prcdata()
        gene_data <- df %>% filter(base__gene == input$Main_gene)
      
        available_scores <- c()
        
        # Check each score for the required condition
        if (any(!is.na(gene_data$VEP_varity_r))) {
          if (any(gene_data$clinvar == "P/LP" & !is.na(gene_data$VEP_varity_r)) &&
              any(gene_data$clinvar == "B/LB" & !is.na(gene_data$VEP_varity_r))) {
            available_scores <- c(available_scores, "VARITY")
          }
        }
        
        if (any(!is.na(gene_data$VEP_alphamissense__pathogenicity))) {
          if (any(gene_data$clinvar == "P/LP" & !is.na(gene_data$VEP_alphamissense__pathogenicity)) && 
              any(gene_data$clinvar == "B/LB" & !is.na(gene_data$VEP_alphamissense__pathogenicity))) {
            available_scores <- c(available_scores, "AlphaMissense")
          }
        }
        
        if (any(!is.na(gene_data$VEP_revel__score))) {
          if (any(gene_data$clinvar == "P/LP" & !is.na(gene_data$VEP_revel__score)) && 
              any(gene_data$clinvar == "B/LB" & !is.na(gene_data$VEP_revel__score))) {
            available_scores <- c(available_scores, "REVEL")
          }
        }
        
        # Update the checkbox group input with valid predictors
        updateCheckboxGroupInput(session, "Main_scores", choices = available_scores, selected = available_scores)
      })
      
      # For LLR VUS plots
      full_data <- reactive({
        req(input$Main_gene)
        req(rv$full_df)
        #df_full <- read.csv("full.csv", stringsAsFactors = FALSE)
        #df_full$clinvar[df_full$clinvar == "N/A"] <- NA # Treat as actual NA
        df_filtered <- rv$full_df %>% filter(base__gene == input$Main_gene & !is.na(clinvar))
        return(df_filtered)
      })
      
      # Variant selection pop-up
      showVariantSelectionModal <- function(df, gene_selected) {
        
        gene_display_df <- df %>% 
          dplyr::filter(base__gene == gene_selected) %>%
          dplyr::select(base__gene, base__achange, clinvar, revstat, stars) %>%
          dplyr::mutate(is_duplicate = base__achange %in% base__achange[duplicated(base__achange)]) 
        
        # Check for duplicates
        has_duplicates <- any(gene_display_df$is_duplicate)
        
        # Sort to place duplicates on top
        gene_display_df <- gene_display_df %>% dplyr::arrange(desc(is_duplicate))
        
        # Build message text if duplicates found
        duplicate_warning <- if (has_duplicates) {
          tags$div(
            style = "background-color: #FFF3CD; border: 1px solid #FFEEBA; padding: 10px; margin-bottom: 10px;",
            tags$strong("Warning: "),
            "Some variants (shown at top of this list) are listed multiple times in ClinVar under the same hgvs_pro with different classification. You may want to de-select them."
          )
        } else {
          NULL
        }
        
        # Show one single modal
        showModal(
          modalDialog(
            title = paste("Select Variants for Gene:", gene_selected),
            # Place the warning at the top if duplicates exist
            duplicate_warning,
            tags$p("Would you like to de-select some variants? 
             (Click on the variants you do not want to include in the PRC...)"),
            DTOutput("variant_table"),
            footer = tagList(
              modalButton("Cancel"),
              actionButton("confirm_selection", "OK")
            ),
            easyClose = FALSE
          )
        )
        
        # Render the interactive table
        output$variant_table <- renderDT({
          datatable(
            gene_display_df,
            selection = list(target = "row", selected = 1:nrow(gene_display_df)), 
            options = list(
              pageLength = 10,
              scrollX = TRUE,
              columnDefs = list(
                list(targets = 6, visible = FALSE)  # '6' = zero-based index for the 7th column, to hide the is_duplicate column
              )
            )
          )
        }, server = TRUE)
      }
      
      # Generate PRC with user selection
      observeEvent(input$Main_plotButton, {
        gene_s <- input$Main_gene
        prcfiltered <- prcdata() %>%
          filter(base__gene == gene_s)
        if (!is.null(df)) {
          showVariantSelectionModal(prcfiltered, gene_s)
        }
      })
      
      # After confirmation from the modal
      observeEvent(input$confirm_selection, {
        # (a) Check if the large data is done loading:
        if (rv$loadingFull || is.null(rv$full_df)) {
          # If it's not ready, show “Still loading” and do NOT proceed
          showModal(modalDialog(
            title = "Still Loading Large Data",
            "Please wait – the large dataset is still loading in the background.
             Try again once it's finished.",
            easyClose = TRUE
          ))
          return()
        }
        
        # Wrap everything in a withProgress block
        withProgress(message = "Generating plots...", value = 0, {
          
          # Keep the existing “Please wait” modal
          showModal(modalDialog(
            title = "Please wait",
            "Generating plots...",
            footer = NULL,
            easyClose = FALSE
          ))
          
          # Step 1: Select rows
          #incProgress(0.1, detail = "Selecting rows")
          selected_rows <- input$variant_table_rows_selected
          gene_s <- input$Main_gene
          prcfiltered <- prcdata() %>%
            filter(base__gene == gene_s)
          
          if (!is.null(prcfiltered) && length(selected_rows) > 0) {
            # Step 2: Prepare selected data
            #incProgress(0.3, detail = "Preparing selected data")
            selected_df <- prcfiltered[selected_rows, ]
            
            # Check that there is at least one P/LP and one B/LB
            if (!(any(selected_df$clinvar == "P/LP") && any(selected_df$clinvar == "B/LB"))) {
              removeModal() # remove "Please wait" modal if it's open
              showModal(modalDialog(
                title = "Error",
                "You must select at least one P/LP and one B/LB variant.",
                easyClose = TRUE,
                footer = modalButton("Close")
              ))
              return()
            }
            
            gene_s <- input$Main_gene
            exclude_common_variants <- input$Main_common_variant_filter
            selected_scores <- input$Main_scores
            
            # Renaming for consistency
            names(selected_df)[names(selected_df) == "VEP_varity_r"] <- "VARITY"
            names(selected_df)[names(selected_df) == "VEP_alphamissense__pathogenicity"] <- "AlphaMissense"
            names(selected_df)[names(selected_df) == "VEP_revel__score"] <- "REVEL"
            names(selected_df)[names(selected_df) == "gnomad__af"] <- "gnomAD_AF"
            
            # Common variant filter
            if (exclude_common_variants) {
              selected_df <- selected_df[is.na(selected_df$gnomAD_AF) | selected_df$gnomAD_AF <= 0.005, ]
            }
            selected_df <- selected_df[order(selected_df$clinvar), ]
            
            # Filter out any rows with an NA in any `selected_scores`
            selected_df <- selected_df[!(rowSums(is.na(selected_df[selected_scores])) > 0 | 
                                           rowSums(selected_df[selected_scores] == "NA", na.rm = TRUE) > 0), ]
            
            selected_variants(selected_df) 
            
            # Convert label to T/F
            selected_df <- selected_df %>% mutate(clinvar = ifelse(clinvar == "P/LP", TRUE, FALSE))
            
            # Count # of P and B
            P_org <- sum(selected_df$clinvar == TRUE)
            B_org <- sum(selected_df$clinvar == FALSE)
            
            tryCatch({
              #incProgress(0.5, detail = "Generating PRC plot")
              yrobj <- yr2(truth = selected_df[["clinvar"]], scores = selected_df[selected_scores], high = rep(TRUE, length(selected_scores)))
              
              # Added for threshold calculation
              thresh_ranges <- calculate_thresh_range_flipper(yrobj, x = 0.9, balanced = TRUE, high = rep(TRUE, length(selected_scores)))
              
              # Store thresholds for rendering
              threshold_data(thresh_ranges)
              
              # end threshold
              
              plot_data(list(
                yrobj = yrobj,
                lty_styles = c("dashed", "solid", "dashed")[1:length(selected_scores)],
                col_styles = c("purple", "cadetblue2", "orange")[1:length(selected_scores)],
                gene_s = gene_s,
                selected_scores = selected_scores,
                B_org = B_org,
                P_org = P_org,
                common_variant_filter = exclude_common_variants,
                prcfiltered = selected_df  # Save the filtered data for metadata
              ))
              output$Main_ErrorText <- renderText("")
              
            }, error = function(e) {
              plot_data(NULL)
              output$Main_ErrorText <- renderText("Error generating the plot.")
            })
            
            output$Main_PRCPlot <- renderPlot({
              plot_info <- plot_data()
              if (!is.null(plot_info)) {
                tryCatch({
                  draw.prc(plot_info$yrobj, lty = plot_info$lty_styles, col = plot_info$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", ")))
                  abline(h = 90, lty = "dashed")
                  legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
                }, error = function(e) {
                  showModal(modalDialog(
                    title = 'Error',
                    'Not enough data - must have at least one pathogenic and benign. Try selecting less VEPs or unselecting the common variant filter.
                  You can also add your own annotations (see Advanced Tab or download VUS).',
                    easyClose = TRUE,
                    footer = NULL
                  ))
                })
              }
            }, width = 600, height = 600, res = 72)
            
            # Render the threshold table
            output$Main_thresholdTableUI <- renderUI({
              req(threshold_data())
              
              tagList(
                h5("VEP threshold to achieve 90% Balanced Precision"),  # Add the title only if the table exists
                tableOutput("Main_thresholdTable")  
              )
            })
            
            output$Main_thresholdTable <- renderTable({
              req(threshold_data())
              
              df <- threshold_data()
              
              # Round numeric columns to 3 decimals if numeric
              df <- as.data.frame(lapply(df, function(x) {
                if (is.numeric(x)) {
                  format(round(x, 3), nsmall = 3)
                } else {
                  x
                }
              }))
              df
            }, rownames = FALSE)
            
            # LLR Plot
            llr_tabs <- list()
            llr_scores <- intersect(selected_scores, c("VARITY", "REVEL", "AlphaMissense"))
            
            #incProgress(0.5, detail = "Reading full LLR data (this might take a while)")
            
            # Retrieve filtered full data
            full_filtered <- full_data()

            # Renaming for consistency
            names(full_filtered)[names(full_filtered) == "VEP_varity_r"] <- "VARITY"
            names(full_filtered)[names(full_filtered) == "VEP_alphamissense__pathogenicity"] <- "AlphaMissense"
            names(full_filtered)[names(full_filtered) == "VEP_revel__score"] <- "REVEL"
            names(full_filtered)[names(full_filtered) == "gnomad__af"] <- "gnomAD_AF"
            
            for (score_idx in seq_along(llr_scores)) {
              score <- llr_scores[score_idx]
              incProgress(1/length(llr_scores), detail = paste("Generating LLR plot for:", score))
              
              posScores <- na.omit(setNames(
                selected_df[selected_df$clinvar == TRUE, score],
                selected_df[selected_df$clinvar == TRUE, "base__achange"]
              ))
              negScores <- na.omit(setNames(
                selected_df[selected_df$clinvar == FALSE, score],
                selected_df[selected_df$clinvar == FALSE, "base__achange"]
              ))
              
              if (length(posScores) > 0 & length(negScores) > 0) {
                llrObj <- buildLLR.kernel(posScores, negScores, outlierSuppression=0.001) # Change outlier suppression
                
                full_filtered_copy <- full_filtered
                full_filtered_copy$llr <- llrObj$llr(full_filtered_copy[[score]]) # full_filtered_copy has llr values
                
                # Define breakpoints and labels for the categories
                breaks <- c(-Inf, -1.27, -0.32, 0.32, 0.64, 1.27, 2.54, Inf)
                labels <- c(
                  "benign_strong",   # (-Inf, -1.27]
                  "benign_support",  # (-1.27, -0.32]
                  "none",            # (-0.32, 0.32]
                  "patho_support",   # (0.32, 0.64]
                  "patho_moderate",  # (0.64, 1.27]
                  "patho_strong",    # (1.27, 2.54]
                  "patho_vstrong"    # (2.54, Inf)
                )
                
                # Categorize llr based on the defined bins
                full_filtered_copy$category <- cut(full_filtered_copy$llr, breaks=breaks, labels=labels, include.lowest=TRUE)
                
                full_filtered_copy <- full_filtered_copy %>% filter(!is.na(category)) # filter NA - due to score and llr being NA
                
                # Define thresholds and search range for crossings
                llrTs <- llrThresholds(optiLLR(0.1))
                x_range <- range(c(posScores, negScores))
                
                # Find where the LLR function crosses the thresholds
                crossings_df <- findLLRcrossings(llrObj$llr, llrTs, x_range)
                
                # Create the stacked bar plot of clinvar vs category
                # Assign colors as requested:
                # dark blue for benign_strong, red for patho_vstrong, grey for none
                # other categories lighter shades in between
                category_colors <- c(
                  "benign_strong"   = "dodgerblue",
                  "benign_support"  = "lightblue",
                  "none"            = "grey",
                  "patho_support"   = "#FFC0CB",   # light pink
                  "patho_moderate"  = "#FF9999",   # medium pink
                  "patho_strong"    = "#FF6666",   # darker pink
                  "patho_vstrong"   = "red"
                )
                
                plot_stack_id <- paste0("Main_Stacked_", score)
                
                # Preserve current iteration's variables, or else it would be same plot
                local({
                  # Make local copies of the variables to avoid referencing loop variables
                  score_copy <- score
                  posScores_copy <- posScores
                  negScores_copy <- negScores
                  llrObj_copy <- llrObj
                  full_filtered_copy_local <- full_filtered_copy
                  crossings_df_copy <- crossings_df
                  
                  plot_id <- paste0("Main_LLRPlot_", score_copy)
                  table_id <- paste0("Main_LLRCrossingsTable_", score_copy)
                  
                  # -----------------------------
                  # 1) DEFINE DOWNLOAD BUTTON IDS
                  # -----------------------------
                  download_llr_png_id <- paste0("download_", score_copy, "_llr_png")
                  download_bar_svg_id <- paste0("download_", score_copy, "_bar_svg")
                  download_both_pdf_id <- paste0("download_", score_copy, "_both_pdf")
                  download_csv_id      <- paste0("download_", score_copy, "_csv")
                  
                  output[[plot_id]] <- renderPlot({
                    tryCatch({
                      drawDensityLLR_fixedRange(
                        full_filtered_copy_local[[score_copy]], 
                        llrObj_copy$llr, 
                        llrObj_copy$posDens, 
                        llrObj_copy$negDens, 
                        posScores_copy, 
                        negScores_copy
                      )
                    }, error = function(e) {
                      showModal(modalDialog(
                        title = 'Error',
                        'Not enough data',
                        easyClose = TRUE,
                        footer = NULL
                      ))
                    })
                  }, width = 500, height = 600, res = 72)
                  
                  # Render the stacked bar chart
                  output[[plot_stack_id]] <- renderPlot({
                    # Make sure we have factor ordering if needed
                    full_filtered_copy_local$clinvar <- factor(full_filtered_copy_local$clinvar, 
                                                               levels = c("P/LP", "B/LB", "VUS", "Conflicting"))
                    ggplot(full_filtered_copy_local, aes(x=clinvar, fill=category)) +
                      geom_bar(position=position_stack(reverse=TRUE)) + # reverse to put pathogenic at top like the LLR
                      scale_fill_manual(values=category_colors) +
                      theme_minimal() +
                      labs(x="ClinVar Category", y="Variant Count", fill="LLR-Derived Evidence Category") +
                      theme(axis.text.x = element_text(angle=45, hjust=1))
                  }, width = 300, height = 500, res = 72)
                  
                  # Render the table of LLR threshold crossings
                  output[[table_id]] <- renderTable({
                    crossings_df_copy
                  })
                  
                  # -----------------------------------
                  # 3) DEFINE DOWNLOAD HANDLERS (NEW)
                  # -----------------------------------
                  
                  # 3a) Download LLR (PNG)
                  output[[download_llr_png_id]] <- downloadHandler(
                    filename = function() { paste0(score_copy, "_LLRPlot.png") },
                    content = function(file) {
                      # Use a PNG device and re-draw the same LLR plot
                      png(file, width = 1500, height = 1800, res = 300)
                      drawDensityLLR_fixedRange(
                        full_filtered_copy_local[[score_copy]],
                        llrObj_copy$llr,
                        llrObj_copy$posDens,
                        llrObj_copy$negDens,
                        posScores_copy,
                        negScores_copy
                      )
                      dev.off()
                    }
                  )
                  
                  # 3b) Download bar plot (SVG)
                  output[[download_bar_svg_id]] <- downloadHandler(
                    filename = function() {
                      paste0(score_copy, "_BarPlot.svg")
                    },
                    content = function(file) {
                      # Open the SVG device
                      svg(file, width = 5, height = 5)
                      
                      # Construct the ggplot as usual
                      p <- ggplot(full_filtered_copy_local, aes(x=clinvar, fill=category)) +
                        geom_bar(position=position_stack(reverse=TRUE)) +
                        scale_fill_manual(values=category_colors) +
                        theme_minimal() +
                        labs(x="ClinVar Category", y="Variant Count", fill="LLR-Derived Evidence Category") +
                        theme(axis.text.x = element_text(angle=45, hjust=1))
                      
                      # Print the ggplot to render onto the SVG device
                      print(p)
                      
                      # Close the device
                      dev.off()
                    }
                  )
                  
                  # 3c) Download both (PDF)
                  output[[download_both_pdf_id]] <- downloadHandler(
                    filename = function() { paste0(score_copy, "_LLR_and_Bar.pdf") },
                    content = function(file) {
                      # Create a 2-page PDF: first page = LLR plot, second page = bar plot
                      pdf(file, width = 6, height = 6)
                      
                      # Page 1: LLR
                      drawDensityLLR_fixedRange(
                        full_filtered_copy_local[[score_copy]],
                        llrObj_copy$llr,
                        llrObj_copy$posDens,
                        llrObj_copy$negDens,
                        posScores_copy,
                        negScores_copy
                      )
                      
                      # Page 2: threshold table
                      plot.new() 
                      gridExtra::grid.table(crossings_df_copy)
                      
                      # Page 3: bar plot
                      #plot.new()  
                      p <- ggplot(full_filtered_copy_local, aes(x=clinvar, fill=category)) +
                        geom_bar(position=position_stack(reverse=TRUE)) +
                        scale_fill_manual(values=category_colors) +
                        theme_minimal() +
                        labs(x="ClinVar Category", y="Count", fill="LLR Category") +
                        theme(axis.text.x = element_text(angle=45, hjust=1))
                      print(p)
                      
                      dev.off()
                    }
                  )
                  
                  # 3d) Download CSV of variants w/ LLR + category
                  output[[download_csv_id]] <- downloadHandler(
                    filename = function() { paste0(score_copy, "_LLR_data.csv") },
                    content = function(file) {
                      # This CSV includes all variants (including VUS/conflicting) 
                      # that remain in full_filtered_copy_local, along with LLR and category
                      write.csv(full_filtered_copy_local, file, row.names = FALSE)
                    }
                  )
                  
                  # -----------------------------------
                  # 4) BUILD THE UI FOR THIS TAB
                  # -----------------------------------
                  llr_tabs[[length(llr_tabs) + 1]] <<- tabPanel(
                    score_copy,
                    fluidRow(
                      column(
                        5,
                        # LLR plot
                        plotOutput(plot_id, width = "500px", height = "600px")
                      ),
                      column(
                        4, offset = 3,
                        # Bar plot
                        plotOutput(plot_stack_id, width = "300px", height = "500px"),
                        
                        # Download buttons stacked vertically
                        tags$div(
                          style = "margin-top: 40px; text-align: left;",
                          downloadButton(download_llr_png_id,   "Download LLR (PNG)"),
                          tags$br(),
                          downloadButton(download_bar_svg_id,   "Download bar plot (SVG)"),
                          tags$br(),
                          downloadButton(download_both_pdf_id,  "Download both (PDF)"),
                          tags$br(),
                          downloadButton(download_csv_id,       "Download all variants with LLRs (CSV)")
                        )
                      )
                    ),
                    tableOutput(table_id)
                  )
                })
              }
            }
            
            if (length(llr_tabs) > 0) {
              output$Main_LLRTabs <- renderUI({
                do.call(tabsetPanel, c(id="llr_tabs", llr_tabs))
              })
            } else {
              output$Main_LLRTabs <- renderUI({
                "No LLR plots available"
              })
            }
            
            # Final step: completed
            incProgress(1.0, detail = "Done")
            removeModal()
            
            if (P_org < 11 || B_org < 11) {
              showModal(modalDialog(
                title = "Warning",
                "There are fewer than 11 Pathogenic/Likely Pathogenic or Benign/Likely Benign variants. 
           Please use extra caution when interpreting these plots.",
                easyClose = TRUE,
                footer = modalButton("Close")
              ))
            }
            
          } else {
            removeModal() # In case wait modal is still visible
            showModal(modalDialog(
              title = "Error",
              "Please select at least two variants.",
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
          }
        })
      })
      
      output$Main_PRC_Download_Buttons <- renderUI({
        req(plot_data())  # Ensure the plot data exists before showing buttons
        
        tags$div(
          style = "margin-top: 20px;",
          helpText(HTML("<strong>PRC Plot Download Options:</strong>")),
          downloadButton("Main_downloadPlotPNG", "Download PRC Plot as PNG"),
          downloadButton("Main_downloadPlotPDF", "Download PRC Plot and Metadata as PDF"),
          downloadButton("Main_downloadCSV", "Download Variants Used as CSV"),
          div(
            style = "display: inline-flex; align-items: center;",
            downloadButton("Main_downloadCSV_VUS", HTML("Download CSV with VUS")),
            actionLink("Main_helpButton_VUS", label = NULL, icon = icon("question-circle"), style = "margin-left: 5px;")
          )
        )
      })
      
      # Download PNG logic for Main App
      output$Main_downloadPlotPNG <- downloadHandler(
        filename = function() {
          paste("PRC_plot_", input$Main_gene, ".png", sep = "")
        },
        content = function(file) {
          plot_info <- plot_data()
          if (!is.null(plot_info)) {
            png(file, width = 6, height = 6, units = "in", res = 72)
            draw.prc(plot_info$yrobj, lty = plot_info$lty_styles, col = plot_info$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", ")))
            abline(h = 90, lty = "dashed")
            legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
            dev.off()
          }
        }
      )
      
      # Download PDF with metadata logic for Main App
      output$Main_downloadPlotPDF <- downloadHandler(
        filename = function() {
          paste("PRC_Report_", input$Main_gene, ".pdf", sep = "")
        },
        content = function(file) {
          plot_info <- plot_data()
          if (!is.null(plot_info)) {
            rmarkdown::render(input = "report_template.Rmd",
                              output_file = file,
                              params = list(
                                gene_s = plot_info$gene_s,
                                selected_scores = plot_info$selected_scores,
                                common_variant_filter = plot_info$common_variant_filter,
                                B_org = plot_info$B_org,
                                P_org = plot_info$P_org,
                                prcfiltered = plot_info$prcfiltered
                              ),
                              envir = new.env(parent = globalenv()))
          }
        }
      )
      
      # Download variants used as CSV
      output$Main_downloadCSV <- downloadHandler(
        filename = function() {
          paste("PRC_data_", input$Main_gene, ".csv", sep = "")
        },
        
        content = function(file) {
          selected_df <- selected_variants()
          if (!is.null(selected_df)) {
            write.csv(selected_df, file, row.names = FALSE)
          }
        }
      )
      
      # Download CSV with VUS
      output$Main_downloadCSV_VUS <- downloadHandler(
        filename = function() {
          paste("PRC_data_VUS_", input$Main_gene, ".csv", sep = "")
        },
        content = function(file) {
          df_VUS <- read.csv("full.csv", stringsAsFactors = FALSE)
          gene_s <- input$Main_gene
          
          prcfiltered <- df_VUS %>%
            filter(base__gene == gene_s) %>%
            arrange(clinvar) # Order alphabetical
          
          write.csv(prcfiltered, file, row.names = FALSE)
        }
      )
      
      # VUS explanation help button
      observeEvent(input$Main_helpButton_VUS, {
        showModal(modalDialog(
          title = "What does this data include?",
          HTML("In addition to P/LP or B/LB variants, this download includes VUS and Conflicting variants. 
          VUS stands for Variants of Uncertain Significance, and Conflicting means conflicting interpretations of pathogenicity. 
          These variants are not included in the PRC generation as they are not classified as Pathogenic or Benign. 
          However, by downloading this CSV, you can override the VUS or Conflicting annotations and reupload it as a custom reference set in Advanced Mode."),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      }, ignoreInit = TRUE)
      
      output$Main_download_buttons <- renderUI({
        req(plot_data())  # Only show if plot_data has been generated
        
        tagList(
          actionButton("Main_helpButton", "Plot Explanation"),
          helpText(HTML("<span style='color:black;'><strong>Download Options: </strong></span>")),
          downloadButton("Main_downloadAllPlotsZIP", "Download All Plots (ZIP)"),
          downloadButton("Main_downloadAllCSVsZIP", "Download All CSVs (ZIP)")
        )
      })
      
      ### New Download Handlers for ZIP Files ###

      ### --- DOWNLOAD ALL PLOTS AS ZIP ----------------------------------------------
      output$Main_downloadAllPlotsZIP <- downloadHandler(
        filename = function() {
          paste0("AllPlots_", input$Main_gene, ".zip")
        },
        content = function(zipfile) {
          req(plot_data())
          
          # 1) Create a temporary folder for saving images
          tmpdir <- file.path(tempdir(), paste0("plotszip_", Sys.Date()))
          if (!dir.exists(tmpdir)) dir.create(tmpdir)
          
          # 2) PRC Plot as PNG (same as before)
          plot_info <- plot_data()
          prc_png_filename <- file.path(tmpdir, paste0("PRC_plot_", plot_info$gene_s, ".png"))
          png(prc_png_filename, width = 6, height = 6, units = "in", res = 72)
          draw.prc(plot_info$yrobj,
                   lty = plot_info$lty_styles,
                   col = plot_info$col_styles,
                   lwd = 2,
                   balanced = TRUE,
                   main = paste0(plot_info$gene_s, " PRCs for ",
                                 paste(plot_info$selected_scores, collapse = ", ")))
          abline(h = 90, lty = "dashed")
          legend("left",
                 legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org),
                            paste("# of Benign and Likely Benign:", plot_info$B_org)),
                 pch = 15, bty = "n")
          dev.off()
          
          # 3) LLR Plots for each predictor, plus the bar plots as PNG
          llr_scores <- intersect(plot_info$selected_scores, c("VARITY", "REVEL", "AlphaMissense"))
          selected_df <- plot_info$prcfiltered      # The data used for PRC
          full_filtered <- full_data()             # The full data (includes VUS/conflicting)
          
          # Rename to match your "score" columns
          names(full_filtered)[names(full_filtered) == "VEP_varity_r"] <- "VARITY"
          names(full_filtered)[names(full_filtered) == "VEP_alphamissense__pathogenicity"] <- "AlphaMissense"
          names(full_filtered)[names(full_filtered) == "VEP_revel__score"] <- "REVEL"
          names(full_filtered)[names(full_filtered) == "gnomad__af"] <- "gnomAD_AF"
          
          # Use the same color scheme as in your original bar plot code
          category_colors <- c(
            "benign_strong"   = "dodgerblue",
            "benign_support"  = "lightblue",
            "none"            = "grey",
            "patho_support"   = "#FFC0CB",
            "patho_moderate"  = "#FF9999",
            "patho_strong"    = "#FF6666",
            "patho_vstrong"   = "red"
          )
          
          for (score in llr_scores) {
            posScores <- na.omit(setNames(
              selected_df[selected_df$clinvar == TRUE, score],
              selected_df[selected_df$clinvar == TRUE, "base__achange"]
            ))
            negScores <- na.omit(setNames(
              selected_df[selected_df$clinvar == FALSE, score],
              selected_df[selected_df$clinvar == FALSE, "base__achange"]
            ))
            
            if (length(posScores) > 0 && length(negScores) > 0) {
              # Build the LLR object
              llrObj <- buildLLR.kernel(posScores, negScores)
              
              # Assign LLR and categories to full data for bar plot
              full_copy <- full_filtered
              full_copy$llr <- llrObj$llr(full_copy[[score]])
              
              breaks <- c(-Inf, -1.27, -0.32, 0.32, 0.64, 1.27, 2.54, Inf)
              labels <- c("benign_strong", "benign_support", "none",
                          "patho_support", "patho_moderate", "patho_strong", "patho_vstrong")
              full_copy$category <- cut(full_copy$llr, breaks, labels, include.lowest = TRUE)
              
              # --- 3a) LLR Density Plot (PNG) ---
              llr_png_path <- file.path(tmpdir, paste0("LLR_", score, ".png"))
              png(llr_png_path, width = 500, height = 600, res = 72)
              drawDensityLLR_fixedRange(
                full_copy[[score]],
                llrObj$llr,
                llrObj$posDens,
                llrObj$negDens,
                posScores,
                negScores
              )
              dev.off()
              
              # --- 3b) Bar Plot as PNG (mirror the original bar plot code) ---
              bar_png_path <- file.path(tmpdir, paste0("Bar_", score, ".png"))
              png(bar_png_path, width = 300, height = 500, res = 72)
              
              full_copy$clinvar <- factor(full_copy$clinvar,
                                          levels = c("P/LP", "B/LB", "VUS", "Conflicting"))
              
              bar_plot <- ggplot(full_copy, aes(x = clinvar, fill = category)) +
                geom_bar(position = position_stack(reverse = TRUE)) +
                scale_fill_manual(values = category_colors) +
                theme_minimal() +
                labs(x = "ClinVar Category", y = "Variant Count",
                     fill = "LLR-Derived Evidence Category") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
              
              print(bar_plot)
              
              dev.off()
            }
          }
          
          # 4) Zip all images
          old_wd <- setwd(tmpdir)
          on.exit(setwd(old_wd), add = TRUE)
          zip::zipr(zipfile = zipfile, files = list.files(tmpdir))
        }
      )
      
      ### --- DOWNLOAD ALL CSVS AS ZIP -----------------------------------------------
      output$Main_downloadAllCSVsZIP <- downloadHandler(
        filename = function() {
          paste0("AllCSVs_", input$Main_gene, ".zip")
        },
        content = function(zipfile) {
          req(plot_data())
          
          tmpdir <- file.path(tempdir(), paste0("csvzip_", Sys.Date()))
          if (!dir.exists(tmpdir)) dir.create(tmpdir)
          
          # 1) Main CSV of selected variants (already existing in your code)
          selected_df <- selected_variants()
          main_csv_path <- file.path(tmpdir, paste0("PRC_data_", input$Main_gene, ".csv"))
          if (!is.null(selected_df)) {
            write.csv(selected_df, main_csv_path, row.names = FALSE)
          }
          
          # 2) CSV with VUS (full.csv for the chosen gene) (already existing in your code)
          df_VUS <- read.csv("full.csv", stringsAsFactors = FALSE)
          prcfiltered_VUS <- df_VUS %>%
            filter(base__gene == input$Main_gene) %>%
            arrange(clinvar)
          vus_csv_path <- file.path(tmpdir, paste0("PRC_data_VUS_", input$Main_gene, ".csv"))
          write.csv(prcfiltered_VUS, vus_csv_path, row.names = FALSE)
          
          # 3) Single CSV with LLR results for all selected predictors
          plot_info <- plot_data()
          llr_scores <- intersect(plot_info$selected_scores, c("VARITY", "REVEL", "AlphaMissense"))
          
          # (a) Get the final data used for PRC (for posScores/negScores)
          selected_df <- plot_info$prcfiltered
          
          # (b) Pull the full data (includes VUS/conflicting)
          full_filtered <- full_data()
          
          # (c) Rename columns for your typical "score" references:
          names(full_filtered)[names(full_filtered) == "VEP_varity_r"] <- "VARITY"
          names(full_filtered)[names(full_filtered) == "VEP_alphamissense__pathogenicity"] <- "AlphaMissense"
          names(full_filtered)[names(full_filtered) == "VEP_revel__score"] <- "REVEL"
          names(full_filtered)[names(full_filtered) == "gnomad__af"] <- "gnomAD_AF"
          
          # We will create new columns for each predictor in "full_copy".
          # (d) Make a working copy so we don't mutate full_data() globally:
          all_llr_df <- full_filtered
          
          for (score in llr_scores) {
            
            # Prepare positives and negatives from your final selected_df
            posScores <- na.omit(
              setNames(selected_df[selected_df$clinvar == TRUE, score],
                       selected_df[selected_df$clinvar == TRUE, "base__achange"])
            )
            negScores <- na.omit(
              setNames(selected_df[selected_df$clinvar == FALSE, score],
                       selected_df[selected_df$clinvar == FALSE, "base__achange"])
            )
            
            # If at least one P/LP and one B/LB to build LLR:
            if (length(posScores) > 0 && length(negScores) > 0) {
              llrObj <- buildLLR.kernel(posScores, negScores)
              
              # Create the LLR column named PREDICTOR_llr
              llr_colname <- paste0(score, "_llr")
              all_llr_df[[llr_colname]] <- llrObj$llr(all_llr_df[[score]])
              
              # Create the category column named PREDICTOR_category
              cat_colname <- paste0(score, "_category")
              breaks <- c(-Inf, -1.27, -0.32, 0.32, 0.64, 1.27, 2.54, Inf)
              labels <- c("benign_strong", "benign_support", "none",
                          "patho_support", "patho_moderate", "patho_strong", "patho_vstrong")
              
              # Cut on the LLR values we just created
              all_llr_df[[cat_colname]] <- cut(all_llr_df[[llr_colname]], 
                                               breaks = breaks, 
                                               labels = labels, 
                                               include.lowest = TRUE)
            } else {
              # If not enough data to build LLR, fill with NA or skip, your choice
              llr_colname <- paste0(score, "_llr")
              cat_colname <- paste0(score, "_category")
              all_llr_df[[llr_colname]] <- NA
              all_llr_df[[cat_colname]] <- NA
            }
          }
          
          # (e) Write out a single CSV for all LLR-enabled predictors:
          #     E.g., "All_LLR_data.csv"
          llr_csv_path <- file.path(tmpdir, "All_LLR_data.csv")
          write.csv(all_llr_df, llr_csv_path, row.names = FALSE)
          
          # 4) Zip everything
          old_wd <- setwd(tmpdir)
          on.exit(setwd(old_wd), add = TRUE)
          zip::zipr(zipfile = zipfile, files = list.files(tmpdir))
        }
      )
      
      
      # Plot explanation help button
      observeEvent(input$Main_helpButton, {
        showModal(modalDialog(
          title = "Plot Explanation",
          HTML("<p><b>What is balanced precision?</b></p>
            <p>Balanced precision is useful in situations where the class distribution of the test set is imbalanced. In other words, it is the precision that would have been expected had the proportion of 'P/LP' examples been balanced (equal to 50%).</p>"),
          p("Definition ", 
            a("here", href = "https://doi.org/10.1016/j.ajhg.2021.08.012", target = "_blank")),
          HTML("<p><b>What are R90BP and AUBPRC?</b></p>
          <p><b>R90BP:</b> The recall achieved at a stringent (90%) balanced precision threshold </p>
          <p><b>AUBPRC:</b> The area under the BPRC curve</p>"),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      })
      
    } else {
      # Logic for Advanced
      rv <- reactiveValues(
        state = reactiveValues(gnomad_filter_inserted = FALSE)
      )
      
      # Define reactive values
      rv$prcdata_own <- reactiveVal(NULL)
      rv$prcdata_fetch <- reactiveVal(NULL)
      rv$plot_data <- reactiveVal(NULL)
      rv$variant_data_df <- reactiveVal(NULL)
      rv$threshold_data <- reactiveVal(NULL) # for storing PRC thresholds - REFACTOR: include as part of prc data?
      
      # Observe 'input$input_type' to reset reactive values and inputs when user switches option
      observeEvent(input$input_type, {
        rv$plot_data(NULL)
        rv$prcdata_fetch(NULL)
        rv$prcdata_own(NULL)
        rv$gnomad_filter_inserted <- FALSE
        
        updateCheckboxGroupInput(session, "scores", choices = NULL, selected = NULL)
        removeUI(selector = "#gnomad_filter_wrapper", immediate = TRUE)
      })
      
      # Logic for 'own' input_type
      observeEvent(input$file_full, {
        if (input$input_type == "own") {
          df <- tryCatch({
            req(input$file_full)
            df <- read.csv(input$file_full$datapath, stringsAsFactors = FALSE)
            
            # Check for necessary columns
            if (!all(c("base__gene", "clinvar") %in% colnames(df))) {
              showModal(modalDialog(
                title = "Error",
                "The uploaded CSV must contain the columns 'base__gene' and 'clinvar'.",
                easyClose = TRUE,
                footer = modalButton("Close")
              ))
              return(NULL)
            }
            
            if (!(any(df$clinvar == "P/LP") && any(df$clinvar == "B/LB"))) {
              showModal(modalDialog(
                title = "Error",
                "You must have at least one P/LP and one B/LB variant.",
                easyClose = TRUE,
                footer = modalButton("Close")
              ))
              return(NULL)
            }
            
            # Consistency for mandatory columns
            colnames(df)[colnames(df) == "gnomad__af"] <- "gnomAD_AF"
            colnames(df)[colnames(df) == "clinvar"] <- "clinvar"
            
            df
          }, error = function(e) {
            showModal(modalDialog(
              title = "Error",
              "There was an error reading the uploaded file. Please ensure it is a valid CSV file.",
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
            NULL
          })
          rv$prcdata_own(df)  # Set the reactiveVal
        }
      })
      
      # Guide on how to format user-inputted csv
      observeEvent(input$upload_guide, {
        showModal(modalDialog(
          title = "Reference Set Format Information",
          HTML("Please ensure your own reference set is a CSV file.<br><br>
            Mandatory columns:<br>
            <b>base__gene:</b> Gene name<br>
            <b>clinvar:</b> Variant interpretation (B/LB or P/LP)<br><br>
            
            Optional columns: <br>
            <b>gnomad_af:</b> GnomAD allele frequency - for filtering out common variants<br><br>
            
            For any predictors or custom scores that you would like to include, put 'VEP_' before the name of the column. For example:<br>
            <b>VEP_alphamissense:</b> AlphaMissense score<br>
            <b>VEP_custom:</b> Additional scores you would like to include<br><br>"),
          easyClose = FALSE,
          footer = tagList(
            downloadButton("download_template_own", "Download CSV Template"),
            modalButton("Close")
          )
        ))
      }, ignoreInit = TRUE)
      
      # Download template for own reference set
      output$download_template_own <- downloadHandler(
        filename = function() {
          "reference_set_template.csv"
        },
        content = function(file) {
          file.copy("own_template.csv", file)
        }
      )
      
      # Download OC reference template
      output$download_template_oc <- downloadHandler(
        filename = function() {
          "reference_set_template.csv"
        },
        content = function(file) {
          file.copy("oc_template.csv", file)
        }
      )
      
      # Logic for 'fetch' input_type
      rv$fetching_in_progress <- FALSE  # Track fetch state
      observeEvent(input$fetchButton, {
        if (input$input_type == "fetch") {
          if (rv$fetching_in_progress) return()  # Prevent duplicate execution
          rv$fetching_in_progress <- TRUE  # Lock fetching process
          
          # Disable the fetch button while fetching
          disable("fetchButton")
          # Initialize empty list to store results
          all_results <- list()
          
          tryCatch({
            req(input$file_fetch)
            df <- read.csv(input$file_fetch$datapath, stringsAsFactors = FALSE)
            
            # Ensure 'chrom' column in df always has the "chr" prefix for consistency
            df <- df %>%
              mutate(chrom = ifelse(grepl("^chr", as.character(chrom)), as.character(chrom), paste0("chr", chrom)))
            
            # Show progress while fetching data
            withProgress(message = 'Fetching Variant Data', value = 0, {
              for (i in 1:nrow(df)) {
                
                chrom <- df$chrom[i]
                pos <- df$pos[i]
                ref <- df$ref_base[i]
                alt <- df$alt_base[i]
                
                # Construct API URL for chrom-based input
                api_url <- paste0(
                  "https://run.opencravat.org/submit/annotate?",
                  "chrom=", chrom,
                  "&pos=", pos,
                  "&ref_base=", ref,
                  "&alt_base=", alt,
                  "&annotators=", paste(sub("varity_er", "varity_r", tolower(input$Fetch_scores)), collapse = ",") 
                )
                
                response <- GET(api_url)
                
                result <- fromJSON(content(response, "text"), flatten = TRUE)
                
                # Extract gene and achange into corresponding columns
                gene <- ifelse(!is.null(result$crx$hugo), result$crx$hugo, NA)
                achange <- ifelse(!is.null(result$crx$achange), result$crx$achange, NA)
                
                # Initialize result_list for required fields
                result_list <- list(
                  chrom = chrom,
                  pos = pos,
                  ref_base = ref,
                  alt_base = alt
                )
                
                # Add gene and achange columns only if they don't already exist in df
                if (!"base__gene" %in% colnames(df)) {
                  result_list[["base__gene"]] <- gene
                }
                if (!"base__achange" %in% colnames(df)) {
                  result_list[["base__achange"]] <- achange
                }
                
                # Mapping possible scores from API result to data
                score_mappings <- list(
                  "gnomad" = list(
                    column_name = "gnomAD_AF",
                    value = if (!is.null(result$gnomad$af)) result$gnomad$af else NA
                  ),
                  "varity_r" = list(
                    column_name = "VEP_VARITY_R",
                    value = if (!is.null(result$varity_r$varity_r)) result$varity_r$varity_r else NA
                  ),
                  "varity_er" = list(
                    column_name = "VEP_VARITY_ER",
                    value = if (!is.null(result$varity_r$varity_er)) result$varity_r$varity_er else NA
                  ),
                  "revel" = list(
                    column_name = "VEP_REVEL",
                    value = if (!is.null(result$revel$score)) result$revel$score else NA
                  ),
                  "alphamissense" = list(
                    column_name = "VEP_AlphaMissense",
                    value = if (!is.null(result$alphamissense$am_pathogenicity)) result$alphamissense$am_pathogenicity else NA
                  ),
                  "provean" = list(
                    column_name = "VEP_Provean",
                    value = if (!is.null(result$provean$rankscore)) result$provean$rankscore else NA
                  ),
                  "sift" = list(
                    column_name = "VEP_SIFT",
                    value = if (!is.null(result$sift$rankscore)) result$sift$rankscore else NA
                  ),
                  "polyphen2" = list(
                    column_name = "VEP_PolyPhen2",
                    value = if (!is.null(result$polyphen2$hvar_rank)) result$polyphen2$hvar_rank else NA
                  ),
                  # Convert clinvar to P or B, ignore VUS or conflicting
                  "clinvar" = list(
                    column_name = "clinvar",
                    value = if (!is.null(result$clinvar$sig)) {
                      clinvar_sig <- tolower(result$clinvar$sig)
                      if (grepl("benign", clinvar_sig)) {
                        "B/LB"
                      } else if (grepl("conflicting", clinvar_sig)) {
                        "Conflicting"
                      } else if (grepl("pathogenic", clinvar_sig)) {
                        "P/LP"
                      } else if (grepl("uncertain", clinvar_sig)) {
                        "VUS"
                      } else {
                        NA
                      }
                    } else {
                      NA
                    }
                  )
                )
                
                for (score in input$Fetch_scores) {
                  if (score %in% names(score_mappings)) {
                    mapping <- score_mappings[[score]]
                    if (mapping$column_name %in% colnames(df)) {
                      # Update column if it already exists
                      df[[mapping$column_name]][i] <- mapping$value
                    } else {
                      # Add as new column
                      result_list[[mapping$column_name]] <- mapping$value
                    }
                  }
                }
                
                result_df <- as.data.frame(result_list, stringsAsFactors = FALSE)
                all_results[[i]] <- result_df
                
                # Show progress bar
                incProgress(1 / nrow(df))
                showNotification(paste(i, "/", nrow(df), "variants fetched"), duration = 3, type = "message")
              }
            })
            
            # Combine all results into a data frame
            variant_data_df <- do.call(rbind, all_results)
            
            variant_data_df <- dplyr::left_join(df, variant_data_df, by = c("chrom", "pos", "ref_base", "alt_base"))
            
            output$errorText <- renderText("")
            rv$prcdata_fetch(variant_data_df)  # Set the reactiveVal
            rv$variant_data_df <- variant_data_df

          }, error = function(e) {
            showModal(modalDialog(
              title = "Error",
              "An error occurred while fetching data.",
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
            rv$prcdata_fetch(NULL)
          })
        }
        # Re-enable the fetch button & reset flag
        enable("fetchButton")
        rv$fetching_in_progress <- FALSE
      })
      
      
      # Define 'prcdata' reactive depending on where it comes from (own or fetched)
      prcdata <- reactive({
        if (input$input_type == "own") {
          rv$prcdata_own()
        } else if (input$input_type == "fetch") {
          rv$prcdata_fetch()
        } else {
          NULL
        }
      })
      
      ### Common variant filter insertion
      # Ensure the checkbox is reset when switching tabs
      observeEvent(input$navbar, {
        if (input$navbar == "Advanced") {
          rv$gnomad_filter_inserted <- FALSE  # Reset flag when entering Advanced mode
          removeUI(selector = "#gnomad_filter_wrapper", immediate = TRUE)
        }
      })
      
      # Ensure the checkbox is removed when switching input types (own <-> fetch)
      observeEvent(input$input_type, {
        rv$gnomad_filter_inserted <- FALSE  # Reset when switching input source
        removeUI(selector = "#gnomad_filter_wrapper", immediate = TRUE)
      })
      
      # Observe when prcdata changes to add gnomAD checkbox
      observeEvent(prcdata(), {
        df <- prcdata()
        if (is.null(df)) return()
        
        # Check for gnomAD columns, case-insensitive
        colnames_lower <- tolower(colnames(df))
        gnomad_columns <- grep("gnomad", colnames_lower, value = TRUE)
        
        # If gnomAD columns are found and the checkbox is not yet added
        if (length(gnomad_columns) > 0 && isFALSE(rv$gnomad_filter_inserted)) {
          # Remove any existing checkbox to prevent duplicates
          removeUI(selector = "#gnomad_filter_wrapper", immediate = TRUE)
          
          # Insert the checkbox
          insertUI(
            selector = "#scores",
            where = "beforeBegin",
            ui = div(id = "gnomad_filter_wrapper",
                     checkboxInput("common_variant_filter", 
                                   "Exclude Common Variants (gnomAD AF > 0.005)", 
                                   value = TRUE)
            )
          )
          rv$gnomad_filter_inserted <- TRUE
          
        } else if (length(gnomad_columns) == 0 && isTRUE(rv$gnomad_filter_inserted)) {
          # Remove the checkbox if no gnomAD columns are found and it is still there
          removeUI(selector = "#gnomad_filter_wrapper", immediate = TRUE)
          rv$gnomad_filter_inserted <- FALSE
        }
        
        # Remove "VEP_" prefix from predictor column names for display
        predictor_columns <- colnames(df)[grepl("^VEP_", colnames(df))]
        predictor_columns <- gsub("^VEP_", "", predictor_columns)
        # Update checkboxGroupInput with trimmed predictor columns
        updateCheckboxGroupInput(session, "scores", choices = predictor_columns, selected = predictor_columns)
      })
      
      # Observe 'plotButton'
      observeEvent(input$plotButton, {
        # Wrap everything in a withProgress block
        withProgress(message = "Generating plots...", value = 0, {
          df <- prcdata()
          if (is.null(df) || nrow(df) == 0) {
            output$errorText <- renderText("Not enough rows to generate the PRC plot.")
            return()
          }
          
          gene_s <- ifelse("base__gene" %in% colnames(df), df$base__gene[1], "Custom Gene")
          exclude_common_variants <- input$common_variant_filter
          selected_scores <- input$scores
          
          if (isTRUE(exclude_common_variants) && "gnomAD_AF" %in% colnames(df)) { # Added isTRUE
            df <- df[is.na(df$gnomAD_AF) | df$gnomAD_AF <= 0.005, ]
          }
          
          # Remove "VEP_" prefix from predictor column names
          colnames(df) <- gsub("^VEP_", "", colnames(df))
          
          df <- df[order(df$clinvar), ]
          
          prcfiltered <- df %>%
            mutate(clinvar = ifelse(clinvar == "P/LP", TRUE, ifelse(clinvar == "B/LB", FALSE, NA)))
          
          prcfiltered <- prcfiltered[!is.na(prcfiltered$clinvar), ]
          
          # Filter out any rows with an NA in any `selected_scores`
          prcfiltered <- prcfiltered[rowSums(is.na(prcfiltered[selected_scores])) == 0, ]
          
          P_org <- sum(prcfiltered$clinvar == TRUE)
          B_org <- sum(prcfiltered$clinvar == FALSE)
          
          tryCatch({
            yrobj <- yr2(truth = prcfiltered[["clinvar"]], scores = prcfiltered[selected_scores], high = c(rep(TRUE, length(selected_scores) - 1), FALSE)) # high = rep(TRUE, length(selected_scores)
            
            # Added for threshold calculation
            thresh_ranges <- calculate_thresh_range_flipper(yrobj, x = 0.9, balanced = TRUE, high = c(rep(TRUE, length(selected_scores) - 1), FALSE))
            
            # Store thresholds for rendering
            threshold_data(thresh_ranges)
            
            # Extra colors for lines more than 3
            # Exclude white and near-white colors
            available_colors <- grDevices::colors()[!grepl("white|ivory|seashell|snow|honeydew|azure|aliceblue|mintcream|ghostwhite", grDevices::colors(), ignore.case = TRUE)]
            
            # Generate random colors for additional predictors if more than three
            num_scores <- min(length(selected_scores), 3)  # Only take a maximum of three predictors for specific styles
            extra_colors <- if (length(selected_scores) > 3) {
              sample(available_colors, length(selected_scores) - 3)
            } else {
              NULL
            }
            
            # Assign styles for the first three predictors, followed by additional random colors if needed
            rv$plot_data(list(
              yrobj = yrobj,
              lty_styles = c("dashed", "solid", "dashed")[1:num_scores],
              col_styles = c("purple", "cadetblue2", "orange")[1:num_scores],
              gene_s = gene_s,
              selected_scores = selected_scores,
              B_org = B_org,
              P_org = P_org,
              common_variant_filter = exclude_common_variants,
              prcfiltered = prcfiltered,
              extra_colors = extra_colors  # Store extra colors for additional predictors if any
            ))
            
            output$errorText <- renderText("")
            
          }, error = function(e) {
            rv$plot_data(NULL)
            output$errorText <- renderText("Error generating plot.")
          })
          
          output$PRCPlot <- renderPlot({
            plot_info <- rv$plot_data()
            if (!is.null(plot_info)) {
              tryCatch({
                colors_to_use <- c(plot_info$col_styles, plot_info$extra_colors)
                draw.prc(
                  plot_info$yrobj,
                  lty = c(plot_info$lty_styles, rep(c("solid", "dashed"), length(plot_info$extra_colors))),
                  col = colors_to_use,
                  lwd = 2,
                  balanced = TRUE,
                  main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", "))
                )
                abline(h = 90, lty = "dashed")
                legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
              }, error = function(e) {
                showModal(modalDialog(
                  title = 'Error',
                  'Not enough data - must have at least one pathogenic and benign.',
                  easyClose = TRUE,
                  footer = NULL
                ))
              })
            }
          }, width = 600, height = 600, res = 72)
          
          # Render the threshold table
          output$thresholdTableUI <- renderUI({
            req(threshold_data())
            
            tagList(
              h5("VEP threshold to achieve 90% Balanced Precision"),  # Add the title only if the table exists
              tableOutput("thresholdTable")  
            )
          })
          
          output$thresholdTable <- renderTable({
            req(threshold_data())
            
            df <- threshold_data()
            
            # Round numeric columns to 3 decimals if numeric
            df <- as.data.frame(lapply(df, function(x) {
              if (is.numeric(x)) {
                format(round(x, 3), nsmall = 3)
              } else {
                x
              }
            }))
            df
          }, rownames = FALSE)
          
          ##############################################################
          ## ADDED LLR SECTION                                        ##
          ##############################################################
          
          # 1) Initialize list for LLR tabs
          llr_tabs <- list()
          llr_scores <- paste0("VEP_", selected_scores)
          
          incProgress(0.5, detail = "Reading full LLR data (this might take a while)")
          # Retrieve filtered full data
          full_filtered <- prcdata()

          for (score_idx in seq_along(llr_scores)) {
            score <- llr_scores[score_idx]
            incProgress(0.5+0.5/length(llr_scores), detail = paste("Generating LLR plot for:", score))

            posScores <- na.omit(setNames(
              full_filtered[full_filtered$clinvar == "P/LP", score],
              full_filtered[full_filtered$clinvar == "P/LP", "base__achange"]
            ))
            negScores <- na.omit(setNames(
              full_filtered[full_filtered$clinvar == "B/LB", score],
              full_filtered[full_filtered$clinvar == "B/LB", "base__achange"]
            ))

            if (length(posScores) > 0 & length(negScores) > 0) {
              llrObj <- buildLLR.kernel(posScores, negScores, outlierSuppression=0.001) # Change outlier suppression
              
              full_filtered_copy <- full_filtered
              full_filtered_copy$llr <- llrObj$llr(full_filtered_copy[[score]]) # full_filtered_copy has llr values
              
              # Define breakpoints and labels for the categories
              breaks <- c(-Inf, -1.27, -0.32, 0.32, 0.64, 1.27, 2.54, Inf)
              labels <- c(
                "benign_strong",   # (-Inf, -1.27]
                "benign_support",  # (-1.27, -0.32]
                "none",            # (-0.32, 0.32]
                "patho_support",   # (0.32, 0.64]
                "patho_moderate",  # (0.64, 1.27]
                "patho_strong",    # (1.27, 2.54]
                "patho_vstrong"    # (2.54, Inf)
              )
              
              # Categorize llr based on the defined bins
              full_filtered_copy$category <- cut(full_filtered_copy$llr, breaks=breaks, labels=labels, include.lowest=TRUE)
              
              full_filtered_copy <- full_filtered_copy %>% filter(!is.na(category)) # filter NA - due to score and llr being NA
              
              # Define thresholds and search range for crossings
              llrTs <- llrThresholds(optiLLR(0.1))
              x_range <- range(c(posScores, negScores))
              
              # Find where the LLR function crosses the thresholds
              crossings_df <- findLLRcrossings(llrObj$llr, llrTs, x_range)
              
              # Create the stacked bar plot of clinvar vs category
              category_colors <- c(
                "benign_strong"   = "dodgerblue",
                "benign_support"  = "lightblue",
                "none"            = "grey",
                "patho_support"   = "#FFC0CB",   # light pink
                "patho_moderate"  = "#FF9999",   # medium pink
                "patho_strong"    = "#FF6666",   # darker pink
                "patho_vstrong"   = "red"
              )
              
              plot_stack_id <- paste0("Main_Stacked_", score)
              
              # Preserve current iteration's variables, or else it would be same plot
              local({
                # Make local copies of the variables to avoid referencing loop variables
                score_copy <- score
                posScores_copy <- posScores
                negScores_copy <- negScores
                llrObj_copy <- llrObj
                full_filtered_copy_local <- full_filtered_copy
                crossings_df_copy <- crossings_df
                
                plot_id <- paste0("Main_LLRPlot_", score_copy)
                table_id <- paste0("Main_LLRCrossingsTable_", score_copy)
                
                # -----------------------------
                # 1) DEFINE DOWNLOAD BUTTON IDS
                # -----------------------------
                download_llr_png_id <- paste0("download_", score_copy, "_llr_png")
                download_bar_svg_id <- paste0("download_", score_copy, "_bar_svg")
                download_both_pdf_id <- paste0("download_", score_copy, "_both_pdf")
                download_csv_id      <- paste0("download_", score_copy, "_csv")
                
                output[[plot_id]] <- renderPlot({
                  tryCatch({
                    drawDensityLLR_fixedRange(
                      full_filtered_copy_local[[score_copy]], 
                      llrObj_copy$llr, 
                      llrObj_copy$posDens, 
                      llrObj_copy$negDens, 
                      posScores_copy, 
                      negScores_copy
                    )
                  }, error = function(e) {
                    showModal(modalDialog(
                      title = 'Error',
                      'Not enough data',
                      easyClose = TRUE,
                      footer = NULL
                    ))
                  })
                }, width = 500, height = 600, res = 72)
                
                # Render the stacked bar chart
                output[[plot_stack_id]] <- renderPlot({
                  # Factor ordering if needed
                  full_filtered_copy_local$clinvar <- factor(full_filtered_copy_local$clinvar, 
                                                             levels = c("P/LP", "B/LB", "VUS", "Conflicting"))
                  ggplot(full_filtered_copy_local, aes(x=clinvar, fill=category)) +
                    geom_bar(position=position_stack(reverse=TRUE)) + # reverse to put pathogenic at top like the LLR
                    scale_fill_manual(values=category_colors) +
                    theme_minimal() +
                    labs(x="ClinVar Category", y="Variant Count", fill="LLR-Derived Evidence Category") +
                    theme(axis.text.x = element_text(angle=45, hjust=1))
                }, width = 300, height = 500, res = 72)
                
                # Render the table of LLR threshold crossings
                output[[table_id]] <- renderTable({
                  crossings_df_copy
                })
                # jumptag
                
                # -----------------------------------
                # 3) DEFINE DOWNLOAD HANDLERS
                # -----------------------------------
                
                # 3a) Download LLR (PNG)
                output[[download_llr_png_id]] <- downloadHandler(
                  filename = function() { paste0(score_copy, "_LLRPlot.png") },
                  content = function(file) {
                    # Use a PNG device and re-draw the same LLR plot
                    png(file, width = 1500, height = 1800, res = 300)
                    drawDensityLLR_fixedRange(
                      full_filtered_copy_local[[score_copy]],
                      llrObj_copy$llr,
                      llrObj_copy$posDens,
                      llrObj_copy$negDens,
                      posScores_copy,
                      negScores_copy
                    )
                    dev.off()
                  }
                )
                
                # 3b) Download bar plot (SVG)
                output[[download_bar_svg_id]] <- downloadHandler(
                  filename = function() {
                    paste0(score_copy, "_BarPlot.svg")
                  },
                  content = function(file) {
                    # Open the SVG device
                    svg(file, width = 5, height = 5)
                    
                    # Construct the ggplot as usual
                    p <- ggplot(full_filtered_copy_local, aes(x=clinvar, fill=category)) +
                      geom_bar(position=position_stack(reverse=TRUE)) +
                      scale_fill_manual(values=category_colors) +
                      theme_minimal() +
                      labs(x="ClinVar Category", y="Variant Count", fill="LLR-Derived Evidence Category") +
                      theme(axis.text.x = element_text(angle=45, hjust=1))
                    
                    # Print the ggplot to render onto the SVG device
                    print(p)
                    
                    # Close the device
                    dev.off()
                  }
                )
                
                # 3c) Download both (PDF)
                output[[download_both_pdf_id]] <- downloadHandler(
                  filename = function() { paste0(score_copy, "_LLR_and_Bar.pdf") },
                  content = function(file) {
                    # Create a 2-page PDF: first page = LLR plot, second page = bar plot
                    pdf(file, width = 6, height = 6)
                    
                    # Page 1: LLR
                    drawDensityLLR_fixedRange(
                      full_filtered_copy_local[[score_copy]],
                      llrObj_copy$llr,
                      llrObj_copy$posDens,
                      llrObj_copy$negDens,
                      posScores_copy,
                      negScores_copy
                    )
                    
                    # Page 2: threshold table
                    plot.new() 
                    gridExtra::grid.table(crossings_df_copy)
                    
                    # Page 3: bar plot
                    #plot.new()  
                    p <- ggplot(full_filtered_copy_local, aes(x=clinvar, fill=category)) +
                      geom_bar(position=position_stack(reverse=TRUE)) +
                      scale_fill_manual(values=category_colors) +
                      theme_minimal() +
                      labs(x="ClinVar Category", y="Count", fill="LLR Category") +
                      theme(axis.text.x = element_text(angle=45, hjust=1))
                    print(p)
                    
                    dev.off()
                  }
                )
                
                # 3d) Download CSV of variants w/ LLR + category
                output[[download_csv_id]] <- downloadHandler(
                  filename = function() { paste0(score_copy, "_LLR_data.csv") },
                  content = function(file) {
                    # This CSV includes all variants (including VUS/conflicting) 
                    # that remain in full_filtered_copy_local, along with LLR and category
                    write.csv(full_filtered_copy_local, file, row.names = FALSE)
                  }
                )
                
                # -----------------------------------
                # 4) BUILD THE UI FOR THIS TAB
                # -----------------------------------
                llr_tabs[[length(llr_tabs) + 1]] <<- tabPanel(
                  score_copy,
                  fluidRow(
                    column(
                      5,
                      # LLR plot
                      plotOutput(plot_id, width = "500px", height = "600px")
                    ),
                    column(
                      4, offset = 3,
                      # Bar plot
                      plotOutput(plot_stack_id, width = "300px", height = "500px"),
                      
                      # Download buttons stacked vertically
                      tags$div(
                        style = "margin-top: 40px; text-align: left;",
                        downloadButton(download_llr_png_id,   "Download LLR (PNG)"),
                        tags$br(),
                        downloadButton(download_bar_svg_id,   "Download bar plot (SVG)"),
                        tags$br(),
                        downloadButton(download_both_pdf_id,  "Download both (PDF)"),
                        tags$br(),
                        downloadButton(download_csv_id,       "Download all variants with LLRs (CSV)")
                      )
                    )
                  ),
                  tableOutput(table_id) # jumptag
                )
              })
            }
          }
          
          if (length(llr_tabs) > 0) {
            output$LLRTabs <- renderUI({
              do.call(tabsetPanel, c(id="llr_tabs", llr_tabs))
            })
          } else {
            output$LLRTabs <- renderUI({
              "No LLR plots available"
            })
          }
          
          # Final step: completed
          incProgress(1.0, detail = "Done")
          removeModal()
          
          if (P_org < 11 || B_org < 11) {
            showModal(modalDialog(
              title = "Warning",
              "There are fewer than 11 Pathogenic/Likely Pathogenic or Benign/Likely Benign variants. 
             Please use extra caution when interpreting these plots.",
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
          }
        })
      })
      
      # PRC Download logic 
      output$PRC_Download_Buttons <- renderUI({
        req(rv$plot_data())  # Ensure the plot data exists before showing buttons
        
        tags$div(
          style = "margin-top: 20px;",
          helpText(HTML("<strong>PRC Plot Download Options:</strong>")),
          downloadButton("downloadPlotPNG", "Download PRC Plot as PNG"),
          downloadButton("downloadPlotPDF", "Download PRC Plot and Metadata as PDF"),
          downloadButton("downloadCSV", "Download Variants Used as CSV"),
          div(
            style = "display: inline-flex; align-items: center;",
            downloadButton("downloadCSV_VUS", HTML("Download CSV with VUS")),
            actionLink("helpButton_VUS", label = NULL, icon = icon("question-circle"), style = "margin-left: 5px;")
          )
        )
      })
      
      # Download PNG logic
      output$downloadPlotPNG <- downloadHandler(
        filename = function() {
          paste("PRC_plot_", Sys.Date(), ".png", sep = "")
        },
        content = function(file) {
          plot_info <- rv$plot_data()
          if (!is.null(plot_info)) {
            num_scores <- length(plot_info$selected_scores)
            
            available_colors <- grDevices::colors()[!grepl("white|ivory|seashell|snow|honeydew|azure|aliceblue|mintcream|ghostwhite", grDevices::colors(), ignore.case = TRUE)]
            
            col_styles <- if (num_scores <= 3) {
              c("purple", "cadetblue2", "orange")[1:num_scores]
            } else {
              sample(available_colors, num_scores)
            }
            
            lty_styles <- rep(c("solid", "dashed"), length.out = num_scores)
            
            png(file, width = 6, height = 6, units = "in", res = 72)
            draw.prc(
              plot_info$yrobj,
              lty = lty_styles,
              col = col_styles,
              lwd = 2,
              balanced = TRUE,
              main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", "))
            )
            abline(h = 90, lty = "dashed")
            legend(
              "left",
              legend = c(
                paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org),
                paste("# of Benign and Likely Benign:", plot_info$B_org)
              ),
              pch = 15,
              bty = "n"
            )
            dev.off()
          }
        }
      )
      
      # Download PDF with metadata logic
      output$downloadPlotPDF <- downloadHandler(
        filename = function() {
          paste("PRC_Report_", Sys.Date(), ".pdf", sep = "")
        },
        content = function(file) {
          plot_info <- rv$plot_data()
          if (!is.null(plot_info)) {
            rmarkdown::render(input = "report_template_custom.Rmd",
                              output_file = file,
                              params = list(
                                gene_s = plot_info$gene_s,
                                selected_scores = plot_info$selected_scores,
                                common_variant_filter = plot_info$common_variant_filter,
                                B_org = plot_info$B_org,
                                P_org = plot_info$P_org,
                                prcfiltered = plot_info$prcfiltered
                              ),
                              envir = new.env(parent = globalenv()))
          }
        }
      )
      
      output$downloadCSV <- downloadHandler(
        filename = function() {
          paste("PRC_data_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          plot_info <- rv$plot_data()
          if (!is.null(plot_info) && !is.null(plot_info$prcfiltered) && nrow(plot_info$prcfiltered) > 0) {
            prcfiltered <- plot_info$prcfiltered
            prcfiltered <- prcfiltered %>%
              mutate(clinvar = ifelse(clinvar == TRUE, "P/LP", "B/LB"))
            write.csv(prcfiltered, file, row.names = FALSE)
          }
        }
      )
      
      output$downloadCSV_VUS <- downloadHandler(
        filename = function() {
          paste("PRC_data_VUS", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          req(rv$variant_data_df)
          write.csv(rv$variant_data_df, file, row.names = FALSE)
        }
      )
      
      # VUS explanation help button
      observeEvent(input$helpButton_VUS, {
        showModal(modalDialog(
          title = "What does this data include?",
          HTML("In addition to P/LP or B/LB variants, this download includes VUS and Conflicting variants. 
          VUS stands for Variants of Uncertain Significance, and Conflicting means conflicting interpretations of pathogenicity. 
          These variants are not included in the PRC generation as they are not classified as Pathogenic or Benign. 
          However, by downloading this CSV, you can override the VUS or Conflicting annotations and reupload it as a custom reference set in Advanced Mode."),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      }, ignoreInit = TRUE)
      
      output$download_buttons <- renderUI({
        req(rv$plot_data())  # Only show if plot_data has been generated
        
        tagList(
          actionButton("helpButton", "Plot Explanation"),
          helpText(HTML("<span style='color:black;'><strong>Download Options: </strong></span>")),
          downloadButton("downloadAllPlotsZIP", "Download All Plots (ZIP)"),
          downloadButton("downloadAllCSVsZIP", "Download All CSVs (ZIP)")
        )
      })
      
      ### New Download Handlers for ZIP Files ###
      ### ---------------------- HELPER FUNCTION TO CREATE ZIP FILES ---------------------- ###
      generate_zip <- function(zipfile, file_list, tmpdir) {
        old_wd <- setwd(tmpdir)
        on.exit(setwd(old_wd), add = TRUE)
        zip::zipr(zipfile = zipfile, files = file_list)
      }
      
      ### ---------------------- DOWNLOAD ALL PLOTS AS ZIP -------------------------------- ###
      output$downloadAllPlotsZIP <- downloadHandler(
        filename = function() {
          paste0("AllPlots_", Sys.Date(), ".zip")
        },
        content = function(zipfile) {
          req(rv$plot_data())
          
          tmpdir <- file.path(tempdir(), paste0("plotszip_", Sys.Date()))
          if (!dir.exists(tmpdir)) dir.create(tmpdir)
          
          plot_info <- rv$plot_data()
          file_list <- c()
          
          # 1) Save PRC Plot
          prc_png <- file.path(tmpdir, paste0("PRC_plot_", plot_info$gene_s, ".png"))
          png(prc_png, width = 6, height = 6, units = "in", res = 72)
          draw.prc(plot_info$yrobj, lty = plot_info$lty_styles, col = plot_info$col_styles, 
                   lwd = 2, balanced = TRUE, 
                   main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", ")))
          abline(h = 90, lty = "dashed")
          legend("left", legend = c(paste("# of Pathogenic:", plot_info$P_org), paste("# of Benign:", plot_info$B_org)), 
                 pch = 15, bty = "n")
          dev.off()
          
          if (file.exists(prc_png)) file_list <- c(file_list, prc_png)
          
          # 2) Save LLR & Bar Plots
          llr_scores <- paste0("VEP_", plot_info$selected_scores)
          full_filtered <- prcdata()
          
          # Define category colors (FIX: Ensure it exists in function scope)
          category_colors <- c(
            "benign_strong"   = "dodgerblue",
            "benign_support"  = "lightblue",
            "none"            = "grey",
            "patho_support"   = "#FFC0CB",
            "patho_moderate"  = "#FF9999",
            "patho_strong"    = "#FF6666",
            "patho_vstrong"   = "red"
          )
          
          for (score in llr_scores) {
            posScores <- na.omit(full_filtered[full_filtered$clinvar == "P/LP", score])
            negScores <- na.omit(full_filtered[full_filtered$clinvar == "B/LB", score])
            
            if (length(posScores) > 0 && length(negScores) > 0) {
              llrObj <- buildLLR.kernel(posScores, negScores)
              
              # Assign LLR values
              full_filtered$llr <- llrObj$llr(full_filtered[[score]])
              
              # ✅ FIX: Generate the category column
              breaks <- c(-Inf, -1.27, -0.32, 0.32, 0.64, 1.27, 2.54, Inf)
              labels <- c("benign_strong", "benign_support", "none",
                          "patho_support", "patho_moderate", "patho_strong", "patho_vstrong")
              
              full_filtered$category <- cut(full_filtered$llr, breaks, labels, include.lowest = TRUE)
              
              # Generate LLR density plot
              llr_png_path <- file.path(tmpdir, paste0("LLR_", score, ".png"))
              png(llr_png_path, width = 500, height = 600, res = 72)
              drawDensityLLR_fixedRange(full_filtered[[score]], llrObj$llr, llrObj$posDens, llrObj$negDens, posScores, negScores)
              dev.off()
              if (file.exists(llr_png_path)) file_list <- c(file_list, llr_png_path)
              
              # Generate Bar plot
              bar_png_path <- file.path(tmpdir, paste0("Bar_", score, ".png"))
              png(bar_png_path, width = 300, height = 500, res = 72)
              
              full_filtered$clinvar <- factor(full_filtered$clinvar,
                                              levels = c("P/LP", "B/LB", "VUS", "Conflicting"))
              
              bar_plot <- ggplot(full_filtered, aes(x = clinvar, fill = category)) +
                geom_bar(position = position_stack(reverse = TRUE)) +
                scale_fill_manual(values = category_colors) +
                theme_minimal() +
                labs(x = "ClinVar Category", y = "Variant Count",
                     fill = "LLR-Derived Evidence Category") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
              
              print(bar_plot)
              dev.off()
              
              if (file.exists(bar_png_path)) file_list <- c(file_list, bar_png_path)
            }
          }
          # 3) Ensure at least one file exists before zipping
          if (length(file_list) == 0) {
            showModal(modalDialog(
              title = "Error", "No plots were generated. Check inputs and try again.", easyClose = TRUE
            ))
            return(NULL)
          }
          
          # 4) Zip all files
          generate_zip(zipfile, file_list, tmpdir)
        }
      )
      
      ### ---------------------- DOWNLOAD ALL CSVS AS ZIP -------------------------------- ###
      output$downloadAllCSVsZIP <- downloadHandler(
        filename = function() {
          paste0("AllCSVs_", Sys.Date(), ".zip")
        },
        content = function(zipfile) {
          req(rv$plot_data())
          
          tmpdir <- file.path(tempdir(), paste0("csvzip_", Sys.Date()))
          if (!dir.exists(tmpdir)) dir.create(tmpdir)
          
          plot_info <- rv$plot_data()
          file_list <- c()
          
          # 1) Save PRC Data CSV
          prc_csv <- file.path(tmpdir, "PRC_data.csv")
          selected_df <- plot_info$prcfiltered
          
          if (!is.null(selected_df) && nrow(selected_df) > 0) {
            write.csv(selected_df, prc_csv, row.names = FALSE)
            if (file.exists(prc_csv)) file_list <- c(file_list, prc_csv)
          } else {
            print("WARNING: PRC data CSV is empty, skipping...")
          }
          
          # 2) Save VUS Data CSV 
          all_df <- prcdata() 
          
          vus_csv <- file.path(tmpdir, "AllData_includingVUS.csv")
          if (!is.null(all_df) && nrow(all_df) > 0) {
            write.csv(all_df, vus_csv, row.names = FALSE)
            if (file.exists(vus_csv)) file_list <- c(file_list, vus_csv)
          } else {
            print("WARNING: VUS data CSV is empty, skipping...")
          }
          
          # 3) Save LLR CSV
          llr_scores <- paste0("VEP_", plot_info$selected_scores)
          full_filtered <- prcdata()
          
          # Ensure full_filtered is not NULL
          if (!is.null(full_filtered) && nrow(full_filtered) > 0) {
            
            all_llr_df <- full_filtered  # Make a copy - NOTE: change name
            
            for (score in llr_scores) {
              # Ensure the score exists in the dataframe
              if (score %in% names(all_llr_df) && score %in% names(all_llr_df)) {
                posScores <- na.omit(setNames(all_llr_df[all_llr_df$clinvar == "P/LP", score],
                                              all_llr_df[all_llr_df$clinvar == "P/LP", "base__achange"]))
                negScores <- na.omit(setNames(all_llr_df[all_llr_df$clinvar == "B/LB", score],
                                              all_llr_df[all_llr_df$clinvar == "B/LB", "base__achange"]))
                
                # If at least one P/LP and one B/LB to build LLR:
                if (length(posScores) > 0 && length(negScores) > 0) {
                  llrObj <- buildLLR.kernel(posScores, negScores)
                  
                  # Create the LLR column named PREDICTOR_llr
                  llr_colname <- paste0(score, "_llr")
                  all_llr_df[[llr_colname]] <- llrObj$llr(all_llr_df[[score]])
                  
                  # Create the category column named PREDICTOR_category
                  cat_colname <- paste0(score, "_category")
                  breaks <- c(-Inf, -1.27, -0.32, 0.32, 0.64, 1.27, 2.54, Inf)
                  labels <- c("benign_strong", "benign_support", "none",
                              "patho_support", "patho_moderate", "patho_strong", "patho_vstrong")
                  
                  all_llr_df[[cat_colname]] <- cut(all_llr_df[[llr_colname]], breaks = breaks,
                                                   labels = labels, include.lowest = TRUE)
                }
              }
            }
            
            # Write out LLR CSV if data exists
            if (ncol(all_llr_df) > ncol(full_filtered)) {  # Check if new columns were added
              llr_csv <- file.path(tmpdir, "All_LLR_data.csv")
              write.csv(all_llr_df, llr_csv, row.names = FALSE)
              if (file.exists(llr_csv)) file_list <- c(file_list, llr_csv)
            }
          }
          
          # 4) Ensure at least one file exists before zipping
          if (length(file_list) == 0) {
            showModal(modalDialog(
              title = "Error",
              "No CSVs were generated. Check inputs and try again.",
              easyClose = TRUE
            ))
            return(NULL)
          }
          
          # 5) Zip all CSVs
          generate_zip(zipfile, file_list, tmpdir)
        }
      )
      
      
      
      
      # ObserveEvent for helpButton
      observeEvent(input$helpButton, {
        showModal(modalDialog(
          title = "Plot Explanation",
          HTML("<p><b>What is balanced precision?</b></p>
      <p>Balanced precision is useful in situations where the class distribution is imbalanced. In other words, it is the precision that would have been expected had the proportion of positive examples been balanced (equal to 50%).</p>"),
          p("Definition ", 
            a("here", href = "https://doi.org/10.1016/j.ajhg.2021.08.012", target = "_blank")),
          HTML("<p><b>What are R90BP and AUBPRC?</b></p>
      <p><b>R90BP:</b> The recall achieved at a stringent (90%) balanced precision threshold </p>
      <p><b>AUBPRC:</b> The area under the BPRC curve</p>"),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      })
    }
  })
}
