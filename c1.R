# Copyright 2014 Andrew A. Hill
#
#  This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

require(flowCore)
require(flowViz)

# ceiling on predicted class 1 and 2 events, based upon prior probabilites
# of classes 1 & 2 in entire training set (n=202 FCS)
ACCEPT.FRAC.MAX <- (10168 + 24303) / 56221534 * 2;

myAshz <- function(x, a=1, b=1, c=0) {
  as.vector(scale(asinh(a+b*x)+c))
}

ashz.transform <- function(transformId, a=1, b=1, c=0) {
  t = new("transform", .Data = function(x) myAshz(x, a, b, c))
  t@transformationId <- transformId
  t
}

fcsFileForNumber <- function(i) {
  file.path(fcs.root, paste(i,"fcs",sep="."))
}

flowFrameForNumber <- function(i, ...) {
  ff <- read.FCS(fcsFileForNumber(i), ...)
  tl <- transformList( colnames(ff), myAshz)
  ff <- transform(ff, tl)
}

fixRanges <- function(ff) {
  min.max <- t(summary(ff)[c("Min.", "Max."),])
  pp <- parameters(ff)
  pdp <- pData(pp)
  pdp[,c("minRange", "maxRange")] <- min.max
  pData(pp) <- pdp
  parameters(ff) <- pp
  ff
}

plotFrameEventClasses <- function(ff, cls, mychans, ...) {
  ff <- ff[,mychans]
  fmat <- exprs(ff)
  cn <- colnames(fmat)
  lcn <- length(cn)
  for (j in 1:(lcn-1)) {
    for(k in (j+1):lcn) {
      cc <- 0
      plot(ff[cls==cc,c(j,k)], pch=".", ...) # xlim=range(fmat[,j]), ylim=range(fmat[,k]), ...)
      for (cc in 1:2) {
        lines(fmat[cls==cc,c(j,k)], col=cc+1,
              type="p", pch=cc)
      }
    }
  }
}


wrapLegend <- function(grps) {
  tt <- table(grps)
  legend(x="bottomright", legend=paste(names(tt), tt), pch=c(".",1,2), col=1:3)
}

wf1.core <- function(ff, grps=rep(0,nrow(ff)), PCT.MAX=8, pe.bwFac=5, plot=FALSE) {
  # ff = unlabelled flow frame
  # returns: predicted classes of events in ff (0, 1, or 2)
  
  events.orig <- 1:nrow(ff)
  
  # applyi SSC-FSC curv2filt (major population)
  f.sc <- curv2Filter("SSC", "FSC", bwFac=12, filterId="ssc.fsc.curv2")
  system.time(f.sc.result <- filter(ff, f.sc))
  ff.1 <- split(ff,f.sc.result)
  ff.1.main <- ff.1[["area 1"]]

  if (plot) {
    plotFrameEventClasses( ff, grps, c("FSC", "SSC"),
                        xlim=c(-3,3), ylim=c(-3,3), main=identifier(ff))
    wrapLegend(grps)
  }
  i.filt <- f.sc.result[["area 1"]]@subSet
  if (plot) {
    plotFrameEventClasses( ff[i.filt,], grps[i.filt], c("FSC", "SSC"),
                          xlim=c(-3,3), ylim=c(-3,3), main=identifier(ff))
    wrapLegend(grps[i.filt])
  }
   
  # apply filtres.pe to exclude class 0
  ff.next <- fixRanges(ff[i.filt,])
  grps.next <- grps[i.filt]
  events.next <- events.orig[i.filt]
  f.pe.fitc <- curv1Filter("PE", filterId="pe.fitc.curv1", bwFac=pe.bwFac)
  system.time( f.pe.fitc.result <- filter(ff.next, f.pe.fitc) )

  # select high PE region 
  i.region.select <- filt.select.by.size(f.pe.fitc.result, PCT.MAX)
  bound.list <- f.pe.fitc.result@filterDetails$pe.fitc.curv1$boundaries
  pe.cutoff <- max( bound.list[[i.region.select]][2], quantile(exprs(ff.next)[,"PE"],1-ACCEPT.FRAC.MAX))

  i.next2 <- exprs(ff.next)[,"PE"]>pe.cutoff
  ff.next2 <- fixRanges(ff.next[i.next2,])
  grps.next2 <- grps.next[i.next2]
  events.next2 <- events.next[i.next2]
  print(table(grps.next2))
  if (plot) {
    plotFrameEventClasses( ff.next, grps.next, c("FITC", "PE"),
                        xlim=c(-3,3), ylim=c(-3,3),
                        main=identifier(ff.next))
    abline(h=unlist(bound.list), col=1, lty=2)
    abline(h=pe.cutoff, col=2)
    wrapLegend(grps.next2)
  }
  
  # apply FITC filter to separate group 1 from group 2
  # f.fitc <- curv1Filter("FITC", filterId="fitc.curv1", bwFac=10)
  fitc.cut <- 0
  f.fitc <- rectangleGate( filterId="fitc.rect", "FITC"=c(fitc.cut,Inf))
  system.time( f.fitc.result <- filter( ff.next2, f.fitc) )
  print(toTable(summary(f.fitc.result)))
  
  i.next4 <- exprs(ff.next2)[,"FITC"]>fitc.cut
  ff.next4 <- fixRanges(ff.next2[i.next4,])
  grps.next4 <- grps.next2[i.next4]
  events.next4 <- events.next2[i.next4]
  
  if (plot) {
    plotFrameEventClasses( ff.next2, grps.next2, c("FITC", "PE"),
                        xlim=c(-3,3), ylim=c(-3,3),
                        main=identifier(ff.next4), smooth=FALSE)
    abline(v=0)
    wrapLegend(grps.next4)
  }
  grps.pred <- rep(0,length(events.orig))
  grps.pred[events.next4] <- 1 # place >FITC in the predicted class 1
  grps.pred[setdiff(events.next2,events.next4)] <- 2 # remaining events >PE cut but <FITC cut in class 2
  grps.pred
}

wf1 <- function(i, PCT.MAX=8, pe.bwFac=5, plot=FALSE) {
  ff <- flowFrameForNumber(i)
  grps <- c()
  if (i>202) {
    # this is the validation set
    grps <-  factor(rep(0,nrow(ff),levels=c(0,1,2)))
  } else {
    # this is the training set
    grps <- fcsEventClassesForNumber(i)
  }
  grps.pred <- wf1.core( ff, grps=grps, PCT.MAX, pe.bwFac, plot=plot)
  pred.tab <- table(grps,factor(grps.pred,levels=0:2))
  ret <- list(pred.tab=pred.tab, grps.pred=grps.pred)
}

wf1.flowSet <- function(ii, PCT.MAX=8, pe.bwFac=1.2, out.root="test-pred-1", plot=FALSE) {
  ntry <- 1
  orig.root <- out.root
  while (file.exists(out.root) & ntry<100) {
    out.root <- paste(orig.root, sprintf("%02d",ntry),sep="")
  }
  if (ntry>=100) {
    stop("your output directory ", out.root, " already exists, will not overwrite")
  }
  dir.create(out.root)
  message("created output directory: " , out.root)
  
  res <- array(NA, dim=c(3,3,length(ii)),
               dimnames=list(truth=as.character(0:2),
                 predicted=as.character(0:2), as.character(ii)))
  for (i in 1:length(ii)) {
    message("working on ", paste(ii[i],"fcs",sep="."), " (", i, " of ", length(ii), ")")
    ret <- wf1(ii[i], PCT.MAX, pe.bwFac, plot=plot)
    res[,,i] <- ret$pred.tab
    cat(as.numeric(ret$grps.pred), file=file.path(out.root, paste(ii[i],"csv",sep=".")), sep="\n")
  }
  res
}

wrap.final.wf1.flowSet <- function(PCT.MAX=8, pe.bwFac=1.2, out.root="validate-pred-final", plot=FALSE) {
  ii <- 203:405
  ret <- wf1.flowSet(ii, PCT.MAX, pe.bwFac, out.root, plot=FALSE)
}

filt.select.by.size <- function(filt.res, pct.cut=10) {
  tab <- toTable(summary(filt.res))
  print(tab)
  pct <- tab[-1,"percent"]
  # last region accounting for >= pct.cut% of events
  if (any(pct>=pct.cut)) {
    i.region.select <- max(which(pct>=pct.cut))
  } else {
    stop(paste("failed to find any region with pct >= ", pct.cut))
  }
  message(sprintf("selected pe region %d of %d", i.region.select, length(pct)))
  i.region.select
}
