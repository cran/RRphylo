#' @title Producing simulated phenotypes with trends
#'
#' @description The function \code{setBM} is wrapper around \pkg{phytools} \code{fastBM} function, which generates BM simulated phenotypes with or without a trend.
#' @usage setBM(tree, nY = 1, s2 = 1, a = 0, type = c("", "brown","trend", "drift"),
#' tr = 10, t.shift = 0.5, trend.type = c("linear", "stepwise"))
#' @param tree a phylogenetic tree.
#' @param nY the number of phenotypes to simulate.
#' @param s2 value of the Browian rate to use in the simulations.
#' @param a the phenotype at the tree root.
#' @param type the type of phenotype to simulate. With the option \code{"brown"} the phenotype will have no trend in the mean or in the rate of evolution (actually the residuals of the phenotype versus time regression, established via Breusch-Pagan test). A variation in the phenotypic mean over time (a phenotypic trend) is obtained by selection the option \code{"drift"}. A trend in the rate of evolution should produce an increased variance in the residuals over time, which can be inspected under the option \code{"trend"}, by means of Breusch-Pagan test for residuals heteroscedasticity.
#' @param trend.type two kinds of heteroscedastic residuals are generated under the \code{"trend"} type. The option \code{"linear"} produces a linear increase in heteroscedasticity, whereas the \code{"stepwise"} option produces an increase after a specified point in time.
#' @param t.shift the relative time distance from the tree root where the stepwise increase in the rate of evolution is indicated to apply.
#' @param tr the intensity of the trend in homoscedasticity with the \code{"stepwise"} option is controlled by the \code{'tr'} argument. The scalar \code{'tr'} is the multiplier of the branches extending after the shift point as indicated by \code{'t.shift'}.
#' @details Note that \code{setBM} differs from \code{fastBM} in that the produced phenotypes are checked for the existence of a temporal trend either in the phenotypes themselves, or in the residuals of the phenotype versus age regression. The user may specify whether she wants trendless data (option \code{"brown"}), phenotypes trending in time (option \code{"drift"}), or phenotypes whose variance increases over time, consistently with the possible existence of a trend in the rate of evolution (option \code{"trend"}). In the latter case, the user may indicate the intensity of the change in variance, and whether it should occur after a given proportion of the tree height (hence a given point back in time). Trees in \code{setBM} are treated as non ultrametric. If an ultrametric tree is fed to the function, \code{setBM} alters slightly the leaf lengths multiplying randomly half of the leaves by 1 * 10-3.
#' @return Either an object of class \code{'array'} containing a single phenotype or an object of class \code{'matrix'} of \emph{n} phenotypes as columns, where \emph{n} is indicated as \code{nY} = \emph{n}.
#' @export
#' @importFrom phytools fastBM
#' @importFrom ape vcv
#' @importFrom lmtest bptest
#' @importFrom stats coefficients
#' @author Pasquale Raia, Silvia Castiglione, Carmela Serio, Alessandro Mondanaro, Marina Melchionna, Mirko Di Febbraro, Antonio Profico, Francesco Carotenuto
#' @examples
#' data("DataOrnithodirans")
#' DataOrnithodirans$treedino->treedino
#'
#' setBM(tree=treedino, nY= 1, type="brown")
#' setBM(tree=treedino, nY= 1, type="drift")
#' setBM(tree=treedino, nY= 1, type="trend", trend.type="linear")



setBM<-function (tree, nY = 1, s2 = 1, a = 0, type = c("", "brown",
                                                       "trend", "drift"), tr = 10, t.shift = 0.5,
                 trend.type = c("linear",
                                "stepwise"))
{
  #require(lmtest)
  if (type == "")
    stop("argument 'type' must be defined")
  switch(type, brown = {
    i = 1
    times <- diag(vcv(tree))
    if(var(times)<0.0001) {
      makeFossil(tree,ex=1.001)->tree
      times <- diag(vcv(tree))
    }
    yy <- list()
    while (length(yy) < nY) {
      y <- fastBM(tree, sig2 = s2, a = a)
      res <- bptest(y ~ times)[4]
      if (res > 0.1 & coefficients(summary(lm(y ~ times)))[8] >
          0.1) {
        yy[[i]] <- y
      }
      yy <- Filter(Negate(is.null), yy)
      i = i + 1
    }
    yy <- do.call(rbind, yy)
    yy <- t(yy)
    if (nY == 1) {
      nam <- rownames(yy)
      yy <- array(yy)
      names(yy) <- nam
    }
  }, trend = {
    i = 1
    times <- diag(vcv(tree))
    if(var(times)<0.0001) {
      makeFossil(tree,ex=1.001)->tree
      times <- diag(vcv(tree))
    }
    yy <- list()
    while (length(yy) < nY) {
      y <- fastBM(tree, sig2 = s2, a = a)
      res <- bptest(y ~ times)[4]
      if (res > 0.05 & coefficients(summary(lm(y ~ times)))[8] >
          0.05) {
        match.arg(trend.type)
        if (trend.type == "linear") {
          y <- (diag(vcv(tree))^2)/min(diag(vcv(tree))) *
            y
        } else {
          y[match(names(which(times > t.shift * max(nodeHeights(tree)))),
                  names(y))] <- y[match(names(which(times >
                                                      t.shift * max(nodeHeights(tree)))), names(y))] *
            tr
        }
        yy[[i]] <- y
      }
      yy <- Filter(Negate(is.null), yy)
      i = i + 1
    }
    yy <- do.call(rbind, yy)
    yy <- t(yy)
    if (nY == 1) {
      nam <- rownames(yy)
      yy <- array(yy)
      names(yy) <- nam
    }
  }, drift = {
    i = 1
    times <- diag(vcv(tree))
    if(var(times)<0.0001) {
      makeFossil(tree,ex=1.001)->tree
      times <- diag(vcv(tree))
    }
    yy <- list()
    while (length(yy) < nY) {
      y <- fastBM(tree, sig2 = s2, a = a)
      res <- bptest(y ~ times)[4]
      if (res > 0.05 & coefficients(summary(lm(y ~ times)))[8] <
          0.001) {
        yy[[i]] <- y
      }
      yy <- Filter(Negate(is.null), yy)
      i = i + 1
    }
    yy <- do.call(rbind, yy)
    yy <- t(yy)
    if (nY == 1) {
      nam <- rownames(yy)
      yy <- array(yy)
      names(yy) <- nam
    }
  })
  return(yy)
}
