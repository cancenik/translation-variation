"map" <- function(x, ...)
{
  UseMethod("map")
}


### Map data to a trained supersom network using any of the layers, or any
### subset of them. Should work for som, xyf and bdk as well...

### Aug 16, 2013
### Modified by Alan Boyle

### 21/2/07:
### - check mapping using subsets of maps: the whatmap argument should
### be checked when it is
###
### o) absent ---------------------------> OK
### i) numerical ------------------------> OK
### ii) a vector of string names; -------> OK
### both cases in
### A) correct usage and ----------------> OK
### B) incorrect usage. -----------------> OK

### FIXME: inefficient (but correct) if there are zeros in x$weights;
### whatmap is a more elegant way to exclude layers.

"map.kohonen" <- function(x, newdata, whatmap=NULL, weights,
                          scale.distances = (nmaps > 1), ...)
{
  codes <- x$codes
  if (is.matrix(codes)) codes <- list(codes)

  if (missing(newdata) & !is.null(x$data))
    newdata <- x$data
  if (is.matrix(newdata)) newdata <- list(newdata)
 
  if (is.null(whatmap) && !is.null(x$whatmap)) {
    whatmap <- x$whatmap
  } else {
    whatmap <- 1
  }
  
  nmaps <- length(whatmap)

  nd <- nrow(newdata[[1]])
  ncodes <- nrow(codes[[ whatmap[1] ]])

  distances <- matrix(0, nd*2, nmaps)
  ## first calculate all distance matrices: every column contains the
  ## distances of the objects to all unit in that particular layer.
  for (i in seq(along = whatmap)) {
    np <- ncol(codes[[ whatmap[i] ]])
    dt <- newdata[[ whatmap[i] ]]
    cd <- codes[[ whatmap[i] ]]
    
    distances[,i] <- .C("mapKohonen",
                        as.double(dt),
                        as.double(cd),
                        as.integer(ncodes),
                        as.integer(nd),
                        as.integer(np),
                        dists = double(nd * 2),
                        NAOK = TRUE,
                        PACKAGE = "kohonen2")$dists
  }
  
#APB - Here we have moved a large part of the proecssing to C and no
# longer need a nd x codes matrix but just a nd x 2 matrix

#  if (scale.distances)
#    distances <- sweep(distances, 2,
#                       apply(distances, 2, max, na.rm=TRUE), FUN="/")
  
  ## next determine overall closest
  if (nmaps > 1) {
#    if (missing(weights)) {
#      weights <- x$weights[whatmap] / sum(x$weights[whatmap])
#    } else {
#      if (abs(sum(weights)) < .Machine$double.eps) {
#        warning("sum of weights equals zero! Unscaled weights are used...")
#      } else {
#        weights <- weights / sum(weights)
#      }
#    }
#    overall.distances <- matrix(rowMeans(sweep(distances, 2, weights,
#                                               FUN="*"),
#                                         na.rm=TRUE),
#                                nd, ncodes, byrow=TRUE)

    ## The next bit seems simple but stumbles on NAs
    ##    overall.distances <- 
    ##      matrix(distances %*% weights, nd, ncodes, byrow=TRUE)
  } else {
    weights <- 1
		overall.distances <- matrix(distances, nd, 2, byrow=TRUE)
  }

  NArows <- which(apply(overall.distances, 1, function(x) all(is.na(x))))
  if (length(NArows) > 0) {
    classif <- mindists <- rep(NA, nrow(overall.distances))
    mindists[-NArows] <- overall.distances[-NArows,]
    classif[-NArows] <- as.integer(overall.distances[-NArows,])
  } else {
		mindists <- overall.distances[,1]
		classif <- as.integer(overall.distances[,2])
  }

	#fix the 0 based classif issue
	classif <- classif + 1
  
  list(unit.classif = classif,
       distances = mindists,
       whatmap = whatmap,
       weights = weights,
       scale.distances = scale.distances)
}

