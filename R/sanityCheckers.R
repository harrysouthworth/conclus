# Perform sanity checks on the connectivity matrix
checkM <- function(x){
  if (!isSymmetric(x)){
    stop("M is not symmetric!")
  } else if (any(is.na(x))){
    stop("M contains NAs!")
  } else if (max(x) != 1){
    stop(paste("maximum of M is", max(x), "should be 1"))
  } else if (min(x) < 0){
    stop("minimum of M is < 0!")
  } else if (any(diag(x) == 0)){
    stop("at least one diagonal of M is 0 - not enough subsamples were run")
  } else if (any(diag(x) != 1)){
    stop("diagonals of M should all be 1")
  }
}

checkClusterFun <- function(fun){
  arguments <- names(formals(fun))
  if (length(arguments) != 2 || !all(arguments == c("x", "k"))){
    stop("cluster should be a function with precisely 2 arguments, 'x' and 'k'")
  }

  test <- fun(daisy(pluton), k=2)

  if (length(class(test)) != 1 | !("integer" %in% class(test)) | length(test) != nrow(pluton)){
    stop("cluster should return an integer vector with length equal to the number of rows in the data")
  }
  invisible()
}

checkDiss <- function(x){
  if (!("dist" %in% class(x))){
    stop("diss should be a distance or dissimilarity matrix")
  }
  x <- as.matrix(x)
  if (!isSymmetric(x)){
    stop("diss must be a symmetric after coercion to matrix")
  }
  if (length(rownames(x)) != nrow(x)){
    stop("diss should have unique rownames after coercion to matarix")
  }
}

checkK <- function(x){
  if (is.null(x)){
    stop("You must specify K")
  } else if (x < 2){
    stop("K must be at least 2")
  }
}

checkSubsample <- function(x){
  if (!is.numeric(x) | x <= 0 | x > 1){
    stop("s should be a number in (0, 1)")
  }
}
