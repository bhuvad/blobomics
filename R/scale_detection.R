checkMask <- function(mask, spe, plot = FALSE) {
  coords = SpatialExperiment::spatialCoords(spe)
  xmx = max(coords[, 1])
  ymx = max(coords[, 2])
  
  if (ncol(mask) != xmx | nrow(mask) != ymx) {
    stop("mask size does NOT match image size")
  }
  if (any(is.na(mask))) {
    stop("'mask' cannot contain missing values")
  }
  if (all(apply(mask, 2, is.logical))) {
    mask = apply(mask, 2, as.numeric)
  }
  if (any(!as.vector(mask) %in% 0:1)) {
    stop("'mask' should be a binary matrix")
  }

  return(mask)
}

checkFeatures <- function(spe, features) {
  if (is.character(features)) {
    # feature is a column name
    if (length(features) == 1) {
      cdata = SummarizedExperiment::colData(spe)
      if (!features %in% colnames(cdata)) {
        stop("'features' not found in colData(spe)")
      } else {
        features = cdata[, features]
      }
    } else {
      # feature contains rownames
      features = rownames(spe) %in% features
    }
  } else if (is.logical(features)) {
    if (!length(features) %in% c(1, nrow(spe))) {
      stop("'features' should have the same length as nrow(spe)")
    }
  } else {
    stop("'features' should be logical or character")
  }
  return(features)
}

getMaxOctaveScales <- function(octave = 1) {
  max(3 + 2^(octave - 1) * 1:4 * 6)
}

getOctaves <- function(spe) {
  # checks
  stopifnot(is(spe, "SpatialExperiment"))
  stopifnot(length(unique(spe$sample_id)) == 1)

  coords = SpatialExperiment::spatialCoords(spe)
  maxscale = min(max(coords[, 1]), max(coords[, 2]))

  o = 1
  while (getMaxOctaveScales(o) < maxscale) {
    o = o + 1
  }

  return(o - 1)
}

detectIP_intl <- function(spe, assay = "counts", mask = NULL, features, thresh, ipassay = "ipoints", octaves = getOctaves(spe)) {
  # checks
  stopifnot(length(unique(spe$sample_id)) == 1)
  stopifnot(octaves < getOctaves(spe) + 1)

  # check coordinates
  coords = SpatialExperiment::spatialCoords(spe)
  xmx = max(coords[, 1])
  ymx = max(coords[, 2])
  if (any(coords - apply(coords, 2, as.integer) > 0)) {
    stop("'spatialCoords' need to be integers for this function")
  }

  # convert assay id to name
  if (is.numeric(assay)) {
    # convert assay id to assay name
    assay = SummarizedExperiment::assayNames(spe)[assay]
  }
  ipassay = paste(ipassay, assay, sep = "_")

  # check mask
  if (is.null(mask)) {
    mask = rep_len(1, nrow(coords))
  } else {
    mask = checkMask(mask, spe, plot = FALSE)
    # transpose
    mask = mask[coords[, 2:1]]
  }

  # get interest points
  emat = t(as.matrix(SummarizedExperiment::assay(spe[features, ], assay)))
  df = .Call(`_spaceblobs_fastHessian`, emat, mask, coords[, 1], coords[, 2], octaves, thresh)
  df = df[df$x <=xmx & df$y <= ymx, ]

  # resolve row indices
  irow = df$gene
  df = df[, -5]
  ixmap = seq_len(nrow(spe))[features]
  irow = ixmap[irow]
  irow = factor(irow, levels = seq_len(nrow(spe)))

  # resolve column indices
  icol = matrix(NA_integer_, xmx, ymx)
  icol[coords] = seq_len(ncol(spe))
  icol = factor(icol[cbind(round(df$x), round(df$y))], seq_len(ncol(spe)))

  # convert to BumpyDataFrameMatrix and add assay
  mat = BumpyMatrix::splitAsBumpyMatrix(df, row = irow, column = icol, sparse = TRUE)
  rownames(mat) = rownames(spe)
  colnames(mat) = colnames(spe)
  SummarizedExperiment::assay(spe, ipassay) = mat

  return(spe)
}

#' Detect multi-scale interest points using blob detection
#'
#' Treating spatial omics datasets as multi-channel images, this function builds a Gaussian scale-space pyramid and uses it to infer the scale of each gene at each locus. It then identifies local maxima. Local maxima above a threshold are points of interest that can be used to estimate local scales using the [detectScale()] function. The function uses a fast approximation adapted from the SURF (speeded-up robust features) method.
#'
#' @param spe a SpatialExperiment object, containing gene/protein or any such omic measurement. NOTE that spatial coordinates need to be integers. This can be achieved by binning data into coarser bins. Binning will speed up computation without significantly affecting interest point detection and is therefore advised.
#' @param assay a character or numeric, specifying the assay to use (default is the first assay).
#' @param masks a named list of binary or logical matrices, that specifies the mask to use for each sample. The defauls (NULL) is to use mask out areas without any expression measurement. Otherwise, the list must be named using the unique [sample_id]'s and should cover all loci (i.e., nrow should be equal to the maximum y coordinate and ncol should equal the maximum x coordinate).
#' @param features a character or logical vector, specifying the features to be used. If character, these should be rownames of the SpatialExperiment object or the name of a column in colData(spe) containing logical values.
#' @param thresh a numeric, indicating the response cutoff for each interest point.
#' @param ipassay a character, specifying the prefix of the assay where detected interest points are stored.
#' @param octaves a numeric, stating the maximum number of octaves that should be computed (default NULL, max octave is automatically computed for each sample).
#' @param BPPARAM an optional [BiocParallelParam] instance determining the parallel back-end to be used during evaluation.
#'
#' @return a SpatialExperiment object containing an additional assay storing the interest points identified at each location. This assay is stored as a BumpyMatrix where each cell of the matrix (each location, gene/protein pair) stores a dataframe with 4 columns, scale, x, y, response. Each coordinate (x,y) represents the sub-pixel resolution of the interest point.
#' 
#' @examples
#' # ADD_EXAMPLES_HERE
#' 
#' @export 
detectInterestPoints <- function(spe, assay = assayNames(spe)[1], masks = NULL, features = TRUE, thresh = 0, ipassay = "ipoints", octaves = NULL, BPPARAM = BiocParallel::bpparam()) {
  stopifnot(is(spe, "SpatialExperiment"))
  stopifnot(ncol(SpatialExperiment::spatialCoords(spe)) == 2)
  stopifnot(thresh >= 0)
  stopifnot(is.null(octaves) | octaves > 0)
  stopifnot(is.character(ipassay))

  # check features
  features = checkFeatures(spe, features)
  if (!any(features)) {
    message("no features to analyse")
    return(spe)
  }
  
  # get samples
  all_samples = unique(spe$sample_id)

  # check masks
  if (!is.null(masks)) {
    stopifnot(is.list(masks))
    stopifnot(length(masks) == length(all_samples))
    stopifnot(all(all_samples %in% names(masks)))
  } else {
    masks = rep(list(NULL), length(all_samples))
    names(masks) = all_samples
  }

  # split by sample, analyse, merge
  spe = lapply(all_samples, \(sample_id) spe[, spe$sample_id == sample_id])
  spe = BiocParallel::bpmapply(\(x, mask) {
    # compute max octaves if not provided
    if (is.null(octaves))
      octaves = getOctaves(x)
  
    detectIP_intl(x, assay, mask, features, thresh, ipassay, octaves)
  }, spe, masks, BPPARAM = BPPARAM, SIMPLIFY = FALSE) |>
    do.call(what = SummarizedExperiment::cbind)
  
  return(spe)
}

#' Estimate local length-scales using interest points
#'
#' This function estimates the local length-scale for each location by combining scale information from interest points detected using each gene.
#'
#' @inheritParams detectInterestPoints
#' 
#' @return a numeric vector containing the estimated scales.
#' @examples
#' # ADD_EXAMPLES_HERE
#' 
#' @export
detectScale <- function(spe, ipassay = "ipoints", masks = NULL, smooth = TRUE, thresh = 0, BPPARAM = BiocParallel::bpparam()) {
  stopifnot(is(spe, "SpatialExperiment"))
  stopifnot(ncol(SpatialExperiment::spatialCoords(spe)) == 2)
  stopifnot(thresh >= 0)

  # extract required data
  expected_cols = c("scale", "x", "y", "response")
  all_samples = unique(spe$sample_id)

  # convert assay id to name
  if (is.numeric(ipassay)) {
    # convert assay id to assay name
    ipassay = SummarizedExperiment::assayNames(spe)[ipassay]
  } else {
    assay_names = SummarizedExperiment::assayNames(spe)
    if (!ipassay %in% assay_names) {
      assay_names = grep(paste0("^", ipassay), assay_names, value = TRUE)
      if (length(assay_names) == 1) {
        ipassay = assay_names
      } else {
        stop("'ipassay' not found in the 'spe' object")
      }
    }
  }
  stopifnot(is(SummarizedExperiment::assay(spe, ipassay), "BumpyMatrix"))
  
  # check masks
  if (!is.null(masks)) {
    stopifnot(is.list(masks))
    stopifnot(length(masks) == length(all_samples))
    stopifnot(all(all_samples %in% names(masks)))
  } else {
    masks = rep(list(NULL), length(all_samples))
    names(masks) = all_samples
  }

  # if colnames absent, create temp
  if (is.null(colnames(spe))) {
    colnames(spe) = paste0("X", seq_len(ncol(spe)))
  }
  cnames = colnames(spe)

  # estimate scale per sample
  spe = lapply(all_samples, \(sample_id) spe[, spe$sample_id == sample_id])
  spe_scales = BiocParallel::bpmapply(\(x, mask) {
    # extract required data
    coords = SpatialExperiment::spatialCoords(x)
    ip = SummarizedExperiment::assay(x, ipassay)
    
    # check mask
    if (is.null(mask)) {
      mask = rep_len(1, nrow(coords))
    } else {
      mask = checkMask(mask, x, plot = FALSE)
      # transpose
      mask = mask[coords[, 2:1]]
    }

    # check df
    if (!all(expected_cols %in% colnames(ip[1, 1][[1]]))) {
      stop("BumpyMatrix is not in the correct format (should be the output of 'detectInterestPoints()'")
    }

    keep = ip[, , "response"] > thresh
    if (sum(sum(keep)) == 0) {
      # if no interest points found
      sc = rep_len(NA_real_, sum(is_samp))
    } else {
      # filter and compute average scale per location
      sc = colSums(sum(ip[, , "scale"][keep]))
      sc = sc / colSums(as.matrix(lengths(ip[, , "scale"][keep])))
      sc = sc / 2
      if (smooth) {
        # diffuse scales and return
        df = data.frame("x" = coords[, 1], "y" = coords[, 2], "scale" = sc)
        df = df[!is.nan(sc), ]
        sc = .Call(`_spaceblobs_smoothScales`, df, mask, coords[, 1], coords[, 2])
        sc = sc[coords]
      } else {
        sc[is.nan(sc)] = 3 / 2
      }
    }
    sc[mask == 0] = NA_real_

    names(sc) = colnames(x)
    return(sc)
  }, spe, masks, BPPARAM = BPPARAM, SIMPLIFY = FALSE) |>
    unlist()
  
  # aggregate results
  spe_scales = spe_scales[cnames]
  names(spe_scales) = NULL
  
  return(spe_scales)
}
