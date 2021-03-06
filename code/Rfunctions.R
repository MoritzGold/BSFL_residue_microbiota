# Function for plot

top_taxa_abundance <- function(physeq, numberOfTaxa = 9, raw = FALSE) {
  ## Args:
  ## - physeq: phyloseq class object
  ## - raw: (Required). Defaults to 'FALSE', logical. Should abundances be transformed to
  ##        frequencies before sorting otus from most to less abundant. 
  ## - numberOfTaxa: number of (top) taxa to keep 
  ##
  ## Returns:
  ## - x: Data frame with numberOfTaxa rows (one per final otu) and column 'Abundance' standing
  ##      for  total counts (if raw = TRUE) or average frequency
  ##      in samples of physeq. Rownames of x are either otus' names (if taxaRank = NULL) or names
  ##      corresponding to taxonomic rank 'TaxaRank'. All NA ranks are assigned to 'Unassigned'.
  stopifnot(!is.null(tax_table(physeq, FALSE)))
  otutab <- otu_table(physeq)
  if ( !taxa_are_rows(otutab) ) {otutab = t(otutab)}
  otutab <- as(otutab, "matrix")
  if (raw) {
    Abundance <- rowSums(otutab)
  } else {
    otutab <- apply(otutab, 2, function(x) x / sum(x))
    Abundance <- rowMeans(otutab)
  }
  ## Get top taxa
  mdf <- data.frame(OTU = names(Abundance), Abundance = Abundance)
  ## Add taxonomic information
  tax <- as(tax_table(physeq), "matrix")
  tax <- data.frame(OTU = rownames(tax), tax)
  mdf <- merge(mdf, tax, by.x = "OTU")
  ## Keep only numberOfTaxa top taxa
  topTaxa <- names(sort(Abundance, decreasing = TRUE))[1:numberOfTaxa]
  mdf <- mdf[ match(topTaxa, mdf$OTU), ]
  mdf$OTU <- factor(mdf$OTU, levels = unique(mdf$OTU))
  return(mdf)
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  ## Multiple plot function
  ##
  ## ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  ## - cols:   Number of columns in layout
  ## - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  ##
  ## If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  ## then plot 1 will go in the upper left, 2 will go in the upper right, and
  ## 3 will go all the way across the bottom.
  ##
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


require(parallel)
options(mc.cores= 2)

## Rarefaction curve, ggplot style
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed. 
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data 
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}


phylodiv <- function(physeq) {
  ## Args:
  ## - physeq: phyloseq class object, from which phylogeny and abundance data are extracted
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  phy <- phy_tree(physeq)
  
  ## Construct incidence matrix of the tree
  incidence <- incidenceMatrix(phy)
  
  ## Order incidence matrix according to community tables
  incidence <- incidence[colnames(x), ]
  
  ## Create community phylogeny matrix by multiplying (community x edge matrix)
  ## where cpm_{ij} gives the abundance of OTUs originating from branch j in community i. 
  cpm <- x %*% incidence
  ## Convert to incidence matrix (0/1) and multiply by edge length to obtain PD per community.
  cpm[cpm > 0] <- 1
  pd <-  cpm %*% phy$edge.length
  
  ## Add sample data information
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$pd <- as.vector(pd)
    pd <- sdf
  }
  
  return (pd)
}


ggpdrare <- function(physeq, step = 10, label = NULL, color = NULL,
                     log = TRUE,
                     replace = FALSE, se = TRUE, plot = TRUE, parallel = FALSE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - log:   (Otional). Default 'TRUE'. Logical value. Should sample size
  ##          be represented using a log10 scale?
  ## - replace: If TRUE, population are treated as of infinite size, with probabilities of occurence
  ##            of a taxa computed from the (finite size) community data
  ## - se  : Logical, should standard error be computed in addition to expected pd
  ## - plot:  Logical, should the graphic be plotted. 
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  phy <- phy_tree(physeq)
  
  ## Construct incidence matrix of the tree
  incidence <- incidenceMatrix(phy)
  
  nedges <- nrow(phy$edge)
  ## Order incidence matrix according to community tables and create
  ## community phylogenetic matrix
  incidence <- incidence[colnames(x), ]
  cpm <- x %*% incidence
  if (se) {
    cat("Preliminary computations for se, may take some time\n")
    cpm.var <- array(NA, dim = c(nedges, nedges, nrow(x)))
    cpm.var.fun <- function(i) {
      union.clade <- incidence[, i] + incidence
      union.clade[union.clade > 0] <- 1
      union.clade <- t(x %*% union.clade)
      ## union.clade[s, j] is the number of individuals in subtrees
      ## generated by cutting branches i (from outer loop) and j
      ## in sample s
      return(union.clade) 
    }
    for (i in seq_len(nedges)) {
      if (i %% 100 == 0) {
        cat(paste("Cutting edge", i, "out of", nedges), sep = "\n")
      }
      cpm.var[i, , ] <- cpm.var.fun(i)
    }
    ## Deprecated code, need to work on a better parallel version
    ## if (parallel) {
    ##     cpm.var <- mclapply(seq_len(nedges), cpm.var.fun, mc.preschedule = TRUE)
    ## } else {
    ##     cpm.var <- lapply(seq_len(nedges), cpm.var.fun)
    ## }
    ## cpm.var <- do.call(rbind, cpm.var)
    ## dim(cpm.var) <- c(nedges, nedges, nrow(x))
    dimnames(cpm.var) <- list(phy$edge[, 2], phy$edge[, 2], rownames(x))
  }
  
  ## Compute overall Phylogenetic Diversity
  pd <-  (0 + (cpm > 0) ) %*% phy$edge.length
  
  
  ## Transform community matrices to frequency data
  tot <- rowSums(x)
  nr <- nrow(x)
  ## Rarefy phylogenetic diversity for one sample (i)
  pdrare <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    ## Simplify matrices and tree to remove unnecessary computations. 
    edges.to.keep <- cpm[i, ] > 0
    branch.lengths <- phy$edge.length[edges.to.keep]
    cpm.i <- cpm[i, edges.to.keep]
    if (se) {
      cpm.var.i <- cpm.var[ edges.to.keep, edges.to.keep, i]
    }
    ## sequence of sample sizes
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    ## Mean and variance of pd for different sample sizes
    ## Start with mean
    if (replace) {
      ## Expected cpm
      cpm.rare <- 1 - t(outer((1 - cpm.i/tot[i]), n, "^"))
    } else {
      ## use lchoose instead of choose for numeric stability
      cpm.rare <-  outer(tot[i] - cpm.i, n, lchoose)
      cpm.rare <- sweep(cpm.rare, 2, lchoose(tot[i], n), FUN = "-")
      cpm.rare <- t(1 - exp(cpm.rare))
    }
    pd.rare <- as.vector(cpm.rare %*% branch.lengths)
    ## Continue with se, if necessary
    if (se) {
      cat(paste("Compute se for sample", rownames(x)[i], ", may take some time"), sep = "\n")
      ## Variance of cpm, computed via a loop to optimize memory use
      centering <-  (1 - cpm.rare) %*% branch.lengths
      pd.rare.var <- rep(NA, length(n))
      for (index in seq_along(n)) {
        size <- n[index]
        if (replace) {
          cpm.var.rare <- (1 - cpm.var.i/tot[i])^size
        } else {
          ## use lchoose instead of choose for numeric stability
          cpm.var.rare <- lchoose(tot[i] - cpm.var.i, size) - lchoose(tot[i], size)
          cpm.var.rare <- exp(cpm.var.rare)
        }
        pd.var <- t(branch.lengths) %*% cpm.var.rare %*% branch.lengths - centering[index]^2
        pd.rare.var[index] <- pd.var
      }
      pd.rare <- data.frame(pd.rare = pd.rare, se = sqrt(pd.rare.var))
    }
    return(data.frame(pd.rare, Size = n, Sample = rownames(x)[i]))
  }
  
  if (parallel) {
    out <- mclapply(seq_len(nr), pdrare, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), pdrare)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data 
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y =  pd, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data, aes_string(x = "Size", y = "pd.rare", group = "Sample", color = color))
  if (log) {
    p <- p + scale_x_log10()
  }
  p <- p + labs(x = "Sample Size (# reads)", y = "Phylogenetic Diversity")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) {
    p <- p + geom_ribbon(aes_string(ymin = "pd.rare - se", ymax = "pd.rare + se",
                                    color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}