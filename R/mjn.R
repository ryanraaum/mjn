options(stringsAsFactors=F)

library(methods)
library(igraph)
library(shape)

my.mode <- function(x, n=3) {
  x <- sort.default(x)
  y <- x[-1L] != x[-n]
  i <- c(which(y), n)
  lengths <- diff(c(0L, i))
  j <- which(lengths == max(lengths))
  x[i[j]]
}

my.unique <- function (x) { x[!duplicated.default(x)] }

mandist <- function (x)
{
    N <- nrow(x)
    d <- .C("R_distance", x = as.double(x), nr = N, nc = ncol(x), 
        d = double(N * (N - 1)/2), diag = as.integer(FALSE), 
        method = as.integer(3), p = 2.0, DUP = FALSE, 
        NAOK = TRUE, PACKAGE = "stats")$d
    attr(d, "Size") <- N
    class(d) <- "dist"
    return(d)
}

allunique.f <- function(x) { sum(!duplicated.default(x)) == 3 }

median.vectors <- function(data) {
  if (!nrow(data) == 3) stop("data matrix must have exactly 3 rows")
  varnames <- colnames(data)
  allunique <- apply(data, 2, allunique.f)
  allunique.cols <- names(which(allunique))
  majorities <- apply(data[,!allunique], 2, my.mode)
  majnames <- names(majorities)
  env <- new.env()
  for (n in majnames) { assign(n, majorities[n], envir=env) }
  for (n in allunique.cols) { assign(n, my.unique(data[,n]), envir=env) }
  return(as.matrix(expand.grid(as.list(env), KEEP.OUT.ATTRS=F)[,varnames])) # order gets rearranged in the process
}

# igraph indices are 0-based, so need to subtract 1
#  supposedly R-igraph is moving to 1-based at some point, 
#  so will need to watch this
i2v <- function(x) { x - 1}

test.all.edge.lengths <- function(x, g, value) {
  all(E(g, path=x)$weight < value) 
}

connected.f <- function(edge, g, delta, epsilon) {
  v.ids <- i2v(edge[1:2])
  if (edge.connectivity(g, v.ids[1], v.ids[2]) == 0) { return(FALSE) }
  if (are.connected(g, v.ids[1], v.ids[2])) { return(TRUE) }
  p <- get.shortest.paths(g, v.ids[1], v.ids[2])
  return(all(sapply(p, test.all.edge.lengths, g, delta-epsilon)))
}

add.nodes <- function(g, edges) {
  node.ids <- unlist(t(edges[,c('n1','n2')]))
  node.weights <- edges[,'d']
  add.edges(g, i2v(node.ids), attr=list(weight=node.weights))
}

delta.step.components <- function(edges, epsilon=0) {
  num.nodes <- sum(!duplicated.default(as.vector(edges[,c('n1','n2')])))
  g <- graph.empty(n=num.nodes, directed=F)
  deltas <- sort.default(my.unique(edges[,'d']))
  feasible <- rep(F, nrow(edges))
  i <- 1
  while (i) {
    delta <- deltas[i]
    already.connected <- apply(edges, 1, connected.f, g, delta, epsilon)
    newly.feasible <- edges[,'d'] <= delta & !feasible & !already.connected
    feasible <- feasible | newly.feasible
    if (sum(newly.feasible)) {
      new.edges <- edges[newly.feasible,,drop=F]
      g <- add.nodes(g, new.edges)
    }
    if (is.connected(g)) { i = FALSE } else { i = i + 1 }
  }
  if (epsilon > 0) {
    for (i in which(deltas > delta & deltas <= (delta + epsilon))) { 
      deltaplus <- deltas[i]
      already.connected <- apply(edges, 1, connected.f, g, deltaplus, epsilon)
      newly.feasible <- edges[,'d'] <= deltaplus & !feasible & !already.connected
      feasible <- feasible | newly.feasible
      if (sum(newly.feasible)) {
        new.edges <- edges[newly.feasible,,drop=F]
        g <- add.nodes(g, new.edges)
      }
    }
  }
  return(list(delta=delta, feasible=feasible, f.num=sum(feasible), g=g, epsilon=epsilon))
}

connection.cost <- function(triplet, mvec) {
  curr <- rbind(mvec, triplet)
  curr.d <- mandist(curr)
  sum(curr.d[1:3])
}

new.sequence.types <- function(data, edges, dsc, epsilon=0, lambda=NULL, prev.combos=NULL) {
  feasible.m <- edges[dsc$feasible,c('n1','n2')]
  combos <- combn(1:dsc$f.num, 2)
  if (length(prev.combos) > 0) {
    dups <- duplicated(cbind(prev.combos, combos), MARGIN=2)[-(1:ncol(prev.combos))]
    combos <- combos[,!dups,drop=F]
  }
  mvecs <- data[F,]
  c.costs <- c()
  if (ncol(combos) > 0) {
    for (i in 1:ncol(combos)) {
      curr.m <- feasible.m[combos[,i],]
      nodes <- my.unique(as.vector(curr.m))
      if (length(nodes) == 3) {
        triplet <- data[nodes,]
        curr.mvecs <- median.vectors(triplet)
        # remove vectors that already exist in the data
        dups <- duplicated(rbind(data,mvecs,curr.mvecs))[-(1:(nrow(data)+nrow(mvecs)))]
        curr.mvecs <- curr.mvecs[!dups,,drop=F]
        d <- mandist(rbind(triplet, curr.mvecs))
        nmvecs <- nrow(curr.mvecs)
        if (nmvecs > 0) {
          new.costs <- rowSums(as.matrix(d)[4:(nmvecs+3),1:3,drop=F])
          c.costs <- c(c.costs, new.costs)
          mvecs <- rbind(mvecs, curr.mvecs)
        }
      }
    }
  }
  keep <- F
  if (length(c.costs) > 0) { 
    if (is.null(lambda)) { lambda=min(c.costs) }
    if (length(c.costs)>0) {
      keep <- c.costs <= (lambda + epsilon)
    }
  }
  prev.combos <- cbind(prev.combos, combos)
  if (sum(keep)) {
    return(list(n=sum(keep), costs=c.costs[keep], median.vectors=mvecs[keep,], 
                lambda=lambda, prev.combos=prev.combos))
  } else {
    return(list(n=0, costs=c(), median.vectors=data[F,], lambda=lambda))
  }
}

identify.obsoletes <- function(g, indices) {
  todrop <- c()
  for (i in indices) {
    vindex <- i - 1
    if (length(V(g)[nei(vindex)]) <= 2) {
      todrop <- c(todrop, i)
    }
  }
  todrop
}

mjn <- function(data, epsilon=0, inferred.prefix='i') {
  if (is.null(row.names(data))) { rownames(data) <- 1:nrow(data) }
  starting.rownames <- rownames(data)
  sampled <- 1:nrow(data)
  rownames(data) <- sampled
  x <- list(n=1, lambda=NULL, prev.combos=combn(1:2,2)[,F])
  while (x$n > 0) {
    distances <- mandist(data)
    edges <- t(rbind(combn(as.integer(rownames(data)),2), as.vector(distances)))
    colnames(edges) <- c("n1", "n2", "d")
    dsc <- delta.step.components(edges, epsilon)

# I think for the set epsilon that the new sequence types code already ensures
# that there are no obsolete sequence types
#     obsoletes <- identify.obsoletes(dsc$g, setdiff(1:nrow(data), sampled))
#     if (length(obsoletes) > 0) {
#       data <- data[-obsoletes,]
#       row.names(data) <- 1:nrow(data)
#     } else {
    x <- new.sequence.types(data, edges, dsc, epsilon, x$lambda, x$prev.combos)
    data <- rbind(data,x$median.vectors)
    rownames(data) <- 1:nrow(data)
#     }
  }
  repeat {
    distances <- mandist(data)
    edges <- t(rbind(combn(as.integer(rownames(data)),2), as.vector(distances)))
    colnames(edges) <- c("n1", "n2", "d")
    dsc <- delta.step.components(edges, 0)
    obsoletes <- identify.obsoletes(dsc$g, setdiff(1:nrow(data), sampled))
    if (length(obsoletes) > 0) {
      data <- data[-obsoletes,]
      row.names(data) <- 1:nrow(data)
    } else {
      break
    }
  }
  num.inferred.nodes <- nrow(data) - length(starting.rownames)
  inferred <- c(rep(F, length(starting.rownames)), rep(T, num.inferred.nodes))
  if (num.inferred.nodes > 0) {
    starting.rownames <- c(starting.rownames, 
                           paste(inferred.prefix, 1:num.inferred.nodes, sep=''))
  }
  rownames(data) <- starting.rownames
  g=dsc$g
  V(g)$name <- starting.rownames
  V(g)$inferred <- inferred
  ret <- list(g=g, data=data)
  class(ret) <- 'mjn'
  ret
}

map.colors <- function(color.map, ids, default.color) {
  col <- color.map[ids]
  col[is.na(col)] <- default.color
  col
}

# cmap is a named vector of colors; names are vertex names
# count is a vector of haplotype counts in vertex order
# cprops is a list, named by vertex name, with proportions of colors
plot.mjn <- function(x, vsize=10, vlabel=NULL, 
                     count=1,
                     props=NULL,
                     default.color="grey",
                     cmap=NULL,
                     layout=NULL,
                     inferred.shape='circle',
                     inferred.size=0.5,
                     elabel=NULL) {
  g <- x$g
  count <- count[V(g)$name]
  count[is.na(count)] <- 1
  V(g)$size <- sqrt(count*vsize^2)
  V(g)[V(g)$inferred]$size <- V(g)[V(g)$inferred]$size * inferred.size
  V(g)$shape <- "circle"
  V(g)[V(g)$inferred]$shape <- inferred.shape
  vcolor <- default.color
  if (!is.null(cmap)) {
    vcolor <- map.colors(cmap, V(g)$name, default.color)
  }
  if (is.null(vlabel)) { vlabel <- V(g)$name }
  if (is.null(elabel)) {
    E(g)$label <- round(E(g)$weight)
  }
  if (is.null(layout)) { 
    wmax <- max(E(g)$weight)
    layout <- layout.fruchterman.reingold(g, weights=wmax/E(g)$weight)
  }
  plot(g, 
       layout=layout, 
       vertex.label=vlabel, 
       vertex.color=vcolor)
  
  if (!is.null(props)) {
    selector <- V(g)$name %in% names(props)
    if (any(selector)) {
      coords <- layout.norm(layout,-1,1,-1,1)[selector,]
      rownames(coords) <- V(g)[selector]$name
      for (nam in V(g)[selector]$name) {
        props <- lprops[nam][[1]]
        cols <- names(props)
        from <- 0
        radius <- V(g)[nam]$size/200
        mid <- coords[nam,]
        for (i in 1:length(props)) {
          to <- from + props[i] * 2*pi
          if (to > pi) { to <- to - 2*pi }
          xyvals <- rbind(mid, getellipse(radius, mid=mid, from=from, to=to), mid)
          polygon(xyvals, col=cols[i], border=NA)
          from <- to
        }
      }
    }
  }
}

mjn.merge <- function(mjn.list) {
  combined.data <- do.call("rbind", lapply(mjn.list, function(x) { x$data } ))
  graphs <- lapply(mjn.list, function(x) { x$g } )
  G <- graph.disjoint.union(graphs)
  node.starting.index <- 0
  edge.starting.index <- 0
  for (i in 1:length(graphs)) {
    g <- graphs[[i]]
    num.nodes <- length(V(g))
    if (num.nodes > 0) {
      node.indices <- node.starting.index:(node.starting.index + num.nodes - 1)
      node.starting.index <- max(node.indices) + 1
      V(G)[node.indices]$name <- V(g)$name
      V(G)[node.indices]$inferred <- V(g)$inferred
    }
    num.edges <- length(E(g))
    if (num.edges > 0) {
      edge.indices <- edge.starting.index:(edge.starting.index + num.edges - 1)
      edge.starting.index <- max(edge.indices) + 1
      E(G)[edge.indices]$weight <- E(g)$weight
    }
  }
  ret <- list(g=G, data=combined.data)
  class(ret) <- 'mjn'
  ret
}

