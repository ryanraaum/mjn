options(stringsAsFactors=F)

library(igraph)

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


M <- as.matrix(data.frame(c1=c(0,1,1,0),c2=c(0,3,0,3),c3=c(0,0,2,2),c4=c(0,0,1,0),c5=c(0,0,0,2)))
#row.names(df) <- c('A', 'B', 'C', 'D')

M2 <- as.matrix(data.frame(c1=c(1,1,1,0),c2=c(0,3,1,3),c3=c(0,0,2,2),c4=c(0,0,1,0),c5=c(0,0,0,2)))
#row.names(df2) <- c('A', 'B', 'C', 'D')

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
  dups <- duplicated(cbind(prev.combos, combos), MARGIN=2)[-(1:ncol(prev.combos))]
  combos <- combos[,!dups,drop=F]
  mvecs <- data[F,]
  c.costs <- c()
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
  keep <- F
  if (length(c.costs) > 0) { 
    if (is.null(lambda)) { lambda=min(c.costs) }
    if (length(c.costs)>0) {
      keep <- c.costs <= (lambda + epsilon)
    }
  }
  prev.combos <- cbind(prev.combos, combos)
  if (sum(keep)) {
    return(list(n=sum(keep), costs=c.costs[keep], median.vectors=mvecs[keep,], lambda=lambda, prev.combos=prev.combos))
  } else {
    return(list(n=0, costs=c(), median.vectors=data[F,], lambda=lambda))
  }
}


mjn <- function(data, epsilon=0) {
  starting.rownames <- rownames(data)
  sampled <- 1:nrow(data)
  rownames(data) <- sampled
  x <- list(n=1, lambda=NULL, prev.combos=combn(1:2,2)[,F])
  while (x$n > 0) {
    distances <- mandist(data)
    edges <- t(rbind(combn(as.integer(rownames(data)),2), as.vector(distances)))
    colnames(edges) <- c("n1", "n2", "d")
    dsc <- delta.step.components(edges, epsilon)
    
    # check for obsolete sequence types
    # if obsolete exist, remove them and redo delta step componenets
    inferred <- setdiff(1:nrow(data), sampled)
    todrop <- c()
    for (i in inferred) {
      vindex <- inferred - 1
      if (length(V(dsc$g)[nei(vindex)]) <= 2) {
        todrop <- c(todrop, i)
      }
    }
    if (length(todrop) > 0) {
      data <- data[-todrop,]
      row.names(data) <- 1:nrow(data)
    } else {
      x <- new.sequence.types(data, edges, dsc, epsilon, x$lambda, x$prev.combos)
      data <- rbind(data,x$median.vectors)
      rownames(data) <- 1:nrow(data)
    }
  }
  distances <- mandist(data)
  edges <- t(rbind(combn(as.integer(rownames(data)),2), as.vector(distances)))
  colnames(edges) <- c("n1", "n2", "d")
  final.dsc <- delta.step.components(edges, 0)
  list(g=final.dsc$g, data=data)
}

LABS = c(LETTERS[1:4], LETTERS[21:26])

x <- mjn(M,2)
g <- x$g
plot(g, layout=layout.fruchterman.reingold(g, weights=1/E(g)$weight), edge.label=E(g)$weight, vertex.label=LABS[1:length(V(g))])


