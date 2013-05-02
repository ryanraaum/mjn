M <- as.matrix(data.frame(c1=c(0,1,1,0),c2=c(0,3,0,3),c3=c(0,0,2,2),c4=c(0,0,1,0),c5=c(0,0,0,2)))
rownames(M) <- 1:nrow(M)
Mdist <- mandist(M)
Medges <- t(rbind(combn(as.integer(rownames(M)),2), as.vector(Mdist)))
Mdeltas <- Medges[,3]
Mes <- Medges[,1:2] - 1
colnames(Medges) <- c('n1', 'n2', 'd')

M2 <- as.matrix(data.frame(c1=c(1,1,1,0),c2=c(0,3,1,3),c3=c(0,0,2,2),c4=c(0,0,1,0),c5=c(0,0,0,2)))
rownames(M2) <- 1:nrow(M2)
M2dist <- mandist(M2)
M2edges <- t(rbind(combn(as.integer(rownames(M2)),2), as.vector(M2dist)))
M2deltas <- M2edges[,3]
M2es <- M2edges[,1:2] - 1
colnames(M2edges) <- c('n1', 'n2', 'd')


LABS = c(LETTERS[1:4], LETTERS[21:26])

x <- mjn(M,2)
g <- x$g
plot(g, layout=layout.fruchterman.reingold(g, weights=1/E(g)$weight), edge.label=E(g)$weight, vertex.label=V(g)$name)
