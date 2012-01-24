M <- as.matrix(data.frame(c1=c(0,1,1,0),c2=c(0,3,0,3),c3=c(0,0,2,2),c4=c(0,0,1,0),c5=c(0,0,0,2)))
#row.names(df) <- c('A', 'B', 'C', 'D')

M2 <- as.matrix(data.frame(c1=c(1,1,1,0),c2=c(0,3,1,3),c3=c(0,0,2,2),c4=c(0,0,1,0),c5=c(0,0,0,2)))
#row.names(df2) <- c('A', 'B', 'C', 'D')

LABS = c(LETTERS[1:4], LETTERS[21:26])

x <- mjn(M,2)
g <- x$g
plot(g, layout=layout.fruchterman.reingold(g, weights=1/E(g)$weight), edge.label=E(g)$weight, vertex.label=V(g)$name)
