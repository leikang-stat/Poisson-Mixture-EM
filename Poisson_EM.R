

em_Poisson <- function(X, tol=.Machine$double.eps){ ##x takes input of data.frame
  
  N <- nrow(X) ##number of observations
  x <- X$x
  error <- Inf
  iter <- 1
  
  ##initial guess, random starts
  mu0 <- sample(1:100, 1)
  mu1 <- sample(1:100, 1)
  phi <- runif(1,0,1)

  while(error > tol ){
    
    ##E-step
    X$q_x1 <- phi*dpois(x,mu0)
    X$q_x2 <- (1-phi)*dpois(x,mu1)
    X$P1_x <- X$q_x1/(X$q_x1+X$q_x2) ##P=1|X
    X$P2_x <- X$q_x2/(X$q_x1+X$q_x2) ##P=2|X
    Q <- sum(log(X$q_x1)*X$P1_x)+sum(log(X$q_x2)*X$P2_x)
    
    ##M-step/update parameters
    mu0_k <- sum(x*X$P1_x)/sum(X$P1_x)
    mu1_k <- sum(x*X$P2_x)/sum(X$P2_x)
    phi_k <- sum(X$P1_x)/N
    
    ##compare Q
    X$q_x1_k <- phi_k*dpois(x,mu0_k)
    X$q_x2_k <- (1-phi_k)*dpois(x,mu1_k)
    X$P1_x_k <- X$q_x1/(X$q_x1+X$q_x2) 
    X$P2_x_k <- X$q_x2/(X$q_x1+X$q_x2) 
    Q_k <- sum(log(X$q_x1_k)*X$P1_x_k)+sum(log(X$q_x2_k)*X$P2_x_k)
    
    ##stop criterion
    error <- Q_k-Q
    iter <- iter+1
    mu0 <- mu0_k
    mu1 <- mu1_k
    phi <- phi_k
  }
  theta<-c(mu0,mu1,phi)
  
  ##cluster assingment based on posterior probability, normalizing term is canceled out
  X$cluster[phi*dpois(x,mu0)>=(1-phi)*dpois(x,mu1)] <- 1
  X$cluster[phi*dpois(x,mu0)<(1-phi)*dpois(x,mu1)] <- 2

  return(list(X,theta))
}


library(ggplot2)

##generate simulated data
set.seed(12345)
n <- 1000

##mixture of two poissons
mu0_true <- 5
mu1_true <- 7
phi_true <- 0.4
label1 <- rep(1,n * phi_true)
label2 <- rep(2,n * (1-phi_true))
X1 <- cbind(rpois(n*phi_true,mu0_true),label1)
X2 <- cbind(rpois(n*(1-phi_true),mu1_true),label2)
X <- data.frame(rbind(X1,X2))
names(X)[names(X)=="label1"] <- "label"
names(X)[names(X)=="V1"] <- "x"

options(digits=22)
.Machine$double.eps<-2.220446e-20

####EXAMPLE: EM-based Clustering
a <- em_Poisson(X) ##a is a list containing X with cluster assignment and parameter estimates


####calculate similarity measures, like jaccard, correlation, matching
clust_1<-as.integer(as.character(a[[1]]$label))
clust_2<-as.integer(as.character(a[[1]]$cluster))

inter.dim <- dim(X)[1]
C_1 <- matrix(clust_1, nr = inter.dim, nc = inter.dim) == matrix(clust_1, nr = inter.dim, nc = inter.dim, byrow = TRUE)
C_2 <- matrix(clust_2, nr = inter.dim, nc = inter.dim) == matrix(clust_2, nr = inter.dim, nc = inter.dim, byrow = TRUE)
diag(C_1) <- 0
diag(C_2) <- 0

##three matching metrics, see https://github.com/leikang-stat/Stability-based-clustering for details
jaccard <- sum(C_1 * C_2)/(sum(C_1) + sum(C_2) - sum(C_1 * C_2))
matching <- (sum(C_1 * C_2)+sum((1-C_1) * (1-C_2)))/(sum(C_1 * C_2)+sum((1-C_1) * (1-C_2))+sum((1-C_1)*C_2)+sum((1-C_2)*C_1))
corr <- sum(C_1 * C_2)/sqrt(sum(C_1)*sum(C_2))
print(jaccard)
print(matching)
print(corr)

####plot the clustering results
m <- ggplot(X, aes(x = x,color=as.factor(label)))
m1 <- m + geom_histogram(binwidth = 1)+ggtitle("True Cluster Result")

m2 <- ggplot(a[[1]], aes(x = x,color=as.factor(cluster)))
m2 <- m2 + geom_histogram(binwidth = 1)+ggtitle("Estimated Cluster Result")

multiplot(m1,m2,cols=1)

##multi ploting function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
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