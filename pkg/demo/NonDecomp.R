##############################################################
#### (I) Sample from G-Wishart Distribution for given graphs #####
##############################################################

###############################################################
##### e.g. (1) Sparse Circle graph:  Dobra Lenkoski and Abel (2011, JASA) Example      #######
###############################################################
p <- 10;  bG <- 103;
A <- toeplitz(c( 1,0.5,rep(0,p-2)));
A[1,p]<- 0.4; A[p,1] <- 0.4;

D   <- diag(p) + 100*solve(A);
adj <- 1 * ( abs( A ) > 0.001 );
#DG  <- MLE_GGM( D, adj, 200, 0.0001, 0 ); 


##########################################################
##### e.g. (2) dense graph      ###########################
##########################################################
p     <- 10 
alpha <- 0.3 
J0    <- 0.5*diag(p) 
B0    <- 1*(matrix(runif(p*p),p,p)<alpha)
for( i in 1:p ){ 
    for( j in 1:i ){ 
        if ( B0[i,j] && i!=j ){
            J0[i,j] <- 0.5
        }
    }
} #end end end
J0    <- J0 + t(J0)
tmp   <- eigen( J0 ) w <- tmp$values#w = eig(J0) 
delta <- ( p * min( w ) - max( w ) ) / ( 1 - p ) 
J0    <- J0 + delta * diag( p )
D     <- diag(p) + 100*solve(J0)
DG    <- D
adj   <- 1*(abs(J0)>1e-4) 
bG    <- 103

burnin <- 1000 nmc <- 1000
### algorithm (i) Edgewise Gibbs algorithm  ################
C      <- diag(p) # Initial value

resBIPSpair       <- GWishart_BIPS_pairwise(bG,DG,adj,C,burnin,nmc)
C_save_edgewise   <-resBIPSpair[[1]]
Sig_save_edgewise <-resBIPSpair[[2]]
###  Maximum clique algorithm  ################
C            <- diag( p ) # Initial value
resBIPSmax   <- GWishart_BIPS_maximumClique( bG,DG,adj,C,burnin,nmc)
C_save_maxC  <- resBIPSmax[[1]]
Sig_save_maxC<- resBIPSmax[[2]]
 
#########################################################################
#### (II) Graphical Model Deterimation Without Approximating Normalizing 
##### of G-Wishart Distribution 
#########################################################################

##############################################
###### e.g. 1:  a 6-node circle graph example 
##############################################

p       <- 6 
tedge   <- p * ( p - 1 ) / 2 
beta    <- 0.5 
indmx   <- matrix( 1:( p ^ 2 ), p, p ) 

b_prior <- 3 
D_prior <- diag( p ) 
n       <- 3 * p 
b_post  <- b_prior+n
A       <- toeplitz( c( 1, 0.5, matrix( 0, 1, p - 2 ) ) ) 
A[1,p]  <- 0.4 
A[p,1]  <- 0.4
S       <- solve( A ) * n 
D_post  <- D_prior + S
adjTrue <- 1*( abs( A ) > 0.001 )  

burnin  <- 300 
nmc     <- 100     
C       <- diag( p ) 

resPAS_DMH <- GWishart_PAS_DMH(b_prior,D_prior,n,S,C,beta,burnin,nmc)
C_save     <- resPAS_DMH$C
Sig_save   <- resPAS_DMH$Sig
adj_save   <- resPAS_DMH$adj

