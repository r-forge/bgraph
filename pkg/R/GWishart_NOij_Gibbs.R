#   Start from C\(i,j)  
#
#     Sample C from Gwishart distribution with  density:
#          p(C) \propto |C|^{(bG-2)/2} exp(-1/2 tr(C DG))
#     where     
#      (1)  bG : d.f.
#      (2)  DG: location
#      (3)  adj: adjacency matrix 
#     C: initial partial covariance matrix;
#     burnin, nmc : number of MCMC burnins and saved samples

GWishart_NOij_Gibbs <- function( bG, DG, adj, C, i, j, edgeij, burnin, nmc )
{
  p <- dim( DG )[1] 
  ## If no edge
  if( edgeij == 0 ){          
    C[i,j]  <- 0  
    C[j,i]  <- 0
    l       <- solve( DG[j,j,drop=FALSE] )
    A       <- rwish( bG, l )
    C_12    <- C[ j,-j,drop=FALSE] 
    C_22    <- C[-j,-j,drop=FALSE]
    invC_22 <- solve( C_22 ) 
    C[j,j]  <- A + C_12 %*% invC_22 %*% t( C_12 )
  } else {    
    reorder   <- 1:p
    reorder   <- c(reorder[-c( i, j )] ,i,j)
    C_reorder <- C[reorder,reorder,drop=FALSE]
    R         <- chol(C_reorder)
    b         <- R[( p - 1 ),( p - 1 )]
    m_post    <- -b * DG[i,j]/DG[j,j]               
    sig_post  <- 1/sqrt(DG[j,j]) 
    adj[i,j]  <- 1
    adj[j,i]  <- 1
    R[p-1,p]  <- rnorm(1) * sig_post + m_post 
    R[p,p]    <- sqrt(rgamma(1, shape=bG/2, scale = 2/DG[j,j,drop=FALSE]) )  #gamrnd(bG/2,2/DG[j,j])
    dR        <- dim(R)
    C_updated <- t( R[,(dR[2]-1):dR[2],drop=FALSE]) %*% R[,dR[2],drop=FALSE]
    C[i,j]    <- C_updated[1]
    C[j,i]    <- C_updated[1]
    C[j,j]    <- C_updated[2]       
  }
  Sig = solve(C)

  IsolatedNodeId <- which( colSums( adj ) == 1 ) # isolated node, i.e. nodes that are not connected to any other nodes
  INsize         <- length( IsolatedNodeId )
  if( ( burnin + nmc ) > 0 ){
    for( iter in 1:( burnin + nmc ) ) {           
      ###  Sample isolated nodes
      if( INsize > 0 ){
        for( i in 1:INsize ){    
           cliqueid <- IsolatedNodeId[i]    
           K_c      <- rwish( bG,solve( DG[cliqueid,cliqueid,drop=FALSE] ) )
             C[cliqueid,cliqueid] <- K_c 
           Sig[cliqueid,cliqueid] <- 1/K_c
        }
      }
      for( i in 1:( p - 1) ){
        for( j in ( i + 1 ):p ){    
          if( adj[i,j] == 1 ){     
            cliqueid   <- c(i,j)
            cliquesize <- 2
            l          <- solve( DG[cliqueid,cliqueid,drop=FALSE] )
            l          <- ( l + t( l ) )/2
            A          <-  rwish(bG+cliquesize-1,l)
            C_12       <-   C[ cliqueid,-cliqueid,drop=FALSE] 
            Sig_12     <- Sig[ cliqueid,-cliqueid,drop=FALSE]     
            Sig_22     <- Sig[-cliqueid,-cliqueid,drop=FALSE]
            invSig_11  <- solve( Sig[cliqueid,cliqueid,drop=FALSE] )
            invSig_11  <- (invSig_11 + t(invSig_11))/2
            invC_22    <- Sig_22 - t( Sig_12 ) %*% invSig_11 %*% Sig_12
            K_c        <- A + C_12 %*% invC_22 %*% t( C_12 )
            K_c        <- ( K_c + t( K_c ) ) / 2 # Numerical stable
            Delta      <- solve( C[cliqueid,cliqueid,drop=FALSE] - K_c )
            C[cliqueid,cliqueid] <- K_c
             
            ## Rank 2 update Sig 
            Sig_bb <- Sig[cliqueid,cliqueid,drop=FALSE]        
            aa     <- solve( Delta - Sig_bb )
            aa     <- ( aa + t( aa ) ) / 2
            Sig    <- Sig + Sig[,cliqueid,drop=FALSE] %*% aa %*% t(Sig[,cliqueid,drop=FALSE])
          }
        }
      }
    }
  }
  return( list( C = C, Sig = Sig ) )
}
