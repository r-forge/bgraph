#     BIPS algorithm for sampling from G-Wishart distribution with  density:
#          p(C) \propto |C|^{(bG-2)/2} exp(-1/2 tr(C DG))
#     where     
#      (1)  bG : d.f.
#      (2)  DG: location
#      (3)  adj: adjacency matrix 
#     C: initial precision matrix;
#     burnin, nmc : number of MCMC burnins and saved samples
#     Reference: Wang and Li (2011) "Efficient Gaussian Graphical Model Determination without 
#      Approximating Normalizing Constant of the G-Wishart distribution "

GWishart_BIPS_pairwise <- function( bG, DG, adj, C, burnin, nmc )
{
  p        <- dim( DG )[1] 
  C_save   <- array( 0, c( p, p, nmc ) )
  Sig_save <- C_save
  C        <- C * adj 
  Sig      <- solve( C )

  IsolatedNodeId = which( colSums( adj ) == 1 )# isolated node
  INsize         = length(IsolatedNodeId)

  if( ( burnin + nmc ) > 0 ){
    for( iter in 1:( burnin + nmc ) ){    
      if( ( iter %% 1000 ) == 0 ){
          cat( "iter = ", iter, " \n"  )
      }#end

      ###  Sample isolated nodes
      if( INsize > 0 ){
        for( i in 1:INsize ){   
          cliqueid               <- IsolatedNodeId[i]    
          K_c                    <- rwish( bG,solve( DG[cliqueid,cliqueid] ) )
            C[cliqueid,cliqueid] <- K_c
          Sig[cliqueid,cliqueid] <- 1 / K_c  
        }#end
      }

      for( i in 1:( p - 1 ) ){
        for( j in ( i + 1 ):p ){
          if( adj[i,j] == 1 ){
            cliqueid   <- c( i, j )
            cliquesize <- 2
            A          <- rwish( bG + cliquesize - 1, solve( DG[cliqueid,cliqueid,drop=FALSE] ) ) 
            C_12       <-   C[ cliqueid,-cliqueid,drop=FALSE] 
            Sig_12     <- Sig[ cliqueid,-cliqueid,drop=FALSE]         
            Sig_22     <- Sig[-cliqueid,-cliqueid,drop=FALSE] 
            invSig_11  <- solve( Sig[cliqueid,cliqueid,drop=FALSE] )
            invSig_11  <- ( invSig_11 + t( invSig_11 ) )/2
            invC_22    <- Sig_22 - t( Sig_12 ) %*% invSig_11 %*% Sig_12
            K_c        <- A + C_12 %*% invC_22 %*% t( C_12 )
            K_c        <- ( K_c + t( K_c ) ) / 2 # Numerical stable
            Delta      <- solve( C[cliqueid,cliqueid,drop=FALSE] - K_c )
            C[cliqueid,cliqueid] <- K_c
            
            ## Rank 2 update Sig 
            Sig_bb <- Sig[cliqueid,cliqueid,drop=FALSE]
            aa     <- solve( Delta - Sig_bb )
            aa     <- ( aa + t( aa ) ) / 2
            Sig    <- Sig + Sig[,cliqueid,drop=FALSE] %*% aa %*% t( Sig[,cliqueid,drop=FALSE] )
          } #end
        } #end
      } #end

      if( iter > burnin ){
        Sig_save[,,( iter - burnin )] = Sig 
          C_save[,,( iter - burnin )] = C
      }#end
    }#end
  }
  return( list( C = C_save, Sig = Sig_save ) )
}
