#     Sample C from Gwishart distribution with  density:
#          p(C) \propto |C|^{(bG-2)/2} exp(-1/2 tr(C DG))
#     where     
#      (1)  bG : d.f.
#      (2)  DG: location
#     C: initial partial covariance matrix;
#     burnin, nmc : number of MCMC burnins and saved samples
GWishart_PAS_DMH <- function( b_prior, D_prior, n, S, C, beta, burnin, nmc )
{
  p        <- dim( D_prior )[1] 
  b_post   <- b_prior + n 
  D_post   <- D_prior + S 
  C_save   <- array( 0, c( p, p, nmc ) )
  Sig_save <- C_save
  adj_save <- C_save
  adj      <- 1 * ( abs( C ) > 1e-5 )
  nedge    <- ( colSums( adj - p ) ) / 2 
  C        <- C * adj 
  if( ( burnin + nmc ) > 0 ){
    for( iter in 1:( burnin + nmc )  ){    
      if( (iter %% 1000 ) == 0 & iter > burnin ){
        cat( 'iter = ', iter, ' nedge = ', nedge, '\n' ) 
        apply( adj_save[ , ,1:( iter - burnin - 1 )], c( 1, 2 ), mean )        
      }
      ## Sample off-diagonal elements     
      for( i in  1:( p - 1 )  ){  
        for( j in ( i + 1 ):p ){
          nedge      <- ( colSums( adj - p ) ) / 2 
          w          <- log_H( b_prior, D_prior, n, S, C, i, j ) 
                        + log( 1 - beta ) - log( beta )
          w          <- 1 / ( exp( w ) + 1 ) 
          current_ij <- adj[i,j]
          propose_ij <- 1 * ( runif( 1 ) < w )
          if( propose_ij != current_ij ){
            resN1ij  <- GWishart_NOij_Gibbs( b_prior, D_prior, adj, C, i, j, propose_ij, 0, 1 )
            C_prop   <- resN1ij[[1]]
            Sig_prop <- resN1ij[[2]]
            r2       <- log_GWishart_NOij_pdf( b_prior, D_prior, C_prop, i, j, current_ij ) 
                      - log_GWishart_NOij_pdf( b_prior, D_prior, C_prop, i, j, propose_ij )
            if( log( runif( 1 ) ) < r2 ){
              adj[i,j]   <- propose_ij
              adj[j,i]   <- propose_ij
              current_ij <- propose_ij
            }
          } # Second MH
          resN0ij <- GWishart_NOij_Gibbs( b_post, D_post, adj, C, i, j, current_ij, 0, 0 )
          C       <- resN0ij[[1]]
          Sig     <- resN0ij[[2]]
        }
      }
      #########  Update C and Sigma given graph
      resBIPS <- GWishart_BIPS_maximumClique( b_post, D_post, adj, C, 0, 1 )
      C       <- resBIPS[[1]][,,1]
      Sig     <- resBIPS[[2]][,,1]
      if( iter > burnin ){
        Sig_save[,,( iter - burnin )] <- Sig 
          C_save[,,( iter - burnin )] <- C
        adj_save[,,( iter - burnin )] <- adj
      }
    }
  }
  return( list( C = C_save, Sig = Sig_save, adj = adj_save ) )
}
