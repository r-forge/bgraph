 #     Sample C from Gwishart distribution with  density:
 #          p(C) \propto |C|^{(bG-2)/2} exp(-1/2 tr(C DG))
 #     where     
 #      (1)  bG : d.f.
 #      (2)  DG: location
 #      (3)  adj: adjacency matrix 
 #     C: initial partial covariance matrix;
 #     burnin, nmc : number of MCMC burnins and saved samples
 GWishart_BIPS_maximumClique <- function( bG, DG, adj, C, burnin, nmc )
 {
 
 adj0       <- adj - diag( diag( adj ) ) # Create adjacency matrix with diagonal elements zero
 p          <- dim( DG )[1] 
 cliqueList <- maximal.cliques( graph.adjacency( adj0 ) )
 numcliques <- length( cliqueList )
 C_save     <- array( 0, c( p, p, nmc ) )
 Sig_save   <- C_save  
 C          <- C * adj

  if( !is.null( numcliques ) ){
    if( ( burnin + nmc ) > 0 ){
      for( iter in 1:( burnin + nmc ) ){ 
        for( i in 1:numcliques ){
          cliqueid   <- cliqueList[[i]] + 1 #which(cliqueMatrix[[i]]==1);
          cliquesize <- length( cliqueid )
          A          <- rwish( bG + cliquesize - 1, solve( DG[cliqueid,cliqueid,drop=FALSE] ) )
          C_12       <- C[ cliqueid,-cliqueid,drop = FALSE]
          C_22       <- C[-cliqueid,-cliqueid,drop = FALSE] 
          C121       <- C_12 %*% solve( C_22 ) %*% t( C_12 ) 
          C121       <- ( C121 + t( C121 ) ) / 2
          K_c        <- A + C121
          C[cliqueid,cliqueid] = ( K_c + t( K_c ) ) / 2
        } #end for
        if( iter > burnin ){
            Sig_save[,,( iter - burnin )] = solve( C ) 
              C_save[,,( iter - burnin )] = C
        } #end if
      } #end for
    } #end if
  } #end if
  return( list( C= C_save, Sig = Sig_save ) )
}

