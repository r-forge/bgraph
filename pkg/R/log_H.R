log_H <- function( b_prior, D_prior, n, S, C, i, j )
{
    #########  (i,j) = 0
    e       <- j
    C0      <- C
    C0[i,j] <- 0
    C0[j,i] <- 0     
    C_12    <- C0[ e,-e,drop=FALSE]
    C_22    <- C0[-e,-e,drop=FALSE]
    invC_22 <- solve(C_22)
    C0_ij   <- matrix( c( C[i,i], 0, 0, C_12 %*% invC_22 %*% t( C_12 ) ), 2, 2, byrow = TRUE )            
    #########  (i,j) = 1
    e       <- c( i, j )
    C_12    <- C[ e,-e,drop=FALSE]  
    C_22    <- C[-e,-e,drop=FALSE] 
    invC_22 <- solve(C_22)
    Ce      <- C_12 %*% invC_22 %*% t( C_12 )
    A       <- C[e,e,drop=FALSE] - Ce  
    a11     <- A[1,1,drop=FALSE]
    C1_ij   <- Ce 
    b_post  <- b_prior + n
    D_post  <- D_prior + S
    h       <- ( - log_iwishart_InvA_const( b_post, D_post[j,j,drop=FALSE] )
                 - log_J( b_post, D_post[e,e,drop=FALSE], a11 ) 
                 + ( n + b_prior - 2 ) / 2 * ( log( a11 ) ) 
                 - sum( diag( D_post[e,e,drop=FALSE] %*% ( C0_ij - C1_ij ) ) )/2 
               )
    return( h )
}
