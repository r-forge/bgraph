rand_truncated_normal_rej_below <- function( mu, sig, below ){
# Generate one sample from N(mu,sig^2)1_{x>below}
# Ref: 
#     Proposition 2.3  C. P. Robert 'Simulation of truncated normal variables'

    mu_neg <- (below - mu)/sig;

    if(mu_neg < 0){
        
        x <- rnorm(1);
        while (x < mu_neg){
            x <- rnorm(1);
        }#end    
        x <- x * sig + mu;
    } else {
        
        alpha <- ( mu_neg + sqrt( mu_neg^2+4))/2; # Optimal alpha
        
        x <- rexp(1, 1/alpha) + mu_neg; #exprnd(1/alpha) 
            
        while (log(runif( 1 ) ) > (-(x-alpha)^2/2) ){
            x <- rexp(1, 1/alpha) + mu_neg;
        } #end
        
        x <- x * sig + mu;
        
    }#end
    return(x)
}





