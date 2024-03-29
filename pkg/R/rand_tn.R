
# Sample from truncated normal N(mu,sig^2)1_{lower<r<upper}
rand_tn <- function(mu,sig,lower,upper){

    z = ( c( lower, upper ) - mu)/sig;
    p = pnorm( z );
    if( p[2] < 1e-30 ){
        r <- upper;
    } else {
        u <- ( p[2] - p[1] ) * runif( 1 ) + p[1];
        u <- min( max( 1e-30, u ), 1-1e-16 );
        r <- qnorm( u );
        r <- r * sig + mu;
    } #end
    return(r)
}


# function [r] = rand_tn(mu,sig,lower,upper)
# % Sample from truncated normal N(mu,sig^2)1_{lower<r<upper}
# z=([lower,upper] - mu)/sig;
# p = normcdf(z);

# if(p(2)<1e-30)
#     r = upper;
# else
#     u = (p(2)-p(1))*rand(1)+p(1);
#     u = min(max(1e-30,u),1-1e-16);
#     r = norminv(u);
#     r = r*sig+mu;
# end
