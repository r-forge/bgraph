
#  Generates the normalizing constant "cons."  for a Inv-Wishart(df,S).
# Nonsingular  pdf is p(K) = cons. |K|^{-(df+2|p|)/2} exp(-trace(inv(K) S)/2) 
log_iwishart_InvA_const <- function( df, S )
{
  S   = as.matrix(S)
  p   = dim( S )[1];
  iwc = ( df + p - 1 ) / 2 * ( log( det( S ) ) - p * log( 2 ) ) - log_multi_gamma( p, ( df + p - 1 ) / 2 );
  return( iwc );
}


#function iwc=log_iwishart_InvA_const(df,S)
#%  Generates the normalizing constant "cons."  for a Inv-Wishart(df,S).
#% Nonsingular  pdf is p(K) = cons. |K|^{-(df+2|p|)/2} exp(-trace(inv(K) S)/2) 

#[p,p]=size(S);
#iwc=(df+p-1)/2*(log(det(S))-p*log(2))-log_multi_gamma(p,(df+p-1)/2);
