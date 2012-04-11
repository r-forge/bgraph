# log normalizing constant approximation using importance sampling for a
# G-Wishart(b,D,adj):
#
#  p(K) = I(b,D)^{-1} |K|^{(d-2)/2} exp(-trace(K D)/2) using Monte Carlo sample of size N;
#
# Reference: Atay-Kayis and Massam 2005 Biometrika
#
# INPUT:   b,D,adj
# OUTPUT:  log(I(b,D))

log_GWishart_ud_const_mc <- function(b,D,adj,N){

    T  = chol( solve( D ) );
    T1 = T %*% diag( 1/diag( T ) );
    p  = dim( D )[1];

    indmx   = matrix( 1:(p^2), p, p ); # Uppermatrix index
    upperind=indmx[upper.tri(indmx)];      
    A       = rep( 0, p,p );
    A[upper.tri(A)]= adj[upper.tri(adj)] #A(upperind) = adj(upperind);
    nu      = rowSums( A );

    logJeach = matrix( 0, 1, N ) #zeros(1,N);

    ### D is non-diagnoal then computational cost O(p^4)
    if( sum( abs( c( T1 ) ) ) - sum( diag( T1 ) ) != 0 )   { 
            
        for( iter in 1:N ){
            Psi              = diag(sqrt(rchisq(1,b+nu)));
            Psi[find(A==1)]  = matrix(rnorm(sum(nu)), sum(nu),1) #randn(sum(nu),1);
            
            logJeach[iter]   = 0 ;
            ###### if i = 1 ############
            for( j in 2:p ){
                if( A[1,j] == 0 ){
                    Psi[1,j]       = -sum( t( Psi[1,1:(j-1)]) * T1[1:(j-1),j]);   
                    logJeach[iter] = logJeach[iter] - Psi[1,j]^2/2;
                } #end
            } #end


            ##### if i>1 ###################
            for( i in 2:p ){
                for( j in ( i + 1 ):p ){
                    if( A[i,j] == 0 ){
                        Psi[i,j] = -sum(t(Psi[i,i:(j-1)])*T1[i:(j-1),j]); 
                        for( r in 1:( i-1 ) ){
                            Psi[i,j] = Psi[i,j] - (1/Psi[i,i])*(Psi[r,i]+ sum(Psi[r,r:(i-1)]*t(T1[r:(i-1),i])))*(Psi[r,j]+ sum(Psi[r,r:(j-1)]*t(T1[r:(j-1),j])));
                        } #end
                        logJeach[iter] = logJeach[iter] - Psi[i,j]^2/2;   
                        } #end
                    }#end
            }#end

        }#end # for iter = 1:N

    ##########################  ELSE IF D IS diagonal MATRIX then simplified
    ##########################  calculation at O(p^2)
    }else{  ########

      for( iter in 1:N ){
        Psi = diag(sqrt(rchisq(1,b+nu)));
        Psi[which(A==1)]  = matrix(rnorm(sum(nu)), sum(nu),1) #randn(sum(nu),1);
        
        logJeach[iter] = 0 ;
        ##### if i>1 ###################
        for( i in 2:p ){
            for( j in ( i + 1 ):p ){
                if(A[i,j]==0){       
                    Psi[i,j]       <-  - 1/Psi[i,i] * sum(Psi[1:(i-1),i]*Psi[1:(i-1),j]);   
                    logJeach[iter] <- logJeach[iter] - Psi[i,j]^2/2;   
                   
                    if (is.nan(logJeach[iter])){
                        stop("error, NA");
                    }#end
                }#end   
            }#end
        }#end

      }#end # for iter = 1:N

    }#end # END IF T is diagonal
     
    bi     = rowSums( adj );
    logC   = sum( nu/2*log(2*pi)+(b+nu)/2*log(2)+lgamma((b+nu)/2)+(b+bi-1)*log(diag(T)));

    offset = max( logJeach );    
    logJmc = log( mean( exp( ( logJeach - offset ) ) ) ) + offset;
    log_c  = logJmc + logC;

    return( log_c )

}



# function log_c = log_GWishart_ud_const_mc(b,D,adj,N)
# % log normalizing constant approximation using importance sampling for a
# % G-Wishart(b,D,adj):
# %
# %  p(K) = I(b,D)^{-1} |K|^{(d-2)/2} exp(-trace(K D)/2) using Monte Carlo sample of size N;
# %
# % Reference: Atay-Kayis and Massam 2005 Biometrika
# %
# % INPUT:   b,D,adj
# % OUTPUT:  log(I(b,D))


# T = chol(inv(D));

# T1 = T*diag(1./diag(T));

# p=size(D,1);

# indmx= reshape([1:p^2],p,p); % Uppermatrix index
# upperind=indmx(triu(indmx,1)>0);      

# A = zeros(p);
# A(upperind) = adj(upperind);
# nu = sum(A,2);


    

# logJeach = zeros(1,N);


# %%% D is non-diagnoal then computational cost O(p^4)
# if sum(abs(T1(:))) - sum(diag(T1)) ~= 0    
    
# for iter = 1:N
   

#     Psi = diag(sqrt(chi2rnd(b+nu)));
#     Psi(find(A==1))  = randn(sum(nu),1);
    
#     logJeach(iter) = 0 ;
# %%%%%% if i = 1 %%%%%%%%%%%%
# for j=2:p
#     if(A(1,j)==0)
#     Psi(1,j) = -sum(Psi(1,1:j-1)'.*T1(1:j-1,j));   
#     logJeach(iter) = logJeach(iter) - Psi(1,j)^2/2;
#     end
# end


# %%%%% if i>1 %%%%%%%%%%%%%%%%%%%
# for i=2:p
#     for j=i+1:p
  
#     if(A(i,j)==0)
#     Psi(i,j) = -sum(Psi(i,i:j-1)'.*T1(i:j-1,j));
   
     
#     for r = 1: i-1
       
#     Psi(i,j) = Psi(i,j) - 1/Psi(i,i)*(Psi(r,i)+ sum(Psi(r,r:i-1).*T1(r:i-1,i)'))* ...
#         (Psi(r,j)+ sum(Psi(r,r:j-1).*T1(r:j-1,j)'));
       
#     end
#     logJeach(iter) = logJeach(iter) - Psi(i,j)^2/2;   
#     end
#     end
# end

# end % for iter = 1:N


# %%%%%%%%%%%%%%%%%%%%%%%%%%  ELSE IF D IS diagonal MATRIX then simplified
# %%%%%%%%%%%%%%%%%%%%%%%%%%  calculation at O(p^2)
# else  %%%%%%%%

# for iter = 1:N
   

#     Psi = diag(sqrt(chi2rnd(b+nu)));
#     Psi(find(A==1))  = randn(sum(nu),1);
    
#     logJeach(iter) = 0 ;


# %%%%% if i>1 %%%%%%%%%%%%%%%%%%%
# for i=2:p
#     for j=i+1:p
  
#     if(A(i,j)==0)       
#     Psi(i,j) =  - 1/Psi(i,i)* sum(Psi(1:i-1,i).*Psi(1:i-1,j));   
#     logJeach(iter) = logJeach(iter) - Psi(i,j)^2/2;   
   
#     if isnan(logJeach(iter))
#         error('NA')
#     end
   
#     end
   
   
   
#     end
# end

# end % for iter = 1:N


# end % END IF T is diagonal
 

#      bi = sum(adj,2);

#      logC =sum(nu/2*log(2*pi)+(b+nu)/2*log(2)+gammaln((b+nu)/2)+(b+bi-1).*log(diag(T)));
         
#      offset = max(logJeach);    
#      logJmc = log(mean(exp((logJeach-offset))))+offset;
    
    
#      log_c = logJmc + logC;

