# modified based on original code by HAO WANG (?) written in matlab

makedecompgraph <- function( A, nodesIDs = NULL, nodenames = NULL ){

    Adj = full( A ); p = dim( Adj )[1];
    if ( is.null( nodenames ) ){ nodenames = as.character( 1:p ) }
    if ( !is.null( nodesIDs) ){ order = nodesIDs } else { order = 1:p }
    i = 1
    while( i < p ){
        a = max( colSums( Adj[1:i,(i + 1):p, drop=FALSE] ) );
        b = which ( colSums( Adj[1:i, ( i + 1 ):p, drop=FALSE] ) == a )[1]
        order[c( i + 1, b + i )] = order[c( b + i, i + 1 )];
        i=i+1
        Adj = full( A ); Adj=Adj[order,order];
    }
    numcliq = 1; 
    C = vector( "list", 1 )
    C[[1]]$ID = c( 1 ); 
    i=2;
    while( i <= p ){
        if( sum( Adj[i,C[[numcliq]]$ID] ) == length( C[[numcliq]]$ID ) ){
          C[[numcliq]]$ID = c( C[[numcliq]]$ID, i );
        } else {
            C[[numcliq]]$dim = length( C[[numcliq]]$ID );
            numcliq = numcliq + 1;
            C = append( C, vector( "list", 1 ) )
            C[[numcliq]]$ID = union( i, which( Adj[i, 1:i] == 1 ) );
        }
        i=i+1
    }
    C[[numcliq]]$dim = length( C[[numcliq]]$ID );
    for ( i in 1:numcliq ) {
        C[[i]]$ID = sort( order[C[[i]]$ID] );
        C[[i]]$names = nodenames[C[[i]]$ID];
    }
    # separators
    S = vector( "list", numcliq )
    UN = C[[1]]$ID; 
    S[[1]]=list( ID = c( ), dim = 0, names = c( ) )
    if( numcliq >1 ){
        for ( i in 2:numcliq ){
            S[[i]]$ID = intersect( UN, C[[i]]$ID )
            S[[i]]$dim = length( S[[i]]$ID )
            S[[i]]$names = nodenames[S[[i]]$ID] 
            UN = union( UN, C[[i]]$ID );
         }
    }
    C[[1]]$names = nodenames[C[[1]]$ID];  #S[[1]]$names = ; 
    G = vector( "list", 2 )
    G[[1]] = C; G[[2]] = S;
    names( G ) <- c("cliques","separators")
    class( G ) <- c("list", "dgraph")
    return( G )
}

full<- function( A ){
    A = A / A
    A[is.na( A )]<-0
    return( A )
}

#test
if(0){
setwd("D:\\workspace\\GraphMCMC")

a=read.csv("example.dat",header=F)
gg=makedecompgraph(a)
dhiw(a,g,10,diag(8))

AA <- diag(5)
AA[2,3]=AA[3,2]=1
AA[1,3]=AA[3,1]=1
AA[2,1]=AA[1,2]=1
AA[4,3]=AA[3,4]=1
AA[5,3]=AA[3,5]=1
AA[5,4]=AA[4,5]=1
bb=makedecompgraph(AA)
}