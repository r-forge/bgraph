#   MAXIMALCLIQUES Find maximal cliques using the Bron-Kerbosch algorithm
#   Given a graphs boolean adjacency matrix, A, find all maximal cliques 
#   on A using the Bron-Kerbosch algorithm in a recursive manner.  The 
#   graph is required to be undirected and must contain no self-edges.
#
#   This function can be used to compute the maximal independent sets 
#   (maximal matchings) of a graph by providing the adjacency matrix of the 
#   corresponding conflict graph (complement of the conflict graph).
#
#   V_STR is an optional input string with the version of the Bron-Kerbosch 
#   algorithm to be used (either 'v1' or 'v2').  Version 2 is faster (and 
#   default), and version 1 is included for posterity.
#
#   MC is the output matrix that contains the maximal cliques in its 
#   columns.
#
#   Ref: Bron, Coen and Kerbosch, Joep, "Algorithm 457: finding all cliques
#   of an undirected graph", Communications of the ACM, vol. 16, no. 9, 
#   pp: 575–577, September 1973.
#
#   Ref: Cazals, F. and Karande, C., "A note on the problem of reporting 
#   maximal cliques", Theoretical Computer Science (Elsevier), vol. 407,
#   no. 1-3, pp: 564-568, November 2008.
#
#   Jeffrey Wildman (c) 2011
#   jeffrey.wildman@gmail.com
maximalCliquesBK<-function( A, v_str ){
# first, some input checking

if( dim(A)[1] != dim(A)[2]){
    stop(error='MATLAB:maximalCliques', 'Adjacency matrix is not square.');
} else if ( !all(all( (A==1) | (A==0)))){
    stop(error='MATLAB:maximalCliques', 'Adjacency matrix is not boolean (zero-one valued).')
} else if ( !all(all(A==t(A)))){
    stop(error='MATLAB:maximalCliques', 'Adjacency matrix is not undirected (symmetric).')
} else if ( sum(diag(abs(A))) != 0){
    stop(error='MATLAB:maximalCliques', 'Adjacency matrix contains self-edges (check your diagonal).');
}#end

if( missing(v_str)) #if( !exist('v_str','var'))
    v_str = 'v2';
#end

if( !(v_str=='v1') && !(v_str=='v2')){
    warning('MATLAB:maximalCliques', 'Version not recognized, defaulting to v2.');
    v_str = 'v2';
}#end

# second, set up some variables

n  = dim(A)[2];      # number of vertices
MC = NULL;             # storage for maximal cliques
R  = NULL;             # currently growing clique
P  = 1:n;            # prospective nodes connected to all nodes in R
X  = NULL;             # nodes already processed


# third, run the algorithm!

if( v_str=='v1'){
    MC <- BKv1(R,P,X,MC,A,n);
}else{
    MC <- BKv2(R,P,X,MC,A,n);
} # end
print(MC)
# end # maximalCliques
return(MC)
}


    # version 1 of the Bron-Kerbosch algo 
    BKv1 <- function( R, P, X, MC, A, n){
            
            if( missing(P) && missing(X) ){
                # report R as a maximal clique
                newMC    = matrix(0,1,n);
                newMC[R] = 1;                   # newMC contains ones at indices equal to the values in R   
                MC       = c(MC, newMC); #MC = [MC newMC.'];
            }else{
                for( u in P ){
                    P    = setdiff(union(P,u),intersect(P,u)) #setxor(P,u);
                    Rnew = c(R, u);
                    Nu   = which(A[u,]!=0);
                    Pnew = intersect(P,Nu);
                    Xnew = intersect(X,Nu);
                    MC   = BKv1(Rnew, Pnew, Xnew, MC, A, n );
                    X    = c(X, u);
                } #end
            }#end
        return(MC)
    }#end # BKv1


    # version 2 of the Bron-Kerbosch algo
    BKv2 <- function( R, P, X, MC, A, n ){
    
        if (missing(P) && missing(X)){
            # report R as a maximal clique
            newMC = matrix(0,1,n);
            newMC[R] = 1;                   # newMC contains ones at indices equal to the values in R   
            MC = c(MC, newMC);
        }else{
            # choose pivot
            ppivots = union(P,X);           # potential pivots
            binP    = matrix(0,1,n);
            binP[P] = 1;                    # binP contains ones at indices equal to the values in P          
            # rows of A(ppivots,:) contain ones at the neighbors of ppivots
            pcounts = A[ppivots, ]%*%t(binP);     # cardinalities of the sets of neighbors of each ppivots intersected with P
            ind     = which.max(pcounts);       #[~,ind] = max(pcounts);
            u_p     = ppivots[ind];             # select one of the ppivots with the largest count
            
            for( u in intersect(which(A[u_p,]==0),P)){   # all prospective nodes who are not neighbors of the pivot
                P = setdiff(union(P,u),intersect(P,u)) #P = setxor(P,u);
                Rnew = c(R, u);
                Nu   = which(A[u,]!=0);
                Pnew = intersect(P,Nu);
                Xnew = intersect(X,Nu);
                MC   = BKv2(Rnew, Pnew, Xnew, MC, A, n);
                X = c(X, u);
            }#end
        }#end
        return(MC)
    } #end # BKv2

# function [ MC ] = maximalCliquesBK(A,v_str)
# %MAXIMALCLIQUES Find maximal cliques using the Bron-Kerbosch algorithm
# %   Given a graph's boolean adjacency matrix, A, find all maximal cliques 
# %   on A using the Bron-Kerbosch algorithm in a recursive manner.  The 
# %   graph is required to be undirected and must contain no self-edges.
# %
# %   This function can be used to compute the maximal independent sets 
# %   (maximal matchings) of a graph by providing the adjacency matrix of the 
# %   corresponding conflict graph (complement of the conflict graph).
# %
# %   V_STR is an optional input string with the version of the Bron-Kerbosch 
# %   algorithm to be used (either 'v1' or 'v2').  Version 2 is faster (and 
# %   default), and version 1 is included for posterity.
# %
# %   MC is the output matrix that contains the maximal cliques in its 
# %   columns.
# %
# %   Ref: Bron, Coen and Kerbosch, Joep, "Algorithm 457: finding all cliques
# %   of an undirected graph", Communications of the ACM, vol. 16, no. 9, 
# %   pp: 575–577, September 1973.
# %
# %   Ref: Cazals, F. and Karande, C., "A note on the problem of reporting 
# %   maximal cliques", Theoretical Computer Science (Elsevier), vol. 407,
# %   no. 1-3, pp: 564-568, November 2008.
# %
# %   Jeffrey Wildman (c) 2011
# %   jeffrey.wildman@gmail.com


# % first, some input checking

# if size(A,1) ~= size(A,2)
#     error('MATLAB:maximalCliques', 'Adjacency matrix is not square.');
# elseif ~all(all((A==1) | (A==0)))
#     error('MATLAB:maximalCliques', 'Adjacency matrix is not boolean (zero-one valued).')
# elseif ~all(all(A==A.'))
#     error('MATLAB:maximalCliques', 'Adjacency matrix is not undirected (symmetric).')
# elseif trace(abs(A)) ~= 0
#     error('MATLAB:maximalCliques', 'Adjacency matrix contains self-edges (check your diagonal).');
# end
    
# if ~exist('v_str','var')
#     v_str = 'v2';
# end

# if ~strcmp(v_str,'v1') && ~strcmp(v_str,'v2')
#     warning('MATLAB:maximalCliques', 'Version not recognized, defaulting to v2.');
#     v_str = 'v2';
# end


# % second, set up some variables

# n = size(A,2);      % number of vertices
# MC = [];            % storage for maximal cliques
# R = [];             % currently growing clique
# P = 1:n;            % prospective nodes connected to all nodes in R
# X = [];             % nodes already processed


# % third, run the algorithm!

# if strcmp(v_str,'v1')
#     BKv1(R,P,X);
# else
#     BKv2(R,P,X);
# end
    

#     % version 1 of the Bron-Kerbosch algo 
#     function [] = BKv1 ( R, P, X )
        
#         if isempty(P) && isempty(X)
#             % report R as a maximal clique
#             newMC = zeros(1,n);
#             newMC(R) = 1;                   % newMC contains ones at indices equal to the values in R   
#             MC = [MC newMC.'];
#         else
#             for u = P
#                 P = setxor(P,u);
#                 Rnew = [R u];
#                 Nu = find(A(u,:));
#                 Pnew = intersect(P,Nu);
#                 Xnew = intersect(X,Nu);
#                 BKv1(Rnew, Pnew, Xnew);
#                 X = [X u];
#             end
#         end
        
#     end % BKv1


#     % version 2 of the Bron-Kerbosch algo
#     function [] = BKv2 ( R, P, X )

#         if (isempty(P) && isempty(X))
#             % report R as a maximal clique
#             newMC = zeros(1,n);
#             newMC(R) = 1;                   % newMC contains ones at indices equal to the values in R   
#             MC = [MC newMC.'];
#         else
#             % choose pivot
#             ppivots = union(P,X);           % potential pivots
#             binP = zeros(1,n);
#             binP(P) = 1;                    % binP contains ones at indices equal to the values in P          
#             % rows of A(ppivots,:) contain ones at the neighbors of ppivots
#             pcounts = A(ppivots,:)*binP.';  % cardinalities of the sets of neighbors of each ppivots intersected with P
#             [~,ind] = max(pcounts);
#             u_p = ppivots(ind);             % select one of the ppivots with the largest count
            
#             for u = intersect(find(~A(u_p,:)),P)   % all prospective nodes who are not neighbors of the pivot
#                 P = setxor(P,u);
#                 Rnew = [R u];
#                 Nu = find(A(u,:));
#                 Pnew = intersect(P,Nu);
#                 Xnew = intersect(X,Nu);
#                 BKv2(Rnew, Pnew, Xnew);
#                 X = [X u];
#             end
#         end
        
#     end % BKv2
       

# end % maximalCliques

