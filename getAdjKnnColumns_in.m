function Adjknn = getAdjKnnColumns_in( SimMatrix,  k_neighbors_vec , symmetrized, keepweight )
% Input: 
% SimMatrix  similarity matrix   
% k_neighbors  vector with number of neighbors of each node   
% symmetrized  1 or 0  
% keepweight   keep similarity as weight of edges 
% Output: Adjknn   matrix  
    if isempty(SimMatrix) 
        % sort_dim = 2;
        SimMatrix =rand(10);   % for testing only 
        warning('test test test ');
    end 
    SZ =size( SimMatrix); 
    % 
    if isscalar(k_neighbors_vec)%确定是否输入为标量
        k_neighbors_vec = k_neighbors_vec(1)*ones( SZ(1), 1 ); 
    end
    if isempty(symmetrized) 
        symmetrized = true;
    end
    if isempty(keepweight) 
        keepweight = false;
    end
     
    if any( k_neighbors_vec>SZ(1)-1 )
        k_neighbors_vec(  k_neighbors_vec>SZ(1)   ) = SZ(1)-1;
        warning( ['There is k_neighbors:','>', num2str(SZ(1)-1),' the maximal number of neighbors'] );
    end
    
	% %     
    SimMatrix(   sub2ind( SZ, 1:SZ(1),1:SZ(1) )      ) = -inf;   %   
    % SimMatrix(   ( eye( SZ ) )==1       ) = -inf; 
    [~,II] = sort( SimMatrix ,2, 'descend' );  
    Adjknn = zeros( SZ );
    for ii=1: SZ(1)
        knn = II(ii,1: k_neighbors_vec(ii) ); 
        if keepweight
            Adjknn(ii,  knn ) = SimMatrix(ii,  knn );    
        else
            Adjknn(ii,  knn ) = 1; 
        end
    end 
    Adjknn(sub2ind( SZ, 1:SZ(1),1:SZ(1) )) = 0; 
    if symmetrized
        [i,j,v] = find( Adjknn ); 
        ind = sub2ind( SZ, i ,j );
        Adjknn = Adjknn' ; 
        Adjknn(ind) = v; 
        % %     Adjknn = Adjknn';
        % %     Adjknn(ind) = 
        % Adjknn = full
    end 
end
