function WAdj = getNormalizedMatrix(Adj, NormalizationType, SetIsolatedNodeSelfLoop )
    if ~exist('Adj','var')  
        Adj =rand(5); dim=1;SetIsolatedNodeSelfLoop = true;  
        NormalizationType = 'col' ;
        % NormalizationType = 'laplacian normalization' ;
        istest = 1; 
        warning('Test Test Test Test Test Test Test ');
    end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% Adj  adjecent matrix  
% % NormalizationType:  
% % 'probability normalization' for 'Random Walk' RWR, RWRH  RWRM  RWRMH and more   
% % 'laplacian normalization'for prince and more....
% SetIsolatedNodeSelfLoop    set isolated node
% >= Matlab 2016
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
%     if ~issparse(Adj)
%         Adj = sparse( Adj );
%     end   
    if ischar(NormalizationType)
    %         NormalizationType =  (NormalizationType);
        switch  lower( NormalizationType )
            case lower( { 'column','col',  ...
                    'ProbabilityNormalizationColumn','ProbabilityNormalizationCol',...
                    'ProbabilityColumnNormalization','ProbabilityColNormalization',...
                    'NormalizationColumn','NormalizationCol' , ...
                    'ColumnNormalization','ColNormalization'   })
                NormalizationName = 'ProbabilityNormalization' ;  %  'Random Walk'  
                dim =1;
            case lower({ 'row' ,'ProbabilityNormalizationRow' ,'NormalizationRow' ,'ProbabilityRowNormalization' ,'RowNormalization'   })
                NormalizationName = 'ProbabilityNormalization' ;  %  'Random Walk'  
                dim =2;
            case lower('LaplacianNormalization')
                NormalizationName = NormalizationType; 
            case lower('LaplacianNormalizationMeanDegree')
                NormalizationName = NormalizationType; 
            case lower('ColNorm2')
                NormalizationName = NormalizationType; 
            case lower('RowNorm2')
                NormalizationName = NormalizationType; 
            case lower({'none', 'None', 'NONE'})
                % NormalizationName = 'None'; 
                WAdj = Adj; 
                return; 
            otherwise
                error(['There is no type of normalization: ',char( string(NormalizationType) )] );
        end
        
    elseif isnumeric(  NormalizationType   ) 
        NormalizationName =  ( 'ProbabilityNormalization' ) ;  %  'Random Walk'  
        dim = NormalizationType; 
        
    elseif isempty( NormalizationType )
        WAdj = Adj; 
        return;  
        
    else; error('There is no defintion of NormalizationType')
    end 
    % NormalizationName = lower( NormalizationName );
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     %
    switch lower( NormalizationName )
        case lower( 'ProbabilityNormalization' )
            degrees = sum(Adj,dim);
            if any( degrees~=1)
                WAdj = Adj./ ( degrees+eps  );           
                % % WAdj = Adj./ repmat( degrees +eps,[size(Adj,1),1]); 
            else
                WAdj = Adj; 
            end
            % 
            if SetIsolatedNodeSelfLoop  && size(Adj,1)==size(Adj,2) 
                ii = find( ~degrees ); 
                idx = sub2ind( size(Adj), ii,ii ); 
                WAdj(idx) = 1;  % set to be 1 for isolated nodes, 
            end
            
        case lower( 'LaplacianNormalization')
            deg_rowvec = ( sum(Adj,1) ).^0.5;  
            deg_colvec = ( sum(Adj,2) ).^0.5;   
            WAdj = (Adj./(deg_colvec+eps))./(deg_rowvec+eps) ;    
            % 
            if SetIsolatedNodeSelfLoop && size(Adj,1)==size(Adj,2)
                ii = find( ~sum(Adj,2) ) ; 
                % size(  WAdj )
                % size(  Adj )
                WAdj( sub2ind( size(Adj), ii,ii ) ) = 1;  % set to be 1 for isolated nodes, 
            end
            
        case lower( 'LaplacianNormalizationMeanDegree')
            n_node = length( Adj );
            km = sum( Adj(:) )./ n_node;  
            WAdj = Adj./( (km.^0.5)*(km.^0.5)  +eps) ;    
            % 
            if SetIsolatedNodeSelfLoop  && size(Adj,1)==size(Adj,2)
                ii = find( ~sum(Adj,2) ); 
                WAdj( sub2ind( size(Adj), ii,ii ) ) = 1;  % set to be 1 for isolated nodes, 
            end
            
        case lower( {'ColNorm2'} )   
            WAdj = Adj./ ( sqrt(sum( Adj.^2 ,1 )) +eps ); 
            
        case lower( {'RowNorm2'} )    
            WAdj = Adj./ ( sqrt(sum( Adj.^2 ,2 )) +eps ); 
            
        case lower( {'None','none'} )
            WAdj = Adj;   
        otherwise
            error(['NormalizationName is wrong: ',char(string(NormalizationName) )   ]);
    end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    if exist('istest','var') && istest 
        WAdj(1:5,1:5)
        sum(WAdj,1) 
        sum(WAdj,2)
    end
end