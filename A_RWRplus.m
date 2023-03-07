function  [Pt, WAdj ]= A_RWRplus(Adj, r_restart, P0, N_max_iter, Eps_min_change, IsNormalized,  NormalizationType)
% A_RWRplus is a generalization of RWR algorithm.  
% Including various propagation algorihtms with initial regularization: classical RWR, Label propagation,and so on
% Including a Solver_IterationPropagation, which can be used directly when IsNormalized is TRUE.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% need getNormalizedMatrix  
% Input 
% Adj
% r_restart
% P0   It should be normalized, though it is neccesary. 
% N_max_iter
% Eps_min_change
% IsNormalized   Whether Adj has been normalized: True or False  
% NormalizationType: including two types of methods by different Normalization Types
% (1) random walk with restart {'ProbabilityNormalizationColumn','ProbabilityNormalizationCol','col','column'}
% (2) 'LaplacianNormalization'   similar to PRINCE(PLOS Computational Biology, 2010, 6: e1000641.)               
% it is equivalent to PRINCE if assigned extended P0 with % disease simialrity and logistic function
% (3)label propagation with memory restart { 'row' ,'ProbabilityNormalizationRow'} 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Ouput
% Pt  
% WAdj   normalized Adj 
    if ~exist('N_max_iter','var') || isempty(N_max_iter) || (isnumeric( N_max_iter) && N_max_iter<=1 ) 
        N_max_iter =100; 
    elseif ~isnumeric( N_max_iter)  
        error('N_max_iter should be isnumeric!!!!!') ;
    end
    %
    if ~exist('Eps_min_change','var') || isempty(Eps_min_change) 
        Eps_min_change =10^-6; 
    elseif isnumeric( Eps_min_change) && Eps_min_change>=1 
        warning('The Eps_min_change is nomenaning. Reset Eps_min_change to be 10^-6.'); 
        Eps_min_change =10^-6;  
    elseif ~isnumeric( Eps_min_change)  
        error('Eps_min_change should be isnumeric!!!!!') ;
    end
    
    if ~exist('IsNormalized','var') || isempty(IsNormalized) 
        IsNormalized = false;  % Adj has been normalized for fast run.   
    end
    
    if ~exist('NormalizationType','var') || isempty(NormalizationType) 
        NormalizationType = 'ProbabilityNormalizationColumn'; %%for 'Random Walk' RWR, RWRH  RWRM  RWRMH and more   
    end        
    % % % 鍙栨秷绋?枏鐭╅樀鐨勮浆鎹紝鍦ㄨ皟鐢ㄧ殑澶栭儴鏍规嵁鎯呭喌鎸囧畾鐭╅樀褰㈠紡,闄ら潪鏋佺鎯呭喌 % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
    P0  = full(P0); 
    % 
    % AdjIsSparse = isparse(Adj); 
    AdjDensity =nnz(Adj)/numel(Adj); 
    if  (size(P0,2)==1 && AdjDensity<0.3 ) || (size(P0,2)>1 && AdjDensity<0.05 )
        Adj = sparse(Adj);          
    elseif (size(P0,2)==1 && AdjDensity>0.3 ) || (size(P0,2)>1 && AdjDensity>0.05 )
        Adj = full(Adj); 
    else 
        % no operation 
    end
    %
    if IsNormalized 
        WAdj = Adj; 
    else
        % WAdj = getNormalizedMatrix(Adj, 'col', true );
        switch NormalizationType
            case {'ProbabilityNormalizationColumn','ProbabilityNormalizationCol','col','column'}
                % random walk with restart
                WAdj = getNormalizedMatrix(Adj, 'col', true );
                %%%P0 = P0./(sum(P0,1)+eps);    % total probability is 1. 
                
            case 'LaplacianNormalization'  
                % propagation similar to PRINCE(PLOS Computational Biology, 2010, 6: e1000641.)
                % it is equivalent to PRINCE if assigned extended P0 with
                % disease simialrity and logistic function
                % A_PRINCEplus is better. 
                WAdj = getNormalizedMatrix(Adj, 'LaplacianNormalization', true );   
                
            case { 'row' ,'ProbabilityNormalizationRow'} 
                % label propagation with memory restart  
                WAdj = getNormalizedMatrix(Adj, 'row', true );  
                
            otherwise
                error(['NormalizationType is wrong: ',char( string(NormalizationType) )]); 
        end        
    end   
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    % % Solver_IterationPropagation
    % % It can be used directly when IsNormalized is TRUE.  
    Pt = P0;
    for T = 1: N_max_iter
        Pt1 = (1-r_restart)*WAdj*Pt + r_restart*P0;
        if all( sum( abs( Pt1-Pt )) < Eps_min_change )
            break;
        end
        Pt = Pt1;
    end
    Pt = full(Pt); 
end