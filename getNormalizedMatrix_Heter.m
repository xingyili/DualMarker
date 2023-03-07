function [ combMatrix , IsNormalized ] = getNormalizedMatrix_Heter(AdjGfG,AdjGfD,AdjDfD, pro_jump,  NormalizationType, isdebug) 
% % % % % % % % % % % % % % % % % % % % 
if ~exist('AdjGfG','var')|| isempty (AdjGfG)
    N_gene=30;
    N_disease = 10; 
    AdjGfG = rand(N_gene,N_gene)>0.2 ; 
    AdjGfD = rand(N_gene,N_disease)>0.2 ; 
    AdjDfD = rand(N_disease,N_disease)>0.2 ; 
    P0_G = zeros(N_gene,1); P0_G(1:3)=1; P0_G = P0_G./sum(P0_G);
    P0_D = zeros(N_disease,1); P0_D(1:2)=1; P0_D = P0_D./sum(P0_D);
    % 
    pro_jump = 0.5;  
    warning('TestTestTestTestTestTestTestTestTestTestTestTestTestTestTest'); 
    isdebug =true; 
    istest =true; 
    NormalizationType = 'ProbabilityNormalizationColumn'; %%for 'Random Walk' RWR, RWRH  RWRM  RWRMH and more   
%     NormalizationType = 'ProbabilityNormalizationRow'; %%for label propagation    
    % NormalizationType = 'LaplacianNormalization'; %%  for label propagation, prince and more....    
%     NormalizationType = 'Weight'; %% Weighting  ....    
%     NormalizationType = 'None'; %%  without normalization ....    
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Input % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% AdjGfG: associatins from (f) genes (G) to Genes (G)  
% AdjGfD: associatins from Diseases (D) to Genes (G) GfD
% AdjDfD  associatins from Diseases (D) to Disease (G)  
% pro_jump: jumping Probability from first layer to second layer or weighting the effect of second layer on the first layer.    
% NormalizationType = 'ProbabilityNormalizationColumn'; %%for 'Random Walk' RWR, RWRH  RWRM  RWRMH and more   
% NormalizationType = 'ProbabilityNormalizationRow'; %%for label propagation    
% NormalizationType = 'LaplacianNormalization'; %%  for label propagation, prince and more....    
% NormalizationType = 'Weight'; %% Weighting  ....    
% NormalizationType = 'None'; %%  without normalization ....    
% Ouput % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% combMatrix is matrix after normalization.  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    % global   Global_Var_RWRH 
    if ~exist('pro_jump','var') || isempty (pro_jump)
        pro_jump = 0.5; 
    elseif pro_jump>1 || pro_jump <0
        error('pro_jump is wrong. it should be between 0 and 1');
    end      
% %     if ~exist('NormalizationType','var') || isempty (NormalizationType)
% %         NormalizationType = 'ProbabilityNormalizationColumn' ; 
% %     elseif ~ismember(NormalizationType,{ 'LaplacianNormalization', 'column','col',  'ProbabilityNormalizationColumn','ProbabilityNormalizationCol'     } )
% %         error(['NormalizationType is wrong: ',char(string(NormalizationType)) ]);
% %     end   
   if ~exist('isdebug','var') || isempty (isdebug)
        isdebug = false;  
    end   
    %  
    [N_gene, N_disease] = size( AdjGfD );
    if isempty( AdjDfD )
       AdjDfD = speye(N_disease);  
       warning('AdjDfD is empty.');
    end 
    %
    if ~exist('NormalizationType','var') || isempty(NormalizationType)
        NormalizationType = 'None'; 
    end
    %  
    IsNormalized = true; 
    switch lower( NormalizationType )
        case lower( {'None'} )
            combMatrix = [ AdjGfG, AdjGfD; AdjGfD', AdjDfD    ] ;
            IsNormalized = false;
            
        case lower( {'Weight'} )  
            combMatrix = [ (1-pro_jump).*AdjGfG, pro_jump.*AdjGfD; pro_jump.*AdjGfD', (1-pro_jump).*AdjDfD    ] ;
            IsNormalized = false;
            
        case lower( {'col','ProbabilityNormalizationColumn','NormalizationColumn', 'Column'} )   %姒傜巼瑙ｉ噴  锛岀‘淇濆垪鍜屼负1 
            idxDis_WithDiseaseGene =  sum( AdjGfD, 1)~=0;   % mark diseases with disease-genes
            idxGene_WithDisease    = (sum( AdjGfD, 2)~=0)';   % mark genes that are associated with diseases
            % WAdj = getNormalizedMatrix(Adj, NormalizationType, SetIsolatedNodeSelfLoop )
            M_GfG = getNormalizedMatrix(AdjGfG   , NormalizationType, true  ); 
            M_DfD = getNormalizedMatrix(AdjDfD   , NormalizationType, true  ); 
            M_GfD = getNormalizedMatrix(AdjGfD   , NormalizationType, false );  % probabilities from disease space to gene space 从疾病空间到基因空间的概率
            M_DfG = getNormalizedMatrix(AdjGfD'  , NormalizationType, false );  % probabilities from gene space to disease space
            %
            M_GfG(:,idxGene_WithDisease)       = (1-pro_jump).*M_GfG(:,idxGene_WithDisease); 
            M_DfD(:,idxDis_WithDiseaseGene )   = (1-pro_jump).*M_DfD(:,idxDis_WithDiseaseGene ) ; 
            M_GfD                           = pro_jump.*M_GfD; % Disease-columns without disease-genes is all zeros. So no use idxDis_WithDiseaseGene 没有疾病基因的疾病列都是零。所以不用idxDis_WithDiseaseGene
            M_DfG                           = pro_jump.*M_DfG; % Gene-columns without diseases is all zeros. So no use idxGene_WithDisease 基因行没有疾病的都是零
            %
            combMatrix = [ M_GfG, M_GfD; M_DfG, M_DfD    ] ; 
            
        case lower( {'row','ProbabilityNormalizationRow','NormalizationRow'} )
            idxDis_WithDiseaseGene = (sum( AdjGfD, 1)~=0);   % mark diseases with disease-genes
            idxGene_WithDisease    = (sum( AdjGfD, 2)~=0);   % mark genes that are associated with diseases
            % WAdj = getNormalizedMatrix(Adj, NormalizationType, SetIsolatedNodeSelfLoop )
            M_GfG = getNormalizedMatrix(AdjGfG   , NormalizationType, true ); 
            M_DfD = getNormalizedMatrix(AdjDfD   , NormalizationType, true ); 
            M_GfD = getNormalizedMatrix(AdjGfD   , NormalizationType, false );  % probabilities from disease space to gene space 
            M_DfG = getNormalizedMatrix(AdjGfD'  , NormalizationType, false );  % probabilities from gene space to disease space
            %
            M_GfG(idxGene_WithDisease,:)       = (1-pro_jump).*M_GfG(idxGene_WithDisease,:) ; 
            M_DfD(idxDis_WithDiseaseGene,: )   = (1-pro_jump).*M_DfD(idxDis_WithDiseaseGene,: ) ; 
            M_GfD                           = pro_jump.*M_GfD; % Disease-columns without disease-genes is all zeros. So no use idxDis_WithDiseaseGene
            M_DfG                           = pro_jump.*M_DfG; % Gene-columns without diseases is all zeros. So no use idxGene_WithDisease
            %
            combMatrix = [ M_GfG, M_GfD; M_DfG, M_DfD    ] ;            
            
        case lower( {'LaplacianNormalization'} )
            % WAdj = getNormalizedMatrix(Adj, NormalizationType, SetIsolatedNodeSelfLoop )
            M_GfG = getNormalizedMatrix(AdjGfG   , NormalizationType, false );   
            M_DfD = getNormalizedMatrix(AdjDfD   , NormalizationType, false ); 
            M_GfD = getNormalizedMatrix(AdjGfD   , NormalizationType, false );  % probabilities from disease space to gene space 
            M_DfG = getNormalizedMatrix(AdjGfD'  , NormalizationType, false );  % probabilities from gene space to disease space
            %
            combMatrix = [ (1-pro_jump).*M_GfG, pro_jump.*M_GfD; pro_jump.*M_DfG, (1-pro_jump).*M_DfD    ] ;            
            
        otherwise
            error('No definition.');
    end
end