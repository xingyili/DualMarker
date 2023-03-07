%Main function of DualMarker
%Input:
%AdjGfG:The adjacency matrix of gene-gene
%AdjGfP:The adjacency matrix of gene-protein
%AdjGfD:The adjacency matrix of gene-Disease
%AdjGfGO:The adjacency matrix of gene-GO
%GSE2034:The gene expression matrix of GSE2034
%Genelist:genelist
%breastGenediease:Scores of breast cancer-related genes downloaded from DisGeNet database
%Output:
%geneRANK:The ranking of biomarkers

%Adjacency matrix of heterogeneous network
Aunion = sparse( ...
    [AdjGfG,  AdjGfP,  AdjGfD,  AdjGfGO; ...
    AdjGfP', AdjPfP,  AdjPfD,  AdjPfGO; ...
    AdjGfD', AdjPfD', AdjDfD,  AdjDfGO; ...
    AdjGfGO',AdjPfGO',AdjDfGO',AdjGOfGO  ]  );

ind_gene = [1: length(AdjGfG)];
ind_protein  = [length(AdjGfG)+1: 2*length(AdjGfG)];

paras.d       = 128 ;
paras.Ortho   = 1  ;
paras.seed    = 0;

worktype = 'classification';

features = getRandNEemb_in(Aunion,paras, [], [], worktype);
Aunion = [];


features_gene = features(ind_gene,2:end);
features_pro  = features(ind_protein,2:end);
features =[];

% Construct the dual-layer heterogeneous network
AdjGfG = corr(features_gene',features_gene');
AdjPfP = corr(features_pro',features_pro');

pro_jump = 0.7;
NormalizationType = 'ProbabilityNormalizationColumn';
[ M_cos , IsNormalized ] = getNormalizedMatrix_Heter(AdjGfG,AdjGfP,AdjPfP, pro_jump, NormalizationType, []) ; 
M_cos=Network_Enhancement(M_cos);

%% Obtain the breast cancer gene score and gene expression data T-test score
%Score of gene disease associations
BreastGeneDiease = "";
a=0;
for i=1:size(breastGenediease,1)
    if  ismember(breastGenediease(i,1),Genelist )
        a=a+1;
        BreastGeneDiease(a,1) = breastGenediease(i,1);
        BreastGeneDiease(a,2) = breastGenediease(i,2);
        BreastGeneDiease(a,3) = breastGenediease(i,3);
    end
end

Breastgenes = unique(BreastGeneDiease(:, 1));
[~, gene_idx] = ismember(BreastGeneDiease(:, 1), Breastgenes);
Breastscore = accumarray(gene_idx, str2double(BreastGeneDiease(:, 3)), [], @sum);
BreastGenescores = [Breastgenes, num2cell(Breastscore)];

col2 = str2double(BreastGenescores(:, 2));
col2_normalized = (col2 - min(col2)) / (max(col2) - min(col2));
BreastGenescores(:, 2) = col2_normalized;



% T-test
Genes_expression_value_cluster_EXP=[];Label_cluster_EXP=[];
Genes_expression_value_cluster_EXP = string(GSE2034(2:end,2:end));
Label_cluster_EXP = string((GSE2034(1,2:end)));
UniqueNode = Genelist;
Genes_ID = GSE2034(2:end,1);

StatisticScore_exp=[];
[StatisticScore_exp,stats] = Obtain_StatisticScore(Genes_expression_value_cluster_EXP, Label_cluster_EXP, UniqueNode,Genes_ID);
StatisticScore_exp_train_0=[];
a=0;
for i=1:length(AdjGfG)
    if  ismember(StatisticScore_exp(i,3),1 )
        a=a+1;
        StatisticScore_exp_train_0(a,1) = StatisticScore_exp(i,1);
        StatisticScore_exp_train_0(a,2) = StatisticScore_exp(i,2);
        StatisticScore_exp_train_0(a,3) = StatisticScore_exp(i,3);
        StatisticScore_exp_train_0(a,4) = StatisticScore_exp(i,4);
    end
end
StatisticScore_exp_0_nor=[];
StatisticScore_exp_0_nor(:,1) = StatisticScore_exp_train_0(:,1);
StatisticScore_exp_0_nor(:,2) = mapminmax(StatisticScore_exp_train_0(:,2)',0,1)';

%Combin BreastGenescores and expGeneSocre
Combine = vertcat(BreastGenescores, StatisticScore_exp_0_nor);

BRCAgenes = unique(Combine(:, 1));
[~, gene_idx2] = ismember(Combine(:, 1), BRCAgenes);
Breastscore2 = accumarray(gene_idx2, str2double(Combine(:, 2)), [], @sum);
Breastgenevalue = [BRCAgenes, num2cell(Breastscore2)];

col = str2double(Breastgenevalue(:, 2));
col_normalized = (col - min(col)) / (max(col) - min(col));
Breastgenevalue(:, 2) = col_normalized;



geneIniValue1 = zeros(size(UniqueNode,1),2);
Breastgenenode = Breastgenevalue(:,1);
geneIniValue1(:,1) = UniqueNode;
for k = 1:size(Breastgenenode,1)
    [~,~,kUniqueNode] = intersect(Breastgenenode(k,1),UniqueNode);
    geneIniValue1(kUniqueNode,2) = Breastgenevalue(k,2);
end

expIntiVALUE = geneIniValue1;


% Network propagation
if any( strcmpi(NormalizationType,{'col','ProbabilityNormalizationColumn','NormalizationColumn', 'Column'}) )
    eta = 0.5;
    P0=[];
    P0 = [ (1-eta)*expIntiVALUE(:,2); eta*expIntiVALUE(:,2)];
    P0 = P0./( sum(P0,1)+eps )  ;
else
    P0 = [ (1-eta)*expIntiVALUE(:,2); eta*expIntiVALUE(:,2)];
end
% [Pt, WAdj ]= A_RWRplus(Adj, r_restart, P0, N_max_iter, Eps_min_change, IsNormalized,  NormalizationType)
r_restart=0.7;
Pt= [];
Pt = A_RWRplus(M_cos, r_restart, P0 , [],[], IsNormalized);
P_G = Pt(1:length(AdjGfG),:);
P_P = Pt(length(AdjGfG)+1:end,:);
Pt = [P_G,P_P];
Pt(:,3) = max(Pt,[],2);

%Order the gene  after network propagation
geneRANK= [];
index = [];
[geneRANK(:,2),index]=sort( Pt(:,3),'descend');
l = [];
for l=1:length(expIntiVALUE(:,1))
    geneRANK(l,1) = expIntiVALUE(index(l),1);
end





