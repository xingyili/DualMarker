%% 5-fold CV - GSE2034, Top-100/150/200
% Input: expression matrix GSE2034, ranked genes GSE2034geneRANK
GSE2034geneRANK = geneRANK;

%% Buffers
GSE2034_aucScore_list_100 = [];
GSE2034_aucScore_list_150 = [];
GSE2034_aucScore_list_200 = [];

%% Data prep
raw_data = double(GSE2034(2:end,2:end)');  % samples x genes
raw_data = zscore(raw_data);
raw_data = mapstd(raw_data);               % row-wise standardization
labels   = double(GSE2034(1,2:end)');      % binary labels

%% 50 repeats
repeat = 50;
rng(2025);
for rep = 1:repeat
    fprintf('[Repeat %d/%d]\n', rep, repeat);
    indices     = crossvalind('Kfold', labels, 5);
    auc_100_list = zeros(1,5);
    auc_150_list = zeros(1,5);
    auc_200_list = zeros(1,5);

    for fold = 1:5
        %% split
        is_test  = (indices == fold);
        is_train = ~is_test;

        train_data  = raw_data(is_train,:);
        train_label = labels(is_train);
        test_data   = raw_data(is_test,:);
        test_label  = labels(is_test);

        GSE2034_train = [train_label'; train_data'];  % columns = samples
        GSE2034_test  = [test_label';  test_data'];

        %% ---- feature selection ----
        gene_names = string(GSE2034(2:end,1));

        % indices of top-k genes in the expression matrix (+1 for the label row)
        select_features = @(k) find(ismember(gene_names, string(GSE2034geneRANK(1:k,1))));
        idx_100 = select_features(100) + 1;
        idx_150 = select_features(150) + 1;
        idx_200 = select_features(200) + 1;

        train_100 = double(GSE2034_train(idx_100,:));
        train_150 = double(GSE2034_train(idx_150,:));
        train_200 = double(GSE2034_train(idx_200,:));

        test_100 = double(GSE2034_test(idx_100,:));
        test_150 = double(GSE2034_test(idx_150,:));
        test_200 = double(GSE2034_test(idx_200,:));

        train_Y = GSE2034_train(1,:)';
        test_Y  = GSE2034_test(1,:)';

        %% ---- train & predict ----
        % Top-100
        model100      = classRF_train(train_100', train_Y, 500);
        [~, vote100]  = classRF_predict(test_100', model100);
        score100      = vote100(:,2) ./ sum(vote100, 2);
        auc_100_list(fold) = AUC(test_Y, score100');

        % Top-150
        model150      = classRF_train(train_150', train_Y, 500);
        [~, vote150]  = classRF_predict(test_150', model150);
        score150      = vote150(:,2) ./ sum(vote150, 2);
        auc_150_list(fold) = AUC(test_Y, score150');

        % Top-200
        model200      = classRF_train(train_200', train_Y, 500);
        [~, vote200]  = classRF_predict(test_200', model200);
        score200      = vote200(:,2) ./ sum(vote200, 2);
        auc_200_list(fold) = AUC(test_Y, score200');
    end

    % mean AUC across 5 folds for this repeat
    GSE2034_aucScore_list_100(end+1) = mean(auc_100_list);
    GSE2034_aucScore_list_150(end+1) = mean(auc_150_list);
    GSE2034_aucScore_list_200(end+1) = mean(auc_200_list);

    fprintf('[Repeat %d/%d] AUC@100=%.4f | AUC@150=%.4f | AUC@200=%.4f\n', ...
        rep, repeat, mean(auc_100_list), mean(auc_150_list), mean(auc_200_list));
end

%% Final means
fprintf('\nFinal Mean AUC Values:\n');
fprintf('Top-100: %.4f\n', mean(GSE2034_aucScore_list_100));
fprintf('Top-150: %.4f\n', mean(GSE2034_aucScore_list_150));
fprintf('Top-200: %.4f\n', mean(GSE2034_aucScore_list_200));
