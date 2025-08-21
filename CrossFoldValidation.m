%% 5折交当验 - GSE2034, Top-50/100/150
% 输入: GSE2034表达矩阵 GSE2034, 排序后的基因 GSE2034geneRANK
GSE2034geneRANK = geneRANK;
%% 编码设置
GSE2034_aucScore_list_50 = [];
GSE2034_aucScore_list_100 = [];
GSE2034_aucScore_list_150 = [];

%% 数据处理
raw_data = double(GSE2034(2:end,2:end)');
raw_data = zscore(raw_data);
raw_data = mapstd(raw_data);  % 行对比规一化
labels = double(GSE2034(1,2:end)');

%% 50次重复
repeat=50;
rng(2025);
for xuanhuan = 1:repeat
    fprintf('[Repeat %d/%d]\n', xuanhuan, repeat);
    indices = crossvalind('Kfold', labels, 5);
    auc_50_list = zeros(1,5);
    auc_100_list = zeros(1,5);
    auc_150_list = zeros(1,5);

    for fold = 1:5
        %% 分割数据
        is_test = (indices == fold);
        is_train = ~is_test;

        train_data = raw_data(is_train,:);
        train_label = labels(is_train);
        test_data = raw_data(is_test,:);
        test_label = labels(is_test);

        GSE2034_train = [train_label'; train_data'];  % 列是样本
        GSE2034_test  = [test_label';  test_data'];

        %% ---- 特征选择 ----
        gene_names = string(GSE2034(2:end,1));

        % 特征 Top-k 选择函数
        select_features = @(k) find(ismember(gene_names, string(GSE2034geneRANK(1:k,1))));

        idx_50 = select_features(50) + 1;
        idx_100 = select_features(100) + 1;
        idx_150 = select_features(150) + 1;

        train_50  = double(GSE2034_train(idx_50,:));
        train_100 = double(GSE2034_train(idx_100,:));
        train_150 = double(GSE2034_train(idx_150,:));

        test_50  = double(GSE2034_test(idx_50,:));
        test_100 = double(GSE2034_test(idx_100,:));
        test_150 = double(GSE2034_test(idx_150,:));

        train_Y = GSE2034_train(1,:)';
        test_Y  = GSE2034_test(1,:)';

        %% ---- 训练 + 预测 ----
        % Top-50
        model50 = classRF_train(train_50', train_Y, 500);
        [~, vote50] = classRF_predict(test_50', model50);
        score50 = vote50(:,2) ./ sum(vote50, 2);
        auc_50_list(fold) = AUC(test_Y, score50');

        % Top-100
        model100 = classRF_train(train_100', train_Y, 500);
        [~, vote100] = classRF_predict(test_100', model100);
        score100 = vote100(:,2) ./ sum(vote100, 2);
        auc_100_list(fold) = AUC(test_Y, score100');

        % Top-150
        model150 = classRF_train(train_150', train_Y, 500);
        [~, vote150] = classRF_predict(test_150', model150);
        score150 = vote150(:,2) ./ sum(vote150, 2);
        auc_150_list(fold) = AUC(test_Y, score150');
    end

    % 每次 5折 平均 AUC
    GSE2034_aucScore_list_50(end+1) = mean(auc_50_list);
    GSE2034_aucScore_list_100(end+1) = mean(auc_100_list);
    GSE2034_aucScore_list_150(end+1) = mean(auc_150_list);

    fprintf('[Repeat %d/50] AUC@50=%.4f | AUC@100=%.4f | AUC@150=%.4f\n', ...
        xuanhuan, auc_50_list(fold), auc_100_list(fold), auc_150_list(fold));
end

%% 最终平均
fprintf('\nFinal Mean AUC Values:\n');
fprintf('Top-50:  %.4f\n', mean(GSE2034_aucScore_list_50));
fprintf('Top-100: %.4f\n', mean(GSE2034_aucScore_list_100));
fprintf('Top-150: %.4f\n', mean(GSE2034_aucScore_list_150));
