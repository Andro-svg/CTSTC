function idx = improved_kmeans(img, K, Size, lambda_local, lambda_global)
    % img: 输入图像矩阵
    % K: 聚类数目
    % lambda_local: 局部特征的权重
    % lambda_global: 非局部特征的权重
%     lambda_local 是一个在 0 到 1 之间的参数，用于加权局部特征的贡献。当 lambda_local 较大时，算法更注重局部信息，这对于处理图像中的细节和纹理等局部特征很有帮助。
% 
%     lambda_global 也是一个在 0 到 1 之间的参数，用于加权非局部特征的贡献。当 lambda_global 较大时，算法更注重全局信息，这对于处理整体结构和背景等非局部特征很有帮助。
% 这两个参数的和应该小于或等于 1，以确保权重的合理性。您可以根据实际问题和图像特性进行调整，以获得最佳的聚类效果。
% 
% 例如，如果您认为图像中的局部细节对于聚类任务很重要，您可能希望增加 lambda_local。如果您认为全局结构更为重要，您可能会调整 lambda_global。通常，这两个参数的调整需要一些试验和经验。    


    % 转换图像矩阵为特征矩阵，例如将每个像素看作具有三个通道的数据点
%     features = reshape(img, [], size(img, 3));

    % 计算初始聚类中心
    initial_centers = kmeans_init(img, K);

    % 运行 K均值 聚类
    idx = weighted_kmeans(img, K, initial_centers, lambda_local, lambda_global);

    % 将聚类结果映射回图像
    idx = reshape(idx, Size(1), Size(2));
end

function initial_centers = kmeans_init(features, K)
    % 随机选择 K 个数据点作为初始聚类中心
    indices = randperm(size(features, 1), K);
    initial_centers = features(indices, :);
end

function idx = weighted_kmeans(features, K, initial_centers, lambda_local, lambda_global)
    % 加权 K均值 聚类

    % 设置迭代次数和容差
    max_iters = 100;
    tol = 1e-6;

    % 初始化变量
    prev_centers = initial_centers;
    centers = initial_centers;
    idx = zeros(size(features, 1), 1);

    for iter = 1:max_iters
        % 分配数据点到最近的聚类中心
        [~, idx] = pdist2(centers, features, 'euclidean', 'Smallest', 1);

        % 更新聚类中心
        for i = 1:K
            cluster_points = features(idx == i, :);
            local_feature = mean(cluster_points);
            global_feature = mean(features);
            
            % 使用加权平均计算新的聚类中心
            centers(i, :) = lambda_local * local_feature + lambda_global * global_feature + (1 - lambda_local - lambda_global) * prev_centers(i, :);
        end

        % 检查收敛条件
        if norm(centers - prev_centers, 'fro') < tol
            break;
        end

        % 更新前一轮的聚类中心
        prev_centers = centers;
    end
end
