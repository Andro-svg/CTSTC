function [B,T,Noise] = LRSD(D, opts)

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');            tol      = opts.tol;                  end
if isfield(opts, 'max_iter');       max_iter = opts.max_iter;             end
if isfield(opts, 'lambda1');        lambda1  = opts.lambda1;              end
if isfield(opts, 'lambda2');        lambda2  = opts.lambda2;              end
if isfield(opts, 'lambda3');        lambda3  = opts.lambda3;              end
if isfield(opts, 'gamma');          gamma    = opts.gamma;                end
if isfield(opts, 'mu');             mu       = opts.mu;                   end
if isfield(opts, 'max_mu');         max_mu   = opts.max_mu;               end
if isfield(opts, 'K');              K   = opts.K;               end

%% Initialization
Nway = size(D);
Noise = ones(Nway);
B = ones(Nway);
T = ones(Nway);
M1 = zeros(Nway);
M2 = zeros(Nway);
M3 = zeros(Nway);


Y1 = 0.5*ones(Nway);  % for background
Y2 = 0.5*ones(Nway);  % for target
W = ones(size(T)); 

for u=1:Nway(3)
    A(:,:,u) = eye(Nway(1:2));
end
Sigma_obs=A;
Sigma_trans=A;
for u=1:Nway(3)
    Sigma_obs_inv(:,:,u) = inv(Sigma_obs(:,:,u));
end
Sigma_trans_inv = Sigma_obs_inv;

preNumT = numel(T);
 
for iter = 1 : max_iter
    %% update Y1
    Y1_flat = reshape(Y1, [], size(Y1, 3));
    lambda_local = 0.4;
    lambda_global = 0.6;
    idx = improved_kmeans(Y1_flat, K, size(Y1), lambda_local, lambda_global);
    Y1_clusters = reshape(idx, size(Y1, 1), size(B, 2));
    all_cluster = {};

    for i = 1:K
        indices = (Y1_clusters == i);
        Y1_i = Y1_flat(indices, :);
        all_cluster{i} = Y1_i;
    end 
    sbd = sum_between_class_distances(all_cluster, K);

    for i = 1:K
        indices = (Y1_clusters == i);
        swd = sum_within_class_distances(all_cluster{i});
        rank_Y1_i = swd / sbd;
        weight = sum(indices(:)) / numel(Y1_clusters);
        Y1 = Y1 + weight * rank_Y1_i;
    end

    Y1 = Y1 + mu * (B + M2);
    
    %% update B
    tmpB = B;
    B = (D-T-Noise-M1+Y1-M2)/2;
    
    %% update T
    tmpT = T;
    tempT = D-B-Noise+Y2-M1-M3;
    thres = W*lambda1/mu;
    T = real( 0.5* prox_l1(tempT, thres));
    
    W = 1 ./ ((abs(T))+ 0.001);

    %% update Y2
    GCC=[];
    for k=1:Nway(3)
        [Gx, Gy, Gxx, Gxy, Gyx, Gyy] = Hessian_matrix(T(:,:,k), 3);
        GCC(:,:,k)= (Gxx .* Gyy - Gxy .* Gyx) ./ (1 + Gx .* Gx + Gy .* Gy);
    end
    Z = double(abs(GCC) > 0.9*max(max(max(GCC))));

    Y2 = update_Y2(Y2, T, M3, Z, Sigma_obs_inv, Sigma_trans_inv, mu, lambda3, A);

    %% update Noise
    Noise = real((mu*(D-B-T-M1))/(2*lambda2+mu));
    
    %% check the convergence
    currNumT = sum(T(:) > 0); 
    chg =norm(D(:)-B(:)-T(:)-Noise(:))/norm(D(:));
    fprintf('iter = %d   res=%.10f \n', iter, chg);    

    if (chg < tol) || (currNumT == preNumT)
        break;
    end
    if norm(tmpT(:) - T(:))< tol
        break;
    end
    if norm(tmpB(:) - B(:))< tol
        break;
    end
    preNumT = currNumT;
    
    %% update Lagrange multipliers M and penalty parameter mu
    M1 = M1 - (D-B-T-Noise);
    M2 = M2 - (Y1-B);
    M3 = M3 - (Y2-T);
    mu = min(gamma*mu,max_mu);  

end
end

function sum_within = sum_within_class_distances(Z1_i)
    sum_within = norm(Z1_i, 'fro');
end

function sum_between = sum_between_class_distances(Z1_i, K)
    F=[];
    for i = 1:K
        F(i) = norm(Z1_i{i}, 'fro');
    end
    sum_between = var(F);
end


function Y2 = update_Y2(Y2, T, M3, Z, Sigma_obs_inv, Sigma_trans_inv, mu, lambda3, A)
    % Update Y2 using closed-form solution
    n3 = size(T, 3);   
    for t = 1:n3
        rhs = lambda3 * (Sigma_obs_inv(:,:,t) * Z(:,:,t) + Sigma_trans_inv(:,:,t) * A(:,:,t) * Z(:,:,max(t-1,1)) * Y2(:,:,max(t-1,1)))...
            + mu * (T(:,:,t) + M3(:,:,t));
        lhs = mu * eye(size(T(:,:,t))) + lambda3 * (Sigma_obs_inv(:,:,t) + Sigma_trans_inv(:,:,t));
        Y2(:,:,t) = rhs / lhs;

    end
end