function [All_Num,time_per_image] = target_detection(readPath, savePath, tuneopts)

if isfield(tuneopts, 'temporal_step');       temporal_step  = tuneopts.temporal_step;   end
if isfield(tuneopts, 'lambdaL');             lambdaL        = tuneopts.lambdaL;         end
if isfield(tuneopts, 'mu');                  mu              = tuneopts.mu;               end
if isfield(tuneopts, 'per');                  per              = tuneopts.per;               end


filesdir = dir([char(readPath) '/*.jpg']);
if isempty( filesdir )
    filesdir = dir( [char(readPath) '/*.bmp'] );
end
if isempty( filesdir )
    filesdir = dir([char(readPath) '/*.png']);
end

files = { filesdir.name };
files = sort_nat(files);
All_Num = length(files);


t_list = [temporal_step+1 : temporal_step : length(files)-temporal_step, length(files)-temporal_step,length(files)];
iteration = 0;
tempT_before=[];
t1 = clock;

for t = t_list
    iteration = iteration + 1;
    fprintf('%s %d%s : \n','Starting', iteration, '-th iteration');
    spat_temp_ten = [];

    %% Construct tensor
    if t == (temporal_step+1)
        for tt=1:temporal_step  
            img = imread([char(readPath) '/GT/' strtok([files{tt}],'.') '_tar.jpg']);
            if size(img, 3) > 1
                img = rgb2gray( img );
            end
            img = double(img);
            spat_temp_ten(:,:,tt) = img;
        end
    else
        for tt=1:temporal_step  
            spat_temp_ten(:,:,tt) = tempT_before(:,:,1+temporal_step+tt);
        end
    end
        if t == t_list(length(t_list))
            for tt=1:(temporal_step+1)
                img = imread([char(readPath) '/' files{t}]);
                if size(img, 3) > 1
                    img = rgb2gray( img );
                end
                img = DesignedFilter(img);
                spat_temp_ten(:,:,temporal_step+tt) = img;
            end
        else
            for tt=1:(temporal_step+1)
                img = imread([char(readPath) '/' files{tt+t-1}]);
                if size(img, 3) > 1
                    img = rgb2gray( img );
                end                
                img = DesignedFilter(img);
                spat_temp_ten(:,:,temporal_step+tt) = img;
            end
        end

    %% Determine k
    [n1,n2,n3]=size(spat_temp_ten);
    X = reshape(spat_temp_ten, n1*n2,n3);
    [~,Sigma,~]=svd(X,'econ');
    Sigmalst=diag(Sigma);
    k=[];
    for i = 2:length(Sigmalst)
        if Sigmalst(i) > per * Sigmalst(i-1)
            k=i-1;
            break
        end
    end
    if isempty(k)
        k=length(Sigmalst);
    end

    %% The LRSD model
    Nway = size(spat_temp_ten);
    opts=[];
    opts.max_iter = 100;
    opts.tol =1e-4;
    opts.gamma = 0.5;
    opts.frame = temporal_step;
    opts.lambda1 = lambdaL/(sqrt(max(Nway(1),Nway(2))) * Nway(3));
    opts.lambda2 = 80 * opts.lambda1;
    opts.lambda3 = 10 * opts.lambda1;
    opts.mu  = mu;
    opts.max_mu = 1e5;
    opts.K = k;
    tenT=[];
    tenB=[];
    tenN=[];
    [~, tenT, ~] = LRSD(spat_temp_ten, opts);
    
    %% reconstrcut imgNum images
    if t == (temporal_step+1)
        for k=1:temporal_step
            tar = spat_temp_ten(:,:,k);
            E = tar;
            imwrite(mat2gray(E), [savePath '/' strtok([files{k}],'.') '_tar.jpg']);  
        end
        for k=(1+temporal_step):(2*temporal_step+1)
            tar = tenT(:,:,k);
            E = tar;
            imwrite(mat2gray(E), [savePath  '/' strtok([files{k}],'.') '_tar.jpg']);
        end
    elseif t == t_list(length(t_list))
        for k=1:(temporal_step)
            tar = tempT_before(:,:,temporal_step+k) + tenT(:,:,k);
            E = tar;
            imwrite(mat2gray(E), [savePath  '/' strtok([files{t-temporal_step+k-1}],'.') '_tar.jpg']);
        end
    else
        for k=1:(temporal_step)
            tar = tempT_before(:,:,temporal_step+k) + tenT(:,:,k);
            E = tar;
            imwrite(mat2gray(E), [savePath   '/' strtok([files{t-temporal_step+k-1}],'.') '_tar.jpg']);
        end 
    end
    tempT_before = tenT;
end
t2=clock;
disp(['Programming has running:',num2str(etime(t2,t1))]);
time_all = etime(t2,t1);
time_per_image = time_all/(iteration*(2*temporal_step+1));
disp(['Each image consumes time: ', num2str(time_per_image)]);
end