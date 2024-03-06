function filtered_img = DesignedFilter(gray_img)
        % 对灰度图像进行傅里叶变换
        fft_img = fft2(gray_img);
        
        % 将频域表示进行中心化
        fft_img = fftshift(fft_img);

        [x, y] = meshgrid(-128:127, -128:127);        
        sigma = 40;%45; % 方差越大，高斯分布越平缓，滤波器的频率响应越低，可以抑制低频信号
        filter = exp(-(x.^2 + y.^2)/(2*sigma^2));
        filter = 1 - filter/max(filter(:));
        
        % 将频域滤波器应用于频域表示
        filtered_fft_img = fft_img .* filter;

%         % 设计Laplacian增强滤波器
%         H = fspecial('laplacian', 0.5);
%         [M,N]=size(gray_img);
%         % 应用滤波器到傅里叶变换图像上
%         filtered_fft_img = 1 - fft_img .* fft2(H, M, N);

%         [M,N]=size(gray_img);
%         % 设计理想高通滤波器
%         D0 = 8;  % 截止半径
%         H = ones(M, N);
%         [X, Y] = meshgrid(1:N, 1:M);
%         center_x = ceil(N/2);
%         center_y = ceil(M/2);
%         H(sqrt((X-center_x).^2 + (Y-center_y).^2) <= D0) = 0;
%         
%         % 应用滤波器到傅里叶变换图像上
%         filtered_fft_img = fft_img .* H;

        % 将滤波后的频域表示进行反中心化
        filtered_fft_img = ifftshift(filtered_fft_img);
        
        % 对滤波后的频域表示进行傅里叶反变换
        filtered_img = real(ifft2(filtered_fft_img));




%     fft_image = fft2(double(gray_img));
%     获取图像大小
%     [M, N] = size(fft_image);
%     
%     设计高斯滤波器，决定哪些频率成分应该被抑制
%     sigma = 10;  % 标准差
%     [X, Y] = meshgrid(1:N, 1:M);
%     center_x = ceil(N/2);
%     center_y = ceil(M/2);
%     gaussian_filter = exp(-((X-center_x).^2 + (Y-center_y).^2) / (2*sigma^2));
%     
%     将高斯滤波器应用到傅里叶变换图像上
%     fft_image_filtered = fft_image .* gaussian_filter;
%     
%     进行逆傅里叶变换，将图像变换回空域
%     filtered_img = real(ifft2(fft_image_filtered));
end
