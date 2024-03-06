img = imread('E:\2024work\87.bmp');
if size(img, 3) > 1
    img = rgb2gray( img );
end
img = double(img);

% G = fspecial('gaussian', [2 2], 2); % Gaussian kernel
% u = imfilter(img, G, 'symmetric');
[Gx, Gy] = gradient(u);

[Gxx, Gxy] = gradient(Gx);
[Gyx, Gyy] = gradient(Gy);

GCCmap= (Gxx .* Gyy - Gxy .* Gyx) ./ (1 + Gx .* Gx + Gy .* Gy);
imshow(mat2gray(GCCmap))
