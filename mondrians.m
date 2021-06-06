function [ img_matrix ] = mondrians( img_size_pix, n, size_range, plot_flag)
% function [ img_matrix ] = mondrians( img_size_pix, n, size_range)
% [ img_matrix ] = mondrians( [600 400], 30*8, [0.05 0.2], 1)
%
% function computes rgb image matrix filled with mondrians
%
% INPUTS:
%
% IMG_SIZE_PIX - matrix dimentions of the output image
% N            - number of mondrians; for balanced colors should be a multiple of 8
% SIZE_RANGE   - range of sizes on fraction of the screen size
% PLOT_FLAG    - if 1 plots the figure
%
% NOTES:
%
% - each shape is a square, i.e. width = height
% - each shape has a unique size
% - sizes, posx and posiy and cols vectors are randomized with randperm
% - if N is a multiple of 8(=number of unique colors), each color will be
% assigned to the equal number of shapes
%
% Natalia Zaretskaya 16.08.2012


color_range = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1; 0 0 0];

% define vectors
img_matrix = zeros(img_size_pix(1), img_size_pix(2), 3);
sizes = round(linspace(img_size_pix(1)*size_range(1), img_size_pix(1)*size_range(2), n));
posx = round(linspace(1, img_size_pix(1), n));
posy = round(linspace(1, img_size_pix(2), n));
cols = repmat(color_range, n/size(color_range,1), 1);

% randomize vectors
sizes = sizes(randperm(length(sizes)));
posx = posx(randperm(length(posx)));
posy = posy(randperm(length(posy)));
cols = cols(randperm(size(cols,1)), :);

% assign values to the matrix
for i = 1:n
    img_matrix(posx(i):posx(i)+sizes(i)-1, posy(i):posy(i)+sizes(i)-1, :) = repmat(reshape(cols(i,:),1,1,3), sizes(i), sizes(i));
end

% center the image, cut the edges:
new_img_size_pix = size(img_matrix);
x_start = round((new_img_size_pix(1)-img_size_pix(1))./2);
x_end = x_start+img_size_pix(1)-1;
y_start = round((new_img_size_pix(2)-img_size_pix(2))./2);
y_end = y_start+img_size_pix(2)-1;

img_matrix = img_matrix(x_start:x_end, y_start:y_end, :);

if plot_flag
    figure;
    imagesc(img_matrix)
end

end

