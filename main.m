

clear
clc

image = 'image.png';
sigma_s = 4;
sigma_r = 0.1;
iteration = 4;
GaussianPrecision = 0.05; %sigma=4 0.1-17*17   0.05-19*19
                          %sigma=3 0.1-13*13

I = im2double(imread(image));

tic;
result = RollingGuidanceFilter(I, sigma_s, sigma_r, iteration, GaussianPrecision);
toc;

figure,imshow(I);
figure,imshow(result);