close all
clear all

% Some images
load clown % scale = 81
clown = X;
load spine % scale = 64
spine = X;
load flujet % scale = 64
flujet = X;
camera = imread('cameraman.tif'); % scale = 256

% Image to use to test 2D
im = camera;
scale = 256;
im1 = double(im);

% Generate Meta Image and Filter
inImg = GenerateMetaImage(im1);
outImg = RecursiveGaussianBlur(inImg,5,0);

% Plot
figure;
subplot(2,2,1);
imshow(inImg.Data,colormap(gray(scale)));
axis equal tight
title('Prefilter Image');
subplot(2,2,2);
imshow(outImg.Data,colormap(gray(scale)));
axis equal tight
title('Postfilter Image');
subplot(2,2,3);
imshow(inImg.Data-outImg.Data,colormap(gray(scale)));
axis equal tight
title('Difference Image');

colorbar