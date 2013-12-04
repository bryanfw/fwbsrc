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
inImg = vuGenerateMetaImage(im1,[1./256 1./256 ],[0 0]);
outImg = GradientAnisotropicDiffusion(inImg,5,0.125,1.0,0);

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

im1= imread('cameraman.tif'); % scale = 256

% 3D
[a,b] = size(im1);
in3D = zeros(3,a,b);
in3D(1,:,:) = im1;
in3D(2,:,:) = im1;
in3D(3,:,:) = im1;

% Generate Meta Image and Filter
inImg3D = vuGenerateMetaImage(in3D,[1./256 1./256 1./256],[0 0 0]);
outImg3D = GradientAnisotropicDiffusion(inImg3D,5,0.0625,1.0,0);

in2D = zeros(a,b);
out2D = zeros(a,b);
in2D(:,:) = inImg3D.Data(1,:,:);
out2D(:,:) = outImg3D.Data(1,:,:);

% Plot
figure;
subplot(2,2,1);
imshow(in2D,colormap(gray(scale)));
axis equal tight
title('Prefilter Image');
subplot(2,2,2);
imshow(out2D,colormap(gray(scale)));
axis equal tight
title('Postfilter Image');
subplot(2,2,3);
imshow(in2D-out2D,colormap(gray(scale)));
axis equal tight
title('Difference Image');

colorbar

%fid = fopen('C:\Documents and Settings\segalj\Desktop\ImageProcessing\Registration\Rigid\MEXRigidRegistration\RigidRegistrationData.raw');
%im1 = reshape(fread(fid,'float'),128,128,128);
%meta_im1 = GenerateMetaImage(im1,[1 1 1],[0 0 0]);
%fclose(fid);
%dims = size(im1);

%outImg = MEX3DGradientAnisotropicDiffusion(meta_im1, 10, 0.125, 2.0);

%figure;
%subplot(1,2,1);imagesc(im1(:,:,64));title('Prefilter Image (slice 64)');
%subplot(1,2,2);imagesc(outImg.Data(:,:,64));title('Postfilter Image (slice 64)');