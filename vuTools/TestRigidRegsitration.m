%% 2D Registration test

path(path, './release')
close all
clear all

rand('state',sum(100*clock));

im = imread('cameraman.tif');
im1 = double(im);
im1_scale = rand(1)+1;
im1_resize = imresize(im1,im1_scale,'bilinear');
meta_im1 = vuGenerateMetaImage(im1_resize,[1./im1_scale 1./im1_scale],[0 0]);
dims1 = size(im1_resize);

a = 30*rand(1)-15
x = 20*rand(1)-10
y = 20*rand(1)-10

cosa = cos(pi*a/180);
sina = sin(pi*a/180);

in_tform = maketform('affine', [cosa sina 0;-sina cosa 0;x y 1]);
R = makeresampler('linear','bound');
im2 = tformarray(im1,in_tform,R, [1 2], [1 2], size(im1),[],[]);

im2_scale = rand(1);
im2_resize = imresize(im2,im2_scale,'bilinear');
dims2 = size(im2_resize);
meta_im2 = vuGenerateMetaImage(im2_resize,[1./im2_scale 1./im2_scale],[0 0]);

tic
[reg_im, reg_xform] = RigidRegistration(meta_im1,meta_im2,1);
toc

% Display results
figure;
subplot(2,2,1);imagesc(im1_resize);title('Fixed Image');
subplot(2,2,2);imagesc(im2_resize);title('Moving Image');
subplot(2,2,3);imagesc(reg_im.Data);title('ITK Rigid Registration Output Image');
subplot(2,2,4);imagesc(im1_resize-reg_im.Data);title('Difference Image');colorbar;
%% 3D Registration test
close all
clear all

fid = fopen('RigidRegistrationData.raw');
im1 = reshape(fread(fid,'float'),128,128,128);
meta_im1 = vuGenerateMetaImage(im1,[1 1 1],[0 0 0]);
fclose(fid);
dims = size(im1);

rand('state',sum(100*clock));
a = 30*rand(1)-15
b = 30*rand(1)-15
g = 30*rand(1)-15
x = 20*rand(1)-10
y = 20*rand(1)-10
z = 20*rand(1)-10

cosa = cos(pi*a/180); sina = sin(pi*a/180);
cosb = cos(pi*b/180); sinb = sin(pi*b/180);
cosg = cos(pi*g/180); sing = sin(pi*g/180);

rota = eye(4);
rota(2:3,2:3) = [cosa sina;-sina cosa];

rotb = eye(4);
rotb(1:2:3,1:2:3) = [cosb -sinb;sinb cosb];

rotg = eye(4);
rotg(1:2,1:2) = [cosg sing;-sing cosg];

xlate = eye(4);
xlate(4,1:3) = [x y z];

mcenter = eye(4);
mcenter(4,1:3) = -dims/2.;

center = eye(4);
center(4,1:3) = dims/2.;

in_xform = mcenter*xlate*rotg*rotb*rota*center;

in_tform = maketform('affine',in_xform);
R = makeresampler('linear','bound');
im2 = tformarray(im1,in_tform,R, [1 2 3], [1 2 3], size(im1),[],[]);
meta_im2 = vuGenerateMetaImage(im2,[1 1 1],[0 0 0]);

tic
[reg_im, reg_xform] = RigidRegistration(meta_im1,meta_im2,1);
toc

% Display results
figure;
subplot(2,2,1);imagesc(im1(:,:,64));title('Fixed Image (slice 64)');
subplot(2,2,2);imagesc(im2(:,:,64));title('Moving Image (slice 64)');
subplot(2,2,3);imagesc(reg_im.Data(:,:,64));title('ITK Rigid Registration Output Image (slice 64)');
temp = im1-reg_im.Data;
subplot(2,2,4);imagesc(temp(:,:,64));title('Difference Image (slice 64)');colorbar;