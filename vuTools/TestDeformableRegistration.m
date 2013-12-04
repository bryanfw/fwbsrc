%% 2D Registration test

close all
clear all

im = imread('cameraman.tif');
dims = size(im);
fixed_im = vuGenerateMetaImage(single(im),[1 1],[0 0]);

im = imread('cameraman_deformed.tif');
moving_im = vuGenerateMetaImage(single(im),[1 1], [0 0]);

[out_im] = MEX2DDeformableRegistration(fixed_im,moving_im,1);

%% 3D Registration test

close all
clear all

load MRData.mat

fixed_im = vuGenerateMetaImage(ana_data, ana_parms.tags(1,[29 30 23]), [0 0 0]);
moving_im = vuGenerateMetaImage(dti_data(:,:,:,1,1,1,1), dti_parms.tags(1,[29 30 23]), [0 0 0]);

[out_im, rigid_xform] = RigidRegistration(fixed_im,moving_im,1);
[out_im2] = MEX3DResampleImage(moving_im,rigid_xform);
[out_im3] = MEX3DDeformableRegistration(fixed_im,out_im2,1);

%%
resample_xform = struct('Matrix',eye(3),'Offset',zeros(1,3));
fixed_resampled_im = MEX3DResampleImage(fixed_im,resample_xform,moving_im);
ss_fixed_resampled_im = vuGenerateMetaImage(fixed_resampled_im.Data(:,:,18), fixed_resampled_im.Spc(1:2), fixed_resampled_im.Origin(1:2));
ss_fixed_im = vuGenerateMetaImage(fixed_im.Data(:,:,30), fixed_im.Spc(1:2), fixed_im.Origin(1:2));
ss_moving_im = vuGenerateMetaImage(out_im.Data(:,:,30), out_im.Spc(1:2), out_im.Origin(1:2));

%%
tic;
out_im2 = MEX2DDeformableRegistration(ss_fixed_im,ss_moving_im,1);
toc
%%
subplot(1,4,1);imagesc(fixed_im.Data(:,:,30))
subplot(1,4,2);imagesc(out_im2.Data);
subplot(1,4,3);imagesc(ss_moving_im.Data);
subplot(1,4,4);quiver(deform_xform.Data(:,:,1),deform_xform.Data(:,:,2));

%%

fid = fopen('Fixed.raw','w')
fwrite(fid,fixed_im.Data(:,:,31),'float32');
fclose(fid);

fid = fopen('Moving.raw','w')
fwrite(fid,out_im.Data(:,:,31),'float32');
fclose(fid);