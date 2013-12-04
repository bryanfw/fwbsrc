% Set up im and MEANS and STD outside 
tic
corr_im = BiasFieldCorrection(im,MEANS,STD,20,3);
toc
% Plot
scale = 256;
orgIm = im3(:,:,1);
corrIm = corr_im(:,:,1);
figure;
subplot(2,2,1);
imshow(orgIm,colormap(gray(scale)));
axis equal tight
title('Prefilter Image');
subplot(2,2,2);
imshow(corrIm,colormap(gray(scale)));
axis equal tight
title('Postfilter Image');
subplot(2,2,3);
imshow(corrIm - orgIm,colormap(gray(scale)));
axis equal tight
title('Difference Image');

colorbar