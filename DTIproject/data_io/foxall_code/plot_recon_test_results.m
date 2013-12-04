%-------------------------------------------------------------------------------------------------------
% PLOT_RECON_TEST_RESULTS - Make Plots To show Intermediate Reconstruction Stages
%-------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------
% CONSTANTS
%-------------------------------------------------------------------------------------
GRAY_VALUES = 256;
GRAY_RANGE  = 96;
GRAY_CENTER = 128;

%-------------------------------------------------------------------------------------
% PLOT 1 - Kspace Images
%-------------------------------------------------------------------------------------

figure

abs_scale = GRAY_RANGE/max(max(max(abs(kspace))));


for S = 1:Ns

    plot_index = (S-1)*2 + 1;

    subplot(Ns,2,plot_index);
    pic = GRAY_CENTER + abs_scale*squeeze(abs(pspace(S,:,:)));
    image(pic);
    axis('image');
    title(['P Echo K Space - Shot ',int2str(S)]);
    colormap(gray(GRAY_VALUES));  

    plot_index = (S-1)*2 + 2;

    subplot(Ns,2,plot_index);
    pic = GRAY_CENTER + abs_scale*squeeze(abs(nspace(S,:,:)));
    image(pic);
    axis('image');
    title(['N Echo K Space - Shot ',int2str(S)]);
    colormap(gray(GRAY_VALUES));

end

%----------------------------------------------------------------------------------------
% PLOT 2 - Undersampled Images
%----------------------------------------------------------------------------------------

figure

abs_scale = GRAY_RANGE/max(max(max(abs([pimage,nimage]))));
phz_scale = GRAY_RANGE/(2*pi);


for S = 1:Ns

plot_index = (S-1)*4 + 1;

    subplot(Ns,4,plot_index);
    pic = GRAY_CENTER + abs_scale*squeeze(abs(pimage(S,:,:)));
    image(pic);
    axis('image');
    title(['P Echo Image - Shot ',int2str(S)]);
    colormap(gray(GRAY_VALUES));  

plot_index = (S-1)*4 + 2;

    subplot(Ns,4,plot_index);
    pic = GRAY_CENTER + phz_scale*squeeze(angle(pimage(S,:,:)));
    image(pic);
    axis('image');
    title(['P Echo Phase - Shot ',int2str(S)]);
    colormap(gray(GRAY_VALUES));


plot_index = (S-1)*4 + 3;

    subplot(Ns,4,plot_index);
    pic = GRAY_CENTER + abs_scale*squeeze(abs(nimage(S,:,:)));
    image(pic);
    axis('image');
    title(['N Echo Image - Shot ',int2str(S)]);
    colormap(gray(GRAY_VALUES));

plot_index = (S-1)*4 + 4;

    subplot(Ns,4,plot_index);
    pic = GRAY_CENTER + phz_scale*squeeze(angle(nimage(S,:,:)));
    image(pic);
    axis('image');
    title(['N Echo Phase - Shot ',int2str(S)]);
    colormap(gray(GRAY_VALUES));

end


%----------------------------------------------------------------------------------------
% PLOT 3 - First Order Corrected Images
%----------------------------------------------------------------------------------------
figure

abs_scale = GRAY_RANGE/max(max(max(abs([fpc_pimage,fpc_nimage]))));
phz_scale = GRAY_RANGE/(2*pi);


for S = 1:Ns

plot_index = (S-1)*4 + 1;

    subplot(Ns,4,plot_index);
    pic = GRAY_CENTER + abs_scale*squeeze(abs(fpc_pimage(S,:,:)));
    image(pic);
    axis('image');
    title(['1st Order Corrected P Echo Image - Shot ',int2str(S)]);
    colormap(gray(GRAY_VALUES));  

plot_index = (S-1)*4 + 2;

    subplot(Ns,4,plot_index);
    pic = GRAY_CENTER + phz_scale*squeeze(angle(fpc_pimage(S,:,:)));
    image(pic);
    axis('image');
    title(['1st Order Corrected P Echo Phase - Shot ',int2str(S)]);
    colormap(gray(GRAY_VALUES));


plot_index = (S-1)*4 + 3;

    subplot(Ns,4,plot_index);
    pic = GRAY_CENTER + abs_scale*squeeze(abs(fpc_nimage(S,:,:)));
    image(pic);
    axis('image');
    title(['1st Order Corrected N Echo Image - Shot ',int2str(S)]);
    colormap(gray(GRAY_VALUES));

plot_index = (S-1)*4 + 4;

    subplot(Ns,4,plot_index);
    pic = GRAY_CENTER + phz_scale*squeeze(angle(fpc_nimage(S,:,:)));
    image(pic);
    axis('image');
    title(['1st Order Corrected N Echo Phase - Shot ',int2str(S)]);
    colormap(gray(GRAY_VALUES));

end

%----------------------------------------------------------------------------------------
% PLOT 4 - Zero Order Corrected Images
%----------------------------------------------------------------------------------------
figure

abs_scale = GRAY_RANGE/max(max(max(abs(zpc_image))));
phz_scale = GRAY_RANGE/(2*pi);


for S = 1:Ns

plot_index = (S-1)*2 + 1;

    subplot(Ns,2,plot_index);
    pic = GRAY_CENTER + abs_scale*squeeze(abs(zpc_image(S,:,:)));
    image(pic);
    axis('image');
    title(['Zero Order Corrected P+N Image - Shot ',int2str(S)]);
    colormap(gray(GRAY_VALUES));  

plot_index = (S-1)*2 + 2;

    subplot(Ns,2,plot_index);
    pic = GRAY_CENTER + phz_scale*squeeze(angle(zpc_image(S,:,:)));
    image(pic);
    axis('image');
    title(['Zero Order Corrected P+N Phase - Shot ',int2str(S)]);
    colormap(gray(GRAY_VALUES));


end


%----------------------------------------------------------------------------------------
% PLOT 5 - Combined Image
%----------------------------------------------------------------------------------------
figure 

abs_scale = GRAY_RANGE/max(max(max(abs(cpc_image))));
phz_scale = GRAY_RANGE/(2*pi);

subplot(2,1,1)
pic = GRAY_CENTER + abs_scale*abs(cpc_image);
image(pic);
axis('image');
title('Zero Order Corrected Combined Image');
colormap(gray(GRAY_VALUES));

subplot(2,1,2)
pic = GRAY_CENTER + phz_scale*angle(cpc_image);
image(pic);
axis('image');
title('Zero Order Corrected Combined Phase');
colormap(gray(GRAY_VALUES));
