%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% NEED TO RUN cuso_recon before you do this next bit


% that other idea


% bad slices:(note intereave-ed-ness)
% 494
% 496
% 498
% 554
% 556
% 558
[]
f  = cell(length(ksz)-2,1);
[f{:}] = ind2sub(ksz(3:end), [494,496,498,554,556,558,674,676,678,734,736,738,794,796,798,854,856,858,914,916,918]);


% of different sort
[f{:}] = ind2sub(ksz(3:end), [2444,2504,2564,2566,2574,2624,2626,2634,268]);

coil_dim_no = find(strcmp(STDstr.dim_labels,'chan'));


[f{:}] = ind2sub(ksz(3:end), 674)
% 674 is the 14th slice in the 4th direction, 4st coil, 2nd fdirection
% :: 
% ec
% 1st volume 

rawk = k_data_rs_uc(:,:,674);
fixk = k_data_rs(:,:,674);

X_range = list_file_info.X_range-list_file_info.kx_offset__calc+1;
% corrs are just applied to the 1-D FT data, lets figure out what they were
corrs = ft(fixk)./builtin('_paren',ft(rawk),X_range(1):X_range(2),:);


% get parrex
[d i]= loadPARREC('Landman_20131030_04_01_10.00.56_(WIP_DTI_low_basic_SENSE).REC');
parim = d(:,:,14,2,1);

% size(d)
% d=squeeze(d);
% size(d)
% 
% d2 = reshape(d,[64 64 60 2 7]);
% d3 = permute(d2,[4 5 1 2 3]);
% d4 = reshape(d3, [2 7 64*64*60]);
% figure; imagesc(sum(d4,3)); % shows which of those last dims is populated
% 
% figure; slice_slider(abs(d2(:,:,:,2,1))); 
% 
% 
% d5 = squeeze(d4(2,:,:)); % get 2nd things (doesnt include b0)
% d6 = sum(d5(1:6,:));
% mean(d6(:))

figure; imagesc(img_normalize(abs(rawk))); colormap gray; axis image; axis off; title('Raw k-space');
figure; imagesc(img_normalize(abs(ft2(rawk)),[0 100]),[.1 1]); colormap gray; axis image; axis off;  title('reconstruction from raw k-space');
figure; imagesc(img_normalize(abs(fixk))); colormap gray; axis image; axis off; title('Raw k-space after Phillips corrections');
figure; imagesc(img_normalize(abs(ft2(fixk)),[0 100]),[.1 1]); colormap gray; axis image; axis off; title('Recon from data after after Phillips corrections');
figure; imagesc(real(corrs)); axis image; axis off; title('Phillips frequency (FRX) and EPI (PHX) corrections (real part)');
figure; imagesc(parim); colormap gray; axis image; axis off; title('Corresponding slice from corresponding PAR/REC');
figure; imagesc(angle(im_data_rs(:,:,674))); axis image; axis off; title('Phase image');
 

im2fix =rawk;

r1 = 64-4;
r2 = 64+3;
c1 = 32-4;
c2 = 32+3;
szm = size(middle);


% make mask
% subtract im fron median filterd version of itself, revealing lumps
tmpim = abs(im2fix);
middle = tmpim(r1:r2,c1:c2);

tmpim(r1:r2,c1:c2) = repmat(median(middle(:)),szm);
g = medfilt2(tmpim,[4 4]);
f = tmpim-g;
f = f-min(f(:));

% get rid of top 1%
tmp = sort(f(:));
thresh = tmp(round(numel(tmp)*.99));

mask = (f>thresh);
mask(r1:r2,c1:c2)= 0;
imagesc(mask)


% apply mask
im2fix(mask) = 0;

% recon 
im2fix(isnan(im2fix)) = 0; %just incase
fixedim = ft(im2fix); % readout dir
fixedim = fixedim(X_range(1):X_range(2),:);
fixedim = fixedim.*corrs;
fixedim(isnan(fixedim)) = 0;
fixedim = ft(fixedim,2);

figure; imagesc(img_normalize( abs(fixedim))); colormap gray; axis image; axis off; title('Corrected');
figure; imagesc(angle(fixedim));  axis image; axis off; title('Corrected - Phase image'); colorbar;






