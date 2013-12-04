% % MUSEdemo.m
% script: (no args, no output)
% 
% % demonstrate theory section of:
% "A robust multi-shot scan strategy for high-resolution diffusion weighted MRI 
% enabled by multiplexed sensitivity-encoding (MUSE)"
% Nan-kuei Chan, et al 2013
%
% Frederick Bryan

clear; 

% seed rand
rng('shuffle');

% load data; 
load('brain_8ch'); % contains 'im' and 'map'

Ncoils = size(map,3);

im = padarray(phantom(160)',[0 30],'both');

% add some phase
[r,c] = ndgrid(-linspace(0.1*pi,0.2*pi,size(im,1)),linspace(.3*pi,0,size(im,2)) );
im = im .*exp(1i*(r.*c) );
% sprinkle in a tiny amount of noise
im = (real(im)+.01*randn(size(im))) + 1i*(imag(im)+.01*randn(size(im)));
tim = im; % keep a "true copy"

% add the coil sensitivity encoding
im = repmat(im,[1 1 Ncoils]).*map;

% sprinkle in some coil-dependent noise
im = (real(im)+.2*randn(size(im))) + 1i*(imag(im)+.2*randn(size(im)));

% reference recon, full data
% parallel imaging recon
full_recon = 1./sum(map .* conj(map),3) .* sum(conj(map).*im,3); 
% sum of squares approximation
SOS_recon = sqrt(sum(conj(im).*im,3)); % fully sampled

% Pulse sequence stuff
% lets do 2 shot-interleave (2 shots, sense factor = 2)
Sfactor = 2;
Nshots = Sfactor;
PFperc = 9/16; % (56%)

% lets go k-space
k = ft2(im); 
% apply partial fourier
% k((floor(PFperc*size(k,1))+1):end, :, :) = 0;

% create the two shots
Rx = 1;
Ry = Sfactor; 
% generate mask w/ correct Rx,Ry undersampling
mask = zeros(size(k)); 
mask(1:Ry:end,:,:) = 1;
mask(:,1:Rx:end,:) = mask(:,1:Rx:end,:)+1;
mask(mask<max(mask(:)))=0; mask(mask>0)=1;
% apply mask and circshifted-by-1 mask to get the two interleaves
k1 = k .* mask * numel(mask)/sum(mask(:)); % 1st interleave (w/ fixed intensity)
k2 = k .* circshift(mask,[1 0 0]) * numel(mask)/sum(mask(:)); % 2nd interleave

% add random (smooth) phase to both images
% generate random phase
phstd = 4*pi; % noise standard deviation
ph1 = phstd*randn(size(tim)); 
ph2 = phstd*randn(size(tim)); 
% smooth phase
gstd = 10; gsz = 50; % std dev and size of smoothing kernel
g = fspecial('gaussian',gsz,gstd); 
ph1 = imfilter(ph1,g); ph1 = zeros(size(ph1)); 
ph2 = imfilter(ph2,g); ph2 = ph2+ones(size(ph2))*deg2rad(45);
% add phase to existing base phase
ph1 = repmat(ph1,[1 1 size(im,3)]);
ph2 = repmat(ph2,[1 1 size(im,3)]);
k1 = abs(k1).*exp(1i*(angle(k1)+ph1));
k2 = abs(k2).*exp(1i*(angle(k2)+ph2)); 
% figure;subplot(2,1,1); imagesc(ph(:,:,1)); colorbar; axis image; axis off;
% subplot(2,1,2); imagesc(angle(k2(:,:,1))); colorbar; axis image; axis off; 

% Demo/testing - leave commented
% figure;dispimg(abs(SOS_recon)); title('full recond, unshifted')
% figure;dispimg(abs(ift2(ft2(SOS_recon).*ph))); 
%     title('full recon, .5 vox shift');

% reconstruct each of the shots (produces aliased images)
kreg = ift2(k1+k2); % regular case, what the scanner might actually output and 
                    % you try to reconstruct if you were ig'nant
aim1 = ift2(k1); % aka u
aim2 = ift2(k2); % aka v

% sense to make a bad image
imreg = sense(kreg,map,1,1);

% sense to make an unaliased images
im1 = sense(aim1, map, Rx, Ry);
im2 = sense(aim2, map, Rx, Ry);


% total variation denoising each image %(cheating b/c I don't want to implement
%   TV filtering
% tv1 = medfilt2(real(im1)) + 1i* medfilt2(imag(im1));
% tv2 = medfilt2(real(im2)) + 1i* medfilt2(imag(im2));
tvstd = 1/gstd;
gg = fspecial('gaussian',gsz,tvstd); 
tv1 = imfilter(real(im1),gg) + 1i* imfilter(imag(im1),gg);
tv2 = imfilter(real(im2),gg) + 1i* imfilter(imag(im2),gg);

% get phase estimates
e1 = tv1./abs(tv1); e1(isnan(e1))=0; 
e2 = tv2./abs(tv2); e2(isnan(e2))=0;


% actual MUSE work
% create circshifted version of map and phase estimates
circmap = circshift(map,[size(map,1)/2 0 0]);
circe1 =  circshift(e1, [size(map,1)/2 0 0]);
circe2 =  circshift(e2, [size(map,1)/2 0 0]);

% this sets up the system of eqns.  
% THIS NEEDS TO BE GENERALIZED FROM THE 2-SHOT CASE TO THE N-SHOT CASE.
% THAT EFFECTS THE NEGATIVE SIGN ABOUT 4 LINES DOWN FROM HERE. 
% IT CAN BE CALCULATED BY CONSIDERING THE EFFECT OF A 1/2/3/WHATEVER 
%    VOXEL SHIFT FOR EACH CIRCSHIFTED VERSION.
Yim = cat(3,aim1,aim2); 
Xim = cat(4, ...
    cat(3,map.*repmat(e1,[1 1 Ncoils]), ...
          map.*repmat(e2,[1 1 Ncoils])),... 
    cat(3,circmap.*repmat(circe1,[1 1 Ncoils]), ...
          -circmap.*repmat(circe2,[1 1 Ncoils])) );
      

% go through each pixel - this is verrry similar to sense
D = zeros(size(map,1),size(map,2),2);
for row = 1:size(map,1);
    for col = 1:size(map,2);
        if sum(map(row,col,:),3) ~= 0 && e1(row,col) ~= 0 && e2(row,col) ~= 0
            % this is mostly redudundant, except for those few dots on the top
            % of the map that go to zero when mean filtering occurs

            % set up system of eqns

            % get each aliased pixel
            Y = squeeze(Yim(row,col,:));
            
            X = squeeze(Xim(row, col,:,:));
            keep = any(X);
            X = X(:,keep);
            
            temp=((X'*X)^-1)*X'*Y;
%             temp = X\Y; 
            D(row,col,:) = temp(:);
            % read sense function for description of this math

        end
    end
end
   
figure; 
foo = cat(4,abs(full_recon),abs(im1), abs(im2), abs((im1+im2)/2),abs(imreg/2), abs(D(:,:,1)));
montage(foo); 
title(sprintf(['L->R:   truth               senseshot1         senseshot2 \n '...
               'senseshot1+senseshot2   sense(regular)      muse']),...
               'FontName','FixedWidth', 'HorizontalAlignment','Center');
figure; 
c = 3; r = 2; 
subplot(c,r,1); imagesc(abs(full_recon),[0 1]); 
    title('abs(truth)'); colorbar; axis image; axis off;
subplot(c,r,2); imagesc(abs(im1)); 
    title('abs(im1)'); colorbar; axis image; axis off;
subplot(c,r,3); imagesc(abs(im2),[0 1]);  
    title('abs(im2)'); colorbar; axis image; axis off;
subplot(c,r,4); imagesc(abs(im1+im2)/2);  
    title('abs(im1+im2)'); colorbar; axis image; axis off;
subplot(c,r,5); imagesc(abs(imreg/2));  
    title('abs(sense)'); colorbar; axis image; axis off;
subplot(c,r,6); imagesc(abs(D(:,:,1)));  
    title('abs(MUSE)'); colorbar; axis image; axis off;

    
    
    % % add a little phase shift to the 2nd shot
% shift = .25;% .25 pixel shift 
% [phx,phy] = meshgrid((-floor(size(k,2)/2):floor(size(k,2)/2)-1)/size(k,2),...
%     (-floor(size(k,1)/2):1:floor(size(k,1)/2)-1)/size(k,1));
% ph = exp(1i*2*pi*(phx+phy)*shift); 
% k2 = k2 .* repmat(ph,[1 1 size(k2,3)]);
    


% figure; 
% c = 4; r = 2; 
% subplot(c,r,1); imagesc(real(full_recon),[0 1]); 
%     title('real(truth)'); colorbar; axis image;
% subplot(c,r,2); imagesc(imag(full_recon)); 
%     title('imag(truth)'); colorbar; axis image;
% subplot(c,r,3); imagesc(real(im1),[0 1]);  
%     title('real(im1)'); colorbar; axis image;
% subplot(c,r,4); imagesc(imag(im1));  
%     title('imag(im1)'); colorbar; axis image;
% subplot(c,r,5); imagesc(real(im2),[0 1]); 
%     title('real(im2)'); colorbar; axis image;
% subplot(c,r,6); imagesc(imag(im2)); 
%     title('imag(im2)'); colorbar; axis image;
% subplot(c,r,5); imagesc(real(im2),[0 1]);
%     title('real(im2)'); colorbar; axis image;
% subplot(c,r,6); imagesc(imag(im2)); 
%     title('imag(im2)'); colorbar; axis image;
% subplot(c,r,5); imagesc(real(im2),[0 1]); 
%     title('real(im2)'); colorbar; axis image;
% subplot(c,r,6); imagesc(imag(im2)); 
%     title('imag(im2)'); colorbar; axis image;
% subplot(c,r,7); imagesc(real((im2+im1)/2),[0 1]); 
%     title('real(im2+im1)'); colorbar; axis image;
% subplot(c,r,8); imagesc(imag((im2+im1)/2)); 
%     title('imag(im2+im1)'); colorbar; axis image;
    
    
    