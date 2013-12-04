% % MUSEdemo.m
% script: (no args, no output)
% 
% % demonstrate theory section of:
% "A robust multi-shot scan strategy for high-resolution diffusion weighted MRI 
% enabled by multiplexed sensitivity-encoding (MUSE)"
% Nan-kuei Chan, et al 2013
%
% Frederick Bryan

% load data; 
load('brain_8ch'); % contains 'im' and 'map'

Ncoils = size(map,3);

% % create synthetic data
% im = padarray(phantom(160)',[0 30],'both');
% % add some phase
% [r,c] = ndgrid(-linspace(pi,pi/4,size(im,1)),linspace(pi,0,size(im,2)) );
% im = im + 1i*(r.*c); 
% % sprinkle in a tiny amount of noise
% im = (real(im)+.01*randn(size(im))) + 1i*(imag(im)+.01*randn(size(im)));
% tim = im; % keep a "true copy"

% add the coil sensitivity encoding
im = repmat(im,[1 1 Ncoils]).*map;
% sprinkle in some coil-dependent noise
im = (real(im)+.4*randn(size(im))) + 1i*(imag(im)+.4*randn(size(im)));

% reference recon, full data
% parallel imaging recon
full_recon = 1./sum(map .* conj(map),3) .* sum(conj(map).*im,3); 
% sum of squares approximation
SOS_recon = sqrt(sum(conj(im).*im,3)); % fully sampled
figure; 
subplot(2,1,1); dispimg(SOS_recon);
    title('Fully sampled, Sum of Squares Recon');
subplot(2,1,2); dispimg(real(full_recon));
    title('Fully sampled, "correct" recon');
    
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

% % add a little phase shift to the 2nd shot
% shift = .25;% .25 pixel shift 
% [phx,phy] = meshgrid((-floor(size(k,2)/2):floor(size(k,2)/2)-1)/size(k,2),...
%     (-floor(size(k,1)/2):1:floor(size(k,1)/2)-1)/size(k,1));
% ph = exp(1i*2*pi*(phx+phy)*shift); 
% k2 = k2 .* repmat(ph,[1 1 size(k2,3)]);
%
% Demo/testing - leave commented
% figure;dispimg(abs(SOS_recon)); title('full recond, unshifted')
% figure;dispimg(abs(ift2(ft2(SOS_recon).*ph))); 
%     title('full recon, .5 vox shift');

% reconstruct each of the shots (produces aliased images)
aim1 = ift2(k1); % aka u
aim2 = ift2(k2); % aka v

figure; dispimg(abs(aim1(:,:,4))); title('Shot 1, coil 4, showing aliasing');

% sense to make an unaliased images
im1 = sense(aim1, map, Rx, Ry);
im2 = sense(aim2, map, Rx, Ry);

figure; 
c = 4; r = 2; 
subplot(c,r,1); imagesc(real(full_recon),[0 1]); 
    title('real(truth)'); colorbar; axis image;
subplot(c,r,2); imagesc(imag(full_recon)); 
    title('imag(truth)'); colorbar; axis image;
subplot(c,r,3); imagesc(real(im1),[0 1]);  
    title('real(im1)'); colorbar; axis image;
subplot(c,r,4); imagesc(imag(im1));  
    title('imag(im1)'); colorbar; axis image;
subplot(c,r,5); imagesc(real(im2),[0 1]); 
    title('real(im2)'); colorbar; axis image;
subplot(c,r,6); imagesc(imag(im2)); 
    title('imag(im2)'); colorbar; axis image;
subplot(c,r,5); imagesc(real(im2),[0 1]);
    title('real(im2)'); colorbar; axis image;
subplot(c,r,6); imagesc(imag(im2)); 
    title('imag(im2)'); colorbar; axis image;
subplot(c,r,5); imagesc(real(im2),[0 1]); 
    title('real(im2)'); colorbar; axis image;
subplot(c,r,6); imagesc(imag(im2)); 
    title('imag(im2)'); colorbar; axis image;
subplot(c,r,7); imagesc(real((im2+im1)/2),[0 1]); 
    title('real(im2+im1)'); colorbar; axis image;
subplot(c,r,8); imagesc(imag((im2+im1)/2)); 
    title('imag(im2+im1)'); colorbar; axis image;

% total variation denoising each image %(cheating b/c I don't want to implement
%   TV filtering
tv1 = medfilt2(real(im1)) + 1i* medfilt2(imag(im1));
tv2 = medfilt2(real(im2)) + 1i* medfilt2(imag(im2));

% get phase estimates
e1 = tv1./abs(tv1); e1(isnan(e1))=0; % this may not be necessary if using regress
e2 = tv2./abs(tv2); e2(isnan(e2))=0;

% preallocate images
B = zeros(size(im1,1), size(im1,2), 2);
% Bint = B;

% create circshifted version of map and phase estimates
circmap = circshift(map,[size(map,1)/2 0 0]);
circe1 =  circshift(e1, [size(map,1)/2 0 0]);
circe2 =  circshift(e2, [size(map,1)/2 0 0]);

Yim = cat(3,aim1,aim2); 
Xim = cat(4, ...
    cat(3,map.*repmat(e1,[1 1 Ncoils]), ...
          map.*repmat(e2,[1 1 Ncoils])),... 
    cat(3,circmap.*repmat(circe1,[1 1 Ncoils]), ...
          circmap.*repmat(circe2,[1 1 Ncoils])) );

% figure; imagesc(real((im2+im1)/2),[0 1]); axis image; axis off; colormap gray; hold on;
% go through each pixel - this is verrry similar to sense
D = zeros(size(map,1),size(map,2));
for row = 1:size(map,1);
    for col = 1:size(map,2);
%         plot(col,row,'yo'); drawnow;
        if sum(map(row,col,:),3) ~= 0 && e1(row,col) ~= 0 && e2(row,col) ~= 0
            % this is mostly redudundant, except for those few dots on the top
            % of the map that go to zero when mean filtering occurs

            if row ==69 && col ==156; keyboard; end
            % set up system of eqns

            % get each aliased pixel
            Y = squeeze(Yim(row,col,:));
            
            X = squeeze(Xim(row, col,:,:));
            keep = sum(X,1)~=0;
            X = X(:,keep);
%             Y= Y(keep);
            
%             temp=((X'*X)^-1)*X'*Y;
            temp = X\Y;
            D(row,col,:) = temp(:);
            % read sense function for description of this math

%             [B(row,col,:),Bint{row,col}] = regress(Y,X);
            
%         end
%         unplot;
    end
end
    
final_im = abs(D(:,:,1));
    
    
    
    
    
    
    