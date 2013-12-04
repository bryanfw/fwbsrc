close all
clear all

fixed_im = single(imread('cameraman.tif'));
moving_im = single(imread('cameraman.tif'));
%moving_im = imrotate(moving_im,10,'bilinear','crop');

rbf_center = 255.*rand(2,1)';
rbf_radius = 127.5.*rand(1);
rbf_coeff = [0 rbf_radius./2.*rand(1,1)'];
rbf([175 84;0 0;180 90;135 50],rbf_center,rbf_radius)

[XX,YY] = meshgrid(0:size(fixed_im,1)-1,0:size(fixed_im,2)-1);
rbf_weights = rbf([XX(:) YY(:)],rbf_center,rbf_radius);
def_XX = XX(:) + rbf_weights(:).*rbf_coeff(1);
def_YY = YY(:) + rbf_weights(:).*rbf_coeff(2);

def_im = interp2(XX,YY,moving_im,def_XX,def_YY);
def_im(isnan(def_im)) = 1;
moving_im = reshape(def_im,256,[]);

%%
minmax_fixed = [min(fixed_im(:)) max(fixed_im(:))]
minmax_moving = [min(moving_im(:)) max(moving_im(:))]

nBins = 32;

step_fixed = (minmax_fixed(2)-minmax_fixed(1)) ./ (nBins-1);
step_moving = (minmax_moving(2)-minmax_moving(1)) ./ (nBins-1);

fixed_lut = floor((fixed_im-minmax_fixed(1))/step_fixed);
moving_lut = floor((moving_im-minmax_moving(1))/step_moving);

mi(fixed_lut,moving_lut,nBins)

%% MATLAB MI
mutual_info = [];
el_time = [];
for ii = -24:0.5:24
    rbf_coeff(2) = ii;
    def_YY = YY(:) + rbf_weights(:).*rbf_coeff(2);
    def_im = interp2(XX,YY,moving_im,def_XX,def_YY);
    def_im(isnan(def_im)) = minmax_moving(1);
    def_lut = floor((def_im-minmax_moving(1))/step_moving);
    subplot(1,2,1);imagesc(fixed_im);
    subplot(1,2,2);imagesc(reshape(def_im,256,[]));
    drawnow;
    tic
    mutual_info = [mutual_info mi(fixed_lut,def_lut,nBins)];
    el_time = [el_time toc];
end

%% IPP MI
fixed_ipp_im = MEXAllocateAlignedMemory(prod(size(fixed_im)),0);
moving_ipp_im = MEXAllocateAlignedMemory(prod(size(moving_im)),0);
MEXCopyArrayToAlignedMemory(fixed_im,fixed_ipp_im);
MEXCopyArrayToAlignedMemory(moving_im,moving_ipp_im);

in_fixed_ipp_im = MEXAllocateAlignedMemory(size(find(rbf_weights>0)),0);
out_fixed_ipp_im = MEXAllocateAlignedMemory(size(find(rbf_weights==0)),0);

in_moving_ipp_im = MEXAllocateAlignedMemory(size(find(rbf_weights>0)),0);
out_moving_ipp_im = MEXAllocateAlignedMemory(size(find(rbf_weights==0)),0);

MEXCopyArrayToAlignedMemory(fixed_im(rbf_weights>0),in_fixed_ipp_im);
MEXCopyArrayToAlignedMemory(fixed_im(rbf_weights==0),out_fixed_ipp_im);

MEXCopyArrayToAlignedMemory(moving_im(rbf_weights==0),out_moving_ipp_im);

nBins = 32;
values = MEXAllocateAlignedMemory(nBins+1,0);
fixed_levels = MEXAllocateAlignedMemory(nBins+1,0);
moving_levels = MEXAllocateAlignedMemory(nBins+1,0);

% THIS FUNCTION NEEDS TO BE SPLIT UP FOR THE TWO IMAGES %
MEXIPPCalculateImageBins(fixed_ipp_im,values,fixed_levels, ...
                         nBins);
MEXIPPCalculateImageBins(moving_ipp_im,values,moving_levels, ...
                         nBins);
fixed_entropy = MEXIPPCalculateFixedImageEntropy(fixed_ipp_im,values, ...
                                                 fixed_levels,nBins);

sample_in_pdf = MEXAllocateAlignedMemory(nBins*nBins,0);
sample_out_pdf = MEXAllocateAlignedMemory(nBins*nBins,0);

MEXIPPCalculateSampleJointPDF(out_fixed_ipp_im,out_moving_ipp_im,values,fixed_levels, ...
                              moving_levels,sample_out_pdf,nBins);

                          
moving_pdf = MEXAllocateAlignedMemory(nBins,1);
ln_joint_pdf = MEXAllocateAlignedMemory(nBins*nBins,0);
ln_moving_pdf = MEXAllocateAlignedMemory(nBins,1);

mutual_info_ipp = [];
el_time_ipp = [];
for i = -24:0.5:24
    rbf_coeff(2) = i;
    def_YY = YY(:) + rbf_weights(:).*rbf_coeff(2);
    def_im = interp2(XX,YY,moving_im,def_XX(rbf_weights>0),def_YY(rbf_weights>0));
    def_im(isnan(def_im)) = minmax_moving(1);

    
    subplot(1,2,1);imagesc(fixed_im);
    temp_im = moving_im;
    temp_im(rbf_weights>0) = def_im;
    subplot(1,2,2);imagesc(reshape(temp_im,256,[]));
    drawnow;
    tic;
    MEXCopyArrayToAlignedMemory(def_im,in_moving_ipp_im);
    MEXIPPCalculateSampleJointPDF(in_fixed_ipp_im,in_moving_ipp_im, ...
        values,fixed_levels,moving_levels, ...
        sample_in_pdf,nBins);
    mutual_info_ipp = [mutual_info_ipp 
        MEXIPPCalculateLimitedMemoryMutualInformation(fixed_entropy, sample_in_pdf, ...
        sample_out_pdf, moving_pdf, ...
        ln_joint_pdf, ln_moving_pdf, ...
        prod(size(fixed_im)), 32)];
    el_time_ipp = [el_time_ipp toc];                                  
end

MEXFreeAlignedMemory(fixed_ipp_im);
MEXFreeAlignedMemory(moving_ipp_im);
MEXFreeAlignedMemory(in_fixed_ipp_im);
MEXFreeAlignedMemory(out_fixed_ipp_im);
MEXFreeAlignedMemory(in_moving_ipp_im);
MEXFreeAlignedMemory(out_moving_ipp_im);
MEXFreeAlignedMemory(values);
MEXFreeAlignedMemory(fixed_levels);
MEXFreeAlignedMemory(moving_levels);
MEXFreeAlignedMemory(sample_in_pdf);
MEXFreeAlignedMemory(sample_out_pdf);
MEXFreeAlignedMemory(moving_pdf);
MEXFreeAlignedMemory(ln_joint_pdf);
MEXFreeAlignedMemory(ln_moving_pdf);

%% Plot the results
figure(2)
plot(1:length(mutual_info),mutual_info,1:length(mutual_info_ipp),mutual_info_ipp);
legend(['MATLAB: ' num2str(sum(el_time),'%3.2f') ' sec'],['IPP: ' num2str(sum(el_time_ipp),'%3.2f') ' sec'])

disp(['RBF size: ' num2str(sum(rbf_weights>0)./prod(size(rbf_weights)).*100) '% of image.']);
disp(['Speed up: ' num2str((1-sum(el_time_ipp)./sum(el_time)).*100) '%.']);