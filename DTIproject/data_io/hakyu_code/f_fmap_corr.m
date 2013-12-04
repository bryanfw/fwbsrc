
function data = f_fmap_corr(I,F,BW,PEdir,blipPolarity,orientation)
%[f_fmap_corr] corrects B0 field inhomogeneity using B0 fieldmap. This is
%for SENSE and multi-shot data.
%
% USAGE
%   data =
%   f_fmap_corr(I,F,FOV,sense_fac,ovs_fac,acq_vox_p,Ns,Tro,PEdir,blipPolar
%   ity,orientation)
%
% INPUT
%   I   : Input image, 2D or 3D 
%   F   : B0 field map including Gradient Echo (GRE) images, 3D or 4D
%   BW  : Bandwidth along phase-encoding direction, Hz 
%   PEdir   : Phase encoding (PE) direction, Y=1, X=2, Z=3 in [ap,rl,fh]
%             direction. 
%   blipPolarity    : Positive or negative PEdir, 1 or -1.
%   orientation     : Slice orientation, 1:axial, 2:coronal, 3:sagital
%
% NOTE
%   This doesn't use Jacobian determinant multiplication.
%
%
% Last modified
%   2009.09.14. v1. Simplified (and different) version of
%   [correctB0field.m]. 2009.09.15. v2. Modified from v1, [fmap_corr.m],
%   based on [test_susceptibility.m].
%                   Now don't use [fmap_corr.m].
%   2009.11.19. v3. Change the sign of blipPolarity (out of this function)
%   and DeltaR. 2010.02.02. v4. Take pre-calculated bandwidth as input.
%   2010.08.05. Put image size validation.
% 2010.08.09
%   This is generated from [fmap_corr_v4.m].
% 2010.12.08.
%   Modified 'PEdir'.



%% Preparation
[ny,nx,ns] = size(I);

% Check image size.
if numel(size(I))==2
    if numel(size(F))~=3
        error('DIM of I is %d, then DIM of F must be %d', ...
            numel(size(I)),numel(size(I))+1)
    end
end
if numel(size(I))==3
    if numel(size(F))~=4
        error('DIM of I is %d, then DIM of F must be %d', ...
            numel(size(I)),numel(size(I))+1)
    end
end
if ~(numel(size(I))==2 || numel(size(I))==3)
    error('DIM of I must be 2D or 3D')
end

if numel(size(F))==3
    if size(F,3)~=2
        error('DIM of F is %d, then size(F,3) must be %d', ...
            numel(size(F)),2)
    end
end
if numel(size(F))==4
    if size(F,4)~=2
        error('DIM of F is %d, then size(F,4) must be %d', ...
            numel(size(F)),2)
    end
end        



%% Generate B0 inhomogeneity corrected image

% (FOV/R/Ns) / acq_vox_p * ovs_fac = EPI_FACTOR (IMG)

% Calculate voxel displacement map.
if numel(size(F))==3
    gamDeltaB0 = blipPolarity * F(:,:,2);
end
if numel(size(F))==4
    gamDeltaB0 = blipPolarity * F(:,:,:,2);
end
DeltaR = gamDeltaB0 / BW;   % voxel displacement

% Correct B0 inhomogeneity.
Icorr = zeros(size(I));
y_v = -floor(ny/2):(ceil(ny/2)-1);
x_v = -floor(nx/2):(ceil(nx/2)-1);
[x_m,y_m] = meshgrid(x_v,y_v);  % cm

h = waitbar(0,'Please wait ...');
for indSlice = 1:ns
    wbar_s = sprintf('Processing page %d',indSlice);
    waitbar(indSlice/ns,h,wbar_s)
    
    % *** Take PE direction.
%     if orientation==1 && PEdir==1
%         xp_m = x_m;
%         yp_m = y_m+DeltaR(:,:,indSlice);
%     elseif orientation==1 && PEdir==2
%         xp_m = x_m+DeltaR(:,:,indSlice);
%         yp_m = y_m;
%     elseif orientation==2 && PEdir==2
%         xp_m = x_m+DeltaR(:,:,indSlice);
%         yp_m = y_m;
%     elseif orientation==2 && PEdir==3
%         xp_m = x_m;
%         yp_m = y_m+DeltaR(:,:,indSlice);
%     elseif orientation==3 && PEdir==1
%         xp_m = x_m+DeltaR(:,:,indSlice);
%         yp_m = y_m;
%     elseif orientation==3 && PEdir==3
%         xp_m = x_m;
%         yp_m = y_m+DeltaR(:,:,indSlice);
%     else
%         error('f_fmap_corr:main','Unknown orientation and PEdir.')
%     end
    if PEdir==1
        xp_m = x_m;
        yp_m = y_m+DeltaR(:,:,indSlice);
    end
    Icorr(:,:,indSlice) = interp2(x_m,y_m,I(:,:,indSlice),xp_m,yp_m,'cubic');
end
close(h)
Icorr(isnan(Icorr(:))) = 0;

% Output data.
data.Icorr = Icorr;
data.DeltaR = DeltaR;
data.version = '[f_fmap_corr.m]';





%% Subfunctions

%_______________________________________________________________________
function [fx,fy,fz] = sub_gradientForBoundary(I,step_v,flagMask)
% [sub_gradientForBoundary] calculate gradient especially when gradient 
% bumps are around object boundary when mask is used.
% 
% Usage:
%   [fx,fy,fz] = sub_gradientForBoundary(I,step_v,flagMask)
%
% Input:
%   I           :   Input image, 2D or 3D
%   step_v      :   Gradient step, [x,y,z] direction
%   flagMask    :   Flag for use of mask, 1:use, 0:no use
%
% Output:
%   [fx,fy,fz]  :   X, Y and Z directional gradient
%
% Last modified:
%   10/25/07.
% HKJ



% Check input.
nDim = length(size(I));
if nDim > 3 || nDim < 2
    error('f_fmap_corr:sub_gradientForBoundary','Input image must be 2D or 3D.')
end
I(isnan(I)) = 0;
[ny,nx,nz] = size(I);    % nz=1 for 2D image

if length(step_v)<2 || length(step_v)>3
    error('f_fmap_corr:sub_gradientForBoundary','Check step_v.')
end

% Calculate gradient.
if nDim==2
    [fx,fy] = gradient(I,step_v(1),step_v(2));
    fz = 0;
else
    [fx,fy,fz] = gradient(I,step_v(1),step_v(2),step_v(3));
end

% Check if the image is masked.
if flagMask~=1
    return;
end

% Treat boundary.
for indz = 1:nz
    %   *** Y direction.
    for indy = 1:ny
        indNoZero = find(I(:,indy,indz));
        if isempty(indNoZero)
            continue
        end
        if length(indNoZero)==1
            fy(max(indNoZero-1,1):min(indNoZero+1,ny),indy,indz) = 0;
            continue
        end
        if indNoZero(1) > 1
            fy(indNoZero(1)-1,indy,indz) = 0;
        end
        fy(indNoZero(1),indy,indz) = ...
            (I(indNoZero(1)+1,indy,indz)-I(indNoZero(1),indy,indz))/step_v(2);
        if indNoZero(end) < ny        
            fy(indNoZero(end)+1,indy,indz) = 0;
        end
        fy(indNoZero(end),indy,indz) = ...
            (I(indNoZero(end),indy,indz)-I(indNoZero(end)-1,indy,indz))/step_v(2);
    end
    
    %   *** X direction.
    for indx = 1:nx
        indNoZero = find(I(indx,:,indz));
        if isempty(indNoZero)
            continue
        end
        if length(indNoZero)==1
            fx(indx,max(indNoZero-1,1):min(indNoZero+1,nx),indz) = 0;
            continue;
        end
        if indNoZero(1) > 1
            fx(indx,indNoZero(1)-1,indz) = 0;
        end
        fx(indx,indNoZero(1),indz) = ...
            (I(indx,indNoZero(1)+1,indz)-I(indx,indNoZero(1),indz))/step_v(1);
        if indNoZero(end) < nx
            fx(indx,indNoZero(end)+1,indz) = 0;
        end 
        fx(indx,indNoZero(end),indz) = ...
            (I(indx,indNoZero(end),indz)-I(indx,indNoZero(end)-1,indz))/step_v(1);
    end
    
end

if nDim==3
    %   *** Z direction.
    for indx = 1:nx
        for indy = 1:ny
            indNoZero = find(I(indx,indy,:));
            if isempty(indNoZero)
                continue
            end
            if length(indNoZero)==1
                fz(indx,indy,max(indNoZero-1,1):min(indNoZero+1,nx)) = 0;
            else
                if indNoZero(1) > 1
                    fz(indx,indy,indNoZero(1)-1) = 0;
                end
                fz(indx,indy,indNoZero(1)) = ...
                    (I(indx,indy,indNoZero(1)+1)-I(indx,indy,indNoZero(1)))/step_v(3);
                if indNoZero(end) < nz
                    fz(indx,indy,indNoZero(end)+1) = 0;
                end
                fz(indx,indy,indNoZero(end)) = ...
                    (I(indx,indy,indNoZero(end))-I(indx,indy,indNoZero(end)-1))/step_v(3);
            end
        end
    end
 
end


return
            





%% END










