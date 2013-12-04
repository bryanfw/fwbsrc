
function [S,BW] = f_coilsensemap(K_ref_4d,K_body_3d,opts_SMAP)
% [f_coilsensemap] generates coil sensitivity  map.
%
% INPUT
%   K_ref_4d:   Individual coil samples, [y,x,z,coil]
%   K_body_3d:  Bodycoil samples, [y,x,z]
%   opts_SMAP:  Options for generating coil senstivity map
%
% OUTPUT
%   S:  Coil sensitivity map, [y,x,coil] --> Then this function must be
%       repeated for each slice, z
%   BW: Mask used for generating the map
%
%
% Last modified
%   2009.05.01.
%   2009.07.17.     Use [sensitivity_map_v3.m].
%   2009.07.20.     Use flag_polyFit=2 in [sensitivity_map_v3.m].
%   2009.07.30.     Use [sensitivity_map_v4.m] considering Pruessmann's method
%                   in his original paper (1999). This is tested in [test_gen_smap_v2.m].
%   2009.09.17.     Apply flipdim for reversing slice order of REF and BODY image data.
%                   Slices (partitions) are from H (head) to F (foot) direction, so must be reversed.
%   2009.09.26.     Add [sensitivity_map_v5.m]. This version uses smoothing (whole) and 2D polynomial
%                   fitting on rim region.
%   2009.11.18.     Use [sensitivity_map_v6.m]. It doesn't calculate PSI.
%   2009.11.23.     Use [sensitivity_map_v8.m]. It uses Gaussian smoothing and
%                   2-D polynomial fitting only on outer rim regions.
%                   [sensitivity_map_v7.m] is test version of v8.
% 2010.08.09.
%   This is generated from [coilsensemap_v3.m].
% 2010.09.22.
%   Do circshift by -1 along slice direction for I_ref_4d and I_body_3d.
% 2010.09.30.
%   Add [f_sensitivity_map_tps.m] for generating sensitivity map using 2D
%   Thin-Plate Spline.
% 2010.10.06.
%   Add SMAPparams.fit_method for choosing between [local 2D polynomial
%   fit, thin-plate spline].
% 2010.10.29.
%   Do circshift by -1 along slice direction for I_ref_4d and I_body_3d
%   regardless of the number of kz.
% 2012.02.01.
%   Use [ift3.m] and [ft3.m].
%
% Ha-Kyu



%% Preparation

% TEST
% opts_SMAP = SMAPparams;

% Get data size.
[Ny,Nx,Nz,Nc] = size(K_ref_4d);

% Generate image-space data.
I_ref_4d = zeros(size(K_ref_4d));
for ind_coil = 1:Nc
%     I_ref_4d(:,:,:,ind_coil) = ifftshift(ifftshift(ifftn( ...
%         fftshift(fftshift(K_ref_4d(:,:,:,ind_coil),2),1) ),2),1);
    I_ref_4d(:,:,:,ind_coil) = ift3(K_ref_4d(:,:,:,ind_coil));
end
% if mod(size(I_ref_4d,3),2)==0
%     I_ref_4d = circshift(I_ref_4d,[0,0,-1,0]);
% end
I_ref_4d = circshift(I_ref_4d,[0,0,-1,0]);
% I_body_3d = ifftshift(ifftshift(ifftn( fftshift(fftshift(K_body_3d,2),1) ),2),1);
I_body_3d = ift3(K_body_3d);
% if mod(size(I_body_3d,3),2)==0
%     I_body_3d = circshift(I_body_3d,[0,0,-1]);
% end
I_body_3d = circshift(I_body_3d,[0,0,-1]);

% Reverse kz (slice) order. 2009.09.17
% opts_SMAP.flip_slice_order==0 if 3D DWI acquisition. 2012.07.05.
if opts_SMAP.flip_slice_order==1
    I_ref_4d = flipdim(I_ref_4d,3);
    I_body_3d = flipdim(I_body_3d,3);
    fprintf('      Slice order is reversed for REF and BODY.\n')
end

% Get target slice, which is opts_SMAP.kz_REF.
I_ref_3d = squeeze(I_ref_4d(:,:,opts_SMAP.kz_REF,:));
I_body_m = squeeze(I_body_3d(:,:,opts_SMAP.kz_REF));

% Clear unused data.
clear  K_ref_4d  K_body_3d  I_ref_4d  I_body_3d



%% Coil sensitivity map

% Generate smap.
fitorder = opts_SMAP.fitorder;           % order of 2D polynomial fitting
fitsize_out = opts_SMAP.fitsize(1);      % supporting region for fit, OUTSIDE object (rim region)
fitsize_in = opts_SMAP.fitsize(2);       % supporting region for fit, INSIDE object
vox_dilate = opts_SMAP.vox_dilate;       % number of voxels to dilate mask
vox_erode = opts_SMAP.vox_erode;         % number of voxels to erode mask, (default:3)

%%% Supporting region size is for one dimension (y or x) only. 
%%% Each dimension has 2*n+1 voxels.
%%% If fitsize_out=10, then [21x21] region is fitted for the rim region.
if strcmpi(opts_SMAP.fit_method,'l2p')
    % Use local 2D polynomial fitting by Pruessmann.
    %fprintf('      Sensitivity map generated for kz[%d] using [f_sensitivity_map.m].\n',opts_SMAP.kz_REF)
    [I_smap_3d,mask_m] = f_sensitivity_map(I_body_m, I_ref_3d, ...
        fitorder, fitsize_out, fitsize_in, vox_dilate);
elseif strcmpi(opts_SMAP.fit_method,'tps')
    % Use 2D Thin-Plate Spline (TPS).
    %fprintf('      Sensitivity map generated for kz[%d] using [f_sensitivity_map_tps.m].\n',opts_SMAP.kz_REF)
    [I_smap_3d,mask_m] = f_sensitivity_map_tps(I_body_m, I_ref_3d, vox_dilate);
else
    error('f_coilsensemap:main','Unknown ''SMAPparams.fit_method''.')
end

% Clear.
clear  I_ref_3d  I_body_m



%% Output
S = I_smap_3d;
BW = mask_m;

clear  I_smap_3d  mask_m  fitorder  fitsize_out  fitsize_in  vox_dilate  vox_erode
pack


%% END

