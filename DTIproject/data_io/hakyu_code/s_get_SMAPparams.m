
%[s_get_SMAPparams] get SMAP parameters, SMAPparams.
%
% USAGE:
%   s_get_SMAPparams
%
% OUTPUT:
%   SMAPparams: A struct with fields as,
%       fitorder, 2D polynomial fit order
%       fitsize, Supporting region for fit, [outside,inside], of object
%       vox_dilate, Number of voxels to dilate mask
%       vox_erode, Number of voxels to erode mask
%       flip_slice_order, Zero or 1 to reverse slice order for REF and BODY image data.
%           This is because the 3D GRE k-space data has reversed slice order.
%
%
% Last modified
% 2010.08.06.
% 2010.10.06.
%   Add SMAPparams.fit_method for choosing between [local 2D polynomial
%   fit, thin-plate spline fit].
% 2012.05.14.
%   Add gen_bodycoil for selecting bodycoil equivalent between 'body' and
%   'ref' for using bodycoil and reference coil data, respectively.
% 2012.07.05.
%   Don't flip REF (and hence SMAP) slice order when 3D DWI is acquired.
%   See [f_coilsensemap.m] for implementing this.
%
% Ha-Kyu



%% Define SMAPparams

fitorder = 2;       % order of 2D polynomial fitting, (default:2)
fitsize_out = 5;    % supporting region for fit, outside object (rim region)
fitsize_in = 5;     % supporting region for fit, inside object
vox_dilate = 15;    % number of voxels to dilate mask
vox_erode = 1;      % number of voxels to erode mask

SMAPparams.fitorder = fitorder;
SMAPparams.fitsize = [fitsize_out,fitsize_in];
SMAPparams.vox_dilate = vox_dilate;
SMAPparams.vox_erode = vox_erode;
if strcmpi(DWIparams.scan_mode,'3D')
    SMAPparams.flip_slice_order = false;
else
    SMAPparams.flip_slice_order = true;   % reverse slice order for REF and BODY image data
end
SMAPparams.fit_method = smap_fit_method; % [l2p, tps] for [local 2D polynomial fit, thin-plate spline fit]
SMAPparams.gen_bodycoil = smap_gen_bodycoil; % ['body','ref']

clear  fitorder  fitsize_out  fitsize_in  vox_dilate  vox_erode  smap_fit_method




