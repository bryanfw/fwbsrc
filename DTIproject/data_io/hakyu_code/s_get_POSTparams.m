
%[s_get_POSTparams] get data post processing parameters, POSTparams.
%
% USAGE:
%   s_get_POSTparams
%
% OUTPUT:
%   s_get_POSTparams: A struct with fields related with,
%       Image co-registration - 
%           optim, Optimization method for co-registration, 'powell' or 'simplex'
%           method, Co-registration method, 'rigid' or 'affine'
%           cost, Cost function for co-registration, 'nmi','mi','ncc','cor','ecc
%       Tensor calculation - 
%           medfilt, Use of median filtering with [2,2] voxels, 'on' or 'off'
%
%
% Last modified
% 2011.01.28
%   Generated.
% 2011.08.24
%   Neighboring voxel size for medfilt is set as medfilt_kern_siz.
% 2012.02.13.
%   Consider not using co-registration.
%
% Ha-Kyu



%% Define POSTparams

% Check parameters for co-registration and tensor calculation.
if ~(strcmpi(optim_s,'powell') || ...
        strcmpi(optim_s,'simplex'))
    error('s_get_POSTparams:main','Unknown POSTparams.coreg.optim')
end
if ~(strcmpi(cost_s,'nmi') || ...
        strcmpi(cost_s,'cor') || ...
        strcmpi(cost_s,'mi') || ...
        strcmpi(cost_s,'ecc') || ...
        strcmpi(cost_s,'ncc'))
    error('s_get_POSTparams:main','Unknown POSTparams.coreg.cost')
end
if ~(strcmpi(method_s,'affine') || ...
        strcmpi(method_s,'rigid') || ...
        strcmpi(method_s,'none'))
    error('s_get_POSTparams:main','Unknown POSTparams.coreg.method')
end
if ~(strcmpi(medfilt_s,'on') || ...
        strcmpi(medfilt_s,'off'))
    error('s_get_POSTparams:main','Unknown POSTparams.tensor.medfilt')
end

% Image co-registration related.
POSTparams.coreg.optim = optim_s;
POSTparams.coreg.method = method_s;
POSTparams.coreg.cost = cost_s;

% Tensor calculation related.
POSTparams.tensor.medfilt = medfilt_s;
POSTparams.tensor.medfilt_kern_siz = medfilt_kern_siz;

% Clear.
clear  optim_s  method_s  cost_s  medfilt_s  medfilt_kern_siz


















