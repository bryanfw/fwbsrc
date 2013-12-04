
function k = f_gauss_kernel(sz,fwhm)
% [f_gauss_kernel] calculates Gaussian kernel.
%
% Note:
%   This function originates from [smoothkern.m].
%
% Usage:
%   k = f_gauss_kernel(sz,fwhm)
%
% Input:
%   sz      :   Kernel size, [X,Y,Z] in voxels up to 3D (can be double), vector
%   fwhm    :   FWHM in voxels (can be double), scalar
%
% Output:
%   k       :   Kernel, up to 3D
%
% See also smoothkern
%
% Last modified:
%   10/15/07.
% HKJ.




%% Check input

if length(sz) > 3
    error('Kernel size must be 1D, 2D or 3D.')
end
fwhm_half = fwhm/2;
sigma = fwhm_half/sqrt(2*log(2));   % standard deviation of Gaussian
if length(sz)==1
    [x] = -(sz-1)/2:(sz-1)/2;
    k = exp(-(x.*x)/(2*sigma*sigma));
    k = k/sum(k(:));
elseif length(sz)==2
    [x,y] = meshgrid(-(sz(2)-1)/2:(sz(2)-1)/2, -(sz(1)-1)/2:(sz(1)-1)/2);
    k = exp(-(x.*x + y.*y)/(2*sigma*sigma));
    k = k/sum(k(:));
else
    [x,y,z] = meshgrid(-(sz(2)-1)/2:(sz(2)-1)/2, -(sz(1)-1)/2:(sz(1)-1)/2, -(sz(3)-1)/2:(sz(3)-1)/2);
    k = exp(-(x.*x + y.*y + z.*z)/(2*sigma*sigma));
    k = k/sum(k(:));
end


