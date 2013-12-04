
function m = ift3(data,varargin)
%[ift3] does 3-D IFFT.
%
% USAGE:
%   m = ift3(data,varargin)
%
% INPUT:
%   data : 3-D k-space data.
%   varargin{1} : nROWS for output data
%   varargin{2} : nCOLS for output data
%   varargin{3} : nPAGES for output data
%
% See also, ift1 ift2
%
% Last modified
% 2012.02.01.
%   Generate this function from ift2.
%
% Ha-Kyu


if nargin==1
    % do this silly business so that the third dim isn't fftshifted
    m = fftshift(fftshift(ifftn(ifftshift(ifftshift(data,2),1)),2),1); % use this
    %m = ifftshift(ifftshift(ifftn(fftshift(fftshift(data,2),1)),2),1);
elseif nargin==3
    nROWS = varargin{1};
    nCOLS = varargin{2};
    nPAGES = varargin{3};
    m = fftshift(fftshift(ifftn(ifftshift(ifftshift(data,2),1),[nROWS,nCOLS,nPAGES]),2),1); % use this
    %m = ifftshift(ifftshift(ifftn(fftshift(fftshift(data,2),1),[nROWS,nCOLS,nPAGES]),2),1);
elseif nargin~=1 && nargin~=3
    error('Number of input arguments must be 1 or 3.')
end






