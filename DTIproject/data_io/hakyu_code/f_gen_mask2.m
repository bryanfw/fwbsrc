
function mask = f_gen_mask2(data,thresh,varargin)
%[f_gen_mask2] use threshold to generate mask which has the same dimension
% as the input data.
%
% Usage:
%   mask = gen_mask(data,thresh)
%   mask = gen_mask(data,thresh,'erode',5)
%   mask = gen_mask(data,thresh,'dilate',7,'ball')
%
% Input:
%   data        Input data to be thresholded. It uses abs value.
%   thresh      Threshold to be used for generating mask
%   varargin
%       varargin{1}     Text string, 'erode' or 'dilate'
%       varargin{2}     Number of voxels to erode or dilate
%       varargin{3}     Structure element, 'disk', 'line', 'ball' etc
%
% Note:
%   This function is from [gen_mask.m].
%
% See also, f_gen_mask, gen_mask
%
%
% Last modified
% 2011.01.19.


%% Get size
% Get size of the data up to 7-D.
[nd1,nd2,nd3,nd4,nd5,nd6,nd7] = size(data);

%% Get input
if nargin == 2
    str = [];
elseif nargin == 5
    str = varargin{1};
    nvox = varargin{2};
    el = varargin{3};
else
    error('gen_mask:main','Number of input arguments must be 2 or 5')
end

%% Generate mask
mask = abs(data)>thresh;
if isempty(str)
    % do nothing
else
    
    for ind3=1:nd3
        for ind4=1:nd4
            for ind5=1:nd5
                for ind6=1:nd6
                    for ind7=1:nd7
                        if strcmpi(str,'erode')
                            mask(:,:,ind3,ind4,ind5,ind6,ind7) = ...
                                imerode(abs(data(:,:,ind3,ind4,ind5,ind6,ind7)),strel(el,nvox));
                        elseif strcmpi(str,'dilate')
                            mask(:,:,ind3,ind4,ind5,ind6,ind7) = ...
                                imdilate(abs(data(:,:,ind3,ind4,ind5,ind6,ind7)),strel(el,nvox));
                        else
                            error('gen_mask:main','Unknown str')
                        end
                    end
                end
            end
        end
    end
    
end





