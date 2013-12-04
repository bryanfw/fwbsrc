function BLUR_IM = vuRecursiveGaussianBlur(X, varargin)
% vuRecursiveGaussianBlur performs a convolution with an appoximated
% Gaussian Kernal
%
%   SYNTAX:
%       BLUR_IM = vuRecursiveGaussianBlur(X)
%       BLUR_IM = vuRecursiveGaussianBlur(X, options)
%
%   OPTIONS & DEFAULTS:
%       sigma = 5 (5% width of X)
%       verbose = 0
%
%   OUTPUTS:
%       BLUR_IM is the blurred image
%
%   EXAMPLE:
%       im = imread('cameraman.tif');
%       blur_im = vuRecursiveGaussianBlur(im);
%       figure;
%       subplot(1,2,1);imshow(im,colormap(gray(256)));title('Image Before Filtering');
%       subplot(1,2,2);imshow(blur_im,colormap(gray(256)));title('Image After Filtering');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 1
  error('MATLAB:vuRecursiveGaussianBlur:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage structs.')
  else
      error('MATLAB:vuRecursiveGaussianBlur:InvalidStruct', 'The image structure is invalid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 3)
  error('MATLAB:vuRecursiveGaussianBlur:UnknownDims', 'vuRecursiveGaussianBlur can only handle images of 2 or 3 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('sigma',5,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuRecursiveGaussianBlur';
p.parse(varargin{:});

% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);

% Row-major for ITK
X = vuRowMajorColumnMajorSwitch(X);

% Get sigma from percentage to size of image
sigma = (p.Results.sigma/100)*(X.Dims(1)*X.Spc(1));

if isStruct
    if length(X.Dims) == 2
        BLUR_IM = MEX2DRecursiveGaussianBlur(X, sigma, p.Results.verbose);
    else
        BLUR_IM = MEX3DRecursiveGaussianBlur(X, sigma, p.Results.verbose);
    end
else
    if length(X.Dims) == 2
        blur_image = MEX2DRecursiveGaussianBlur(X, sigma, p.Results.verbose);
        BLUR_IM = blur_image.Data;
    else
        blur_image = MEX3DRecursiveGaussianBlur(X, sigma, p.Results.verbose);
        BLUR_IM = blur_image.Data;
    end
end

% Copy Orientation if it exists
if (isfield(X,'Orientation'))
    BLUR_IM.Orientation = X.Orientation;
end
% Copy parmeters if they exist
if (isfield(X,'Parms'))
    BLUR_IM.Parms = X.Parms;
end

% Column-major for Matlab
BLUR_IM = vuRowMajorColumnMajorSwitch(BLUR_IM);