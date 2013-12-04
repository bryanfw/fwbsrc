function GRAD_IM = vuGradientMagnitude(X, varargin)
% vuGradientMagnitude(...) calculate the gradient magnitude for each 
% pixel in X
%
%   SYNTAX:
%       GRAD_IM = vuGradientMagnitude(X)
%       GRAD_IM = vuGradientMagnitude(X, options)       
%
%   OPTIONS & DEFAULTS:
%       verbose = 0
%
%   OUTPUTS:
%       GRAD_IM is the filtered image
%
%   EXAMPLE:
%       im = imread('cameraman.tif');
%       grad_im = vuGradientMagnitude(im,'verbose',1);
%       figure;
%       subplot(1,2,1);imshow(im,colormap(gray(256)));title('Image Before Filtering');
%       subplot(1,2,2);imshow(curv_im,colormap(gray(256)));title('Image After Filtering');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 1
  error('MATLAB:vuGradientMagnitude:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuGradientMagnitude:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 3)
  error('MATLAB:vuGradientMagnitude:UnknownDims', 'vuGradientMagnitude can only handle images of 2 or 3 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuGradientMagnitude';
p.parse(varargin{:});


% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);

% Row-major for ITK
X = vuRowMajorColumnMajorSwitch(X);

if isStruct
    if length(X.Dims) == 2
        GRAD_IM = MEX2DGradientMagnitude(X,p.Results.verbose);
    else
        GRAD_IM = MEX3DGradientMagnitude(X, p.Results.verbose);
    end
else
    if length(X.Dims) == 2
        grad_image = MEX2DGradientMagnitude(X, p.Results.verbose);
        GRAD_IM = grad_image.Data;
    else
        grad_image = MEX3DGradientMagnitude(X, p.Results.verbose);
        GRAD_IM = grad_image.Data;
    end
end

% Copy Orientation if it exists
if (isfield(X,'Orientation'))
    GRAD_IM.Orientation = X.Orientation;
end
% Copy parmeters if they exist
if (isfield(X,'Parms'))
    GRAD_IM.Parms = X.Parms;
end

% Column-major for Matlab
GRAD_IM = vuRowMajorColumnMajorSwitch(GRAD_IM);