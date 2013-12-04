function SIG_IM = vuSigmoidFunction(X, varargin)
% vuSigmoidFunction(...) calculates the sigmoid function with value
%	of alpha and beta, for each pixel in X.  The output image will have
% values between minValue and maxValue.
%
%   SYNTAX:
%       SIG_IM = vuSigmoidFunction(X)
%       SIG_IM = vuSigmoidFunction(X, options)       
%
%   OPTIONS & DEFAULTS:
%				alpha = 1
%				beta = 0
%				minValue = 0
%				maxValue = 1
%       verbose = 0
%
%   OUTPUTS:
%       SIG_IM is the filtered image
%
%   EXAMPLE:
%       im = imread('cameraman.tif');
%       grad_im = vuSigmoidFunction(im,'maxValue',256,'verbose',1);
%       figure;
%       subplot(1,2,1);imshow(im,colormap(gray(256)));title('Image Before Filtering');
%       subplot(1,2,2);imshow(curv_im,colormap(gray(256)));title('Image After Filtering');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 1
  error('MATLAB:vuSigmoidFunction:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuSigmoidFunction:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 3)
  error('MATLAB:vuSigmoidFunction:UnknownDims', 'vuSigmoidFunction can only handle images of 2 or 3 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('alpha',1,@(x) isa(x,'double'));
p.addParamValue('beta',0,@(x) isa(x,'double'));
p.addParamValue('minValue',0,@(x) isa(x,'double'));
p.addParamValue('maxValue',1,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuSigmoidFunction';
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
        SIG_IM = MEX2DSigmoidFunction(X, p.Results.alpha, p.Results.beta, p.Results.minValue, p.Results.maxValue, p.Results.verbose);
    else
        SIG_IM = MEX3DSigmoidFunction(X, p.Results.alpha, p.Results.beta, p.Results.minValue, p.Results.maxValue, p.Results.verbose);
    end
else
    if length(X.Dims) == 2
        sig_image = MEX2DSigmoidFunction(X, p.Results.alpha, p.Results.beta, p.Results.minValue, p.Results.maxValue, p.Results.verbose);
        SIG_IM = sig_image.Data;
    else
        sig_image = MEX3DSigmoidFunction(X, p.Results.alpha, p.Results.beta, p.Results.minValue, p.Results.maxValue, p.Results.verbose);
        SIG_IM = sig_image.Data;
    end
end

% Copy Orientation if it exists
if (isfield(X,'Orientation'))
    SIG_IM.Orientation = X.Orientation;
end
% Copy parmeters if they exist
if (isfield(X,'Parms'))
    SIG_IM.Parms = X.Parms;
end

% Column-major for Matlab
SIG_IM = vuRowMajorColumnMajorSwitch(SIG_IM);