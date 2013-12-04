function RESULT = vuHoughCircleTransform(X, varargin)
% vuHoughCircleTransform(...) performs a hough transform on image X,
% in order to find circles within the image.
%
%   SYNTAX:
%       RESULT = vuHoughCircleTransform(X)
%       RESULT = vuHoughCircleTransform(X, options)       
%
%   OPTIONS & DEFAULTS:
%       numCircles = 1
%       minRadius = 0
%       maxRadius = 1
%       threshold = 0
%       discRadiusRatio = 1
%       gradientSigma = 1
%       verbose = 0
%
%   OUTPUTS:
%       RESULT is the form Circle Center (x,y), and Radius
%
%   EXAMPLE:
%       im = imread('moon.tif');
%       RESULT = vuHoughCircleTransform(im,'maxRadius',300,'verbose',1);
%       
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 1
  error('MATLAB:vuHoughCircleTransform:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuHoughCircleTransform:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 2)
  error('MATLAB:vuHoughCircleTransform:UnknownDims', 'vuHoughCircleTransform can only handle images of 2 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('numCircles',1,@(x) isa(x,'double'));
p.addParamValue('minRadius',0,@(x) isa(x,'double'));
p.addParamValue('maxRadius',10,@(x) isa(x,'double'));
p.addParamValue('threshold',0,@(x) isa(x,'double'));
p.addParamValue('discRadiusRatio',1,@(x) isa(x,'double'));
p.addParamValue('gradientSigma',1,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuHoughCircleTransform';
p.parse(varargin{:});



% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);

% Row-major for ITK
X = vuRowMajorColumnMajorSwitch(X);

RESULT = MEX2DHoughCircleTransform(X, p.Results.numCircles, p.Results.minRadius, p.Results.maxRadius, p.Results.threshold, p.Results.discRadiusRatio, p.Results.gradientSigma, p.Results.verbose);


% Copy Orientation if it exists
if (isfield(X,'Orientation'))
    RESULT.Orientation = X.Orientation;
end
% Copy parmeters if they exist
if (isfield(X,'Parms'))
    RESULT.Parms = X.Parms;
end

% Column-major for Matlab
RESULT = vuRowMajorColumnMajorSwitch(RESULT);