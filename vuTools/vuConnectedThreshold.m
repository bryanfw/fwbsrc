function SEGM_IM = vuConnectedThreshold(X, seed, lower, upper, varargin)
% vuConnectedThreshold(...) performs segmentation on image X,
% by thresholding the pixels connected to 'seed', [x y z], with the lower and upper
% limits
%
%   SYNTAX:
%       SEGM_IM = vuConnectedThreshold(X, seed, lower, upper)
%       SEGM_IM = vuConnectedThreshold(X, seed, lower, upper, options)     
%
%   OPTIONS & DEFAULTS:
%       replacementValue = 1
%       numIterations = 5
%       timeStep = 0.125
%       verbose = 0
%
%   OUTPUTS:
%       SEGM_IM is the segmented image
%
%   EXAMPLE:
%       im = load('Brain4.dat');
%       seed = [165 118];
%       lower = 37;
%       upper = 50;
%       conn_im = vuConnectedThreshold(im,seed,lower,upper);
%       figure;
%       subplot(1,2,1);imshow(im,colormap(gray(256)));title('Image Before Segmenting');
%       subplot(1,2,2);imshow(conn_im,colormap(gray(256)));title('Image After Segmenting');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science


if nargin < 4
  error('MATLAB:vuConnectedThreshold:NotEnoughInputs', 'Not enough input arguments.');
elseif nargin > 8
  warning('MATLAB:vuConnectedThreshold:TooManyInputs', 'Too many arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuConnectedThreshold:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 3)
  error('MATLAB:vuConnectedThreshold:UnknownDims', 'vuConnectedThreshold can only handle images of 2 or 3 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('numIterations',5,@(x) isa(x,'double'));
p.addParamValue('replacementValue',1,@(x) isa(x,'double'));
p.addParamValue('timeStep',0.125,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuConnectedThreshold';
p.parse(varargin{:});

% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);
lower = double(lower);
upper = double(upper);
seed = double(seed);

% Row-major for ITK
X = vuRowMajorColumnMajorSwitch(X);

if isStruct
    if length(X.Dims) == 2
        SEGM_IM = MEX2DConnectedThreshold(X, seed, lower, upper, p.Results.replacementValue, p.Results.numIterations, p.Results.timeStep, p.Results.verbose);
    else
        SEGM_IM = MEX3DConnectedThreshold(X, seed, lower, upper, p.Results.replacementValue, p.Results.numIterations, p.Results.timeStep, p.Results.verbose);
    end
else
    if length(X.Dims) == 2
        segm_image = MEX2DConnectedThreshold(X, seed, lower, upper, p.Results.replacementValue, p.Results.numIterations, p.Results.timeStep, p.Results.verbose);
        SEGM_IM = segm_image.Data;
    else
        segm_image = MEX3DConnectedThreshold(X, seed, lower, upper, p.Results.replacementValue, p.Results.numIterations, p.Results.timeStep, p.Results.verbose);
        SEGM_IM = segm_image.Data;
    end
end

% Copy Orientation if it exists
if (isfield(X,'Orientation'))
    SEGM_IM.Orientation = X.Orientation;
end
% Copy parmeters if they exist
if (isfield(X,'Parms'))
    SEGM_IM.Parms = X.Parms;
end

% Column-major for Matlab
SEGM_IM = vuRowMajorColumnMajorSwitch(SEGM_IM);