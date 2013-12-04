function SEGM_IM = vuVesselnessMeasure(X, seed, lower, upper, varargin)
% vuVesselnessMeasure(...) performs measure on image X,
% to find a vesselness measure through a Hessian matrix
%
%   SYNTAX:
%       SEGM_IM = vuVesselnessMeasure(X)
%       SEGM_IM = vuVesselnessMeasure(X, options)     
%
%   OPTIONS & DEFAULTS:
%       sigma = 1.0
%       alpha1 = 0.5
%       alpha2 = 2.0
%       verbose = 0
%
%   OUTPUTS:
%       SEGM_IM is the segmented image
%
%   EXAMPLE:
%       im = vuOpenImage('Vessels.dat');
%       seg_im = vuVesselnessMeasure(im,seed,lower,upper);
%       figure;
%       subplot(1,2,1);imshow(im,colormap(gray(256)));title('Image Before Segmenting');
%       subplot(1,2,2);imshow(seg_im,colormap(gray(256)));title('Image After Segmenting');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science


if nargin < 1
  error('MATLAB:vuVesselnessMeasure:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuVesselnessMeasure:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 3 || length(X.Dims) > 3)
  error('MATLAB:vuVesselnessMeasure:UnknownDims', 'vuVesselnessMeasure can only handle images of 3 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('sigma',1.0,@(x) isa(x,'double'));
p.addParamValue('alpha1',0.5,@(x) isa(x,'double'));
p.addParamValue('alpha2',2.0,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuVesselnessMeasure';
p.parse(varargin{:});

% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);

% Row-major for ITK
X = vuRowMajorColumnMajorSwitch(X);

if isStruct
    SEGM_IM = MEX3DVesselSegmentation(X, p.Results.sigma, p.Results.alpha1, p.Results.alpha2, p.Results.verbose);
else
    segm_image = MEX3DVesselSegmentation(X, p.Results.sigma, p.Results.alpha1, p.Results.alpha2, p.Results.verbose);
    SEGM_IM = segm_image.Data;
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