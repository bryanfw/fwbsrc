function DIST_IM = vuSignedDistanceMap(X, varargin)
% vuSignedDistanceMap(...) calculate the gradient magnitude for each 
% pixel in X
%
%   SYNTAX:
%       DIST_IM = vuSignedDistanceMap(X)
%       DIST_IM = vuSignedDistanceMap(X, options)       
%
%   OPTIONS & DEFAULTS:
%		insideValue = 0
%		outsideValue = 255
%       verbose = 0
%
%   OUTPUTS:
%       DIST_IM is the distance map image
%
%   EXAMPLE:
%       im = ones(256,256)*255;
%		im(75:175,75:175) = 0;
%       dist_im = vuSignedDistanceMap(im,'verbose',1);
%       figure;
%       subplot(1,2,1);imagesc(im);axis image;title('Binary Image');
%       subplot(1,2,2);imagesc(dist_im);axis image;colormap gray;title('Distance Map Image');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 1
  error('MATLAB:vuSignedDistanceMap:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuSignedDistanceMap:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 3)
  error('MATLAB:vuSignedDistanceMap:UnknownDims', 'vuSignedDistanceMap can only handle images of 2 or 3 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('insideValue',0,@(x) isa(x,'double'));
p.addParamValue('outsideValue',255,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuSignedDistanceMap';
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
        DIST_IM = MEX2DSignedDistanceMap(X,p.Results.insideValue,p.Results.outsideValue,p.Results.verbose);
    else
        DIST_IM = MEX3DSignedDistanceMap(X,p.Results.insideValue,p.Results.outsideValue, p.Results.verbose);
    end
else
    if length(X.Dims) == 2
        dist_image = MEX2DSignedDistanceMap(X,p.Results.insideValue,p.Results.outsideValue, p.Results.verbose);
        DIST_IM = dist_image.Data;
    else
        dist_image = MEX3DSignedDistanceMap(X,p.Results.insideValue,p.Results.outsideValue, p.Results.verbose);
        DIST_IM = dist_image.Data;
    end
end

% Copy Orientation if it exists
if (isfield(X,'Orientation'))
    DIST_IM.Orientation = X.Orientation;
end
% Copy parmeters if they exist
if (isfield(X,'Parms'))
    DIST_IM.Parms = X.Parms;
end

% Column-major for Matlab
DIST_IM = vuRowMajorColumnMajorSwitch(DIST_IM);