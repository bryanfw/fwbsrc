function SEGM_IM = vuWatershedSegmentation(X, varargin)
% vuWatershedSegmentation(...) performs segmentation on image X, by taking
% the gradient image of X, and performing a watershed segmentation algorithm
%
%   SYNTAX:
%       SEGM_IM = vuWatershedSegmentation(X)
%       SEGM_IM = vuWatershedSegmentation(X, options)     
%
%   OPTIONS & DEFAULTS:
%       level = 0.001;
%       threshold = 0.15;
%       verbose = 0
%
%   OUTPUTS:
%       SEGM_IM is the segmented image (RGB format)
%
%   EXAMPLE:
%       im = load('Brain4.dat');
%       conn_im = vuWatershedSegmentation(im);
%       figure;
%       subplot(1,2,1);imshow(im,colormap(gray(256)));title('Image Before Segmenting');
%       subplot(1,2,2);imshow(conn_im,colormap(gray(256)));title('Image After Segmenting');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science


if nargin < 1
  error('MATLAB:vuWatershedSegmentation:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuWatershedSegmentation:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 2)
  error('MATLAB:vuWatershedSegmentation:UnknownDims', 'vuWatershedSegmentation can only handle images of 2 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('level',0.001,@(x) (x>=0)&&(x<=1)&&isa(x,'double'));
p.addParamValue('threshold',0.15,@(x) (x>=0)&&(x<=1)&&isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuWatershedSegmentation';
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
        SEGM_IM = MEX2DWatershed(X, p.Results.level, p.Results.threshold, p.Results.verbose);
        % Switch for Matlab RGB format
        SEGM_IM.Data = permute(SEGM_IM.Data,[3 2 1]);
    end
else
    if length(X.Dims) == 2
        segm_image = MEX2DWatershed(X, p.Results.level, p.Results.threshold, p.Results.verbose);
        % Switch for Matlab RGB format
        segm_image.Data = permute(segm_image.Data,[3 2 1]);
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
