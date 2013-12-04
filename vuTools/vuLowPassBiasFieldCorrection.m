function [CORR_IM, MASK] = vuLowPassBiasFieldCorrection(X, varargin)
% vuLowPassBiasFieldCorrection performs a bias field correction estimation
% using a Gaussian blur as a low-pass filter
%
%   SYNTAX:
%       [CORR_IM, MASK] = vuLowPassBiasFieldCorrection(X)
%       [CORR_IM, MASK] = vuLowPassBiasFieldCorrection(X, options)
%
%       Quick correction of bias field in image X, using a mask with
%       threshold maskThres, and using a sigma (in %) blur
%
%   OPTIONS & DEFAULTS:
%       maskThres = 20
%       sigma = 10 (10% width of X)
%       verbose = 0
%
%   OUTPUTS:
%       CORR_IM is the corrected image
%
%   EXAMPLE:
%       im = load('Brain1.txt');
%       corr_im = vuLowPassBiasFieldCorrection(im);
%       figure;
%       subplot(1,2,1);imshow(im,colormap(gray(256)));title('Image Before Filtering');
%       subplot(1,2,2);imshow(corr_im,colormap(gray(256)));title('Image After Filtering');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science


if nargin < 1
  error('MATLAB:vuLowPassBiasFieldCorrection:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage structs.')
  else
      error('MATLAB:vuLowPassBiasFieldCorrection:InvalidStruct', 'The image structure is invalid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 3)
  error('MATLAB:vuLowPassBiasFieldCorrection:UnknownDims', 'vuLowPassBiasFieldCorrectioon can only handle images of 2 or 3 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('maskThres',20,@(x) isa(x,'double'));
p.addParamValue('sigma',10,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuLowPassBiasFieldCorrection';
p.parse(varargin{:});

% Store max pixel value
maxBefore = max(X.Data(:));

if (length(X.Dims)==2)
    % Generate automatic mask
    mask = zeros(X.Dims(2),X.Dims(1));
    for i = 1:X.Dims(2)
        for j = 1:X.Dims(1)
            if X.Data(i,j) > p.Results.maskThres
                mask(i,j)= 1;
            end
        end
    end
end

if (length(X.Dims)==3)
    % Generate automatic mask
    mask = zeros(X.Dims(2),X.Dims(1),X.Dims(3));
    for i = 1:X.Dims(2)
        for j = 1:X.Dims(1)
            for k = 1:X.Dims(3)
                if X.Data(i,j,k) > p.Results.maskThres
                    mask(i,j,k)= 1;
                end
            end
        end
    end
end

% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);

% Blur image
blur_im = vuRecursiveGaussianBlur(X,'sigma',p.Results.sigma,'verbose',p.Results.verbose);

% Invert each element
blur_im.Data = 1./blur_im.Data;

% Multiple out bias field
X.Data = X.Data.*blur_im.Data;

% Rescale pixel values
maxAfter = max(X.Data(:));
X.Data = X.Data.*(maxBefore/maxAfter);

% Apply mask
X.Data = X.Data.*mask;

% Return
if isStruct
    CORR_IM = X;
else
    CORR_IM = X.Data;
end

% Copy Orientation if it exists
if (isfield(X,'Orientation'))
    CORR_IM.Orientation = X.Orientation;
end
% Copy parmeters if they exist
if (isfield(X,'Parms'))
    CORR_IM.Parms = X.Parms;
end

MASK = mask;
