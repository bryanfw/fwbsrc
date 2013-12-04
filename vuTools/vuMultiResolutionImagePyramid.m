function MULT_IM = vuMultiResolutionImagePyramid(X, varargin)
% vuMultiResolutionImagePyramid(...) creates any array of images with different
% resolutions of image X.
%
%   SYNTAX:
%       MULT_IM = vuMultiResolutionImagePyramid(X)
%       MULT_IM = vuMultiResolutionImagePyramid(X, options)      
%
%   OPTIONS & DEFAULTS:
%       levels = 3
%       verbose = 0
%
%   OUTPUTS:
%       MULT_IM is the array of images
%
%   EXAMPLE:
%       im = imread('cameraman.tif');
%       mult_im = vuMultiResolutionImagePyramid(im);
%       figure;
%       imshow(mult_im(1).Data,colormap(gray(256)));title('Image 1 (Lowest Resolution)');
%       figure;
%       imshow(mult_im(2).Data,colormap(gray(256)));title('Image 2');
%       figure;
%       imshow(mult_im(3).Data,colormap(gray(256)));title('Image 3 (Full Resolution)');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science


if nargin < 1
  error('MATLAB:vuMultiResolutionImagePyramid:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuMultiResolutionImagePyramid:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 3)
  error('MATLAB:vuMultiResolutionImagePyramid:UnknownDims', 'vuMultiResolutionImagePyramid can only handle images of 2 or 3 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('levels',3,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuMultiResolutionImagePyramid';
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
        MULT_IM = MEX2DMultiResolutionImagePyramid(X, p.Results.levels, p.Results.verbose);
    else
        MULT_IM = MEX3DMultiResolutionImagePyramid(X, p.Results.levels, p.Results.verbose);
    end
else
    if length(X.Dims) == 2
        mult_image = MEX2DMultiResolutionImagePyramid(X, p.Results.levels, p.Results.verbose);
        for i = 1:p.Results.levels
            MULT_IM(i).Data = mult_image(i).Data;
        end
    else
        mult_image = MEX3DMultiResolutionImagePyramid(X, p.Results.levels, p.Results.verbose);
        for i = 1:p.Results.levels
            MULT_IM(i).Data = mult_image(i).Data;
        end
    end
end


for i = 1:p.Results.levels
    % Copy Orientation if it exists
    if (isfield(X,'Orientation'))
        MULT_IM(i).Orientation = X.Orientation;
    end
    % Copy parmeters if they exist
    if (isfield(X,'Parms'))
      MULT_IM(i).Parms = X.Parms;
    end
    % Column-major for Matlab
    MULT_IM(i) = vuRowMajorColumnMajorSwitch(MULT_IM(i));
end