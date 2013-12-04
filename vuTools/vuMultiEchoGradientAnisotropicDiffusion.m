function GRAD_IM = vuMultiEchoGradientAnisotropicDiffusion(X, varargin)
% vuMultiEchoGradientAnisotropicDIFFUSION performs gradient anisotropic diffusion
% filtering on an mult-echo image.
%
%   SYNTAX:
%       GRAD_IM = vuMultiEchoGradientAnisotropicDiffusion(X)
%       GRAD_IM = vuMultiEchoGradientAnisotropicDiffusion(X, options)
%
%   OPTIONS & DEFAULTS:
%       numIterations = 5
%       timeStep = 0.125 (for 2D mult-echo images) or 0.0625 (for mult-echo 3D images)
%       condParameter = 3.0
%       verbose = 0
%
%   OUTPUTS:
%       GRAD_IM is the filtered image
%
%   EXAMPLE:
%       im = imread('cameraman.tif');
%				im(:,:,2) = im;
%       grad_im = vuMultiEchoGradientAnisotropicDiffusion(im);
%       figure;
%       subplot(1,2,1);imshow(im(:,:,1),colormap(gray(256)));title('Image Before Filtering');
%       subplot(1,2,2);imshow(grad_im(:,:,1),colormap(gray(256)));title('Image After Filtering');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science


if nargin < 1
  error('MATLAB:vuMultiEchoGradientAnisotropicDiffusion:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuMultiEchoGradientAnisotropicDiffusion:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  isStruct = false;
  disp('X is a matrix')
end

if (length(X.Dims) < 3 || length(X.Dims) > 4)
  error('MATLAB:vuMultiEchoGradientAnisotropicDiffusion:UnknownDims', 'vuMultiEchoGradientAnisotropicDiffusion can only handle mult-echo images of 2 or 3 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('numIterations',5,@(x) isa(x,'double'));
p.addParamValue('timeStep',0.125,@(x) isa(x,'double'));
p.addParamValue('condParameter',3.0,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuMultiEchoGradientAnisotropicDiffusion';
p.parse(varargin{:});

% Change Default timestep for 3D 
if(~isempty(p.UsingDefaults))
    if ((~isempty(regexpi(cell2mat(p.UsingDefaults),'timeStep'))) && (length(X.Dims)==4))
        timeStep = 0.0625;
    else
        timeStep = p.Results.timeStep;
    end
end

% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);

% Row-major for ITK
X = vuRowMajorColumnMajorSwitch(X);

if isStruct
    if length(X.Dims) == 3
        GRAD_IM = MEX2DMultiEchoGradientAnisotropicDiffusion(X, p.Results.numIterations, timeStep, p.Results.condParameter, p.Results.verbose);
    else
        GRAD_IM = MEX3DMultiEchoGradientAnisotropicDiffusion(X, p.Results.numIterations, timeStep, p.Results.condParameter, p.Results.verbose);
    end
else
    if length(X.Dims) == 3
        grad_image = MEX2DMultiEchoGradientAnisotropicDiffusion(X, p.Results.numIterations, timeStep, p.Results.condParameter, p.Results.verbose);
        GRAD_IM = grad_image.Data;
    else
        grad_image = MEX3DMultiEchoGradientAnisotropicDiffusion(X, p.Results.numIterations, timeStep, p.Results.condParameter, p.Results.verbose);
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