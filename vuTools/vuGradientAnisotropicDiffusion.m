function GRAD_IM = vuGradientAnisotropicDiffusion(X, varargin)
% vuGRADIENTANISOTROPICDIFFUSION performs gradient anisotropic diffusion
% filtering on an image.
%
%   SYNTAX:
%       GRAD_IM = vuGradientAnisotropicDiffusion(X)
%       GRAD_IM = vuGradientAnisotropicDiffusion(X, options)
%
%   OPTIONS & DEFAULTS:
%       numIterations = 5
%       timeStep = 0.125 (for 2D images) or 0.0625 (for 3D images)
%       condParameter = 3.0
%       verbose = 0
%
%   OUTPUTS:
%       GRAD_IM is the filtered image
%
%   EXAMPLE:
%       im = imread('cameraman.tif');
%       grad_im = vuGradientAnisotropicDiffusion(im);
%       figure;
%       subplot(1,2,1);imshow(im,colormap(gray(256)));title('Image Before Filtering');
%       subplot(1,2,2);imshow(grad_im,colormap(gray(256)));title('Image After Filtering');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science


if nargin < 1
  error('MATLAB:vuGradientAnisotropicDiffusion:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuGradientAnisotropicDiffusion:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  isStruct = false;
  disp('X is a matrix')
end

if (length(X.Dims) < 2 || length(X.Dims) > 3)
  error('MATLAB:vuGradientAnisotropicDiffusion:UnknownDims', 'vuGradientAnisotropicDiffusion can only handle images of 2 or 3 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('numIterations',5,@(x) isa(x,'double'));
p.addParamValue('timeStep',0.125,@(x) isa(x,'double'));
p.addParamValue('condParameter',3.0,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuGradientAnisotropicDiffusion';
p.parse(varargin{:});

% Change Default timestep for 3D 
if(~isempty(p.UsingDefaults))
    if ((~isempty(regexpi(cell2mat(p.UsingDefaults),'timeStep'))) && (length(X.Dims)==3))
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
    if length(X.Dims) == 2
        GRAD_IM = MEX2DGradientAnisotropicDiffusion(X, p.Results.numIterations, timeStep, p.Results.condParameter, p.Results.verbose);
    else
        GRAD_IM = MEX3DGradientAnisotropicDiffusion(X, p.Results.numIterations, timeStep, p.Results.condParameter, p.Results.verbose);
    end
else
    if length(X.Dims) == 2
        grad_image = MEX2DGradientAnisotropicDiffusion(X, p.Results.numIterations, timeStep, p.Results.condParameter, p.Results.verbose);
        GRAD_IM = grad_image.Data;
    else
        grad_image = MEX3DGradientAnisotropicDiffusion(X, p.Results.numIterations, timeStep, p.Results.condParameter, p.Results.verbose);
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