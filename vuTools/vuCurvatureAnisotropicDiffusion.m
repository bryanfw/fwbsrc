function CURV_IM = vuCurvatureAnisotropicDiffusion(X, varargin)
% vuCurvatureAnisotropicDiffusion(...) performs anisotropic diffusion on image X,
% using a modified curvature diffusion equation.
%
%   SYNTAX:
%       CURV_IM = vuCurvatureAnisotropicDiffusion(X)
%       CURV_IM = vuCurvatureAnisotropicDiffusion(X, options)       
%
%   OPTIONS & DEFAULTS:
%       numIterations = 5
%       timeStep = 0.125 (for 2D images) or 0.0625 (for 3D images)
%       condParameter = 3.0
%       verbose = 0
%
%   OUTPUTS:
%       CURV_IM is the filtered image
%
%   EXAMPLE:
%       im = imread('cameraman.tif');
%       curv_im = vuCurvatureAnisotropicDiffusion(im,'verbose',1);
%       figure;
%       subplot(1,2,1);imshow(im,colormap(gray(256)));title('Image Before Filtering');
%       subplot(1,2,2);imshow(curv_im,colormap(gray(256)));title('Image After Filtering');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 1
  error('MATLAB:vuCurvatureAnisotropicDiffusion:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuCurvatureAnisotropicDiffusion:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 3)
  error('MATLAB:vuCurvatureAnisotropicDiffusion:UnknownDims', 'vuCurvatureAnisotropicDiffusion can only handle images of 2 or 3 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('numIterations',5,@(x) isa(x,'double'));
p.addParamValue('timeStep',0.125,@(x) isa(x,'double'));
p.addParamValue('condParameter',3.0,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuCurvatureAnisotropicDiffusion';
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
        CURV_IM = MEX2DCurvatureAnisotropicDiffusion(X, p.Results.numIterations, timeStep, p.Results.condParameter, p.Results.verbose);
    else
        CURV_IM = MEX3DCurvatureAnisotropicDiffusion(X, p.Results.numIterations, timeStep, p.Results.condParameter, p.Results.verbose);
    end
else
    if length(X.Dims) == 2
        curv_image = MEX2DCurvatureAnisotropicDiffusion(X, p.Results.numIterations, timeStep, p.Results.condParameter, p.Results.verbose);
        CURV_IM = curv_image.Data;
    else
        curv_image = MEX3DCurvatureAnisotropicDiffusion(X, p.Results.numIterations, timeStep, p.Results.condParameter, p.Results.verbose);
        CURV_IM = curv_image.Data;
    end
end

% Copy Orientation if it exists
if (isfield(X,'Orientation'))
    CURV_IM.Orientation = X.Orientation;
end
% Copy parmeters if they exist
if (isfield(X,'Parms'))
    CURV_IM.Parms = X.Parms;
end

% Column-major for Matlab
CURV_IM = vuRowMajorColumnMajorSwitch(CURV_IM);