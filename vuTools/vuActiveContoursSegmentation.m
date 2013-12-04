function SEGM_IM = vuActiveContoursSegmentation(IM, INIT, varargin)
% vuActiveContoursSegmentation segments regions within an image by driving 
% a region-based active contour through a local binary fitting energy.
%   SYNTAX:
%       SEGM_IM = vuActiveContoursSegmentation(IM, INIT)
%       SEGM_IM = vuActiveContoursSegmentation(IM, INIT, options)
%
%   PARAMETERS:
%       IM : Image to be segmented
%       INIT : Initial Contour
%
%   OPTIONS
%       numIter : The number of iterations through the evolution
%       lambda1, lambda2 : Variables controls edge detection
%       nu : Relating to image intensities (usually 0.001*(max(im))^2)
%       tau : Time step
%       sigma : For the Gaussian Kernal
%       mu : 
%       epsilon : Width of delta function
%       discrete : Use discrete Gaussian Kernal Flag
%       verbose : Verbose output flag
%
%   OPTIONS & OUTPUTS
%       numIter = 250
%       lambda1, lambda2 = 1.0 and 1.0
%       nu : 0.001*255*255
%       tau = 0.1
%       sigma = 3.0
%       mu = 1.0
%       epsilon = 1.0
%       discrete = false;
%       verbose = false;
%
%   OUTPUTS:
%       SEGM_IM is the evolved contour.  Contour at 0 represents
%       segmentation
%
%
%   EXAMPLE:
%       im = 
%       init =
%       seg_im = vuActiveContoursSegmentation(im,init);
%       figure;
%       subplot(1,2,1);imshow(im,[]);title('Original Image');
%       contour(double(init),[0 0],'r');
%       subplot(1,2,2);imshow(im,[]);title('Segmented Outline');hold on;
%       contour(double(seg_im),[0 0],'r');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 2
  error('MATLAB:vuActiveContoursSegmentation:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(IM) && isstruct(INIT))
  %Check for meta image structure
  if(isfield(IM,'Data') && isfield(INIT,'Data') && ...
      isfield(IM,'Dims') && isfield(INIT,'Dims') && ...
      isfield(IM,'Spc') && isfield(INIT,'Spc') && ...
      isfield(IM,'Origin') && isfield(INIT,'Origin'))
    disp('IM and INIT are valid MetaImage structs.')
  else
      error('MATLAB:vuActiveContoursSegmentation:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  IM = vuGenerateMetaImage(single(IM),ones(1,length(size(IM))),zeros(1,length(size(IM))));
  INIT = vuGenerateMetaImage(single(INIT),ones(1,length(size(INIT))),zeros(1,length(size(INIT))));
  disp('IM and INIT are matrices.')
  isStruct = false;
end

if ((length(IM.Dims) < 2 || length(IM.Dims) > 3) || (length(INIT.Dims) < 2 || length(INIT.Dims) > 3))
  error('MATLAB:vuActiveContoursSegmentation:UnknownDims', 'vuActiveContoursSegmentation can only handle images of 2 or 3 dimensions.');
end

if length(IM.Dims) ~= length(INIT.Dims)
  error('MATLAB:vuActiveContoursSegmentation:DimsMismatch', 'IM and INIT must have the same number of dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('numIter',250,@(x) isa(x,'double'));
p.addParamValue('lambda1',1,@(x) isa(x,'double'));
p.addParamValue('lambda2',1,@(x) isa(x,'double'));
p.addParamValue('nu',0.001*255*255,@(x) isa(x,'double'));
p.addParamValue('tau',0.1,@(x) isa(x,'double'));
p.addParamValue('sigma',3.0,@(x) isa(x,'double'));
p.addParamValue('mu',1,@(x) isa(x,'double'));
p.addParamValue('epsilon',1,@(x) isa(x,'double'));
p.addParamValue('discrete',0,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuActiveContoursSegmentation';
p.parse(varargin{:})


% Cast meta image elements (double-check)
IM.Data = single(IM.Data);
IM.Dims = double(IM.Dims);
IM.Spc = double(IM.Spc);
IM.Origin = double(IM.Origin);
INIT.Data = single(INIT.Data);
INIT.Dims = double(INIT.Dims);
INIT.Spc = double(INIT.Spc);
INIT.Origin = double(INIT.Origin);

% Row-major for ITK
IM = vuRowMajorColumnMajorSwitch(IM);
INIT = vuRowMajorColumnMajorSwitch(INIT);

if isStruct
    if length(IM.Dims) == 2
      SEGM_IM = MEX2DActiveContours(IM, INIT, p.Results.sigma, p.Results.nu, p.Results.tau, p.Results.mu, p.Results.lambda1, ...
          p.Results.lambda2, p.Results.epsilon, p.Results.numIter, p.Results.discrete, p.Results.verbose);
    else
      SEGM_IM = MEX3DActiveContours(IM, INIT, p.Results.sigma, p.Results.nu, p.Results.tau, p.Results.mu, p.Results.lambda1, ...
          p.Results.lambda2, p.Results.epsilon, p.Results.numIter, p.Results.discrete, p.Results.verbose);
    end
    SEGM_IM = vuGenerateMetaImage(SEGM_IM,IM.Spc,IM.Origin);
else
    if length(IM.Dims) == 2
      SEGM_IM = MEX2DActiveContours(IM, INIT, p.Results.sigma, p.Results.nu, p.Results.tau, p.Results.mu, p.Results.lambda1, ...
          p.Results.lambda2, p.Results.epsilon, p.Results.numIter, p.Results.discrete, p.Results.verbose);
    else
      SEGM_IM = MEX3DActiveContours(IM, INIT, p.Results.sigma, p.Results.nu, p.Results.tau, p.Results.mu, p.Results.lambda1, ...
          p.Results.lambda2, p.Results.epsilon, p.Results.numIter, p.Results.discrete, p.Results.verbose);
    end
end

% Copy Orientation if it exists
if (isfield(IM,'Orientation'))
    SEGM_IM.Orientation = IM.Orientation;
end
% Copy parmeters if they exist
if (isfield(IM,'Parms'))
    SEGM_IM.Parms = IM.Parms;
end

% Column-major for Matlab
SEGM_IM = vuRowMajorColumnMajorSwitch(SEGM_IM);