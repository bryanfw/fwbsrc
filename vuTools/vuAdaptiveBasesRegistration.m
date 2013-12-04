function [deformedIm, XX, YY, ZZ] = vuAdaptiveBasesRegistration(fixedIm,movingIm,varargin)
% vuAdaptiveBasesRegistration(...) performs a deformable
% registration between input images
%
%   SYNTAX:
%       [deformedIm, XX, YY, ZZ] = vuAdaptiveBasesRegistration(fixedIm, movingIm)
%       [deformedIm, XX, YY, ZZ] = vuAdaptiveBasesRegistration(fixedIm, movingIm, options)       
%
%   PARAMETERS:
%       fixedIm - The image movingIm will be registered to
%       movingIm - The image that will be deformed
%
%   OPTIONS:
%       gridResolutions : Cell array.  Each cell is a different image
%       resolution contains that resolution's radial bases grids.  
%       i.e. {[3 5],[8 12]} start at 1/2 the input image, with a 3x3, then
%       5x5 RBF grid, than at full resolution, optimizes a 8x8 and 12x12
%       gradientStep : the step size when determining the gradient map
%       gradientThreshValue : the gradient threshold cut-off (value)
%       gradientThreshPercent : the gradient threshold cut-off (% of RBFs) (0-1)
%       rbfRadiusPercent : the percentage of the rbf function to use (0-1)
%       maxIterations : the maximum iteration of the optimizer
%       initialStep : the first step size guess by Powell optimizer
%       xTol : Tolerance of the optimizer coefficients
%       fTol : Tolerance of the optimizer metric
%       optimizer : Optimizer to use
%       mask : Mask determining where RBFs optimize
%       verbose : Verbose output
%
%   OPTIONS & DEFAULTS:
%       gridResolutions = {[3 5],[5 8],[10]}
%       gradientStep = 0.1  
%       gradientThreshValue = [0 0 0]
%       gradientThreshPercent = [1 1 1]
%       rbfRadiusPercent = [1 1 1]
%       maxIterations = 10
%       initialStep = 1.0 - Powell optimizer only
%       xTol = 0.01
%       linMinXTol = 0.001 - Powell optimizer only
%       fTol = 0.001
%       optimizer = 'powell' - also have 'amoeba','conjugate_gradient'
%       mask = none
%       verbose = 0
%
%   OUTPUTS:
%       deformedIm - The registered image
%       [XX, YY, ZZ] - The deformation fields in x and y and z
%
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 2
  error('MATLAB:vuAdaptiveBasesRegistration:NotEnoughInputs', 'Not enough input arguments.');
end

% Get and parse options
p = inputParser;
p.addParamValue('gridResolutions',{[3 5],[5 8],10},@(x) isa(x,'cell'));
p.addParamValue('gradientStep',0.5,@(x) isa(x,'double'));
p.addParamValue('gradientThreshPercent',[1 1 1],@(x) isa(x,'double'));
p.addParamValue('gradientThreshValue',[0 0 0],@(x) isa(x,'double'));
p.addParamValue('rbfRadiusPercent',[1 1 1],@(x) isa(x,'double'));
p.addParamValue('maxIterations',10,@(x) isa(x,'double'));
p.addParamValue('initialStep',1,@(x) isa(x,'double'));
p.addParamValue('xTol',0.01,@(x) isa(x,'double'));
p.addParamValue('linMinXTol',0.001,@(x) isa(x,'double'));
p.addParamValue('fTol',0.001,@(x) isa(x,'double'));
p.addParamValue('minmax',100,@(x) isa(x,'double'));
p.addParamValue('optimizer','powell',@(x) isa(x,'char'));
p.addParamValue('mask',0,@(x) true);
p.addParamValue('maskPercent',[1 1 1],@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuCurvatureAnisotropicDiffusion';
p.parse(varargin{:});

if(isstruct(fixedIm) && isstruct(movingIm))
  %Check for meta image structure
  if(isfield(fixedIm,'Data') && isfield(movingIm,'Data') && ...
      isfield(fixedIm,'Dims') && isfield(movingIm,'Dims') && ...
      isfield(fixedIm,'Spc') && isfield(movingIm,'Spc') && ...
      isfield(fixedIm,'Origin') && isfield(movingIm,'Origin'))
    disp('fixedIm and movingIm are valid MetaImage structs.')
  else
      error('MATLAB:vuAdaptiveBasesRegistration:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  fixedIm = vuGenerateMetaImage(single(fixedIm),ones(1,length(size(fixedIm))),zeros(1,length(size(fixedIm))));
  movingIm = vuGenerateMetaImage(single(movingIm),ones(1,length(size(movingIm))),zeros(1,length(size(movingIm))));
  disp('fixedIm and movingIm are matrices.')
  isStruct = false;
end

if ((length(fixedIm.Dims) < 2 || length(fixedIm.Dims) > 3) || (length(movingIm.Dims) < 2 || length(movingIm.Dims) > 3))
  error('MATLAB:vuAdaptiveBasesRegistration:UnknownDims', 'vuAdaptiveBasesRegistration can only handle images of 2 or 3 dimensions.');
end

if length(fixedIm.Dims) ~= length(movingIm.Dims)
  error('MATLAB:vuAdaptiveBasesRegistration:DimsMismatch', 'fixedIm and movingIm must have the same number of dimensions.');
end

% Check mask
if(isstruct(p.Results.mask))
  %Check for meta image structure
  if(isfield(p.Results.mask,'Data') &&  isfield(p.Results.mask,'Dims') && ...
          isfield(p.Results.mask,'Spc') && isfield(p.Results.mask,'Origin'))
  else
      error('MATLAB:vuAdaptiveBasesRegistration:InvalidStruct', 'The input mask structure is not valid.');
  end
  isMaskOK = true;
else
  if (size(p.Results.mask,1)~=1)
      if(ndims(p.Results.mask)==ndims(fixedIm.Data))&&(min(size(p.Results.mask)==size(fixedIm.Data))==1)
          mask = vuGenerateMetaImage(p.Results.mask);
      else
          error('MATLAB:vuAdaptiveBasesRegistration:InvalidMask', 'The input mask size does not match input images.');
      end
      isMaskOK = true;
  else
      mask = zeros([1 size(p.Results.gridResolutions,2)]);
      isMaskOK = false;
  end
end
% Copy original movingIm Spc and Origin
movingImOrgSpc = movingIm.Spc;
movingImOrgOrigin = movingIm.Origin;

% Cast meta image elements (double-check)
fixedIm.Data = single(fixedIm.Data);
fixedIm.Dims = double(fixedIm.Dims);
movingIm.Data = single(movingIm.Data);
movingIm.Dims = double(movingIm.Dims);

% Resolutions
numResolutions = size(p.Results.gridResolutions,2);

% Check options
gradientThreshPercent = p.Results.gradientThreshPercent;
rbfRadiusPercent = p.Results.rbfRadiusPercent;
if (max(gradientThreshPercent)>1||min(gradientThreshPercent<0))
    error('MATLAB:vuAdaptiveBasesRegistration:Parm','gradientThreshPercent must be a fraction between 0 and 1')
end
if (max(rbfRadiusPercent)>1||min(rbfRadiusPercent<0))
error('MATLAB:vuAdaptiveBasesRegistration:Parm','rbfRadiusPercent must be a fraction between 0 and 1')
end
if (numel(rbfRadiusPercent)~=numResolutions)
    if(~isempty(p.UsingDefaults))
        if (~isempty(regexpi(cell2mat(p.UsingDefaults),'rbfRadiusPercent')))
            rbfRadiusPercent = ones(1,numResolutions)*0.8;
        else
            error('MATLAB:vuAdaptiveBasesRegistration:Parm','The length of rbfRadiusPercent must equal the number of resolutions')
        end
    end
end
if (numel(gradientThreshPercent)~=numResolutions)
    if(~isempty(p.UsingDefaults))
        if (~isempty(regexpi(cell2mat(p.UsingDefaults),'gradientThreshPercent')))
            gradientThreshPercent = ones(1,numResolutions);
        else
            error('MATLAB:vuAdaptiveBasesRegistration:Parm','The length of gradientThreshPercent must equal the number of resolutions')
        end
    end
end
gradientThreshValue = p.Results.gradientThreshValue;
if (numel(gradientThreshValue)~=numResolutions)
    if(~isempty(p.UsingDefaults))
        if (~isempty(regexpi(cell2mat(p.UsingDefaults),'gradientThreshValue')))
            gradientThreshValue = zeros(1,numResolutions);
        else
            error('MATLAB:vuAdaptiveBasesRegistration:Parm','The length of gradientThreshValue must equal the number of resolutions')
        end
    end
end
maskPercent = p.Results.maskPercent;
if (max(maskPercent)>1||min(maskPercent<0))
    error('MATLAB:vuAdaptiveBasesRegistration:Parm','maskPercent must be a fraction between 0 and 1')
end
if (numel(maskPercent)~=numResolutions)
    if(~isempty(p.UsingDefaults))
        if (~isempty(regexpi(cell2mat(p.UsingDefaults),'maskPercent')))
            maskPercent = ones(1,numResolutions);
        else
            error('MATLAB:vuAdaptiveBasesRegistration:Parm','The length of maskPercent must equal the number of resolutions')
        end
    end
end

% 2D
if length(fixedIm.Dims) == 2
    
    % Set Spacing and Origin to unit values
    fixedIm.Spc = [1 1];
    fixedIm.Origin = [0 0];
    movingIm.Spc = [1 1];
    movingIm.Origin = [0 0];

    % Create an image pyramid
    fixedIm = vuMultiResolutionImagePyramid(fixedIm,'levels',numResolutions);
    movingIm = vuMultiResolutionImagePyramid(movingIm,'levels',numResolutions);
    if (isMaskOK)
        if (isstruct(p.Results.mask))
            mask = vuMultiResolutionImagePyramid(p.Results.mask,'levels',numResolutions);
        else
            mask = vuMultiResolutionImagePyramid(mask,'levels',numResolutions);
        end
    end
    % Create an initial deformation field
    [XX,YY] = meshgrid(0:size(fixedIm(1).Data,2)-1,0:size(fixedIm(1).Data,1)-1);
    XX = vuGenerateMetaImage(XX);
    YY = vuGenerateMetaImage(YY);

    % Loop through each resolution
    for resolution = 1:numResolutions

        % Send in images
        [deformedIm,XX,YY] = vu2DAdaptiveBasesRegistrationWithRadius(vuGenerateMetaImage(fixedIm(resolution).Data),vuGenerateMetaImage(movingIm(resolution).Data),XX,YY,...
            'gridResolutions',p.Results.gridResolutions{resolution}, 'gradientStep',p.Results.gradientStep, ...
            'gradientThreshPercent',gradientThreshPercent(resolution), 'gradientThreshValue',gradientThreshValue(resolution), ...
            'rbfRadiusPercent',rbfRadiusPercent,'maxIterations',p.Results.maxIterations,'initialStep',p.Results.initialStep, ...
            'xTol',p.Results.xTol,'linMinXTol',p.Results.linMinXTol,'fTol',p.Results.fTol,'optimizer',p.Results.optimizer, ...
            'mask',mask(resolution),'maskPercent',maskPercent(resolution),'verbose',p.Results.verbose,'minmax',p.Results.minmax);

        % Upsample deformation fields
        if(resolution < numResolutions)

            % Up scale deformation grids
            scale = fixedIm(resolution+1).Dims./XX.Dims;
            outInfo.Origin = [0 0];
            outInfo.Spc = 1./scale;
            outInfo.Dims = fixedIm(resolution+1).Dims;
            tran.Matrix = eye(2);
            % Fancy stuff to maintain correct numbers
            XX.Data(:,end+1) = XX.Data(:,end)+1;
            XX.Data(end+1,:) = XX.Data(end,:);
            XX = vuGenerateMetaImage(XX.Data.*scale(1));
            YY.Data(:,end+1) = YY.Data(:,end);
            YY.Data(end+1,:) = YY.Data(end,:)+1;
            YY = vuGenerateMetaImage(YY.Data.*scale(2));
            % Resample
            XX = vuResampleImage(XX,tran,'outImageInfo',outInfo);
            YY = vuResampleImage(YY,tran,'outImageInfo',outInfo);
            XX.Spc = [1 1];
            YY.Spc = [1 1];

        % Last Image
        else
            deformedIm = MEX2DDeformImage(vuRowMajorColumnMajorSwitch(movingIm(resolution)),vuRowMajorColumnMajorSwitch(XX), ...
                vuRowMajorColumnMajorSwitch(YY));
            deformedIm = vuRowMajorColumnMajorSwitch(deformedIm);
        end
    end
    ZZ = 0;
    
% 3D
elseif length(fixedIm.Dims) == 3
    
    % Set Spacing and Origin to unit values
    fixedIm.Spc = [1 1 1];
    fixedIm.Origin = [0 0 0];
    movingIm.Spc = [1 1 1];
    movingIm.Origin = [0 0 0];

    % Create an image pyramid
    fixedIm = vuMultiResolutionImagePyramid(fixedIm,'levels',numResolutions);
    movingIm = vuMultiResolutionImagePyramid(movingIm,'levels',numResolutions);
    if (isMaskOK)
        if (isstruct(p.Results.mask))
            mask = vuMultiResolutionImagePyramid(p.Results.mask,'levels',numResolutions);
        else
            mask = vuMultiResolutionImagePyramid(mask,'levels',numResolutions);
        end
    end
    
    % Create an initial deformation field
    [XX,YY,ZZ] = meshgrid(0:size(fixedIm(1).Data,2)-1,0:size(fixedIm(1).Data,1)-1,0:size(fixedIm(1).Data,3)-1);
    XX = vuGenerateMetaImage(XX);
    YY = vuGenerateMetaImage(YY);
    ZZ = vuGenerateMetaImage(ZZ);

    % Loop through each resolution
    for resolution = 1:numResolutions

        % Send in images
        [deformedIm,XX,YY,ZZ] = vu3DAdaptiveBasesRegistrationWithRadius(vuGenerateMetaImage(fixedIm(resolution).Data),vuGenerateMetaImage(movingIm(resolution).Data),XX,YY,ZZ,...
            'gridResolutions',p.Results.gridResolutions{resolution}, 'gradientStep',p.Results.gradientStep, ...
            'gradientThreshPercent',gradientThreshPercent(resolution), 'gradientThreshValue',gradientThreshValue(resolution), ...
            'rbfRadiusPercent',rbfRadiusPercent,'maxIterations',p.Results.maxIterations,'initialStep',p.Results.initialStep, ...
            'xTol',p.Results.xTol,'linMinXTol',p.Results.linMinXTol,'fTol',p.Results.fTol,'optimizer',p.Results.optimizer, ...
            'mask',mask(resolution),'maskPercent',maskPercent(resolution),'verbose',p.Results.verbose,'minmax',p.Results.minmax);

        % Upsample deformation fields
        if(resolution < numResolutions)

            % Up scale deformation grids
            scale = fixedIm(resolution+1).Dims./XX.Dims;
            outInfo.Origin = [0 0 0];
            outInfo.Spc = 1./scale;
            outInfo.Dims = fixedIm(resolution+1).Dims;
            tran.Matrix = eye(3);
            % Fancy stuff to maintain correct numbers
            XX.Data(:,end+1,:) = XX.Data(:,end,:)+1;
            XX.Data(end+1,:,:) = XX.Data(end,:,:);
            XX.Data(:,:,end+1) = XX.Data(:,:,end);
            XX = vuGenerateMetaImage(XX.Data.*scale(1));
            YY.Data(:,end+1,:) = YY.Data(:,end,:);
            YY.Data(end+1,:,:) = YY.Data(end,:,:)+1;
            YY.Data(:,:,end+1) = YY.Data(:,:,end);
            YY = vuGenerateMetaImage(YY.Data.*scale(2));
            ZZ.Data(:,end+1,:) = ZZ.Data(:,end,:);
            ZZ.Data(end+1,:,:) = ZZ.Data(end,:,:);
            ZZ.Data(:,:,end+1) = ZZ.Data(:,:,end)+1;
            ZZ = vuGenerateMetaImage(ZZ.Data.*scale(3));
            % Resample
            XX = vuResampleImage(XX,tran,'outImageInfo',outInfo);
            YY = vuResampleImage(YY,tran,'outImageInfo',outInfo);
            ZZ = vuResampleImage(ZZ,tran,'outImageInfo',outInfo);
            XX.Spc = [1 1 1];
            YY.Spc = [1 1 1];
            ZZ.Spc = [1 1 1];

        % Last Image
        else
            deformedIm = MEX3DDeformImage(vuRowMajorColumnMajorSwitch(movingIm(resolution)),vuRowMajorColumnMajorSwitch(XX), ...
                vuRowMajorColumnMajorSwitch(YY),vuRowMajorColumnMajorSwitch(ZZ));
            deformedIm = vuRowMajorColumnMajorSwitch(deformedIm);
        end
    end
    
else
      error('MATLAB:vuAdaptiveBasesRegistration:UnknownDims', 'vuAdaptiveBasesRegistration can only handle images of 2 or 3 dimensions.');
end

if(~isStruct)
    deformedIm = deformedIm.Data;
    XX = XX.Data;
    YY = YY.Data;
    if length(fixedIm(numResolutions).Dims) == 3
        ZZ = ZZ.Data;
    end
else

    % Copy Original movingIm Spc and Origin
    deformedIm.Spc = movingImOrgSpc;
    deformedIm.Origin = movingImOrgOrigin;

    % Copy Orientation if it exists
    if (isfield(movingIm,'Orientation'))
        deformedIm.Orientation = movingIm(numResolutions).Orientation;
    end
    % Copy parmeters if they exist
    if (isfield(movingIm,'Parms'))
        deformedIm.Parms = movingIm(numResolutions).Parms;
    end
end