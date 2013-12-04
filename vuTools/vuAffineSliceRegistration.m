function [REG_IM, REG_XFORMS] = vuAffineSliceRegistration(X, Y, varargin)
% vuAffineSliceRegistration registers two images slice-by-slice, through a full affine
% transformation (Note: see vuAffineRegistration() for full 3D-3D
% registration)
%
%   SYNTAX:
%       [REG_IM, REG_XFORMS] = vuAffineSliceRegistration(X, Y)
%       [REG_IM, REG_XFORMS] = vuAffineSliceRegistration(X, Y, options)
%
%       Calculates a rotation, translation, scaling, and shearing registration 
%       between images X and Y (slice-by-slice) that optimizes the Mutual
%       Informattion between the two images.  X and Y can be 2D or 3D scalar
%       images (RGB images will not work). 
%
%   OPTIONS
%       May also specified either/both XMask, and YMask for masking
%       registration regions
%
%       initialize : Flag to initialize transform based on image
%       intensity centers
%       scales : scales for the parameters (matrix, center, translation)
%       being optimized
%       percent : Percentage of pixels used to calculate metric 
%
%   OPTIONS & OUTPUTS
%       XMask = 0
%       YMask = 0
%       initialize = 1;
%       scales = [1.0 1.0 1.0 1.0 0.01 0.01 0.01 0.01] 
%       bins = 32;
%       numIterations = 300;
%       minStep = 0.001;
%       maxStep = 0.1;
%       pixelsUsed = 15;
%       verbose = 0;
%
%   OUTPUTS:
%       REG_IM is the rescaled image.
%       REG_XFORMS is an array of the affine transformations of each slice
%
%   EXAMPLE:
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 2
  error('MATLAB:vuAffineSliceRegistration:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X) && isstruct(Y))
  %Check for meta image structure
  if(isfield(X,'Data') && isfield(Y,'Data') && ...
      isfield(X,'Dims') && isfield(Y,'Dims') && ...
      isfield(X,'Spc') && isfield(Y,'Spc') && ...
      isfield(X,'Origin') && isfield(Y,'Origin'))
    disp('X and Y are valid MetaImage structs.')
  else
      error('MATLAB:vuAffineSliceRegistration:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  Y = vuGenerateMetaImage(single(Y),ones(1,length(size(Y))),zeros(1,length(size(Y))));
  disp('X and Y are matrices.')
  isStruct = false;
end

if ((length(X.Dims) < 2 || length(X.Dims) > 3) || (length(Y.Dims) < 2 || length(Y.Dims) > 3))
  error('MATLAB:vuAffineSliceRegistration:UnknownDims', 'vuAffineSliceRegistration can only handle images of 2 or 3 dimensions.');
end

if length(X.Dims) ~= length(Y.Dims)
  error('MATLAB:vuAffineSliceRegistration:DimsMismatch', 'X and Y must have the same number of dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('XMask',0,@(x) true);
p.addParamValue('YMask',0,@(x) true);
p.addParamValue('initialize',1,@(x) isa(x,'double'));
p.addParamValue('bins',32,@(x) isa(x,'double'));
p.addParamValue('numIterations',300,@(x) isa(x,'double'));
p.addParamValue('minStep',0.001,@(x) isa(x,'double'));
p.addParamValue('maxStep',0.1,@(x) isa(x,'double'));
p.addParamValue('scales',[1.0 1.0 1.0 1.0 0.01 0.01 0.01 0.01],@(x) true);
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.addParamValue('percent',15,@(x) isa(x,'double')&&x<=100&&x>0);
p.FunctionName='vuAffineSliceRegistration';
p.parse(varargin{:})

if (length(p.Results.scales)~=8)
    error('MATLAB:vuAffineSliceRegistration:Scales','The scales vector does not have the correct dimensions (8 for 2D)')
end

% Check masks
if (isempty(p.UsingDefaults)||isempty(regexpi(cell2mat(p.UsingDefaults),'XMask')))
  if(isstruct(p.Results.XMask))
    %Check for meta image structure
    if(isfield(p.Results.XMask,'Data') &&  isfield(p.Results.XMask,'Dims') && ...
            isfield(p.Results.XMask,'Spc') && isfield(p.Results.XMask,'Origin'))
      disp('XMask is a valid MetaImage struct.')
    else
         error('MATLAB:vuAffineSliceRegistration:InvalidStruct', 'The input image mask structure is not valid.');
    end
    XMask = p.Results.XMask;
  else
    %Assume matrix inputs with origin and spacing same as X
    XMask = vuGenerateMetaImage(single(p.Results.XMask),X.Spc,X.Origin);
    disp('XMask is a matrix.')
  end
  if max((XMask.Dims ~= size(XMask.Data)))
    error('MATLAB:vuAffineSliceRegistration:DimsMismatch', 'Invalid MetaImage format: XMask Dims does not match size(Data).');
  end

  if length(X.Dims) ~= length(XMask.Dims)
    error('MATLAB:vuAffineSliceRegistration:DimsMismatch', 'X and XMask must have the same number of dimensions.');
  end
end
if (isempty(p.UsingDefaults)||isempty(regexpi(cell2mat(p.UsingDefaults),'YMask')))
  if(isstruct(p.Results.YMask))
    %Check for meta image structure
    if(isfield(p.Results.YMask,'Data') &&  isfield(p.Results.YMask,'Dims') && ...
            isfield(p.Results.YMask,'Spc') && isfield(p.Results.YMask,'Origin'))
      disp('YMask is a valid MetaImage struct.')
    else
         error('MATLAB:vuAffineSliceRegistration:InvalidStruct', 'The input image mask structure is not valid.');
    end
    YMask = p.Results.YMask;
  else
    %Assume matrix inputs with origin and spacing same as Y
    YMask = vuGenerateMetaImage(single(p.Results.YMask),Y.Spc,Y.Origin);
    disp('YMask is a matrix.')
  end
  if max((YMask.Dims ~= size(YMask.Data)))
    error('MATLAB:vuAffineSliceRegistration:DimsMismatch', 'Invalid MetaImage format: YMask Dims does not match size(Data).');
  end

  if length(Y.Dims) ~= length(YMask.Dims)
    error('MATLAB:vuAffineSliceRegistration:DimsMismatch', 'Y and YMask must have the same number of dimensions.');
  end
end

% Cast if masks exist
if (exist('XMask','var'))
    XMask.Data = single(XMask.Data);
    XMask.Dims = double(XMask.Dims);
    XMask.Spc = double(XMask.Spc);
    XMask.Origin = double(XMask.Origin);
    XMask = vuRowMajorColumnMajorSwitch(XMask);
else
    XMask = p.Results.XMask;
end
if (exist('YMask','var'))
    YMask.Data = single(YMask.Data);
    YMask.Dims = double(YMask.Dims);
    YMask.Spc = double(YMask.Spc);
    YMask.Origin = double(YMask.Origin);
    YMask = vuRowMajorColumnMajorSwitch(YMask);
else
    YMask = p.Results.YMask;
end


% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);
Y.Data = single(Y.Data);
Y.Dims = double(Y.Dims);
Y.Spc = double(Y.Spc);
Y.Origin = double(Y.Origin);

% Row-major for ITK
X = vuRowMajorColumnMajorSwitch(X);
Y = vuRowMajorColumnMajorSwitch(Y);

if isStruct
    if length(X.Dims) == 2
      [REG_IM, REG_XFORMS] = MEX2DAffineRegistration(X,Y,XMask,YMask,p.Results.initialize,p.Results.scales,p.Results.maxStep,p.Results.minStep,p.Results.numIterations,p.Results.bins,p.Results.percent/100,p.Results.verbose);
      % Column-major for Matlab
      REG_XFORMS = vuRowMajorColumnMajorSwitch(REG_XFORMS);
    else
        for i = 1:Y.Dims(3)
            XSlice = vuGenerateMetaImage(X.Data(:,:,i),X.Spc(1:2),X.Origin(1:2));
            YSlice = vuGenerateMetaImage(Y.Data(:,:,i),Y.Spc(1:2),Y.Origin(1:2));
            if (isstruct(XMask))
                XMaskSlice = vuGenerateMetaImage(XMask.Data(:,:,i),XMask.Spc(1:2),XMask.Origin(1:2));
            else
                XMaskSlice = 0;
            end
            if (isstruct(YMask))
                YMaskSlice = vuGenerateMetaImage(YMask.Data(:,:,i),YMask.Spc(1:2),YMask.Origin(1:2));
            else
                YMaskSlice = 0;
            end     
            [reg_image, reg_xform] = MEX2DAffineRegistration(XSlice,YSlice,XMaskSlice,YMaskSlice,p.Results.initialize,p.Results.scales,p.Results.maxStep,p.Results.minStep,p.Results.numIterations,p.Results.bins,p.Results.percent/100,p.Results.verbose);
            REG_IM.Data(:,:,i) = reg_image.Data;
            REG_XFORMS(i) = reg_xform;
            % Column-major for Matlab
            REG_XFORMS(i) = vuRowMajorColumnMajorSwitch(REG_XFORMS(i));
        end
        REG_IM.Dims = Y.Dims;
        REG_IM.Spc = Y.Spc;
        REG_IM.Origin = Y.Origin;
    end
else
    if length(X.Dims) == 2
      [reg_image, REG_XFORMS] = MEX2DAffineRegistration(X,Y,XMask,YMask,p.Results.initialize,p.Results.scales,p.Results.maxStep,p.Results.minStep,p.Results.numIterations,p.Results.bins,p.Results.percent/100,p.Results.verbose);
      REG_IM = reg_image.Data;
      % Column-major for Matlab
      REG_XFORMS = vuRowMajorColumnMajorSwitch(REG_XFORMS);
    else
        for i = 1:Y.Dims(3)
            XSlice = vuGenerateMetaImage(X.Data(:,:,i),X.Spc(1:2),X.Origin(1:2));
            YSlice = vuGenerateMetaImage(Y.Data(:,:,i),Y.Spc(1:2),Y.Origin(1:2));
            if (isstruct(XMask))
                XMaskSlice = vuGenerateMetaImage(XMask.Data(:,:,i),XMask.Spc(1:2),XMask.Origin(1:2));
            else
                XMaskSlice = 0;
            end
            if (isstruct(YMask))
                YMaskSlice = vuGenerateMetaImage(YMask.Data(:,:,i),YMask.Spc(1:2),YMask.Origin(1:2));
            else
                YMaskSlice = 0;
            end
            [reg_image, reg_xform] = MEX2DAffineRegistration(XSlice,YSlice,XMaskSlice,YMaskSlice,p.Results.initialize,p.Results.scales,p.Results.maxStep,p.Results.minStep,p.Results.numIterations,p.Results.bins,p.Results.percent/100,p.Results.verbose);
            REG_IM(:,:,i) = reg_image.Data;
            REG_XFORMS(i) = reg_xform;
            % Column-major for Matlab
            REG_XFORMS(i) = vuRowMajorColumnMajorSwitch(REG_XFORMS(i));
        end
    end
end

% Copy Orientation if it exists
if (isfield(X,'Orientation'))
    REG_IM.Orientation = X.Orientation;
end
% Copy parmeters if they exist
if (isfield(X,'Parms'))
    REG_IM.Parms = X.Parms;
end

% Column-major for Matlab
REG_IM = vuRowMajorColumnMajorSwitch(REG_IM);