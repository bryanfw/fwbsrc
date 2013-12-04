function [REG_IM, REG_XFORM, PARMS] = vuRigidRegistration(X, Y, varargin)
% vuRigidRegistration Rigidly register two images using Mutual Information.
%   SYNTAX:
%       [REG_IM, REG_XFORM, PARMS] = vuRigidRegistration(X, Y)
%       [REG_IM, REG_XFORM, PARMS] = vuRigidRegistration(X, Y, options)
%
%       Calculates a rigid registration between images X and Y that optimizes the Mutual
%       Information between the two images.  X and Y can be 2D or 3D scalar
%       images (RGB images will not work). 
%
%   OPTIONS
%       May also specified either/both XMask, and YMask for masking
%       registration regions
%       scales : scales for the parameters (matrix, center, translation)
%       being optimized
%       initialize : Flag to initialize transform based on image
%       intensity centers
%
%   OPTIONS & OUTPUTS
%       XMask = 0
%       YMask = 0
%       initialize = 1
%       scales = [1.0 0.001 0.001 0.001 0.001] (2D)
%       scales = [1.0 1.0 1.0 0.001 0.001 0.001 0.001 0.001 0.001] (3D)
%       verbose = 0
%
%   OUTPUTS:
%       REG_IM is the registered image that provides maximum Mutual
%       Information.
%
%       REG_XFORM is the affine transformation matrix that can be used to
%       transform Y into X.  REG_XFORM is in a format compatible with
%       MAKETFORM (see HELP MAKETFORM).
%
%		PARMS 
%       (2D) is the paramters of the registration.  8 element vector : Angle Z
%       (in radians), Center of Rotation [x y], Translation [x y], 
%		number of iterations (at full resolution), starting MI, ending MI
%       (3D) is the paramters of the registration.  12 element vector : Angles X,Y,Z
%       (in radians), Center of Rotation [x y z], Translation [x y z], 
%		number of iterations (at full resolution), starting MI, ending MI
%
%   EXAMPLE:
%       im = imread('cameraman.tif');
%       tran.Matrix = [cosd(10) -sind(10); sind(10) cosd(10)];
%       tran.Offset = [10 10];
%       im2 = vuResampleImage(im,tran);
%       reg_im = vuRigidRegistration(im,im2);
%       figure;
%       subplot(1,3,1);imshow(im,colormap(gray(256)));title('Original Image');
%       subplot(1,3,2);imshow(im2,colormap(gray(256)));title('Transformed Image');
%       subplot(1,3,3);imshow(reg_im,colormap(gray(256)));title('Registered Image');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 2
  error('MATLAB:vuRigidRegistration:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X) && isstruct(Y))
  %Check for meta image structure
  if(isfield(X,'Data') && isfield(Y,'Data') && ...
      isfield(X,'Dims') && isfield(Y,'Dims') && ...
      isfield(X,'Spc') && isfield(Y,'Spc') && ...
      isfield(X,'Origin') && isfield(Y,'Origin'))
    disp('X and Y are valid MetaImage structs.')
  else
      error('MATLAB:vuRigidRegistration:InvalidStruct', 'The input image structure is not valid.');
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
  error('MATLAB:vuRigidRegistration:UnknownDims', 'vuRigidRegistration can only handle images of 2 or 3 dimensions.');
end

if length(X.Dims) ~= length(Y.Dims)
  error('MATLAB:vuRigidRegistration:DimsMismatch', 'X and Y must have the same number of dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('XMask',0,@(x) true);
p.addParamValue('YMask',0,@(x) true);
p.addParamValue('initialize',1,@(x) isa(x,'double'));
p.addParamValue('scales',[1.0 0.001 0.001 0.001 0.001],@(x) true);
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuRigidRegistration';
p.parse(varargin{:})

% Change defaults scales for 3D
if (~isempty(p.UsingDefaults))
    if ((~isempty(regexpi(cell2mat(p.UsingDefaults),'scales'))) && (length(X.Dims)==2))
        scales = p.Results.scales;
    elseif ((~isempty(regexpi(cell2mat(p.UsingDefaults),'scales'))) && (length(X.Dims)==3))
        scales = [1.0 1.0 1.0 0.001 0.001 0.001 0.001 0.001 0.001];
    else
        scales = p.Results.scales;
    end
else
    scales = p.Results.scales;
end

% Check scales is correct
if ((length(X.Dims)==3) && (length(scales)~=9))
    error('MATLAB:vuRigidRegistration:Scales','The scales vector does not have the correct dimensions (9 for 3D, 5 for 2D)')
elseif ((length(X.Dims)==2) && (length(scales)~=5))
    error('MATLAB:vuRigidRegistration:Scales','The scales vector does not have the correct dimensions (9 for 3D, 5 for 2D)')
end

% Check masks
if (isempty(p.UsingDefaults)||isempty(regexpi(cell2mat(p.UsingDefaults),'XMask')))
  if(isstruct(p.Results.XMask))
    %Check for meta image structure
    if(isfield(p.Results.XMask,'Data') &&  isfield(p.Results.XMask,'Dims') && ...
            isfield(p.Results.XMask,'Spc') && isfield(p.Results.XMask,'Origin'))
      disp('XMask is a valid MetaImage struct.')
    else
         error('MATLAB:vuRigidRegistration:InvalidStruct', 'The input image mask structure is not valid.');
    end
    XMask = p.Results.XMask;
  else
    %Assume matrix inputs with origin and spacing same as X
    XMask = vuGenerateMetaImage(single(p.Results.XMask),X.Spc,X.Origin);
    disp('XMask is a matrix.')
  end

  if length(X.Dims) ~= length(XMask.Dims)
    error('MATLAB:vuRigidRegistration:DimsMismatch', 'X and XMask must have the same number of dimensions.');
  end
end
if (isempty(p.UsingDefaults)||isempty(regexpi(cell2mat(p.UsingDefaults),'YMask')))
  if(isstruct(p.Results.YMask))
    %Check for meta image structure
    if(isfield(p.Results.YMask,'Data') &&  isfield(p.Results.YMask,'Dims') && ...
            isfield(p.Results.YMask,'Spc') && isfield(p.Results.YMask,'Origin'))
      disp('YMask is a valid MetaImage struct.')
    else
         error('MATLAB:vuRigidRegistration:InvalidStruct', 'The input image mask structure is not valid.');
    end
    YMask = p.Results.YMask;
  else
    %Assume matrix inputs with origin and spacing same as Y
    YMask = vuGenerateMetaImage(single(p.Results.YMask),Y.Spc,Y.Origin);
    disp('YMask is a matrix.')
  end

  if length(Y.Dims) ~= length(YMask.Dims)
    error('MATLAB:vuRigidRegistration:DimsMismatch', 'Y and YMask must have the same number of dimensions.');
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
      [REG_IM, REG_XFORM, PARMS] = MEX2DRigidRegistration(X,Y,XMask,YMask,p.Results.initialize,scales,p.Results.verbose);
    else
      [REG_IM, REG_XFORM, PARMS] = MEX3DRigidRegistration(X,Y,XMask,YMask,p.Results.initialize,scales,p.Results.verbose);
    end
else
    if length(X.Dims) == 2
      [reg_image, REG_XFORM, PARMS] = MEX2DRigidRegistration(X,Y,XMask,YMask,p.Results.initialize,scales,p.Results.verbose);
      REG_IM = reg_image.Data;
    else
      [reg_image, REG_XFORM, PARMS] = MEX3DRigidRegistration(X,Y,XMask,YMask,p.Results.initialize,scales,p.Results.verbose);
      REG_IM = reg_image.Data;
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
REG_XFORM = vuRowMajorColumnMajorSwitch(REG_XFORM);