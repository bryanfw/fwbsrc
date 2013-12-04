function TRAN_IM = vuResampleImage(X, transformation, varargin)
% vuResampleImage(...) resamples image X according to the transformation
% specified by transformation (Matrix/Offset OR Centered).
% Note: Transformations are defined from the output image space to the
% input image space.
%
%   TRANSFORM STRUCTURE:
%       tranformation
%               |->Matrix (Rotation matrix)
%               |->Offset (Offset)
%       OR
%
%       tranformation
%               |->Matrix (Rotation matrix)
%               |->Translation (Translation)
%               |->Center (Center of Rotation, Default (0,0) or (0,0,0))
%
%   SYNTAX:
%       TRAN_IM = vuResampleImage(X, transformation)
%       TRAN_IM = vuResampleImage(X, transformation, options)
%
%   OPTIONS:
%       outImageInfo must be struct with : .Dims, .Spc, .Origin
%       interpolator = 'linear', 'nearest', 'winsinc'
%       padding = 'none', 'zero', 'symmetric' . Padding is useful when
%       maintaining exact FOV dimension.
%
%   OPTIONS & DEFAULTS:
%       padding = none
%       interpolator = 'linear'
%       outImageInfo = same as image X
%       verbose = 0
%
%   OUTPUTS:
%       TRAN_IM is the transformed image
%
%   EXAMPLE:
%       im = imread('cameraman.tif');
%       tran.Matrix = [cosd(15) -sind(15);sind(15) cosd(15)];
%       tran.Offset = [10 10];
%       tran_im = vuResampleImage(im,tran);
%       figure;
%       subplot(1,2,1);imshow(im,colormap(gray(256)));title('Image Before Resample);
%       subplot(1,2,1);imshow(tran_im,colormap(gray(256)));title('Image After Resample);
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 2
    error('MATLAB:vuResampleImage:NotEnoughInputs','Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuResampleImage:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 3)
  error('MATLAB:vuResampleImage:UnknownDims', 'vuResampleImage can only handle images of 2 or 3 dimensions.');
end

if(isstruct(transformation))
  %Check for transformation structure
  if(isfield(transformation,'Matrix'))
    if (isfield(transformation,'Offset'))
        disp('Transformation is Matrix/Offset.')
        transType = 0;
        if (isfield(transformation,'Translation') || isfield(transformation,'Center'))
            warning('MATLAB:vuResampleImage:TransformOverstated','The Transform given has both offset and translation/center information.\n  vuResampleImage using Matrix/Offset transformation.')
        end
    elseif(isfield(transformation,'Translation'))
        disp('Transformation is Centered.')
        transType = 1;
        if(~isfield(transformation,'Center'))
            transformation.Center = [0 0 0];
        end
    else
        if(isfield(transformation,'Center'))
            disp('Transformation is Centered.')
            transformation.Translation = [0 0 0];
            transType = 1;
        else
            disp('Transformation is Matrix/Offset.')
            transformation.Offset = [0 0 0];
            transType = 0;
        end
    end
  else
      error('MATLAB:vuResampleImage:InvalidStruct', 'The input transformation structure is not valid.');
  end
end

% Get and parse options
p = inputParser;
p.addParamValue('outImageInfo',struct('Dims',double(X.Dims),'Spc',double(X.Spc),'Origin',double(X.Origin)),@(x) isstruct(x));
p.addParamValue('interpolator','linear',@ischar);
p.addParamValue('padding','none',@ischar);
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuResampleImage';
p.parse(varargin{:});

if (strcmp(p.Results.interpolator,'linear'))
    interpolator = 1;
elseif (strcmp(p.Results.interpolator,'nearest'))
    interpolator = 2;
elseif (strcmp(p.Results.interpolator,'winsinc'))
    interpolator = 3;
else
    error('MATLAB:vuResampleImage:InvalidInterpolator','The interpolator value given is not valid.');
end

if (strcmp(p.Results.padding,'none'))
    % Do nothing!
elseif (strcmp(p.Results.padding,'zero'))
    % Zero pad
    if (length(X.Dims)==2)
        X.Data = padarray(X.Data,[1 1]);
    else
        X.Data = padarray(X.Data,[1 1 1]);
    end
    X.Dims = X.Dims + 2;
    X.Origin = X.Origin - X.Spc;
    
    % Resample to add half voxel to delaunay space
    if (length(X.Dims)==2)
        tran.Matrix = eye(2);
        tran.Offset = [0 0];
        % out image is +1 dim, same spacing, origin moved half voxel back,
        % from original X
        outImg = struct('Dims',double(X.Dims)-1,'Spc',double(X.Spc),'Origin',double(X.Origin)+double(X.Spc)./2);
        X = MEX2DResampleImage(X, tran, 0, outImg, interpolator,0);
    else
        tran.Matrix = eye(3);
        tran.Offset = [0 0 0];
        % out image is +1 dim, same spacing, origin moved half voxel back,
        % from original X
        outImg = struct('Dims',double(X.Dims)-1,'Spc',double(X.Spc),'Origin',double(X.Origin)+double(X.Spc)./2);
        X = MEX3DResampleImage(X, tran, 0, outImg, interpolator,0);
    end
    
elseif (strcmp(p.Results.padding,'symmetric'))
    % Symmetric pad
    if (length(X.Dims)==2)
        X.Data = padarray(X.Data,[1 1],'symmetric');
    else
        X.Data = padarray(X.Data,[1 1 1],'symmetric');
    end
    X.Dims = X.Dims + 2;
    X.Origin = X.Origin - X.Spc;
    
    % Resample to add half voxel to delaunay space
    if (length(X.Dims)==2)
        tran.Matrix = eye(2);
        tran.Offset = [0 0];
        % out image is +1 dim, same spacing, origin moved half voxel back,
        % from original X
        outImg = struct('Dims',double(X.Dims)-1,'Spc',double(X.Spc),'Origin',double(X.Origin)+double(X.Spc)./2);
        X = MEX2DResampleImage(X, tran, 0, outImg, interpolator,0);
    else
        tran.Matrix = eye(3);
        tran.Offset = [0 0 0];
        % out image is +1 dim, same spacing, origin moved half voxel back,
        % from original X
        outImg = struct('Dims',double(X.Dims)-1,'Spc',double(X.Spc),'Origin',double(X.Origin)+double(X.Spc)./2);
        X = MEX3DResampleImage(X, tran, 0, outImg, interpolator,0);
    end
    
else
    error('MATLAB:vuResampleImage:InvalidPadding','The padding value given is not valid.');
end

if(isstruct(p.Results.outImageInfo))
    if(isfield(p.Results.outImageInfo,'Dims') && isfield(p.Results.outImageInfo,'Spc') && ...
            isfield(p.Results.outImageInfo,'Origin'))
        disp('outImageInfo struct is valid.')
    else
        error('MATLAB:vuResampleImage:InvalidStruct','The output image template is not valid.');
    end
else
    error('MATLAB:vuResampleImage:InvalidStruct','The output image template is not valid.');
end

% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);

% Row-major for ITK
X = vuRowMajorColumnMajorSwitch(X);
transformation = vuRowMajorColumnMajorSwitch(transformation);

if isStruct
    if length(X.Dims) == 2
        TRAN_IM = MEX2DResampleImage(X, transformation, transType, p.Results.outImageInfo, interpolator, p.Results.verbose);
    else
        TRAN_IM = MEX3DResampleImage(X, transformation, transType, p.Results.outImageInfo, interpolator, p.Results.verbose);
    end
else
    if length(X.Dims) == 2
        tran_image = MEX2DResampleImage(X, transformation, transType, p.Results.outImageInfo, interpolator, p.Results.verbose);
        TRAN_IM = tran_image.Data;
    else
        tran_image = MEX3DResampleImage(X, transformation, transType, p.Results.outImageInfo, interpolator, p.Results.verbose);
        TRAN_IM = tran_image.Data;
    end
end

% Copy Orientation if it exists
if (isfield(p.Results.outImageInfo,'Orientation'))
    TRAN_IM.Orientation = p.Results.outImageInfo.Orientation;
end
% Copy parmeters if they exist
if (isfield(X,'Parms'))
    TRAN_IM.Parms = X.Parms;
end

% Column-major for Matlab
TRAN_IM = vuRowMajorColumnMajorSwitch(TRAN_IM);