function Mask = vuGenerateImageMask(X, region)
% vuGenerateImageMask Generate a image mask for X based on the input region
%   SYNTAX:
%       img = vuGenerateImageMask generates an image mask for X based on the
%       region input.
%
%       The region input must have the following properties:
%           region.Span - The span of the region (in mm)
%           region.Origin - The origin of the top left corner
%           Optional:
%           region.Orientation - The orientation matrix for the region
%   
%   DEFAULTS:
%       region.Orientation = [1 0; 0 1] (2D Image)
%       region.Orientation = [1 0 0;0 1 0;0 0 1] (3D Image)
%
%   OUTPUTS:
%       Mask : the image mask
%
%   EXAMPLE:
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

    
if nargin < 2
  error('MATLAB:vuGenerateImageMask:NotEnoughInputs', 'Not enough input arguments.');
end
if nargin > 3
    warning('MATLAB:vuGenerateImageMask:TooManyInputs', 'Too many arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuGenerateImageMask:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 3)
  error('MATLAB:vuGenerateImageMask:UnknownDims', 'vuGenerateImageMask can only handle images of 2 or 3 dimensions.');
end

% Check for orienation
if(~isfield(X,'Orientation'))
    if(length(X.Dims) == 2)
        X.Orientation = [1 0; 0 1];
    else
        X.Orientation = [1 0 0; 0 1 0; 0 0 1];
    end
end

% check input region
if (~isstruct(region))
    error('The input image region must be a structure')
else
    if(isfield(region,'Span')&&isfield(region,'Origin'))
        dims = size(region.Span);
        if~(length(dims)==2 || length(dims)==3)
            error('The input image region must be a 2D or 3D image')
        end
        % Fit the regions dimensions to the input image
        if (ndims(X.Data) == 2)
            % Check for fields
            if(~isfield(region,'Origin'))
                region.Origin = [1 1];
            end
            if(~isfield(region,'Orientation'))
                region.Orientation = [1 0;0 1];
            end
            % Resize fields
            region.Span = region.Span(1:2);
            region.Origin = region.Origin(1:2);
            if~(length(region.Orientation)==2)
                warning('Orientation matrix has incorrect dimensions.  Results may not be correct');
            end
            region.Orientation = region.Orientation(1:2,1:2);
        else
            % Check for fields
            if(~isfield(region,'Origin'))
                region.Origin = [1 1 1];
            elseif(~isfield(region,'Orientation'))
                region.Orientation = [1 0 0;0 1 0;0 0 1];
            end
            % Resize if necessary (assuming 2D)
            if(length(region.Span)==2)
                region.Span(3) = 0;
            end
            if(length(region.Origin)==2)
                region.Origin(3) = 0;
            end
            if(length(region.Orientation)==2)
                region.Orientation(:,3) = [0 0];
                region.Orientation(3,:) = [0 0 1];
            end
        end
    else
        error('MATLAB:vuGenerateImageMask:InvalidRegion','The input region structure is not valid.')
    end
end

% 2D
if (ndims(X.Data)==2)
    % Create mask
    mask.Data = ones(2,2);
    mask.Dims = [2 2];
    mask.Spc = [region.Span(1) region.Span(2)];
    mask.Origin = region.Origin;
    mask.Orientation = inv(region.Orientation);
    
    % Transform to magnetic space then image space
    tran.Matrix = mask.Orientation*inv(X.Orientation);
    tran.Center = vuCalculateImageCenter(region);
    
    % Resample to reference image
    ref_mask = vuResampleImage(mask, tran, 'outImageInfo', X, 'interpolator','nearest');
    % Round to one
    ref_mask.Data = ref_mask.Data;
else
    % Create mask
    mask.Data = ones(2,2,2);
    mask.Dims = [2 2 2];
    mask.Spc = [region.Span(1) region.Span(2) region.Span(3)];
    mask.Origin = region.Origin;
    mask.Orientation = inv(region.Orientation);
    
    % Transform to magnetic space then image space
    tran.Matrix = mask.Orientation*inv(X.Orientation);
    tran.Center = vuCalculateImageCenter(region);

    % Resample to reference image
    ref_mask = vuResampleImage(mask, tran, 'outImageInfo', X, 'interpolator','nearest');
    % Round to one
    ref_mask.Data = ref_mask.Data;
    
end

if(isStruct)
    Mask.Data = ref_mask.Data;
    Mask.Dims = X.Dims;
    Mask.Spc = X.Spc;
    Mask.Origin = X.Origin;
    Mask.Orientation = X.Orientation;
else
    Mask = ref_mask.Data;
end