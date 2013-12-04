function center = vuCalculateImageCenter(X)
% vuCalculateImageCenter calculates the center of the input image X, based on
% the origin (top-left corner) and the spacing
%   SYNTAX:
%       Center = vuCalculateImageCenter(X);
%
%   OUTPUTS:
%       Center : 2D or 3D vecter of the center of im
%
%   EXAMPLE:
%       im = imread('cameraman.tif');
%       meta_im = vuGenerateMetaImage(im,[1 1],[0 0]);
%       center = vuCalculateImageCenter(meta_im);
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 1
  error('MATLAB:vuCalculateImageCenter:NotEnoughInputs', 'Not enough input arguments.');
elseif nargin > 1
  warning('MATLAB:vuCalculateImageCenter:TooManyInputs', 'Too many arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  elseif(isfield(X,'Span') && isfield(X,'Origin'))
    disp('X is a valid ROI.')
    center = X.Origin + X.Span/2;
    return
  else
      error('MATLAB:vuCalculateImageCenter:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 2 || length(X.Dims) > 5)
  error('MATLAB:vuCalculateImageCenter:UnknownDims', 'vuCalculateImageCenter can only handle images of 2, 3, 4, or 5 dimensions.');
end

% Return center
center = (X.Origin-X.Spc./2) + ((X.Dims(1:length(X.Spc))./2).*X.Spc);