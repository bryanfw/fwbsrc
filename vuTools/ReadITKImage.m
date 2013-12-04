function IM = ReadITKImage(filename, dims)
% READITKIMAGE reads in image formats supported by ITK.
%   SYNTAX:
%       IM = ReadITKImage(filename) reads the image identified by filename and
%       returns the image in metaimage format
%
%   OUTPUTS:
%       IM is the image in metaimage format
%
%   EXAMPLE:
%       im = ReadImage('cameraman.tif');
if nargin < 2
  error('MATLAB:readimage:NotEnoughInputs', 'Not enough input arguments.');
elseif nargin > 2
  warning('MATLAB:readimage:TooManyInputs', 'Too many arguments.');
end

if dims == 2
    IM = MEX2DITKImageFileReader(filename);
elseif dims == 3
    IM = MEX3DITKImageFileReader(filename);
else
    error('MATLAB:readimage:Dimensions','Dimensions mus be either 2 or 3D.');
end
