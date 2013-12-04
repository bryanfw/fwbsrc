function WriteITKImage(filename, Image)
% WRITEITKIMAGE writes image file formats recognized by ITK
%   SYNTAX:
%      WriteITKImage(filename, Image) writes Image to a file
%
%   OUTPUTS:
%       None
%
%   EXAMPLE:
%       im = imread('cameraman.tif');
%       im1 = double(im);
%       ReadImage('testdata.gipl', im1);
if nargin < 2
  error('MATLAB:readimage:NotEnoughInputs', 'Not enough input arguments.');
elseif nargin > 2
  warning('MATLAB:readimage:TooManyInputs', 'Too many arguments.');
end

if(isstruct(Image))
  %Check for meta image structure
  if(isfield(Image,'Data') &&  isfield(Image,'Dims') && ...
          isfield(Image,'Spc') && isfield(Image,'Origin'))
    disp('Image is a valid MetaImage structs.')
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  Image = vuGenerateMetaImage(single(Image),ones(1,length(size(Image))),zeros(1,length(size(Image))));
  isStruct = false;
  disp('Image is a matrix')
end

if length(Image.Dims) == 2
    MEXImageFileWriter2D(filename, Image);
elseif length(Image.Dims) == 3
    MEXImageFileWriter3D(filename, Image);
else
    error('MATLAB:readimage:Dimensions','Dimensions mus be either 2 or 3D.');
end
