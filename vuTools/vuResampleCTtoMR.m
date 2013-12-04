function RES_IM = vuResampleCTtoMR(CT,MR)
% vuResampleCTtoMR resamples ct images to mr space to initialize
% the image to the correct orientation
%
%   SYNTAX:
%       RES_IM = vuResampleCTtoMR(ctIm,mrIm);
%
%   OPTIONS:
%       none
%
%   OUTPUT:
%       RES_IM is the resampled ct image
%
%   EXAMPLE:
%       mr = vuOpenImage('./somemrimage/02.fid/fid');
%       ct = vuOpenImage(./ctimage/somectimage.log');
%       reCT = vuResampleCTtoMR(ct,mr);
%
% Copyright (c) 2008 Vanderbilt University Institute of Imaging Science

if nargin < 2
  error('MATLAB:vuResampleCTtoMR:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(MR) && isstruct(CT))
  %Check for meta image structure
  if(isfield(MR,'Data') && isfield(CT,'Data') && ...
      isfield(MR,'Dims') && isfield(CT,'Dims') && ...
      isfield(MR,'Spc') && isfield(CT,'Spc') && ...
      isfield(MR,'Origin') && isfield(CT,'Origin'))
    disp('MR and CT are valid MetaImage structs.')
  else
      error('MATLAB:vuResampleCTtoMR:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  error('MATLAB:vuResampleCTtoMR:InvalidStruct', 'The input image structure is not valid.');
end

if ((length(MR.Dims) < 3 || length(MR.Dims) > 3) || (length(CT.Dims) < 3 || length(CT.Dims) > 3))
  error('MATLAB:vuResampleCTtoMR:UnknownDims', 'vuResampleCTtoMR can only handle images of 3 dimensions.');
end

% Put in default position if not already
if (max(abs(vuCalculateImageCenter(CT))))
    defaultCT = CT;
    defaultCT.Origin = -CT.Dims/2.*CT.Spc + CT.Spc/2;
    tran.Matrix = eye(3);
    CT = vuResampleImage(CT,tran,'outImageInfo',defaultCT);
end

% Put spacing in cm (from mm)
CT.Spc = CT.Spc./10;
CT.Origin = CT.Origin./10;

% Flip to MR coordinate system
CT.Data = flipdim(CT.Data,1);

% Set up transform
tran.Matrix = [cosd(180) 0 -sind(180);0 1 0;sind(180) 0 cosd(180)];
tran.Translation = [0 1.7 2.2]+[0.0937 -0.1520 -0.7792];

% Resample CT
RES_IM = vuResampleImage(CT,tran,'outImageInfo',MR);

