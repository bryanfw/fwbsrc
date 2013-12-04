function [RES_IM] = vuResamplePETtoMR(PET,MR)
% vuResamplePETtoMR resamples pet images to mr space using a fixed
% transformation determined though phantom studies
%   SYNTAX:
%       RES_IM = vuResamplePETtoMR(petIm,mrIm);
%
%   OPTIONS:
%       none
%
%   OUTPUT:
%       RES_IM is the resampled pet image
%
%   EXAMPLE:
%       pet = vuOpenImage('./somepetimage.img.hdr');
%       mr = vuOpenImage(./mrimage/02.fid/fid');
%       rePet = vuResamplePETtoMR(pet,mr);
%
% Copyright (c) 2008 Vanderbilt University Institute of Imaging Science

if nargin < 2
  error('MATLAB:vuResamplePETtoMR:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(PET) && isstruct(MR))
  %Check for meta image structure
  if(isfield(PET,'Data') && isfield(MR,'Data') && ...
      isfield(PET,'Dims') && isfield(MR,'Dims') && ...
      isfield(PET,'Spc') && isfield(MR,'Spc') && ...
      isfield(PET,'Origin') && isfield(MR,'Origin'))
    disp('PET and MR are valid MetaImage structs.')
  else
      error('MATLAB:vuResamplePETtoMR:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  error('MATLAB:vuResamplePETtoMR:InvalidStruct', 'The input image structure is not valid.');
end

if ((length(PET.Dims) < 3 || length(PET.Dims) > 3) || (length(MR.Dims) < 3 || length(MR.Dims) > 3))
  error('MATLAB:vuResamplePETtoMR:UnknownDims', 'vuResamplePETtoMR can only handle images of 3 dimensions.');
end

% Change parameters to cm from mm
PET.Spc = PET.Spc./10;
PET.Origin = PET.Origin./10;

% ---Build transform---

% Put pet in default position
defaultBedSetting = [0 67.0452 227.032];
tran.Matrix = eye(3);
tran.Offset = defaultBedSetting - PET.Parms.BedOffset;

% Manual adjustment
tran.Offset = tran.Offset./10 + [0 1 0];

% From registration
tran.Matrix =  [1.0000 0.0004 0.0006;-0.0004 1.0000 0.0003;-0.0006 -0.0003 1.0000];
tran.Offset = tran.Offset+[-0.0905 0.0170 0.3049];

% Resample PET image to MR
RES_IM = vuResampleImage(PET,tran,'outImageInfo',MR);

% Flip y dimension to match MR coordinate system
RES_IM.Data = flipdim(RES_IM.Data,1);


