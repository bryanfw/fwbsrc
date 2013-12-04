function [RES_IM,XFORM] = vuResamplePETtoCT(PET,CT)
% vuResamplePETtoCT resamples pet images to ct space using a fixed
% transformation determined though phantom studies
%   SYNTAX:
%       RES_IM = vuResamplePETtoCT(petIm,ctIm);
%       [RES_IM,XFORM] = vuResamplePETtoCT(petIm,ctIm);
%
%   OPTIONS:
%       none
%
%   OUTPUT:
%       RES_IM is the resampled pet image
%       XFORM is the transform used
%
%   EXAMPLE:
%       pet = vuOpenImage('./somepetimage.img.hdr');
%       ct = vuOpenImage(./ctimage/somectimage.log');
%       rePet = vuResamplePETtoCT(pet,ct);
%
% Copyright (c) 2008 Vanderbilt University Institute of Imaging Science

if nargin < 2
  error('MATLAB:vuResamplePETtoCT:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(PET) && isstruct(CT))
  %Check for meta image structure
  if(isfield(PET,'Data') && isfield(CT,'Data') && ...
      isfield(PET,'Dims') && isfield(CT,'Dims') && ...
      isfield(PET,'Spc') && isfield(CT,'Spc') && ...
      isfield(PET,'Origin') && isfield(CT,'Origin'))
    disp('PET and CT are valid MetaImage structs.')
  else
      error('MATLAB:vuResamplePETtoCT:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  error('MATLAB:vuResamplePETtoCT:InvalidStruct', 'The input image structure is not valid.');
end

if ((length(PET.Dims) < 3 || length(PET.Dims) > 3) || (length(CT.Dims) < 3 || length(CT.Dims) > 3))
  error('MATLAB:vuResamplePETtoCT:UnknownDims', 'vuResamplePETtoCT can only handle images of 3 dimensions.');
end

% ---Build transform---

% Default Setting
defaultBedSetting = [0 67.0452 227.032];

% Transform to default setting
tran.Matrix = [cosd(180) 0 -sind(180);0 1 0;sind(180) 0 cosd(180)];
tran.Offset = defaultBedSetting - PET.Parms.BedOffset;
tran1 = eye(4);
tran1(1:3,1:3) = tran.Matrix;
tran1(1:3,4) = tran.Offset';

% Transform to initialize
tran2 = [1.0000 0 0 0.5400;0 1.0000 0 10.8000;0 0 1.0000 -19.4400;0 0 0 1.0000];

% Transform from registration
tran3 = [0.9998 0.0016 0.0196 -0.0708;-0.0015 1.0000 -0.0064 -0.8061;-0.0196 0.0064 0.9998 2.0868;0 0 0 1.0000];

% Combine Transforms
tranT = tran1*tran2*tran3;

tranF.Matrix = tranT(1:3,1:3);
tranF.Offset = tranT(1:3,4)';

% Transform image
RES_IM = vuResampleImage(PET,tranF,'outImageInfo',CT);
XFORM = tranF;

