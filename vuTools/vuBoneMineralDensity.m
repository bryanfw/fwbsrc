function [BMD STD] = vuBoneMineralDensity(floodImLow,floodImHigh,darkImLow,darkImHigh,imageLow,imageHigh, varargin)
% vuBoneMineralDensity(...) calculates bone mineral and soft tissue
% densities of the given images
%
%   SYNTAX:
%       [BMD STD] = vuBoneMineralDensity(floodImLow,floodImHigh,darkImLow,darkImHigh,imageLow,imageHigh)
%       [BMD STD] = vuBoneMineralDensity(floodImLow,floodImHigh,darkImLow,darkImHigh,imageLow,imageHigh, options)       
%
%   OPTIONS & DEFAULTS:
%       muTissueHigh = 0.308 // mu of tissue at hi kVp 
%       muTissueLow = 0.412 // mu of tissue at lo kVp
%       muBoneHigh = 0.4996 // mu of bone at hi kVp
%       muBoneLow = 0.9929 // mu of bone at lo kVp
%
%   OUTPUTS:
%       BMD is the bone mineral density
%       STD is the soft tissue density
%
%   EXAMPLE:
%       floodImLow = vuOpenImage('Mean_Flood_Lo.mha');
%       floodImHigh = vuOpenImage('Mean_Flood_Hi.mha');
%       darkImLow = vuOpenImage('Mean_Dark_Lo.mha');
%       darkImHigh = vuOpenImage('Mean_Dark_Hi.mha');
%       imageLow = vuOpenImage('Mean_Img_Lo.mha');
%       imageHigh = vuOpenImage('Mean_Img_Hi.mha');
%       [BMD STD] = vuBoneMineralDensity(floodImLow,floodImHigh,darkImLow,darkImHigh,imageLow,imageHigh)
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 6
  error('MATLAB:vuBoneMineralDensity:NotEnoughInputs', 'Not enough input arguments.');
end

if (~isStructure(floodImLow))
    floodImLow = vuGenerateMetaImage(floodImLow);
else
    floodImLow.Data = single(floodImLow.Data);
end
if (~isStructure(floodImHigh))
    floodImHigh = vuGenerateMetaImage(floodImHigh);
else
    floodImHigh.Data = single(floodImHigh.Data);
end
if (~isStructure(darkImLow))
    darkImLow = vuGenerateMetaImage(darkImLow);
else
    darkImLow.Data = single(darkImLow.Data);
end
if (~isStructure(darkImHigh))
    darkImHigh = vuGenerateMetaImage(darkImHigh);
else
    darkImHigh.Data = single(darkImHigh.Data);
end
if (~isStructure(imageLow))
    imageLow = vuGenerateMetaImage(imageLow);
else
    imageLow.Data = single(imageLow.Data);
end
if (~isStructure(imageHigh))
    imageHigh = vuGenerateMetaImage(imageHigh);
    isStruct = 0;
else
    imageHigh.Data = single(imageHigh.Data);
    isStruct = 1;
end

% Get and parse options
p = inputParser;
p.addParamValue('muTissueLow',0.412,@(x) isa(x,'double'));
p.addParamValue('muTissueHigh',0.308,@(x) isa(x,'double'));
p.addParamValue('muBoneLow',0.9929,@(x) isa(x,'double'));
p.addParamValue('muBoneHigh',0.4996,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuBoneMineralDensity';
p.parse(varargin{:});

if isStruct
    [BMD STD] = MEX2DBoneMineralDensity(floodImLow, floodImHigh, darkImLow, darkImHigh, imageLow, imageHigh, ...
        p.Results.muTissueLow,p.Results.muTissueHigh,p.Results.muBoneLow,p.Results.muBoneHigh,p.Results.verbose);
else
    [bmd std] = MEX2DBoneMineralDensity(floodImLow, floodImHigh, darkImLow, darkImHigh, imageLow, imageHigh, ...
        p.Results.muTissueLow,p.Results.muTissueHigh,p.Results.muBoneLow,p.Results.muBoneHigh,p.Results.verbose);
    BMD = bmd.Data;
    STD = std.Data;
end

% Copy Orientation if it exists
if (isfield(imageHigh,'Orientation'))
    BMD.Orientation = imageHigh.Orientation;
    STD.Orientation = imageHigh.Orientation;
end
% Copy parmeters if they exist
if (isfield(imageHigh,'Parms'))
    BMD.Parms = imageHigh.Parms;
    STD.Parms = imageHigh.Parms;
end

return;

% Function to Check the Structure of the image
function isStruct = isStructure(X)

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
      isStruct = true;
      if (length(X.Dims) ~= 2)
          error('MATLAB:vuBoneMineralDensity:UnknownDims', 'vuBoneMineralDensity can only handle images of 2 dimensions.');
      end
  else
      error('MATLAB:vuBoneMineralDensity:InvalidStruct', 'The input image structure is not valid.');
  end
else
    if (ndims(X)~=2)
        error('MATLAB:vuBoneMineralDensity:UnknownDims', 'vuBoneMineralDensity can only handle images of 2 dimensions.');
    end
  isStruct = false;
end



return