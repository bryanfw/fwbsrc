function [CORR_IM, MASK] = vuParametricBiasFieldCorrection(X, MEAN, STD, varargin)
% vuParametricBiasFieldCorrection performs 2D and 3D Parametric Bias Field Correction
% by estimating and correcting the bias field in the image using Legendre
% polynomials
%
%   SYNTAX:
%       [CORR_IM, MASK] = vuParametricBiasFieldCorrection(X, MEAN, STD)
%       [CORR_IM, MASK] = vuParametricBiasFieldCorrection(X, MEAN, STD, options)
%
%       Corrects image X, with average tissue intensities listed
%       in MEAN, and their corresponding standard deviations
%       listed in STD.  Automatic image MASK is generated to mask
%       background.
%
%       Optional items include:
%       MaskThres - set the threshold value to mask the background
%       Degree - set the degree of the correction polynomial
%       Growth - set the growth factor of the optimizer
%       VolMaxIter - set the maximum iterations 
%       LogField - 1=Multiplicative, 0=Additive
%       SlabIden - Identify Slabs (1=on 0=off)
%       InterSliceCorr - Correct inter-slice intensities
%       MultiResSched - defines the multiresolution schedule
%
%   OPTIONS & DEFAULTS:
%       MaskThres = 20
%       Degree = 3
%       Growth = 1.05 
%       VolMaxIter = 20000
%       LogField = 1 (true)
%       SlabIden = 0 (false)
%       InterSliceCorr = 1 (true)
%       MultiResSched = [2 2 2 1 1 1]
%       verbose = 0
%
%   OUTPUTS:
%       CORR_IM is the corrected image
%       MASK is the output mask generated
%
%   EXAMPLE:
%       im = load('Brain1.txt');
%       MEAN = [154 115];
%       STD = [21.8 11.7];
%       [CORR_IM, MASK] = vuParametricBiasFieldCorrection(im,MEAN,STD);
%       figure
%       subplot(1,2,1);imshow(im,colormap(gray(256)));title('Image Before Filtering');
%       subplot(1,2,2);imshow(CORR_IM,colormap(gray(256)));title('Image After Filtering');
%       figure
%       imshow(MASK);title('Mask Image');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 3
  error('MATLAB:vuParametricBiasFieldCorrection:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage structs.')
    if(length(X.Dims) == 2)
        temp = zeros(X.Dims(2),X.Dims(1),2);
        temp(:,:,1) = X.Data;
        X.Data = temp;
        X = vuGenerateMetaImage(single(X.Data),[X.Spc 1],[X.Origin 0]);
        is2D = true;
    else
        is2D = false;
    end
  else
      error('MATLAB:vuParametricBiasFieldCorrection:InvalidImage','This image struction given in invalid.');
  end
  isStruct = true;
else
    if(length(size(X)) == 2)
        s = size(X);
        temp = zeros(s(1),s(2),2);
        temp(:,:,1) = X;
        X = temp;
        X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
        is2D = true;
    else
        %Assume matrix inputs with (0,0) origin and unit spacing
        X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
        is2D = false;
    end
    isStruct = false;
    disp('X is a matrix')
end

if (length(X.Dims) < 3 || length(X.Dims) > 3)
  error('MATLAB:vuParametricBiasFieldCorrection:UnknownDims', 'vuParametricBiasFieldCorrection can only handle images of 2 or 3 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('MaskThres',20,@(x) isa(x,'double'));
p.addParamValue('Degree',3,@(x) isa(x,'double'));
p.addParamValue('Growth',1.05,@(x) (isa(x,'double')&&(x>1)));
p.addParamValue('VolMaxIter',20000,@(x) isa(x,'double'));
p.addParamValue('LogField',1,@(x) isa(x,'double'));
p.addParamValue('SlabIden',0,@(x) isa(x,'double'));
p.addParamValue('InterSliceCorr',1,@(x) isa(x,'double'));
p.addParamValue('MultiResSched',[2 2 2 1 1 1],@(x) mod(length(x),3)==0);
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuParametricBiasFieldCorrection';
p.parse(varargin{:});

% Generate automatic mask
mask = zeros(X.Dims(2),X.Dims(1),X.Dims(3));
for i = 1:X.Dims(2)
    for j = 1:X.Dims(1)
        for k = 1:X.Dims(3)
            if X.Data(i,j,k) > p.Results.MaskThres
               mask(i,j,k)= 1;
            end
        end
    end
end

if(is2D)
    mask(:,:,2)=0;
end

mask = vuGenerateMetaImage(single(mask),X.Spc,X.Origin);

% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);

% Row-major for ITK
X = vuRowMajorColumnMajorSwitch(X);
mask = vuRowMajorColumnMajorSwitch(mask);

if isStruct
    if is2D
        corr_image = MEX3DMRIBiasFieldCorrectionFilter(X,mask,MEAN,STD, p.Results.Degree, p.Results.Growth, p.Results.VolMaxIter, p.Results.LogField, p.Results.SlabIden, p.Results.InterSliceCorr, p.Results.MultiResSched, p.Results.verbose);
        CORR_IM.Data = corr_image.Data(:,:,1);
        CORR_IM.Dims = corr_image.Dims(1:2);
        CORR_IM.Spc = corr_image.Spc(1:2);
        CORR_IM.Origin = corr_image.Origin(1:2);
        MASK.Data = mask.Data(:,:,1);
        MASK.Dims = mask.Dims(1:2);
        MASK.Spc = mask.Spc(1:2);
        MASK.Origin = mask.Origin(1:2);
    else
        CORR_IM = MEX3DMRIBiasFieldCorrectionFilter(X,mask,MEAN,STD, p.Results.Degree, p.Results.Growth, p.Results.VolMaxIter, p.Results.LogField, p.Results.SlabIden, p.Results.InterSliceCorr, p.Results.MultiResSched, p.Results.verbose);
        MASK = mask;
    end
else
    if is2D
        corr_image = MEX3DMRIBiasFieldCorrectionFilter(X,mask,MEAN,STD, p.Results.Degree, p.Results.Growth, p.Results.VolMaxIter, p.Results.LogField, p.Results.SlabIden, p.Results.InterSliceCorr, p.Results.MultiResSched, p.Results.verbose);
        CORR_IM = corr_image.Data(:,:,1);
        MASK = mask.Data(:,:,1);
    else
        corr_image = MEX3DMRIBiasFieldCorrectionFilter(X,mask,MEAN,STD, p.Results.Degree, p.Results.Growth, p.Results.VolMaxIter, p.Results.LogField, p.Results.SlabIden, p.Results.InterSliceCorr, p.Results.MultiResSched, p.Results.verbose);
        CORR_IM = corr_image.Data;
        MASK = mask.Data;
    end
end

% Copy Orientation if it exists
if (isfield(X,'Orientation'))
    CORR_IM.Orientation = X.Orientation;
end
% Copy parmeters if they exist
if (isfield(X,'Parms'))
    CORR_IM.Parms = X.Parms;
end

% Column-major for Matlab
CORR_IM = vuRowMajorColumnMajorSwitch(CORR_IM);
MASK = vuRowMajorColumnMajorSwitch(MASK);