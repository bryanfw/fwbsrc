function [CORR_IM, B, SEG_IM] = vuILICBiasFieldCorrection(X, varargin)
% vuILICBiasFieldCorrection performs a bias field correction on image X
% based on Intergrated Local Intensity Clustering.  Segmented regions may
% also be found through results
%
%  SYNTAX:
%       [CORR_IM, B, SEG_IM] = vuILICBiasFieldCorrection(X)
%       [CORR_IM, B, SEG_IM] = vuILICBiasFieldCorrection(X, options)
%
%       Corrects image X giving corrected image (CORR_IM), Bias Field (B),
%       and Segmented Regions (SEG_IM)
%
%   OPTIONS:
%       sigma : Size of gaussian kernal
%       numIterations : Number of iteration
%       numRegions : Number of Regions in the image (2,3, or 4)
%
%   DEFAULTS:
%       sigma = 4;
%       numIterations = 50;
%       numRegions = 3;
%
%   OUTPUTS:
%       CORR_IM : Corrected Image
%       B : Bias Field Image
%       SEG_IM : Segmented Regions
%
%   EXAMPLE:
%       im = vuOpenImage('BrainImage.PAR');
%       [corrIm, B, segIm] = vuILICBiasFieldCorrection(im);
%       vuOrthogonalOverlay(im,corrIm);
%
% Copyright (c) 2008 - Vanderbilt University Insititute of Imaging Science

if nargin < 1
    error('MATLAB:vuILICBiasFieldCorrection:NotEnoughInputs','Not enough input arguments');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuILICBiasFieldCorrection:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) < 3 || length(X.Dims) > 3)
  error('MATLAB:vuILICBiasFieldCorrection:UnknownDims', 'vuILICBiasFieldCorrection can only handle images of 3 dimensions.');
end
 
% Get and parse options
p = inputParser;
p.addParamValue('sigma',4,@(x) isa(x,'double'));
p.addParamValue('numIterations',50,@(x) isa(x,'int'));
p.addParamValue('numRegions',3,@(x) isa(x,'int'));
p.FunctionName='vuBrainExtraction';
p.parse(varargin{:});

% Check number of regions
if (p.Results.numRegions < 2)||(p.Results.numRegions > 4)
  error('MATLAB:vuILICBiasFieldCorrection:NumberOfRegions', 'vuILICBiasFieldCorrection only allows number of regions of 2, 3, or 4.');
end

% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = single(X.Dims);
X.Spc = single(X.Spc);
X.Origin = single(X.Origin);

% Row-major for ITK
X = vuRowMajorColumnMajorSwitch(X);

% Normailze data to 1-255
maxData = max(X.Data(:));
minData = min(X.Data(:));
X.Data = (X.Data-minData).*(254/(maxData-minData))+1;

% Create inputs
A=max(X.Data(:));
C(1)=0.1; C(2)=0.3; C(3)=0.5;
C=C*A;
lambda = single(0.001);
Ksigma = single(gaussMask(p.Results.sigma,3));
KONE=single(ones(size(X.Data)));
KONE = convn(KONE,Ksigma,'same'); 

% Initialize regions
if (p.Results.numRegions == 2)
    M(:,:,:,1)=single(rand(size(X.Data)));
    M(:,:,:,2) = 1 - M(:,:,:,1);
elseif (p.Results.numRegions == 3)
    M(:,:,:,1)=single(rand(size(X.Data))*0.5);
    M(:,:,:,2)=single(rand(size(X.Data))*0.5);
    M(:,:,:,3) = 1 - M(:,:,:,1) - M(:,:,:,2);
else
    M(:,:,:,1)=single(rand(size(X.Data))*0.33);
    M(:,:,:,2)=single(rand(size(X.Data))*0.33);
    M(:,:,:,3)=single(rand(size(X.Data))*0.33);
    M(:,:,:,4) = 1 - M(:,:,:,1) - M(:,:,:,2) - M(:,:,:,2);
end

B = single(ones(size(X.Data)));



[SEG_IM, B, C2] = EVOL_ILIC_3D_v1_MEX(X.Data, B, C, M, Ksigma, KONE,lambda,single(p.Results.numIterations),single(1));

% Construct outputs
B = vuRowMajorColumnMajorSwitch(B);
SEG_IM = vuRowMajorColumnMajorSwitch(SEG_IM);
X = vuRowMajorColumnMajorSwitch(X);

% Scale data back
X.Data = (X.Data-1)./(254/(maxData-minData)) + minData;

% Correct Bias field
if (isStruct)
    CORR_IM = X;
    CORR_IM.Data = CORR_IM.Data./B;
    B = vuGenerateMetaImage(B,X.Spc,X.Origin);
    SEG_IM = vuGenerateMetaImage(SEG_IM,X.Spc,X.Origin);
else
    CORR_IM = X.Data./B;
end

% Cast meta data
CORR_IM.Dims = double(CORR_IM.Dims);
CORR_IM.Spc = double(CORR_IM.Spc);
CORR_IM.Origin = double(CORR_IM.Origin);
SEG_IM.Dims = double(SEG_IM.Dims);
SEG_IM.Spc = double(SEG_IM.Spc);
SEG_IM.Origin = double(SEG_IM.Origin);
B.Dims = double(B.Dims);
B.Spc = double(B.Spc);
B.Origin = double(B.Origin);

function g = gaussMask(sig, dim)  
% Generate 2D or 3D Gaussian mask as convolution kernel
% input: 
%       p = [p1 p2 p3] specify the size of the mask
%       sig is the std of Gaussian kernel
% Copyright (c) 2004--2007 by Chunming Li
% Author: Chunming Li, 10/15/2004

r=2*ceil(2*sig)+1;
p=r*ones(1,dim);

siz   = (p-1)/2;
std   = sig;

if length(p)==3    
[x,y,z] = meshgrid(-siz(2):siz(2),-siz(1):siz(1),-siz(3):siz(3));
arg   = -(x.*x + y.*y+z.*z)/(2*std*std);
h     = exp(arg);
h(h<eps*max(h(:))) = 0;    
sumh = sum(h(:));
if sumh ~= 0,
g  = h/sumh;
end
elseif length(p)==2    
[x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
arg   = -(x.*x + y.*y)/(2*std*std);    
h     = exp(arg);
h(h<eps*max(h(:))) = 0;    
sumh = sum(h(:));
if sumh ~= 0,
g  = h/sumh;
end    
end