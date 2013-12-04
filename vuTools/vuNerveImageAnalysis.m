function RESULT = vuNerveImageAnalysis(X, Seg_X, varargin)
% vuNerveImageAnalysis(...) performs analysis on image X, to find nerve
% cells, and compute the myelin and axoplasm count per cell.
%
%   SYNTAX:
%       RESULT = vuNerveImageAnalysis(X, Seg_X)   
%
%   OUTPUTS:
%       RESULT is the analysis results
%       format: Label, Original_Seed (2), Inner_Center (2), Inner_Radius,
%       Inner_Pixel_Count, Inner_Pixel_Circle_Ratio, Outer_Center (2),
%       Outer_Radius, Outer_Pixel_Count, Outer_Pixel_Circle_Ratio
%
%   EXAMPLE:
%       im = vuOpenImage('Nerve.tif');
%       diff_im = vuCurvatureAnisotropicDiffusion(im);
%       seg_im = vuWatershedSegmentation(diff_im);
%       results = vuNerveImageAnalysis(diff_im, seg_im);
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science


if nargin < 2
  error('MATLAB:vuNerveImageAnalysis:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X) && isstruct(Seg_X))
  %Check for meta image structure
  if(isfield(X,'Data') && isfield(Seg_X,'Data') && ...
      isfield(X,'Dims') && isfield(Seg_X,'Dims') && ...
      isfield(X,'Spc') && isfield(Seg_X,'Spc') && ...
      isfield(X,'Origin') && isfield(Seg_X,'Origin'))
    disp('X and Seg_X are valid MetaImage structs.')
  else
      error('MATLAB:vuNerveImageAnalysis:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  Seg_X = vuGenerateMetaImage(single(Seg_X),ones(1,length(size(Seg_X))),zeros(1,length(size(Seg_X))));
  disp('X and Seg_X are matrices.')
  isStruct = false;
end

if ((length(X.Dims) < 2 || length(X.Dims) > 2) || (length(Seg_X.Dims) < 2 || length(Seg_X.Dims) > 2))
  error('MATLAB:vuNerveImageAnalysis:UnknownDims', 'vuNerveImageAnalysis can only handle images of 2 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuNerveImageAnalysis';
p.parse(varargin{:});

% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);
Seg_X.Data = single(Seg_X.Data);
Seg_X.Dims = double(Seg_X.Dims);
Seg_X.Spc = double(Seg_X.Spc);
Seg_X.Origin = double(Seg_X.Origin);

% Row-major for ITK
X = vuRowMajorColumnMajorSwitch(X);
Seg_X = vuRowMajorColumnMajorSwitch(Seg_X);

RESULT = MEX2DNerveImageAnalysis(X,Seg_X,p.Results.verbose);
