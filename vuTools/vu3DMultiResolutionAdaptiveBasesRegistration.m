function [deformedIm, XX, YY, ZZ] = vu3DMultiResolutionAdaptiveBasesRegistration(fixedIm,movingIm,varargin)
% vu3DMultiResolutionAdaptiveBasesRegistration(...) performs a deformable
% registration between input images
%
%   SYNTAX:
%       [deformedIm, XX, YY, ZZ] = vu3DMultiResolutionAdaptiveBasesRegistration(fixedIm, movingIm)
%       [deformedIm, XX, YY, ZZ] = vu3DMultiResolutionAdaptiveBasesRegistration(fixedIm, movingIm, options)       
%
%   PARAMETERS:
%       fixedIm - The image movingIm will be registered to
%       movingIm - The image that will be deformed
%
%   OPTIONS & DEFAULTS:
%       gridResolutions = {[3 5],[5 8],[10 20]}
%       gradientStep = 0.1  
%       maxIterations = 2000
%       initialStep = 0.1 - Powell optimizer only
%       xTol = 0.0001
%       linMinXTol = 0.0001 - Powell optimizer only
%       fTol = 0.00001
%       optimizer = 'powell' - also have 'amoeba','conjugate_gradient'
%       verbose = 0
%
%   OUTPUTS:
%       deformedIm - The registered image
%       [XX, YY, ZZ] - The deformation fields in x and y and z
%
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

% Get and parse options
p = inputParser;
p.addParamValue('gridResolutions',{[3 5],[5 8],[10 20]},@(x) isa(x,'cell'));
p.addParamValue('gradientStep',0.1,@(x) isa(x,'double'));
p.addParamValue('maxIterations',2000,@(x) isa(x,'double'));
p.addParamValue('initialStep',0.1,@(x) isa(x,'double'));
p.addParamValue('xTol',0.0001,@(x) isa(x,'double'));
p.addParamValue('linMinXTol',0.0001,@(x) isa(x,'double'));
p.addParamValue('fTol',0.0001,@(x) isa(x,'double'));
p.addParamValue('optimizer','powell',@(x) isa(x,'char'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuCurvatureAnisotropicDiffusion';
p.parse(varargin{:});

% Resolutions
numResolutions = size(p.Results.gridResolutions,2);

% Create an image pyramid
fixedIm = vuMultiResolutionImagePyramid(fixedIm,'levels',numResolutions);
movingIm = vuMultiResolutionImagePyramid(movingIm,'levels',numResolutions);

% Create an initial deformation field
[XX,YY,ZZ] = meshgrid(0:size(fixedIm(1).Data,2)-1,0:size(fixedIm(1).Data,1)-1,0:size(fixedIm(1).Data,3)-1);
XX = vuGenerateMetaImage(XX);
YY = vuGenerateMetaImage(YY);
ZZ = vuGenerateMetaImage(ZZ);

% Loop through each resolution
for resolution = 1:numResolutions
    
    % Send in images
    [deformedIm,XX,YY,ZZ] = vu3DAdaptiveBasesRegistration(vuGenerateMetaImage(fixedIm(resolution).Data),vuGenerateMetaImage(movingIm(resolution).Data),XX,YY,ZZ,...
        'gridResolutions',p.Results.gridResolutions{resolution}, 'gradientStep',p.Results.gradientStep, ...
        'maxIterations',p.Results.maxIterations,'initialStep',p.Results.initialStep, ...
        'xTol',p.Results.xTol,'linMinXTol',p.Results.linMinXTol,'fTol',p.Results.fTol,'optimizer',p.Results.optimizer,'verbose',p.Results.verbose);

    % Upsample deformation fields
    if(resolution < numResolutions)
        
        % Up scale deformation grids
        scale = fixedIm(resolution+1).Dims./XX.Dims;
        outInfo.Origin = [0 0 0];
        outInfo.Spc = 1./scale;
        outInfo.Dims = fixedIm(resolution+1).Dims;
        tran.Matrix = eye(3);
        % Fancy stuff to maintain correct numbers
        XX.Data(:,end+1,:) = XX.Data(:,end,:)+1;
        XX.Data(end+1,:,:) = XX.Data(end,:,:);
        XX.Data(:,:,end+1) = XX.Data(:,:,end);
        XX = vuGenerateMetaImage(XX.Data.*scale(1));
        YY.Data(:,end+1,:) = YY.Data(:,end,:);
        YY.Data(end+1,:,:) = YY.Data(end,:,:)+1;
        YY.Data(:,:,end+1) = YY.Data(:,:,end);
        YY = vuGenerateMetaImage(YY.Data.*scale(2));
        ZZ.Data(:,end+1,:) = ZZ.Data(:,end,:);
        ZZ.Data(end+1,:,:) = ZZ.Data(end,:,:);
        ZZ.Data(:,:,end+1) = ZZ.Data(:,:,end)+1;
        ZZ = vuGenerateMetaImage(ZZ.Data.*scale(3));
        % Resample
        XX = vuResampleImage(XX,tran,'outImageInfo',outInfo);
        YY = vuResampleImage(YY,tran,'outImageInfo',outInfo);
        ZZ = vuResampleImage(ZZ,tran,'outImageInfo',outInfo);
        XX.Spc = [1 1 1];
        YY.Spc = [1 1 1];
        ZZ.Spc = [1 1 1];

    % Last Image
    else
        deformedIm = MEX3DDeformImage(vuRowMajorColumnMajorSwitch(movingIm(resolution)),vuRowMajorColumnMajorSwitch(XX), ...
            vuRowMajorColumnMajorSwitch(YY),vuRowMajorColumnMajorSwitch(ZZ));
        deformedIm = vuRowMajorColumnMajorSwitch(deformedIm);
    end
end
