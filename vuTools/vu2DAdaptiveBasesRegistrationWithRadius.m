function [deformedIm, XX, YY] = vu2DAdaptiveBasesRegistrationWithRadius(fixedIm,movingIm,XX,YY,varargin)
% vu2DAdaptiveBasesRegistrationWithRadius(...) performs a deformable registration
% between input images
%
%   SYNTAX:
%       [deformedIm, XX, YY] = vu2DAdaptiveBasesRegistrationWithRadius(fixedIm, movingIm, XX, YY)
%       [deformedImIm, XX, YY] = vu2DAdaptiveBasesRegistrationWithRadius(fixedIm, movingIm, XX, YY, options)
%
%   PARAMETERS:
%       fixedIm - The image movingIm will be registered to
%       movingIm - The image that will be deformed
%       XX, YY - Initial deformation fields (usually created with meshgrid)
%
%   OPTIONS & DEFAULTS:
%       gridResolutions = [3 5 10]
%       gradientStep = 0.5  
%       gradientThreshValue = Inf
%       gradientThreshPercent = 1
%				rbfRadiusPercent = 1;
%       maxIterations = 10
%       initialStep = 1.0 - Powell optimizer only
%       xTol = 0.001
%       linMinXTol = 0.0001 - Powell optimizer only
%       fTol = 0.00001
%       optimizer = 'powell' - also have 'amoeba','conjugate_gradient'
%       mask = none
%       maskPercent = 1
%       verbose = 0
%
%   OUTPUTS:
%       deformedIm - The registered image
%       [XX, YY] - The deformation fields in x and y
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

% Get size of image
sizeFixed = fixedIm.Dims;

% Get and parse options
p = inputParser;
p.addParamValue('gridResolutions',[3 5 10],@(x) isa(x,'double'));
p.addParamValue('gradientStep',0.5,@(x) isa(x,'double'));
p.addParamValue('gradientThreshPercent',1,@(x) isa(x,'double'));
p.addParamValue('rbfRadiusPercent',1,@(x) isa(x,'double'));
p.addParamValue('gradientThreshValue',0,@(x) isa(x,'double'));
p.addParamValue('maxIterations',10,@(x) isa(x,'double'));
p.addParamValue('initialStep',1.0,@(x) isa(x,'double'));
p.addParamValue('xTol',0.01,@(x) isa(x,'double'));
p.addParamValue('linMinXTol',0.001,@(x) isa(x,'double'));
p.addParamValue('fTol',0.001,@(x) isa(x,'double'));
p.addParamValue('minmax',0.1,@(x) isa(x,'double'));
p.addParamValue('optimizer','powell',@(x) isa(x,'char'));
p.addParamValue('mask',0,@(x) true);
p.addParamValue('maskPercent',1,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vu2DAdaptiveBasesRegistrationWithRadius';
p.parse(varargin{:});

if (strcmpi(p.Results.optimizer,'powell'))
    optimizer_type = 1;
elseif (strcmpi(p.Results.optimizer,'amoeba'))
    optimizer_type = 2;
elseif (strcmpi(p.Results.optimizer,'conjugate_gradient'))
    optimizer_type = 3;
else
    error('Matlab:vu2DAdaptiveBasesRegistrationWithRadius:Optimizer','The optimizer chosen was not reconized');
end

% Flip for C++
fixedIm = vuRowMajorColumnMajorSwitch(fixedIm);
movingIm = vuRowMajorColumnMajorSwitch(movingIm);
XX = vuRowMajorColumnMajorSwitch(XX);
YY = vuRowMajorColumnMajorSwitch(YY);

% Get Ready ...
deformedIm = movingIm;
fixedIPPIm = MEXAllocateAlignedMemory(numel(fixedIm.Data),0);
movingIPPIm = MEXAllocateAlignedMemory(numel(movingIm.Data),0);
MEXCopyArrayToAlignedMemory(fixedIm.Data,fixedIPPIm);
MEXCopyArrayToAlignedMemory(movingIm.Data,movingIPPIm);
nBins = 32;
values = MEXAllocateAlignedMemory(nBins+1,0);
fixedLevels = MEXAllocateAlignedMemory(nBins+1,0);
movingLevels = MEXAllocateAlignedMemory(nBins+1,0);
MEXIPPCalculateImageBins(fixedIPPIm,values,fixedLevels,nBins);
MEXIPPCalculateImageBins(movingIPPIm,values,movingLevels,nBins);
fixedEntropy = MEXIPPCalculateFixedImageEntropy(fixedIPPIm,values,fixedLevels,nBins);
inJointPDF = MEXAllocateAlignedMemory(nBins*nBins,0);
outJointPDF = MEXAllocateAlignedMemory(nBins*nBins,0);
movingPDF = MEXAllocateAlignedMemory(nBins,1);
lnJointPDF = MEXAllocateAlignedMemory(nBins*nBins,0);
lnMovingPDF = MEXAllocateAlignedMemory(nBins,1);

% Free Some Memory
MEXFreeAlignedMemory(fixedIPPIm);
MEXFreeAlignedMemory(movingIPPIm);

% Loop through each level
for level = 1:size(p.Results.gridResolutions,2)

    % Calculate radius
    rbfGridResolution = p.Results.gridResolutions(level);

    % Calculate Grid spacing and RBF radius
    gridSpacing = max(sizeFixed)/rbfGridResolution;
    rbfRadius = round(sqrt(2*gridSpacing^2));
    
    % Create a RBF
    rbfWeights = CalculateRBF(rbfRadius);
    rbfWeights = single(rbfWeights);
    sizeRBFWeights = size(rbfWeights);
    
    % Calculate where we are going to place the initial rbfs
    rbfXLocations = gridSpacing/2:gridSpacing:sizeFixed(1)-gridSpacing/4;
    rbfYLocations = gridSpacing/2:gridSpacing:sizeFixed(2)-gridSpacing/4;
	if (isempty(rbfXLocations))
        rbfXLocations = sizeFixed(1)/2;
    end
    if (isempty(rbfYLocations))
        rbfYLocations = sizeFixed(2)/2;
    end
    gridSpacing(2) = gridSpacing;
    if (sizeFixed(1) >= sizeFixed(2))
        rbfGridResolution(2) = size(rbfYLocations,2);
    else
        rbfGridResolution(2) = rbfGridResolution(1);
        rbfGridResolution(1) = size(rbfXLocations,2);
    end
    % Round and set for indexing start at 0
    rbfXLocations = round(rbfXLocations);
    rbfYLocations = round(rbfYLocations);
    
    % Calculate the gradients
    rbfCoeff = zeros(8,1);
    rbfOrigin = zeros(2,1);

    % Calculate current deformed image
    deformedIm = MEX2DDeformImage(movingIm,XX,YY);
    
    gradientSquareNormGrid = MEX2DRadialBasesGradientMapWithRadius( fixedIm, movingIm, deformedIm, XX, YY, rbfWeights, ...
        rbfRadius, rbfGridResolution, rbfXLocations, rbfYLocations, p.Results.gradientStep, nBins, values, fixedEntropy, ...
        fixedLevels, movingLevels, movingPDF, lnMovingPDF, inJointPDF, outJointPDF, lnJointPDF);
    gradientSquareNormGrid = vuRowMajorColumnMajorSwitch(gradientSquareNormGrid);
    
    % Filter the gradients for optimization
    [rbfXXLocations,rbfYYLocations] = meshgrid(rbfXLocations,rbfYLocations);
    rbfXXIndicesKeep = [];
    rbfYYIndicesKeep = [];
    [rbfXXIndices,rbfYYIndices] = meshgrid(1:rbfGridResolution(1),1:rbfGridResolution(2));
    rbfIndices = 1:rbfGridResolution(1)*rbfGridResolution(2);
    
    % Thresholding
    [sortGrad, sortIdx] = sort(gradientSquareNormGrid(:),'descend');
    idxNum = ceil(p.Results.gradientThreshPercent*numel(sortGrad));
    threshValue = max(sortGrad(idxNum),p.Results.gradientThreshValue);
    
    while(~isempty(gradientSquareNormGrid))
        [sortGrad, sortIdx] = sort(gradientSquareNormGrid(:),'descend');
        if (sortGrad(1) < threshValue)
            break;
        end
        tempIdx = [rbfYYIndices(sortIdx(1)),rbfXXIndices(sortIdx(1))];
        
        %figure out the four to delete
        removeSub = [];
        if(tempIdx(1)-1 > 0), removeSub = cat(1,removeSub,[tempIdx(1)-1,tempIdx(2)]); end
        if(tempIdx(1)+1 <= rbfGridResolution(2)), removeSub = cat(1,removeSub,[tempIdx(1)+1,tempIdx(2)]); end
        if(tempIdx(2)-1 > 0), removeSub = cat(1,removeSub,[tempIdx(1),tempIdx(2)-1]); end
        if(tempIdx(2)+1 <= rbfGridResolution(1)), removeSub = cat(1,removeSub,[tempIdx(1),tempIdx(2)+1]); end

        removeIdx = sub2ind([rbfGridResolution(2),rbfGridResolution(1)],removeSub(:,1),removeSub(:,2));
        removeIdx = [removeIdx; rbfIndices(sortIdx(1))];
        [dummy,dummy,removeIdx] = intersect(removeIdx,rbfIndices);
        rbfYYLocations(removeIdx) = [];
        rbfXXLocations(removeIdx) = [];
        rbfXXIndices(removeIdx) = [];
        rbfYYIndices(removeIdx) = [];
        rbfIndices(removeIdx) = [];
        gradientSquareNormGrid(removeIdx) = [];
        
        rbfXXIndicesKeep = cat(1,rbfXXIndicesKeep,tempIdx(2));
        rbfYYIndicesKeep = cat(1,rbfYYIndicesKeep,tempIdx(1));

    end

    % Setup our locations
    rbfXLocations = rbfXLocations(rbfXXIndicesKeep(:));
    rbfYLocations = rbfYLocations(rbfYYIndicesKeep(:));

		if (max(rbfGridResolution(:))>5) 
			minmax = p.Results.minmax;
		else
			minmax = 1000;
		end

    PowellCoeff = MEX2DRadialBasesOptimizationWithRadius(rbfXLocations,rbfYLocations, ...
                    rbfWeights,p.Results.rbfRadiusPercent, movingIm, fixedIm, deformedIm, ...
                    XX,YY,fixedEntropy,fixedLevels,movingLevels, nBins, ...
                    values,inJointPDF,outJointPDF,movingPDF,lnMovingPDF,lnJointPDF, ...
                    p.Results.maxIterations,p.Results.initialStep,p.Results.xTol, ...
                    p.Results.linMinXTol,p.Results.fTol,optimizer_type,p.Results.verbose, ...
                    vuRowMajorColumnMajorSwitch(p.Results.mask), gridSpacing,p.Results.minmax);
end

% Free our memory
MEXFreeAlignedMemory(values);
MEXFreeAlignedMemory(fixedLevels);
MEXFreeAlignedMemory(movingLevels);
MEXFreeAlignedMemory(inJointPDF);
MEXFreeAlignedMemory(outJointPDF);
MEXFreeAlignedMemory(movingPDF);
MEXFreeAlignedMemory(lnJointPDF);
MEXFreeAlignedMemory(lnMovingPDF); 


deformedIm = vuRowMajorColumnMajorSwitch(deformedIm);
XX = vuRowMajorColumnMajorSwitch(XX);
YY = vuRowMajorColumnMajorSwitch(YY);

function RBF_Weights = CalculateRBF(Radius)
% Function to calculate a RBF with a given radius
Radius = round(Radius);
XX = repmat((1:Radius*2)',[1 Radius*2])-Radius;
YY = repmat(1:Radius*2,[Radius*2 1])-Radius;

RR = sqrt(sum(XX(:).^2 + YY(:).^2,2))./Radius;
RR = reshape(RR,Radius*2,[]);

RBF_Weights = zeros(size(RR));
RBF_Weights(RR<=1) = (1 - RR(RR<=1)).^4.*(3*RR(RR<=1).^3 + 12*RR(RR<=1).^2 + 16*RR(RR<=1) + 4)./4;