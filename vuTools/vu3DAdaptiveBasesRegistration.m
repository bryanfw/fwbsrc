function [deformedIm, XX, YY, ZZ] = vu3DAdaptiveBasesRegistration(fixedIm,movingIm,XX,YY,ZZ,varargin)
% vu3DAdaptiveBasesRegistration(...) performs a deformable registration
% between input images
%
%   SYNTAX:
%       [deformedIm, XX, YY, ZZ] = vu3DAdaptiveBasesRegistration(fixedIm, movingIm, XX, YY, ZZ)
%       [deformedImIm, XX, YY, ZZ] = vu3DAdaptiveBasesRegistration(fixedIm, movingIm, XX, YY, ZZ, options)
%
%   PARAMETERS:
%       fixedIm - The image movingIm will be registered to
%       movingIm - The image that will be deformed
%       XX, YY, ZZ - Initial deformation fields (usually created with meshgrid)
%
%   OPTIONS & DEFAULTS:
%       gridResolutions = [3 5 10]
%       gradientStep = 0.5  
%       gradientThreshValue = 0
%       gradientThreshPercent = 1
%       maxIterations = 2000
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
%       [XX, YY, ZZ] - The deformation fields in x and y and z
%
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

% Get size of image
sizeFixed = fixedIm.Dims;

% Get and parse options
p = inputParser;
p.addParamValue('gridResolutions',[3 5 10],@(x) isa(x,'double'));
p.addParamValue('gradientStep',0.5,@(x) isa(x,'double'));
p.addParamValue('gradientThreshPercent',1,@(x) isa(x,'double'));
p.addParamValue('gradientThreshValue',0,@(x) isa(x,'double'));
p.addParamValue('maxIterations',2000,@(x) isa(x,'double'));
p.addParamValue('initialStep',1.0,@(x) isa(x,'double'));
p.addParamValue('xTol',0.0001,@(x) isa(x,'double'));
p.addParamValue('linMinXTol',0.001,@(x) isa(x,'double'));
p.addParamValue('fTol',0.0001,@(x) isa(x,'double'));
p.addParamValue('optimizer','powell',@(x) isa(x,'char'));
p.addParamValue('mask',0,@(x) true);
p.addParamValue('maskPercent',1,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vu3DAdaptiveBasesRegistration';
p.parse(varargin{:});

if (strcmpi(p.Results.optimizer,'powell'))
    optimizer_type = 1;
elseif (strcmpi(p.Results.optimizer,'amoeba'))
    optimizer_type = 2;
elseif (strcmpi(p.Results.optimizer,'conjugate_gradient'))
    optimizer_type = 3;
else
    error('Matlab:vu3DAdaptiveBasesRegistration:Optimizer','The optimizer chosen was not reconized');
end

% Flip for C++
fixedIm = vuRowMajorColumnMajorSwitch(fixedIm);
movingIm = vuRowMajorColumnMajorSwitch(movingIm);
XX = vuRowMajorColumnMajorSwitch(XX);
YY = vuRowMajorColumnMajorSwitch(YY);
ZZ = vuRowMajorColumnMajorSwitch(ZZ);

% Create a holder for the bounding image
boundingIm = zeros(sizeFixed(2),sizeFixed(1),sizeFixed(3));

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
    rbfRadius = round(sqrt(3*gridSpacing^2));
    
    % Create a RBF
    rbfWeights = Calculate3DRBF(rbfRadius);
    rbfWeights = single(rbfWeights);
    sizeRBFWeights = size(rbfWeights);
    
    % Calculate where we are going to place the initial rbfs
    rbfXLocations = gridSpacing/2:gridSpacing:sizeFixed(1)-gridSpacing/4;
    rbfYLocations = gridSpacing/2:gridSpacing:sizeFixed(2)-gridSpacing/4;
    rbfZLocations = gridSpacing/2:gridSpacing:sizeFixed(3)-gridSpacing/4;
    gridSpacing(2) = gridSpacing;
    gridSpacing(3) = gridSpacing(1);
    if ((sizeFixed(1) >= sizeFixed(2))&&(sizeFixed(1) >= sizeFixed(3)))
        rbfGridResolution(2) = size(rbfYLocations,2);
        rbfGridResolution(3) = size(rbfZLocations,2);
    elseif (sizeFixed(2) >= sizeFixed(3))
        rbfGridResolution(2) = rbfGridResolution(1);
        rbfGridResolution(1) = size(rbfXLocations,2);
        rbfGridResolution(3) = size(rbfZLocations,2);
    else
        rbfGridResolution(3) = rbfGridResolution(1);
        rbfGridResolution(1) = size(rbfXLocations,2);
        rbfGridResolution(2) = size(rbfYLocations,2);  
    end
    % Round and set for indexing start at 0
    rbfXLocations = round(rbfXLocations);
    rbfYLocations = round(rbfYLocations);
    rbfZLocations = round(rbfZLocations);
    
    % Calculate the gradients
    gradientSquareNormGrid = zeros([numel(rbfXLocations) numel(rbfYLocations) numel(rbfZLocations)]);
    rbfCoeff = zeros(24,1);
    rbfOrigin = zeros(3,1);

    % Calculate current deformed image
    deformedIm = MEX3DDeformImage(movingIm,XX,YY,ZZ);

    gradientSquareNormGrid = MEX3DRadialBasesGradientMap( fixedIm, movingIm, deformedIm, XX, YY, ZZ, rbfWeights, ...
        rbfRadius, rbfGridResolution, rbfXLocations, rbfYLocations, rbfZLocations, p.Results.gradientStep, nBins, values, fixedEntropy, ...
        fixedLevels, movingLevels, movingPDF, lnMovingPDF, inJointPDF, outJointPDF, lnJointPDF);

    gradientSquareNormGrid = vuRowMajorColumnMajorSwitch(gradientSquareNormGrid);
    % Filter the gradients for optimization
    [rbfXXLocations,rbfYYLocations,rbfZZLocations] = meshgrid(rbfXLocations,rbfYLocations,rbfZLocations);
    rbfXXIndicesKeep = [];
    rbfYYIndicesKeep = [];
    rbfZZIndicesKeep = [];
    [rbfXXIndices,rbfYYIndices,rbfZZIndices] = meshgrid(1:rbfGridResolution(1),1:rbfGridResolution(2),1:rbfGridResolution(3));
    rbfIndices = 1:rbfGridResolution(1)*rbfGridResolution(2)*rbfGridResolution(3);
    
    % Thresholding
    [sortGrad, sortIdx] = sort(gradientSquareNormGrid(:),'descend');
    idxNum = ceil(p.Results.gradientThreshPercent*numel(sortGrad));
    threshValue = max(sortGrad(idxNum),p.Results.gradientThreshValue);
       
    while(~isempty(gradientSquareNormGrid))
        [sortGrad, sortIdx] = sort(gradientSquareNormGrid(:),'descend');
        if (sortGrad(1) < threshValue)
            break;
        end
        tempIdx = [rbfYYIndices(sortIdx(1)),rbfXXIndices(sortIdx(1)),rbfZZIndices(sortIdx(1))];
        
        %figure out the four to delete
        removeSub = [];
        if(tempIdx(1)-1 > 0), removeSub = cat(1,removeSub,[tempIdx(1)-1,tempIdx(2),tempIdx(3)]); end
        if(tempIdx(1)+1 <= rbfGridResolution(2)), removeSub = cat(1,removeSub,[tempIdx(1)+1,tempIdx(2),tempIdx(3)]); end
        if(tempIdx(2)-1 > 0), removeSub = cat(1,removeSub,[tempIdx(1),tempIdx(2)-1,tempIdx(3)]); end
        if(tempIdx(2)+1 <= rbfGridResolution(1)), removeSub = cat(1,removeSub,[tempIdx(1),tempIdx(2)+1,tempIdx(3)]); end
        if(tempIdx(3)-1 > 0), removeSub = cat(1,removeSub,[tempIdx(1),tempIdx(2),tempIdx(3)-1]); end
        if(tempIdx(3)+1 <= rbfGridResolution(3)), removeSub = cat(1,removeSub,[tempIdx(1),tempIdx(2),tempIdx(3)+1]); end

        removeIdx = sub2ind([rbfGridResolution(2),rbfGridResolution(1),rbfGridResolution(3)],removeSub(:,1),removeSub(:,2),removeSub(:,3));
        removeIdx = [removeIdx; rbfIndices(sortIdx(1))];
        [dummy,dummy,removeIdx] = intersect(removeIdx,rbfIndices);
        rbfYYLocations(removeIdx) = [];
        rbfXXLocations(removeIdx) = [];
        rbfZZLocations(removeIdx) = [];
        rbfXXIndices(removeIdx) = [];
        rbfYYIndices(removeIdx) = [];
        rbfZZIndices(removeIdx) = [];
        rbfIndices(removeIdx) = [];
        gradientSquareNormGrid(removeIdx) = [];
        
        rbfXXIndicesKeep = cat(1,rbfXXIndicesKeep,tempIdx(2));
        rbfYYIndicesKeep = cat(1,rbfYYIndicesKeep,tempIdx(1));
        rbfZZIndicesKeep = cat(1,rbfZZIndicesKeep,tempIdx(3));

    end

    % Find our areas to optimizer
    for i = 1:length(rbfXXIndicesKeep)
        
        % Set up Origin
        rbfOrigin = [rbfXLocations(rbfXXIndicesKeep(i)) rbfYLocations(rbfYYIndicesKeep(i)) rbfZLocations(rbfZZIndicesKeep(i))];

        percentage = 1;
        % Check if significant portion is in mask
        if (isstruct(p.Results.mask))
            maskBoundingBox = round([rbfOrigin(1:3)-rbfRadius rbfOrigin(1:3)+rbfRadius]);
            maskBoundingBox(maskBoundingBox<1) = 1;
            maskBoundingBox(4) = min(sizeFixed(1),maskBoundingBox(4));
            maskBoundingBox(5) = min(sizeFixed(2),maskBoundingBox(5));
            maskBoundingBox(6) = min(sizeFixed(3),maskBoundingBox(6));  
            maskBoundingIm = p.Results.mask.Data(maskBoundingBox(2):maskBoundingBox(5),maskBoundingBox(1):maskBoundingBox(4),maskBoundingBox(3):maskBoundingBox(6));
            percentage = numel(find(maskBoundingIm==1))/numel(maskBoundingIm);
        end
        if (percentage >= p.Results.maskPercent)
            % Place 4 rbfs around original
            rbfLocations = [rbfOrigin(1:3)-gridSpacing/2 rbfOrigin(1:3)+[gridSpacing(1)/2 -1*gridSpacing(2)/2 -1*gridSpacing(3)/2] ...
            rbfOrigin(1:3)-[gridSpacing(1)/2 -1*gridSpacing(2)/2 gridSpacing(3)/2] rbfOrigin(1:3)+[gridSpacing(1)/2 gridSpacing(2)/2 -1*gridSpacing(3)/2] ...
            rbfOrigin(1:3)+[-1*gridSpacing(1)/2 -1*gridSpacing(2)/2 gridSpacing(3)/2] rbfOrigin(1:3)+[gridSpacing(1)/2 -1*gridSpacing(2)/2 gridSpacing(3)/2] ...
            rbfOrigin(1:3)-[gridSpacing(1)/2 -1*gridSpacing(2)/2 -1*gridSpacing(3)/2] rbfOrigin(1:3)+gridSpacing/2];   
            rbfLocations = round(rbfLocations-rbfRadius);

            % Set up bounding box
            boundingBox = [rbfLocations(1:3)+1 rbfLocations(22:24)+(2*rbfRadius)];
            boundingBox(boundingBox<1) = 1;
            boundingBox(4) = min(sizeFixed(1),boundingBox(4));
            boundingBox(5) = min(sizeFixed(2),boundingBox(5));
            boundingBox(6) = min(sizeFixed(3),boundingBox(6));
            boundingIm(:) = 0;  
            boundingIm(boundingBox(2):boundingBox(5),boundingBox(1):boundingBox(4),boundingBox(3):boundingBox(6))=1;
            boundingIm = vuRowMajorColumnMajorSwitch(boundingIm);
            inIdx = find(boundingIm==1);
            outIdx = find(boundingIm==0);
            boundingIm = vuRowMajorColumnMajorSwitch(boundingIm);

            % Calculate inside region
            inDefIm = zeros(boundingBox(5)-boundingBox(2)+1,boundingBox(4)-boundingBox(1)+1,boundingBox(6)-boundingBox(3)+1);
            inDefIm = vuGenerateMetaImage(inDefIm);
            inDefIm = vuRowMajorColumnMajorSwitch(inDefIm);

            % Create some memory for each rbf
            inFixedIPPIm = MEXAllocateAlignedMemory(size(inIdx,1),0);
            inMovingIPPIm = MEXAllocateAlignedMemory(size(inIdx,1),0);
            MEXCopyArrayToAlignedMemory(fixedIm.Data(inIdx),inFixedIPPIm);
            if (size(outIdx,1)>0)
                outMovingIPPIm = MEXAllocateAlignedMemory(size(outIdx,1),0);
                outFixedIPPIm = MEXAllocateAlignedMemory(size(outIdx,1),0);
                MEXCopyArrayToAlignedMemory(fixedIm.Data(outIdx),outFixedIPPIm);
                MEXCopyArrayToAlignedMemory(deformedIm.Data(outIdx),outMovingIPPIm);
                MEXIPPCalculateSampleJointPDF(outFixedIPPIm,outMovingIPPIm,values,fixedLevels,movingLevels,outJointPDF,nBins);
            else
                MEXZeroAlignedMemory(outJointPDF);
            end

            % Initial Guess
            rbfCoeff = zeros(24,1);

            PowellCoeff = MEX3DRadialBasesOptimization(rbfCoeff,rbfWeights,rbfLocations, ...
                            movingIm, inDefIm, inFixedIPPIm, inMovingIPPIm, ...
                            XX,YY,ZZ,fixedEntropy,fixedLevels,movingLevels, nBins, ...
                            values,inJointPDF,outJointPDF,movingPDF,lnMovingPDF,lnJointPDF, ...
                            p.Results.maxIterations,p.Results.initialStep,p.Results.xTol, ...
                            p.Results.linMinXTol,p.Results.fTol,optimizer_type,p.Results.verbose);


            % Calculate current deformed image
            deformedIm = MEX3DDeformImage(movingIm,XX,YY,ZZ); 

            % Free up some memory
            MEXFreeAlignedMemory(inFixedIPPIm);
            MEXFreeAlignedMemory(inMovingIPPIm);
            if (size(outIdx,1)>0)
                MEXFreeAlignedMemory(outMovingIPPIm);
                MEXFreeAlignedMemory(outFixedIPPIm);
            end
        end
    end
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
ZZ = vuRowMajorColumnMajorSwitch(ZZ);

function RBF_Weights = Calculate3DRBF(Radius)
% Function to calculate a RBF with a given radius
Radius = round(Radius);
[XX,YY,ZZ] = meshgrid(1:Radius*2,1:Radius*2,1:Radius*2);
XX = XX - Radius;
YY = YY - Radius;
ZZ = ZZ - Radius;
RR = sqrt(sum(XX(:).^2 + YY(:).^2 + ZZ(:).^2,2))./Radius;
RR = reshape(RR,Radius*2,Radius*2,[]);

RBF_Weights = zeros(size(RR));
RBF_Weights(RR<=1) = (1 - RR(RR<=1)).^4.*(3*RR(RR<=1).^3 + 12*RR(RR<=1).^2 + 16*RR(RR<=1) + 4)./4;