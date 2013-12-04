function RESULT = vuAxonMyelinAnalysis(X, Seg_X, varargin)
% vuAxonMyelinAnalysis(...) performs anisotropic diffusion on image X,
% using a modified curvature diffusion equation.
%
%   SYNTAX:
%       RESULT = vuAxonMyelinAnalysis(X, Seg_X)
%       RESULT = vuAxonMyelinAnalysis(X, Seg_X, options)       
%
%   OPTIONS & DEFAULTS:
%
%   OUTPUTS:
%       RESULT 
%
%   EXAMPLE:
%       im = vuOpenImage('Nerve.tif');
%       diff_im = vuCurvatureAnisotropicDiffusion(im);
%       seg_im = vuWatershedSegmentation(diff_im);
%       results = vuAxonMyelinAnalysis(diff_im,seg_im);     
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 2
  error('MATLAB:vuAxonMyelinAnalysis:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X) && isstruct(Seg_X))
  %Check for meta image structure
  if(isfield(X,'Data') && isfield(Seg_X,'Data') && ...
      isfield(X,'Dims') && isfield(Seg_X,'Dims') && ...
      isfield(X,'Spc') && isfield(Seg_X,'Spc') && ...
      isfield(X,'Origin') && isfield(Seg_X,'Origin'))
    disp('X and Seg_X are valid MetaImage structs.')
  else
      error('MATLAB:vuAxonMyelinAnalysis:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  Seg_X = vuGenerateMetaImage(single(Seg_X),ones(1,length(size(Seg_X))),zeros(1,length(size(Seg_X))));
  disp('X and Seg_X are matrices.')
  isStruct = false;
end

% Convert Segmented Image to Unique indexed image
Seg_X.Data = single(cmunique(double(Seg_X.Data)));
Seg_X.Dims = Seg_X.Dims(1:2);

if ((length(X.Dims) < 2 || length(X.Dims) > 2) || (length(Seg_X.Dims) < 2 || length(Seg_X.Dims) > 2))
  error('MATLAB:vuAxonMyelinAnalysis:UnknownDims', 'vuAxonMyelinAnalysis can only handle images of 2 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('roundnessThresh',0.6,@(x) isa(x,'double'));
p.addParamValue('myelinRatioThresh',0.8,@(x) isa(x,'double'));
p.addParamValue('axoplasmPixelThresh',100,@(x) isa(x,'double'));
p.addParamValue('myelinPixelThresh',100,@(x) isa(x,'double'));
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

% Create image holder for area
area = false(size(Seg_X.Data));

% Create myelin/non-myelin image (approx.)
BW_X = im2bw(X.Data,graythresh(X.Data));
% Eliminate small regions
BW_X = bwareaopen(BW_X,p.Results.myelinPixelThresh);


% Waitbar
if(p.Results.verbose)
    h = waitbar(0,'Please wait...');
end

numRegions = max(Seg_X.Data(:));

% Cycle through each index
for i = 1:numRegions
    
    if (p.Results.verbose)
        waitbar(i/numRegions, h, [['Region ' int2str(i)] ' of ' int2str(numRegions)]);
    end
    
    % Indices of current area
    indices = find(Seg_X.Data==i);

    % Check over pixel threshold
    if (length(indices) < p.Results.axoplasmPixelThresh)
        continue;
    % Check if not myelin
    elseif ((size(find(BW_X(indices)==0),1))/(size(indices,1)) < p.Results.myelinRatioThresh)
        
        % Create BW area image
        area(:) = 0;
        area(indices) = 1;

        % Find boundaries
        [B,L] = bwboundaries(area);
        
        % Check for closed area
        if (max(L(:)) ~= 1)
            continue;
        else
            % Blob analysis on region
            stats = regionprops(L,'Area','Centroid','MajorAxisLength','MinorAxisLength','Eccentricity','Perimeter');
            
            % Calculate our metric
            metric = 4*pi*stats(1).Area/stats(1).Perimeter^2;
            
            % Check roundness
            if (metric < p.Results.roundnessThresh)
                continue;
            else
                % Draw radial lines in 8 directions
                radialLine = 1;
                radialStats = [];
                for y = -1:1
                    for x = -1:1
                        % Skip 0,0
                        if (x==0&&y==0)
                            continue
                        end
                        % Loop until found two boundaries
                        boundaryCount = 0;
                        axoplasmPixelIndex = i;
                        myelinPixelIndex = [];
                        radialPoint = round(stats(1).Centroid);
                        while(boundaryCount < 2)
                            % Looking for myelin
                            if(boundaryCount==0)
                                % Increase radial line
                                radialPoint = radialPoint + [x y];
                                % Check if it is within the image
                                if (min(radialPoint)<=0 || min(Seg_X.Dims-radialPoint)<0)
                                    radialStats(radialLine).MyelinIndex = [];
                                    radialStats(radialLine).InnerBoundaryPoint = [];
                                    radialStats(radialLine).OuterBoundaryPoint = [];
                                    break
                                end
                                % Check for new region
                                if (Seg_X.Data(radialPoint(2),radialPoint(1))~=axoplasmPixelIndex)
                                    % Set myelin index
                                    myelinPixelIndex = Seg_X.Data(radialPoint(2),radialPoint(1));
                                    % Increase boundary count
                                    boundaryCount = boundaryCount + 1;
                                    % Store Data
                                    radialStats(radialLine).MyelinIndex = myelinPixelIndex;
                                    radialStats(radialLine).InnerBoundaryPoint = radialPoint;
                                end
                            else
                                % Increase radial line
                                radialPoint = radialPoint + [x y];
                                % Check if it is within the image
                                if (min(radialPoint)<=0 || min(Seg_X.Dims-radialPoint)<0)
                                    radialStats(radialLine).OuterBoundaryPoint = [];
                                    break
                                end
                                % Check for new region
                                if (Seg_X.Data(radialPoint(2),radialPoint(1))~=myelinPixelIndex)
                                    % Increase boundary count
                                    boundaryCount = boundaryCount + 1;
                                    % Store Data
                                    radialStats(radialLine).OuterBoundaryPoint = radialPoint;
                                end
                            end
                        end
                        radialLine = radialLine + 1;
                    end
                end
                % Add to results
                myelinIndices = zeros(1,8);
                for a = 1:8
                    if (~isempty(radialStats(a).MyelinIndex))
                        myelinIndices(a) = radialStats(a).MyelinIndex;
                    end
                end
                if(~exist('RESULT','var'))
                    RESULT(1) = stats(1);
                    RESULT(1).AxonIndex = axoplasmPixelIndex;
                    RESULT(1).MyelinArea = size(find(Seg_X.Data==median(myelinIndices)),1);
                    RESULT(1).RadialStats = radialStats;
                else
                    index = size(RESULT,2)+1;
                    RESULT(index).Area = stats(1).Area;
                    RESULT(index).Centroid = stats(1).Centroid;
                    RESULT(index).MajorAxisLength = stats(1).MajorAxisLength;
                    RESULT(index).MinorAxisLength = stats(1).MinorAxisLength;
                    RESULT(index).Eccentricity = stats(1).Eccentricity;
                    RESULT(index).Perimeter = stats(1).Perimeter;
                    RESULT(index).AxonIndex = axoplasmPixelIndex;
                    RESULT(index).MyelinArea = size(find(Seg_X.Data==median(myelinIndices)),1);
                    RESULT(index).RadialStats = radialStats;
                end
            end
                
        end
    end
end

if(p.Results.verbose)
    close(h);   
end

