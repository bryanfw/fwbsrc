function [CORR_IM, MASK] = vuDistortionCorrection(X, FieldMap, varargin)
% vuDistortionCorrection performs a distortion correction estimation
% using a field map.  FieldMap should contain valid field map data only, by
% running original fieldmap data through vuBrainExtractor.
%
%   SYNTAX:
%       [CORR_IM, MASK] = vuDistortionCorrection(X, FieldMap)
%       [CORR_IM, MASK] = vuDistortionCorrection(X, FieldMap, options)
%       Correction of distortion field in image X, discribed in FieldMap
%
%   OPTIONS & DEFAULTS:
%       FieldStrength = 3
%       Jacobian = 0 (false)
%       PEDirection = 'Y'
%       WaterFatShift = X.Parms.water_fat_shift
%       verbose = 0
%
%   OUTPUTS:
%       CORR_IM is the corrected image
%
%   EXAMPLE:
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 1
  error('MATLAB:vuDistortionCorrection:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage structs.')
  else
      error('MATLAB:vuDistortionCorrection:InvalidStruct', 'The image structure is invalid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  disp('X is a matrix')
  isStruct = false;
end

if(isstruct(FieldMap))
  %Check for meta image structure
  if(isfield(FieldMap,'Data') &&  isfield(FieldMap,'Dims') && ...
          isfield(FieldMap,'Spc') && isfield(FieldMap,'Origin'))
    disp('FieldMap is a valid MetaImage structs.')
  else
      error('MATLAB:vuDistortionCorrection:InvalidStruct', 'The image structure is invalid.');
  end
  isStruct = true;
else
  tmFieldMap = FieldMap;
  clear FieldMap;
  FieldMap.Data = tmFieldMap;
  FieldMap.Dims = X.Dims;
  FieldMap.Origin = X.Origin;
  FieldMap.Spc = X.Spc;
end

% Get and parse options
p = inputParser;
p.addParamValue('FieldStrength',3,@(x) isa(x,'double'));
p.addParamValue('Jacobian',0,@(x) isa(x,'double'));
p.addParamValue('PEDirection','Y',@ischar);
p.addParamValue('WaterFatShift',0,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuDistortionCorrection';
p.parse(varargin{:});

Wavelength = 440;

% Change Default Water fat shift if it exists 
if(~isempty(p.UsingDefaults))
    if (~isempty(regexpi(cell2mat(p.UsingDefaults),'WaterFatShift')))
        try
            WaterFatShift = X.Parms.water_fat_shift;
        catch
            error('MATLAB:vuDistortionCorrection:WaterFatShift','Must set Water-Fat-Shift under X.Parms.water_fat_shift');
        end
    else
        WaterFatShift = p.Results.WaterFatShift;
    end
end

if (length(X.Dims) < 2 || length(X.Dims) > 4)
  error('MATLAB:vuDistortionCorrection:UnknownDims', 'vuDistortionCorrectionCorrectioon can only handle images of 2, 3, or 4 dimensions.');
end

if (strcmpi(p.Results.PEDirection,'x'))
    PEDir = 1;
elseif (strcmpi(p.Results.PEDirection,'y'))
    PEDir = 2;
else
    error('MATLAB:vuDistortionCorrection:InvalidPEDirection','Phase Encode Direction must be x or y');
end

% Calculate FOV
FOV = X.Dims(1:length(X.Spc)).*X.Spc;

% Calculate BWValue
BWValue = Wavelength/3*p.Results.FieldStrength/WaterFatShift*X.Dims(PEDir);

% Generate mask from fieldmap
mask = zeros(size(FieldMap.Data(:,:,:,1)));
mask(FieldMap.Data(:,:,:,1)>0) = 1;

% Prep for FieldMapCorrection_fast algorithm
if (PEDir == 1)
    X.Data = permute(X.Data,[2 1 3 4]);
    FieldMap.Data = permute(FieldMap.Data,[2 1 3 4]);
    mask = permute(mask,[2 1 3]);
    FOV = permute(FOV,[2 1 3]);
end

% Mask field map
FieldMap.Data(:,:,:,2) = FieldMap.Data(:,:,:,2).*mask;

% Setup loop variables
if (length(X.Dims)<4)
    dynamics = 1;
else
    dynamics = X.Dims(4);
end
if (length(X.Dims)<3)
    slices = 1;
else
    slices = X.Dims(3);
end

% Waitbar
if(p.Results.verbose)
    h = waitbar(0,'Please wait...');
end

corr_im = zeros(X.Dims,'single');

for a = 1:slices
    if (p.Results.verbose)
        waitbar(a/slices, h, [['slice ' int2str(a)] ' complete']);
    end
    for b = 1:dynamics
        corr_im(:,:,a,b) = FieldmapCorrection_fast(X.Data(:,:,a,b),FieldMap.Data(:,:,a,2), BWValue, FOV, p.Results.Jacobian, 1, mask(:,:,a));
    end
end

if(p.Results.verbose)
    close(h);   
end

% Some of the following code is to help with memory problems

if(isStruct)
    % Store all the information (other than Data)
    X.Data = [];
    CORR_IM = X;
    MASK = X;
end

% Free some memory
clear X

% Flip back corr_im
if (PEDir == 1)
    corr_im = permute(corr_im,[2 1 3 4]);
end

if(isStruct)
    CORR_IM.Data = corr_im;
    MASK.Data = mask;
else
    CORR_IM = corr_im;
    MASK = mask;
end

return;

% Modified version of ...
% correction GE with field map
% written by Ning Xu on Jan 25, 2005
% column is phase encoding direction, DeltaB is Hz
% Given field map
% Modified by: Kevin Wilson, on June 25, 2007
function mGE = FieldmapCorrection_fast(GE, DeltaB, Bwp, FOV, bJac, BlipPolarity, twoDspgrmask)


vSize = size(GE);
mGE = zeros(size(GE));
GRADIENTy = Bwp/FOV(2);%Hz/mm
dx = FOV(1)/vSize(1);
dy = FOV(2)/vSize(2);
if length(vSize)==2
    vSize(3)=1;
    [dDeltaBdx, dDeltaBdy] = CompGrad(DeltaB, [dx, dy]);
else
    dz = FOV(3)/vSize(3);
    [dDeltaBdx, dDeltaBdy, dDeltaBdz] = CompGrad(DeltaB, [dx, dy, dz]);
end

if nargin>=7
    DeltaB = DeltaB.*twoDspgrmask;
end
Jacobian = 1 + BlipPolarity*dDeltaBdy/GRADIENTy;
Jacobian(Jacobian<0) = 0.1;
% change DeltaB from Hz to voxel shifts
DeltaB = DeltaB/GRADIENTy/dy;

[X,Y] = meshgrid(1:vSize(1),1:vSize(2));
XI = X;
for k=1:vSize(3)
    YI = Y-DeltaB(:,:,k);
    mGE(:,:,k) = interp2(GE(:,:,k),XI,YI);
    if bJac
        mGE(:,:,k) = mGE(:,:,k).*Jacobian(:,:,k);
    end
end
mGE(isnan(mGE))=0;

if nargin>=7
    mGE1 = zeros(size(GE));
    mGE1(find(twoDspgrmask>0)) = mGE(find(twoDspgrmask>0));
    mGE = mGE1;
end

function [fx,fy,fz] = CompGrad(f, steps, IsMasked)
% function [fx,fy,fz] = CompGrad(f, steps, IsMasked)
% Compute the approximate of gradient. Based on matlab function: gradient.
% Consider outside boundary of given data. mainly used for masked data.
% Input:
%   IsMasked: 1: default, deal with outside boundary specially;
%             0: same as gradient().
%
% Writen by Yong Li, 7/17/2005
% Last modified on 8/19/2005
% x is columns, y is rows
if nargin < 3, IsMasked = 1; end

dsize = size(f);
ndim = length(dsize);
if ndim ==2, dsize(3) =1;end
if nargin < 2 || isempty(steps), steps = ones(ndim,1);end
if ndim == 3
    [fx fy fz]=gradient(f,steps(1),steps(2),steps(3));
else
    [fx fy] = gradient(f,steps(1),steps(2));
end

if ~IsMasked, return; end

for k=1:dsize(3)
    % 1st dim, y direction
    for j = 1:dsize(2)
        ind = find(f(:,j,k));
        if isempty(ind), continue;end
        if length(ind)==1, fy(max(ind-1,1):min(ind+1,dsize(1)),j,k)=0; continue;end
        if ind(1)>1, fy(ind(1)-1,j,k) = 0; end
        fy(ind(1),j,k) = (f(ind(1)+1,j,k)-f(ind(1),j,k))/steps(2);
        if ind(end)<dsize(1), fy(ind(end)+1,j,k) = 0; end
        fy(ind(end),j,k) = (f(ind(end),j,k)-f(ind(end-1),j,k))/steps(2);
    end
    % 2nd dim, x direction
    for i = 1:dsize(1)
        ind = find(f(i,:,k));
        if isempty(ind), continue;end
        if length(ind)==1, fx(i,max(ind-1,1):min(ind+1,dsize(2)),k)=0; continue;end
        if ind(1)>1, fx(i,ind(1)-1,k) = 0; end
        fx(i,ind(1),k) = (f(i,ind(1)+1,k)-f(i,ind(1),k))/steps(1);
        if ind(end)<dsize(2), fx(i,ind(end)+1,k) = 0; end
        fx(i,ind(end),k) = (f(i,ind(end),k)-f(i,ind(end-1),k))/steps(1);
    end
end

if ndim == 3
    % 3rd dim, z direction
    for j = 1:dsize(2)
        for i= 1:dsize(1)
            ind = find(f(i,j,:));
            if isempty(ind), continue;end
            if length(ind)==1,fz(i,j,max(ind-1,1):min(ind+1,dsize(3)))=0; continue;end
            if ind(1)>1, fz(i,j,ind(1)-1)=0;end
            fz(i,j,ind(1))=(f(i,j,ind(1)+1)-f(i,j,ind(1)))/steps(3);
            if ind(end)<dsize(3),fz(i,j,ind(end)+1)=0;end
            fz(i,j,ind(end))=(f(i,j,ind(end))-f(i,j,ind(end)-1))/steps(3);
        end
    end
end
