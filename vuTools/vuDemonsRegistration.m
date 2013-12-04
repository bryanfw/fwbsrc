function [REG_IM, XX, YY, ZZ] = vuDemonsRegistration(X, Y, varargin)
% vuDemonsRegistration registers two images though Demons' algorithm
% registration.
%
%   SYNTAX:
%       [REG_IM] = vuDemonsRegistration(X, Y)
%       [REG_IM] = vuDemonsRegistration(X, Y, options)
%
%       Calculates vector field to deform Y into X. 
%
%   OPTIONS
%
%       numIterations : The number of iterations of the solver
%       sigma : The stardard deviation of the Gaussian kernal used to
%       smooth the vector field
%       verbose : Verbose output
%
%   OPTIONS & OUTPUTS
%       numIterations = 1000;
%       sigma = 1.0;
%       verbose = 0;
%
%   OUTPUTS:
%       REG_IM is the registered image.
%
%   EXAMPLE:
%
% Copyright (c) 2008 - Vanderbilt University Institute of Imaging Science

if nargin < 2
  error('MATLAB:vuDemonsRegistration:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X) && isstruct(Y))
  %Check for meta image structure
  if(isfield(X,'Data') && isfield(Y,'Data') && ...
      isfield(X,'Dims') && isfield(Y,'Dims') && ...
      isfield(X,'Spc') && isfield(Y,'Spc') && ...
      isfield(X,'Origin') && isfield(Y,'Origin'))
    disp('X and Y are valid MetaImage structs.')
  else
      error('MATLAB:vuDemonsRegistration:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  %Assume matrix inputs with (0,0) origin and unit spacing
  X = vuGenerateMetaImage(single(X),ones(1,length(size(X))),zeros(1,length(size(X))));
  Y = vuGenerateMetaImage(single(Y),ones(1,length(size(Y))),zeros(1,length(size(Y))));
  disp('X and Y are matrices.')
  isStruct = false;
end

if ((length(X.Dims) < 2 || length(X.Dims) > 3) || (length(Y.Dims) < 2 || length(Y.Dims) > 3))
  error('MATLAB:vuDemonsRegistration:UnknownDims', 'vuDemonsRegistration can only handle images of 2 or 3 dimensions.');
end

if length(X.Dims) ~= length(Y.Dims)
  error('MATLAB:vuDemonsRegistration:DimsMismatch', 'X and Y must have the same number of dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('numIterations',1000,@(x) isa(x,'double'));
p.addParamValue('sigma',1.0,@(x) isa(x,'double'));
p.addParamValue('verbose',0,@(x) isa(x,'double'));
p.FunctionName='vuDemonsRegistration';
p.parse(varargin{:})

% Cast meta image elements (double-check)
X.Data = single(X.Data);
X.Dims = double(X.Dims);
X.Spc = double(X.Spc);
X.Origin = double(X.Origin);
Y.Data = single(Y.Data);
Y.Dims = double(Y.Dims);
Y.Spc = double(Y.Spc);
Y.Origin = double(Y.Origin);

% Row-major for ITK
X = vuRowMajorColumnMajorSwitch(X);
Y = vuRowMajorColumnMajorSwitch(Y);

if isStruct
    if length(X.Dims) == 2
      [REG_IM,DEF_IM] = MEX2DDemonsRegistration(X,Y,p.Results.numIterations,p.Results.sigma,p.Results.verbose);
			XX = vuGenerateMetaImage(squeeze(DEF_IM.Data(1,:,:)),DEF_IM.Spc,DEF_IM.Origin);
			YY = vuGenerateMetaImage(squeeze(DEF_IM.Data(2,:,:)),DEF_IM.Spc,DEF_IM.Origin);
    else
      [REG_IM,DEF_IM] = MEX3DDemonsRegistration(X,Y,p.Results.numIterations,p.Results.sigma,p.Results.verbose);
			XX = vuGenerateMetaImage(squeeze(DEF_IM.Data(1,:,:,:)),DEF_IM.Spc,DEF_IM.Origin);
			YY = vuGenerateMetaImage(squeeze(DEF_IM.Data(2,:,:,:)),DEF_IM.Spc,DEF_IM.Origin);
			ZZ = vuGenerateMetaImage(squeeze(DEF_IM.Data(3,:,:,:)),DEF_IM.Spc,DEF_IM.Origin);
    end
else
    if length(X.Dims) == 2
      [reg_image,def_image] = MEX2DDemonsRegistration(X,Y,p.Results.numIterations,p.Results.sigma,p.Results.verbose);
      REG_IM = reg_image.Data;
			XX = squeeze(def_image.Data(1,:,:));
			YY = squeeze(def_image.Data(2,:,:));
    else
      [reg_image,def_image] = MEX3DDemonsRegistration(X,Y,p.Results.numIterations,p.Results.sigma,p.Results.verbose);
      REG_IM = reg_image.Data;
			XX = squeeze(def_image.Data(1,:,:,:));
			YY = squeeze(def_image.Data(2,:,:,:));
			ZZ = squeeze(def_image.Data(3,:,:,:));
    end
end

% Copy Orientation if it exists
if (isfield(X,'Orientation'))
    REG_IM.Orientation = X.Orientation;
end
% Copy parmeters if they exist
if (isfield(X,'Parms'))
    REG_IM.Parms = X.Parms;
end

% Column-major for Matlab
REG_IM = vuRowMajorColumnMajorSwitch(REG_IM);
XX = vuRowMajorColumnMajorSwitch(XX);
YY = vuRowMajorColumnMajorSwitch(YY);
if length(X.Dims) == 3
		ZZ = vuRowMajorColumnMajorSwitch(ZZ);
end