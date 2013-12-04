function IM = vuOpenImage(filename, varargin)
% vuOpenImage allows users to select and open image files
%
%   SYNTAX:
%       IM = vuOpenImage;
%       IM = vuOpenImage(filename, options);
%
%   OPTIONS:
%       dimensions : '2','3' : Specify Dimensions of Image  (Note: only
%       needs to be set for 2D mha, nii, and hdr images
%       transform : For fid data -> Transform data to image space (1 or 0)
%       sliceOrientation : Slice orientation for SPAR (1=TRA,2=SAG,3=COR)
%
%   OPTIONS & DEFAULTS:
%       dimensions = 3;
%       transform = 1;
%       sliceOrientation = 3;
%
%   OUTPUT:
%       IM is the opened image
%
%   EXAMPLE:
%       im = vuOpenImage;
%       im2 = vuOpenImage('Brain1.dat');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

% Get and parse options
p = inputParser;
p.addParamValue('dimensions',3,@(x) isa(x,'double'));
p.addParamValue('transform',1,@(x) isa(x,'double'));
p.addParamValue('sliceOrientation',3,@(x) isa(x,'double'));
p.parse(varargin{:});

% Set default pathname
pathname = '';

if nargin < 1 || strcmp(filename,'click')
    % Get filename
    [filename, pathname, filterindex] = uigetfile({'*.dat;*.DAT;*fid;*FID;*.log;*.LOG;*.par;*.PAR;*.fdf;*.FDF;*.spar;*.SPAR;*.mha;*.MHA;*.hdr;*.HDR;*.nii;*.NII;*.nii.gz;*.NII.GZ;*.dcm;*.DCM','Medical Image Files'; ...
        '*.mat;*.MAT;*.dat;*.DAT','Matlab Data';'*.bmp;*.BMP;*.cur;*.CUR;*.gif;*.GIF;*.hdf;*.HDF;*.ico;*.ICP;*.jpg;*.JPG;*.jpeg;*.JPEG;*.pbm;*.PBM;*.pcx;*.PCX;*.pgm;*.PGM;*.png;*.PNG;*.pnm;*.PNM;*.ppm;*.PPM;*.ras;*.RAS;*.tif;*.TIF;*.tiff;*.TIFF;*.xwd;*.XWD', 'Image Files'; ...
        '*.*', 'All Files (*.*)'},'Pick an image file ...');

    % Check if file was selected
    if (~filterindex)
        IM = -1;
        return
    end
end

% full file name
fullFilename = fullfile(pathname,filename);

if(~exist(fullFilename,'file'))
    error('Could not find specified file.')
end


% Check for fid type, else get extension
if (strcmp(filename(end-2:end),'fid'))
    ext = 'fid';
else
    ext = filename(end-3:end);
end

% Find file type
isMed = regexpi('fid .log .par .fdf .mha spar .hdr .nii i.gz .dcm .dat',ext);
isMatlab = regexpi('.mat .dat',ext);
isImage = regexpi('.bmp .cur .gif .hdf .ico .jpg .jpe .pbm .pcx .pgm .png .pnm .ppm .ras .tif .xwd',ext);

% Medical Image go to PAR, or ITK
if(isMed)
    % Check for PAR-REC (reader below)
    if (regexpi('.par',ext))
        disp('Medical Image Read (PAR-REC) ...')
        [im, spc, origin, orient, parms] = readPAR(fullFilename);
        IM.Data = squeeze(single(im));
        IM.Dims = size(IM.Data);
        IM.Spc = spc;
        IM.Origin = origin;
        IM.Orientation = orient;
        IM.Parms = parms;
    elseif (regexpi('.log',ext))
        disp('Medical Image Read (LOG-RAW) ...')
        [im, spc, origin, orient] = readLOG(fullFilename);
        IM.Data = squeeze(single(im));
        IM.Dims = size(IM.Data);
        IM.Spc = spc;
        IM.Origin = origin;
        IM.Orientation = orient;
    elseif (regexpi('.fdf',ext))
        disp('Medical Image Read (FDF) ...')
        [im, spc, origin, orient, parms] = readFDF(fullFilename);
        IM.Data = squeeze(single(im));
        IM.Dims = size(IM.Data);
        IM.Spc = spc;
        IM.Origin = origin;
        IM.Orientation = orient;
        if (isstruct(parms))
            IM.Parms = parms;
        end
    elseif (regexpi('fid',ext))
        disp('Medical Image Read (FID) ...')
        [im, dims, spc, origin, orient, parms] = readFID(fullFilename,p.Results.transform);
        IM.Data = squeeze(single(im));
        IM.Dims = dims;
        IM.Spc = spc;
        IM.Origin = origin;
        IM.Orientation = orient;
        if (isstruct(parms))
            IM.Parms = parms;
        end
    elseif (regexpi('spar',ext))
        disp('Medical Image Read (SPAR) ...')
        [data, span, origin, orient, parms] = readSPAR(fullFilename,p.Results.sliceOrientation);
        IM.Data = data;
        IM.Dims = size(IM.Data);
        IM.Span = span;
        IM.Origin = origin;
        IM.Orientation = orient;
        if (isstruct(parms))
            IM.Parms = parms;
        end
    elseif (regexpi('.nii',ext))
        disp('Medical Image Read (NII) ...')
        [im, dims, spc, origin, orient, parms] = readNII(fullFilename);
        IM.Data = squeeze(single(im));
        IM.Dims = dims;
        IM.Spc = spc;
        IM.Origin = origin;
        IM.Orientation = orient;
        if (isstruct(parms))
            IM.Parms = parms;
        end
    elseif (regexpi('.dat',ext))
        % Test if it is from NeuroScan
        fid = fopen(fullFilename);
        x = fgets(fid);
        fclose(fid);
        if (~isempty(findstr(x,'[Subject]')))
            disp('Neuro Scan Data Read (DAT) ...')
            IM = readDAT(fullFilename);
        else
            disp('Matlab Data Read ...')
            IM = load(fullFilename);
        end
    elseif (regexpi('.hdr',ext))
        % Test if it is from our PET
        fid = fopen(fullFilename);
        x = fgets(fid);
        fclose(fid);
        if (x(1)=='#')
            disp('Medical Image Read (HDR) ...')
            [data, spc, origin, orient, parms] = readHDR(fullFilename);
            IM.Data = data;
            IM.Dims = size(IM.Data);
            IM.Spc = spc;
            IM.Origin = origin;
            IM.Orientation = orient;
            IM.Parms = parms;
        else
            % Analyis format
            dims = p.Results.dimensions;
            if (dims==2)
                disp('2-D Medical Image Read ...')
                IM = MEX2DITKImageFileReader(fullFilename);
            elseif (dims==3)
                disp('3-D Medical Image Read ...')
                IM = MEX3DITKImageFileReader(fullFilename);
            elseif (dims==4)
                disp('4-D Medical Image Read ...')
                IM = MEX4DITKImageFileReader(fullFilename);
            else
                error('MATLAB:vuOpenImage:UnsupportedDimension','vuOpenImage does not support the input image dimension for this type of file')
            end
            try
                % Read and store header
                [Parms,be] = read_hdr_raw(fullFilename);
                Parms.be = be;
                % Adjust Orientation Matrix
                orient = Q2M(double([Parms.quatern_b Parms.quatern_c Parms.quatern_d]));
                IM.Orientation = orient(1:3,1:3);
                IM.Parms = Parms;
            catch
                disp('NifTi header read failed');
            end
        end
    % Let ITK try
    else
        dims = p.Results.dimensions;
        if (dims==2)
            disp('2-D Medical Image Read ...')
            IM = MEX2DITKImageFileReader(fullFilename);
        elseif (dims==3)
            disp('3-D Medical Image Read ...')
            IM = MEX3DITKImageFileReader(fullFilename);
        elseif (dims==4)
            disp('4-D Medical Image Read ...')
            IM = MEX4DITKImageFileReader(fullFilename);
        else
            error('MATLAB:vuOpenImage:UnsupportedDimension','vuOpenImage does not support the input image dimension for this type of file')
        end
        if (~isempty(regexpi('.nii',ext))||~isempty(regexpi('i.gz',ext)))
           try
                % Read and store header
                [Parms,be] = read_hdr_raw(fullFilename);
                Parms.be = be;
                % Adjust Orientation Matrix
                orient = Q2M(double([Parms.quatern_b Parms.quatern_c Parms.quatern_d]));
                IM.Orientation = orient(1:3,1:3);
                IM.Parms = Parms;
            catch
                disp('NifTi header read failed');
            end
        end
    end
    % Column-major for Matlab
    IM = vuRowMajorColumnMajorSwitch(IM);
    
    
% Matlab data load
elseif(isMatlab)
    disp('Matlab Data Read ...')
    IM = load(fullFilename);
    
% Regular Image go to Matlab
elseif(isImage)
    disp('Regular Image Read ...')
    IM = imread(fullFilename);
    
% Try to open deliminated images
else
    disp('Image Not Reconized ... Attempting to Read ...')
    try
        IM = dlmread(fullFilename);
        
        % Column-major for Matlab
        IM = vuRowMajorColumnMajorSwitch(IM);
    catch
        disp('Image Read Unsuccessful ... File format not reconized.')
        IM = -1;
    end
end

function [im, dims, spc, origin, orient, nifti] = readNII(fullFilename)
% NII Reader
% By: Jeff Luci
% Adapted By: Kevin Wilson
% Date: June, 2008
% Copyright (c) 2008 - Vanderbilt University Institute of Imaging Science

fid=fopen(fullFilename, 'rb');
if (fid == -1)
    disp('.NII image open unsuccessful')
    data = -1;
    return;
end


nifti.sizeof_hdr = fread(fid, 1, 'int'); %sizeof_hdr
nifti.data_type =fread(fid, 10, 'char'); %UNUSED data_type
nifti.db_name =fread(fid, 18, 'char'); %UNUSED db_name
nifti.extents =fread(fid, 1, 'int'); %UNUSED extents
nifti.session_error =fread(fid, 1, 'short'); %UNUSED session_error
nifti.regular =fread(fid, 1, 'char'); %UNUSED regular
nifti.dim_info =fread(fid, 1, 'char'); %dim_info
nifti.dim =fread(fid, 8, 'short'); %dim
nifti.intent_p1 =fread(fid, 1, 'float'); %intent parameter 1 ***
nifti.intent_p2 =fread(fid, 1, 'float'); %intent parameter 2 ***
nifti.intent_p3 =fread(fid, 1, 'float'); %intent parameter 3 ***
nifti.intent_code =fread(fid, 1, 'short'); %NIFIT intent code ***
nifti.datatype =fread(fid, 1, 'short'); %datatype
nifti.bitpix =fread(fid, 1, 'short'); %bitpix
nifti.slice_start =fread(fid, 1, 'short'); %first slice number
nifti.pix_dim =fread(fid, 8, 'float'); %grid spacings ***
nifti.vox_offset =fread(fid, 1, 'float'); %vox_offset
nifti.scl_slope =fread(fid, 1, 'float'); %scl_slope
nifti.scl_inter =fread(fid, 1, 'float'); %scl_inter
nifti.slice_end =fread(fid, 1, 'short'); %last slice number
nifti.slice_code =fread(fid, 1, 'char'); %slice timing code *** WHY CHAR?
nifti.xyzt_units =fread(fid, 1, 'char'); %xyzt_units
nifti.cal_max =fread(fid, 1, 'float'); %display max
nifti.cal_min =fread(fid, 1, 'float'); %display min
nifti.slice_duration =fread(fid, 1, 'float'); %slice time
nifti.toffset = fread(fid, 1, 'float'); %toffset
nifti.glmax =fread(fid, 1, 'int'); %unused
nifti.glmin =fread(fid, 1, 'int'); %unused
nifti.descrip =fread(fid, 80, 'char'); %description, use Varian seqfil
nifti.aux_file =fread(fid, 24, 'char'); %alternate filename
nifti.qform_code =fread(fid, 1, 'short'); %qform code ***
nifti.sform_code =fread(fid, 1, 'short'); %sform code ***
nifti.quartern_b =fread(fid, 1, 'float'); %Quarternion b param ***
nifti.quartern_c =fread(fid, 1, 'float'); %Quarternion c param ***
nifti.quartern_d =fread(fid, 1, 'float'); %Quarternion d param ***
nifti.qoffset_x =fread(fid, 1, 'float'); %Quarternion x param ***
nifti.qoffset_y =fread(fid, 1, 'float'); %Quarternion y param ***
nifti.qoffset_z =fread(fid, 1, 'float'); %Quarternion z param ***
nifti.srow_x =fread(fid, 4, 'float'); %1st row affine transform
nifti.srow_y =fread(fid, 4, 'float'); %2nd row affine transform
nifti.srow_z =fread(fid, 4, 'float'); %3rd row affine transform
nifti.intent_name =fread(fid, 16, 'char'); %name of data
nifti.magic =fread(fid, 4, 'char'); %Magic string

tmp = fread(fid,4,'char'); % 4 blank spaces (Assume not a extended header)
b = nifti.quartern_b;
c = nifti.quartern_c;
d = nifti.quartern_d;
a = sqrt(1.0-(b*b+c*c+d*d));
orient = [a*a+b*b-c*c-d*d 2*b*c-2*a*d 2*b*d+2*a*c;2*b*c+2*a*d a*a+c*c-b*b-d*d 2*c*d-2*a*b;2*b*d-2*a*c 2*c*d+2*a*b a*a+d*d-c*c-b*b];
if (nifti.pix_dim(1) < 0)
    orient = orient*[1 0 0;0 1 0;0 0 -1];
end
origin = [nifti.qoffset_x nifti.qoffset_y nifti.qoffset_z];
spc = nifti.pix_dim(2:nifti.dim(1)+1)';
if (length(spc)>3)
    spc = spc(1:3);
end
dims = nifti.dim(2:nifti.dim(1)+1)';

switch nifti.datatype
    case 2
        datatype = 'uint8=>single';
    case 4
        datatype = 'int16=>single';
    case 8
        datatype = 'int32=>single';
    case 16
        datatype = 'float32=>single';
    case 64
        datatype = 'double=>single';
    case 256
        datatype = 'int8=>single';
    case 512
        datatype = 'uint16=>single';
    case 768
        datatype = 'uint32=>single';
    case 1024 
        datatype = 'int64=>single';
    case 1280
        datatype = 'uint64=>single';
    otherwise
        disp('.NII image open unsuccessful : Datatype Unknown')
        data = -1;
        return;
end

im = fread(fid,datatype);
im = reshape(im,dims);

% Rescale
if isfield(nifti,'scl_slope')||isfield(nifti,'scl_inter')
    slope = 1;
    inter = 0;
    if isfield(nifti,'scl_slope')&&nifti.scl_slope~=0
        slope = nifti.scl_slope;
    end
    if isfield(nifti,'scl_inter')
        inter = nifti.scl_inter;
    end
    im = im*slope + inter;
end


fclose(fid);


function [data, spc, origin, orient] = readLOG(fullFilename)
% LOG/RAW Reader
% By: Kevin Wilson
% Date: May 1, 2007
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science
fid = fopen(fullFilename,'r');
if (fid == -1)
    disp('.LOG/RAW image open unsuccessful')
    data = -1;
    return;
end

% get line until end of header
line = fgets(fid);
while (~feof(fid))    
    
    % Get origin
    if (~isempty(findstr(line,'X Number of Voxels')))
        n = findstr(line,':');
        xdim = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'Y Number of Voxels')))
        n = findstr(line,':');
        ydim = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'Z Number of Voxels')))
        n = findstr(line,':');
        zdim = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'X Voxel Size')))
        n = findstr(line,':');
        xspc = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'Y Voxel Size')))
        n = findstr(line,':');
        yspc = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'Z Voxel Size')))
        n = findstr(line,':');
        zspc = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'X Origin')))
        n = findstr(line,':');
        xoff = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'Y Origin')))
        n = findstr(line,':');
        yoff = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'Z Origin')))
        n = findstr(line,':');
        zoff = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'Bed Axis Position')))
        n = findstr(line,':');
        bedAxis = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'Bed Height')))
        n = findstr(line,':');
        bedHeight = sscanf(line(n+1:length(line)),'%f');
    end

    line = fgets(fid);
end

fclose(fid);

% if no dimensions assume 512^3
if (~exist('xdim','var'))
    xdim = 512;
    ydim = 512;
    zdim = 512;
end

% dimensions
dims = [xdim ydim zdim];

% if no spacing assume  0.250
if (~exist('xspc','var'))
    xspc = 0.250;
    yspc = 0.250;
    zspc = 0.250;
end

% Spacing
spc = [xspc yspc zspc]; 

% if no offset assume  0
if (~exist('xoff','var'))
    xoff = 0;
    yoff = 0;
    zoff = 0;
end

% Choose Bed Axis = 333.108, Bed Height = 38.137 as default
bedAxisDefault = 333.108;
bedHeightDefault = 38.137;

% if no bed info assume defaults 
if (~exist('bedAxis','var'))
    bedAxis = bedAxisDefault;
    bedHeight = bedHeightDefault;
end

% Orientation (LOG file record oblique?)
orient = eye(3);

rawFile = strcat(fullFilename(1:length(fullFilename)-3),'raw');

if (length(dims)==3)&&(dims(1)==512 && dims(2)==512 && dims(3) == 512)
    disp('Resampling to 256.256.256');
    pause(1)
    data = MEX3DResampleOnFly(rawFile,0);
    data = vuRowMajorColumnMajorSwitch(data);
    dims = dims./2;
    spc = spc.*2;
else
    fid2 = fopen(rawFile,'r');
    if (fid2 == -1)
        rawFile = strcat(fullFilename(1:length(fullFilename)-3),'RAW');
        fid2 = fopen(rawFile,'r');
        if (fid2 == -1)
            error('Could not find raw file');
        end
    end
    data = fread(fid2,'int16');
    
    % zero pad if dimension don't match
    if (length(data) < prod(dims))
        for i = length(data):prod(dims)
            data(i) = 0;
        end
    end

    data = reshape(data,dims);

    fclose(fid2);
end

% origin
origin = [xoff yoff zoff].*spc - dims.*spc./2;
origin = origin + spc./2;
origin = origin - [0 bedHeightDefault-bedHeight bedAxis-bedAxisDefault];


return;

function [data, spc, origin, orient, parms] = readHDR(fullFilename)
% HDR Reader (VUIIS PET)
% By: Kevin Wilson
% Date: May 11, 2007
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science
fid = fopen(fullFilename,'r','ieee-be');
if fid == -1
    error('MATLAB:vuOpenImage:FileError',sprintf('Can not open file %s',fullname));
end

line = fgets(fid);
line = fgets(fid);

% Original Image Make some assumptions
rotation = [0 0 0];
translation = [0 0 0];
zflip = 'No';
pixelSpc = [];
pixelSpcX = [];
pixelSpcY = [];
pixelSpcZ = [];
totalFrames = 1;																													
zoom = 4;
calibrationFactor = 1;
scaleFactor = 1;
isotopeBranchingFraction = 1;

% Check if it is the original image, or modifications
if (~isempty(findstr(line,'modifications of')))
    dims = [];
    while(~feof(fid))
        line = fgets(fid);
        if (~isempty(findstr(line,'Rotations_deg(x,y,z)')))
            n = findstr(line,':');
            rotation = sscanf(line(n+1:length(line)),'%f, %f, %f')';
        elseif (~isempty(findstr(line,'Translations_mm(x,y,z)')))
            n = findstr(line,':');
            translation = sscanf(line(n+1:length(line)),'%f, %f, %f')';
        elseif (~isempty(findstr(line,'Images flipped about z-axis')))
            n = findstr(line,':');
            zflip = sscanf(line(n+1:length(line)),'%s');
        elseif (~isempty(findstr(line,'Matrix resliced to')))
            n = findstr(line,'Matrix resliced to');
            n = n+18;
            dims = sscanf(line(n+1:length(line)),'%fx%fx%f')';
            n = findstr(line,'pixel size');
            n = n+10;
            pixelSpc = sscanf(line(n+1:length(line)),'%f')';
            n = findstr(line,'slice thickness');
            n = n+15;
            sliceThk = sscanf(line(n+1:length(line)),'%f')';
        elseif (~isempty(findstr(line,'#----------')))
            break;
        end
    end
    % Check if the shortcut worked
    if (isempty(dims))
        % If not try the full version
        frewind(fid);
        while(~feof(fid))
            line = fgets(fid);
            if (~isempty(findstr(line,'total_frames')))
                n = 13;
                totalFrames = sscanf(line(n+1:length(line)),'%f');
            elseif (~isempty(findstr(line,'x_dimension')))
                n = 12;
                xdim = sscanf(line(n+1:length(line)),'%f');
            elseif (~isempty(findstr(line,'y_dimension')))
                n = 12;
                ydim = sscanf(line(n+1:length(line)),'%f');
            elseif (~isempty(findstr(line,'z_dimension')))
                n = 12;
                zdim = sscanf(line(n+1:length(line)),'%f');
            elseif (~isempty(findstr(line,'pixel_size')))
                if (~isempty(findstr(line,'pixel_size_x')))
                    n = 13;
                    pixelSpcX = sscanf(line(n+1:length(line)),'%f');
                elseif (~isempty(findstr(line,'pixel_size_y')))
                    n = 13;
                    pixelSpcY = sscanf(line(n+1:length(line)),'%f');
                elseif (~isempty(findstr(line,'pixel_size_z')))
                    n = 13;
                    pixelSpcZ = sscanf(line(n+1:length(line)),'%f');
                end
                n = 11;
                pixelSpc = sscanf(line(n+1:length(line)),'%f');
                % Assuming spacing for now
                pixelSpc = pixelSpc*10;
                sliceThk = 0.796;
            end
        end     
        dims = [xdim ydim zdim];
    end
else
    while(~feof(fid))
        line = fgets(fid);
        if (~isempty(findstr(line,'total_frames')))
            n = 13;
            totalFrames = sscanf(line(n+1:length(line)),'%f');
        elseif (~isempty(findstr(line,'x_dimension')))
            n = 12;
            xdim = sscanf(line(n+1:length(line)),'%f');
        elseif (~isempty(findstr(line,'y_dimension')))
            n = 12;
            ydim = sscanf(line(n+1:length(line)),'%f');
        elseif (~isempty(findstr(line,'z_dimension')))
            n = 12;
            zdim = sscanf(line(n+1:length(line)),'%f');
        elseif (~isempty(findstr(line,'ending_bed_offset')))
            n = 18;
            translation(3) = sscanf(line(n+1:length(line)),'%f');
        elseif (~isempty(findstr(line,'vertical_bed_offset')))
            n = 20;
            translation(2) = sscanf(line(n+1:length(line)),'%f');
        elseif (~isempty(findstr(line,'isotope_branching_fraction')))
            n = 27;
            isotopeBranchingFraction = sscanf(line(n+1:length(line)),'%f');
        elseif (~isempty(findstr(line,'scale_factor')))
            n = 13;
            scaleFactor = sscanf(line(n+1:length(line)),'%f');
        elseif (~isempty(findstr(line,'calibration_factor')))
            n = 19;
            calibrationFactor = sscanf(line(n+1:length(line)),'%f');
        elseif (~isempty(findstr(line,'zoom'))&&findstr(line,'zoom')==1)
            n = 5;
            zoom = sscanf(line(n+1:length(line)),'%f');  
        elseif (~isempty(findstr(line,'pixel_size')))
            if (~isempty(findstr(line,'pixel_size_x')))
                n = 13;
                pixelSpcX = sscanf(line(n+1:length(line)),'%f');
            elseif (~isempty(findstr(line,'pixel_size_y')))
                n = 13;
                pixelSpcY = sscanf(line(n+1:length(line)),'%f');
            elseif (~isempty(findstr(line,'pixel_size_z')))
                n = 13;
                pixelSpcZ = sscanf(line(n+1:length(line)),'%f');
            end
            n = 11;
            pixelSpc = sscanf(line(n+1:length(line)),'%f');
            % Assuming spacing for now
            pixelSpc = pixelSpc*10;
            sliceThk = 0.796;
        end
    end     
    dims = [xdim ydim zdim];
end

% Spacing
if (~isempty(pixelSpcX)&&~isempty(pixelSpcY)&&~isempty(pixelSpcZ))
    spc = [pixelSpcX pixelSpcY pixelSpcZ];
elseif (~isempty(pixelSpc))
    spc = [pixelSpc pixelSpc sliceThk]; 
else
    spc = [0.474519 0.474519 0.796];
end

if (spc(1)==0)
    if (zoom ==1)
        spc=[1.8981 1.8981 0.796];
    elseif (zoom == 2)
        spc=[0.949039 0.949039 0.796];
    elseif (zoom ==4)
        spc=[0.474519 0.474519 0.796];
    end
end

% parms
parms.BedOffset = translation*10;
parms.Zoom = zoom;
parms.CalibrationFactor = calibrationFactor;
parms.ScaleFactor = scaleFactor;
parms.IsotopeBranchingFraction = isotopeBranchingFraction;

% origin
origin = -1*dims.*spc/2;
origin = origin + spc./2;

% Orientation
xrot = [1 0 0;0 cosd(rotation(1)) -sind(rotation(1));0 sind(rotation(1)) cosd(rotation(1))];
yrot = [cosd(rotation(2)) 0 sind(rotation(2));0 1 0;-sind(rotation(2)) 0 cosd(rotation(2))];
zrot = [cosd(rotation(3)) -sind(rotation(3)) 0;sind(rotation(3)) cosd(rotation(3)) 0;0 0 1];
orient = yrot*xrot*zrot;
if (strcmp(zflip,'Yes'))
    orient = orient*[1 0 0; 0 1 0; 0 0 -1];
end
% Flip matrix, because it will be flipped back above in
% vuRowMajortoColumnMajorSwitch
orient = permute(orient,[2 1]);

imgFile = fullFilename(1:length(fullFilename)-4);
fid2 = fopen(imgFile,'r');
if (fid2 == -1)
    error('Could not find img file');
end

data = fread(fid2,'float=>single');

if (totalFrames > 1)
    dims(4) = totalFrames;
end

% zero pad if dimension don't match
if (length(data) < prod(dims))
    data(length(data)+1:prod(dims)) = 0;
end

data = reshape(data,dims);

fclose(fid2);
return;

function [im, spc, origin, orient, parms] = readFDF(fullFilename)
% FDF Reader
% By: Kevin Wilson
% Date: March 28, 2007
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science


% Check Endian type
ch=' ';
header='         ';
fid = fopen(fullFilename, 'r');
if fid == -1
        str = sprintf('Can not open file %s',fullFilename);
        error(str);
end
while(uint8(ch)~=0)
    ch=fread(fid, 1, 'uint8=>char');
    header=[header ch];
    if(strcmp(header(numel(header)-8:numel(header)), 'bigendian'))
        str = fread(fid, 4, 'uint8=>char');
        bigEndian = str2num(str(4));
        ch=uint8(10);
    end
end
fclose(fid);

if(~exist('bigEndian', 'var') || bigEndian==1)
    byteOrder='ieee-be';
else
    byteOrder='ieee-le';
end

fid = fopen(fullFilename,'r',byteOrder);
if (fid == -1)
    disp('FDF image open unsuccessful')
    im = -1;
    return;
end
% get line until end of header
line = fgets(fid);
while isempty(findstr('checksum',line))

    % Get binary data type
    if (~isempty(findstr(line,'storage')))
        clist = findstr(line,'"');
        dtype = sscanf(line(clist(1)+1:clist(2)-1),'%s');
    end
    
    % Get size of each number
    if (~isempty(findstr(line,'bits')))
        n = findstr(line,'=');
        bits = sscanf(line(n+1:length(line)),'%d');
    end
    
    % Get matrix size
    if (~isempty(findstr(line,'matrix')))
        n = findstr(line,'{');
        dims = sscanf(line(n+1:length(line)),'%d, %d, %d')';
        if(length(dims)==2)
            dims(3) = 1;
        end
    end
      
    % Get matrix size
    if (~isempty(findstr(line,'roi[]')))
        n = findstr(line,'{');
        span = sscanf(line(n+1:length(line)),'%f, %f, %f')';
    end
    
    % Get matrix size
    if (~isempty(findstr(line,'location')))
        n = findstr(line,'{');
        origin = sscanf(line(n+1:length(line)),'%f, %f, %f')';
    end
    
    % Get Orientation
    if (~isempty(findstr(line,'orientation[]')))
        n = findstr(line,'{');
        orient = sscanf(line(n+1:length(line)),'%f, %f, %f, %f, %f, %f, %f, %f, %f')';
        % Note: This is Column major, will be flipped above in
        % vuRowMajorToColumnMajorSwitch
        orient = reshape(orient,[3 3]);
    end

    line = fgets(fid);
end

% Put together percision
if (strcmp(dtype(1:3),'flo'))
    precision = sprintf('float%d',bits);
elseif (strcmp(dtype(1:3),'int'))
    precision = sprintf('int%d',bits);
end    

% Skip past NULL character that separates header and data
k = fread(fid,1,'uchar');
while k ~= 0
    k = fread(fid,1,'uchar');
end

im = fread(fid,precision);

% zero pad if dimension don't match
if (length(im) < prod(dims))
    for i = length(im):prod(dims)
        im(i) = 0;
    end
end

im = reshape(im,dims);
spc = (span*10)./dims;

% Origin in top-left
origin = origin*10 - ((span*10)./2);
origin = origin + spc./2;

% Read procpar (if exists)
path = fileparts(fullFilename);
procparFile = fullfile(path,'procpar');
if (exist(procparFile,'file'))
    pp = parsepp(procparFile);
    parms = pp;
    % Take of phase shifts from origin (PE shifts not implemented on data from
    % the MR.)
    origin = origin + [(pp.ppe*10) 0 0];
else
    parms = 0;
end

fclose(fid);

return;


function [data, span, origin, orient, par] = readSPAR(fullFilename,sliceOri)
% SPAR Reader
% By: Kevin Wilson
% Date: March 29, 2007
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science
fid = fopen(fullFilename,'r');
if (fid == -1)
    disp('SPAR image open unsuccessful')
    data = -1;
    return;
end

% get line until end of header
line = fgets(fid);
while (~feof(fid))    
    
    % Get origin
    if (~isempty(findstr(line,'ap_off_center')))
        n = findstr(line,':');
        apoff = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'lr_off_center')))
        n = findstr(line,':');
        lroff = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'cc_off_center')))
        n = findstr(line,':');
        ccoff = sscanf(line(n+1:length(line)),'%f');
    end
    
    % Get span
    if (~isempty(findstr(line,'ap_size')))
        n = findstr(line,':');
        apspan = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'lr_size')))
        n = findstr(line,':');
        lrspan = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'cc_size')))
        n = findstr(line,':');
        ccspan = sscanf(line(n+1:length(line)),'%f');
    end
    
    % Get orientation
    if (~isempty(findstr(line,'ap_angulation')))
        n = findstr(line,':');
        apori = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'lr_angulation')))
        n = findstr(line,':');
        lrori = sscanf(line(n+1:length(line)),'%f');
    end
    if (~isempty(findstr(line,'cc_angulation')))
        n = findstr(line,':');
        ccori = sscanf(line(n+1:length(line)),'%f');
    end
    
    % Get dims
    if (~isempty(findstr(line,'rows :')))
        n = findstr(line,':');
        rows = sscanf(line(n+1:length(line)),'%d');
    end
    if (~isempty(findstr(line,'samples :')))
        n = findstr(line,':');
        cols = sscanf(line(n+1:length(line)),'%d');
    end
    
    % Get parameters
    if (~isempty(findstr(line,'synthesizer_frequency :')))
        n = findstr(line,':');
        synthesizer_frequency = sscanf(line(n+1:length(line)),'%d');
    end
    if (~isempty(findstr(line,'sample_frequency :')))
        n = findstr(line,':');
        sample_frequency = sscanf(line(n+1:length(line)),'%d'); 
    end

    line = fgets(fid);
end

fclose(fid);

% parameters
par.synthesizer_frequency = synthesizer_frequency;
par.sample_frequency = sample_frequency;

% dimensions
dims = [rows cols];

% x=rl, y=ap, z=fh
if (sliceOri == 1)
    xoffC = lroff;
    yoffC = apoff;
    zoffC = ccoff;
    xori = lrori;
    yori = apori;
    zori = ccori;
    fov = [lrspan apspan ccspan];
% x=ap, y=hf, z=lr
elseif (sliceOri == 2)
    xoffC = apoff;
    yoffC = -ccoff;
    zoffC = -lroff;
    xori = apori;
    yori = -ccori;
    zori = -lrori;
    fov = [apspan ccspan lrspan];
% x=rl, y=hf, z=ap
elseif (sliceOri == 3)
    xoffC = lroff;
    yoffC = -ccoff;
    zoffC = apoff;
    xori = lrori;
    yori = -ccori;
    zori = apori;
    fov = [lrspan ccspan apspan];
else
    error('Unreconized Slice Orientation');
end

% Span
span = fov; 

% origin
origin = [xoffC yoffC zoffC]-(span./2);

% Orientation
xrot = [1 0 0;0 cosd(xori) -sind(xori);0 sind(xori) cosd(xori)];
yrot = [cosd(yori) 0 sind(yori);0 1 0;-sind(yori) 0 cosd(yori)];
zrot = [cosd(zori) -sind(zori) 0;sind(zori) cosd(zori) 0;0 0 1];
orient = xrot*yrot*zrot;
% Flip matrix, because it will be flipped back above in
% vuRowMajortoColumnMajorSwitch
orient = permute(orient,[2 1]);

sdatFile = strcat(fullFilename(1:length(fullFilename)-4),'sdat');
fid2 = fopen(sdatFile,'r','vaxd');
if (fid2 == -1)
    sdatFile = strcat(fullFilename(1:length(fullFilename)-4),'SDAT');
    fid2 = fopen(sdatFile,'r','vaxd');
    if (fid2 == -1)
        warning('Could not find sdat file');
    end
end
try
    data = fread(fid2, 'float');
    data = data([1:2:length(data)]) + i*data([2:2:length(data)]);
    data = reshape(data,dims);
    fclose(fid2);
catch
    data = [];
end
return;

function [im, spacing, origin, orient, parms] = readPAR(fullFilename)
% PAR-REC Reader
% Original Date: 16-nov-2004
% Author: Benoit Desjardins, MD, PhD
%         Department of Radiology
%         University of Michigan
%
% JCG 1/27/05
% Modified for 1, input of parfile name directly
%              2, Data format output selection
%              3, Platform independent file conversion
%
%
% Edited 2/14/05 to correctly parse multiple echoes
%
% Adapted for vuOpenImage
% By: Kevin Wilson
% Date: April 13, 2007
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

% Get data
[parms,ver] = read_parfile(fullFilename);
flt_data = read_recfile([fullFilename(1:(end-4)) '.REC'],parms,ver,'FP');

% Format for vuOpenImage
im = squeeze(flt_data);
if(strcmp(ver,'V3'))
    sliceOri = parms.tags(1,20);
else
    sliceOri = parms.tags(1,26);
end

if(strcmp(ver,'V3'))
   if(parms.scan_percentage ~= 100)
       warning('MATLAB:vuOpenImage:ScanPercent', 'The scan percentage of the image is not 100 percent.  Check the spacing of the image');
   end
end

% x=rl, y=ap, z=fh
if (sliceOri == 1)
    xoffC = parms.offcenter(3);
    yoffC = parms.offcenter(1);
    zoffC = parms.offcenter(2);
    xori = parms.angulation(3);
    yori = parms.angulation(1);
    zori = parms.angulation(2);
    fov = [parms.fov(3) parms.fov(1) parms.fov(2)];
% x=ap, y=hf, z=lr
elseif (sliceOri == 2)
    xoffC = parms.offcenter(1);
    yoffC = -parms.offcenter(2);
    zoffC = -parms.offcenter(3);
    xori = parms.angulation(1);
    yori = -parms.angulation(2);
    zori = -parms.angulation(3);
    fov = [parms.fov(1) parms.fov(2) parms.fov(3)];
% x=rl, y=hf, z=ap
elseif (sliceOri == 3)
    xoffC = parms.offcenter(3);
    yoffC = -parms.offcenter(2);
    zoffC = parms.offcenter(1);
    xori = parms.angulation(3);
    yori = -parms.angulation(2);
    zori = parms.angulation(1);
    fov = [parms.fov(3) parms.fov(2) parms.fov(1)];
else
    error('Unreconized Slice Orientation');
end

% Calculate outputs
dims = size(im);
if (length(dims) == 2)
		dims = [dims 1];
end
spacing = fov./dims(1:3);

% Reconstruction spacing
if(~strcmp(ver,'V3'))
    spacing(1:2) = [parms.tags(1,29) parms.tags(1,30)];
end

% Top-left
origin = [xoffC yoffC zoffC] - (spacing.*dims(1:3))/2;
origin = origin + spacing/2;

% Orientation
xrot = [1 0 0;0 cosd(xori) -sind(xori);0 sind(xori) cosd(xori)];
yrot = [cosd(yori) 0 sind(yori);0 1 0;-sind(yori) 0 cosd(yori)];
zrot = [cosd(zori) -sind(zori) 0;sind(zori) cosd(zori) 0;0 0 1];
orient = xrot*yrot*zrot;
% Flip matrix, because it will be flipped back above in
% vuRowMajortoColumnMajorSwitch
orient = permute(orient,[2 1]);

% save original par filename
parms.filename = fullFilename;
parms.version = ver;
parms.SliceOrientation = sliceOri;
return


%==========================================================================
function [parms,ver] = read_parfile(file_name) 

% read all the lines of the PAR file
nlines  = 0;
fid = fopen(file_name,'r');
if (fid < 1), error(['.PAR file ', file_name, ' not found.']); end;
while ~feof(fid)
    curline = fgetl(fid);
    if ~isempty(curline)
        nlines = nlines + 1;
        lines(nlines) = cellstr(curline); % allows variable line size
    end
end
fclose(fid);

% identify the header lines
NG = 0;
for L = 1:size(lines,2)
    curline = char(lines(L));
    if (size(curline,2) > 0 && curline(1) == '.')
       NG = NG + 1;
       geninfo(NG) = lines(L); 
    end
end
if (NG < 1), error('.PAR file has invalid format'); end;

% figure out if V3 or V4 PAR files
test_key = '.    Image pixel size [8 or 16 bits]    :';
if strmatch(test_key,geninfo); % only V3 has that key in the headers
    ver = 'V3';
    template=get_template_v3;
elseif strmatch('# CLINICAL TRYOUT             Research image export tool     V4.1',lines);
    ver = 'V4.1';
    template=get_template_v4;
elseif strmatch('# CLINICAL TRYOUT             Research image export tool     V4.2',lines);
    % Treat v4.2 same as v4.1
    ver = 'V4.1';
    template=get_template_v4;
else
    ver = 'V4';
    template=get_template_v4;
end;

% parse the header information
for S=1:size(template,1)
    line_key = char(template(S,1));
    value_type = char(template(S,2));
    field_name = char(template(S,3));
    L = strmatch(line_key,geninfo);

    if ~isempty(L)
        curline = char(geninfo(L));
        value_start = 1 + strfind(curline,':');
        value_end = size(curline,2);
    else
        value_type = ':-( VALUE NOT FOUND )-:';
    end 

    switch (value_type)
    case { 'float scalar' 'int   scalar' 'float vector' 'int   vector'}
         parms.(field_name) = str2num(curline(value_start:value_end));
    case { 'char  scalar' 'char  vector' }
         parms.(field_name) = strtrim(curline(value_start:value_end));
    otherwise
         parms.(field_name) = '';
    end

end

% parse the tags for each line of data
nimg  = 0;
mlines = char(lines);
imglines = mlines;
imglines(:) = '.';
for L = 1:size(mlines,1)
   curline = mlines(L,:);
   firstc = curline(1);
   if (size(curline,2) > 0 && firstc ~= '.' && firstc ~= '#' && firstc ~= '*')
      nimg = nimg + 1;
      imglines(nimg,:) = curline;
   end
end
imglines = imglines(1:nimg,:);
for L = 1:size(imglines,1)
   tags = str2num(imglines(L,:));
   if L==1
      parms.tags = zeros(size(imglines,1),length(tags));
   end
   parms.tags(L,:) = tags;
end
if (nimg < 1), error('Missing scan information in .PAR file'); end;

return;

%==========================================================================

function [image_data,dims] = read_recfile(recfile_name,parms,ver,out_fmt)

types_list = unique(parms.tags(:,5)'); % to handle multiple types
slice_list = unique(parms.tags(:,1)');
phase_list = unique(parms.tags(:,4)');
echo_list = unique(parms.tags(:,2)');
dynamic_list = unique(parms.tags(:,3)');
if(strcmp(ver,'V4.1'))
    diffb_list = unique(parms.tags(:,42)');
    grad_list = unique(parms.tags(:,43)');
    % if dwi reduce ndiffb by one (anatomic will be append on end of first
    % dwi series)
    if (size(diffb_list,2)>1)
        diffb_list = diffb_list(2:size(diffb_list,2));
        grad_list(end+1) = grad_list(end)+1;
    end
else
    diffb_list = 1;
    grad_list = 1;
end
scan_tag_size = size(parms.tags);
nimg = scan_tag_size(1);
nslice = size(slice_list,2); 
nphase = size(phase_list,2);
necho = size(echo_list,2);
ndyn = size(dynamic_list,2);
if(strcmp(ver,'V4.1'))
    ndiffb = size(diffb_list,2);
    ngrad = size(grad_list,2);
else
    ndiffb = 1;
    ngrad = 1;
end
ntype = size(types_list,2);
% no mix# indicated in the tags themselves

if ( isfield(parms,'recon_resolution') )
    nline = parms.recon_resolution(1);
    stride = parms.recon_resolution(2);
else
    nline = parms.tags(1,10);
    stride = parms.tags(1,11);
end

switch(ver)
case {'V3'}, pixel_bits = parms.pixel_bits;
case {'V4'}, pixel_bits = parms.tags(1,8);  % assume same for all imgs
case {'V4.1'}, pixel_bits = parms.tags(1,8);  % assume same for all imgs
end

switch (pixel_bits)
    case { 8 }, read_type = 'int8';
    case { 16 }, read_type = 'short';
    otherwise, read_type = 'uchar';
end

% read the REC file
fid = fopen(recfile_name,'r','l');
[binary_1D,read_size] = fread(fid,inf,[read_type '=>single']);

fclose(fid);
if (read_size ~= nimg*nline*stride)
    disp(sprintf('Expected %d int.  Found %d int',nimg*nline*stride,read_size));
    if (read_size > nimg*nline*stride)
        error('.REC file has more data than expected from .PAR file')
    else
        error('.REC file has less data than expected from .PAR file')
    end
else
    %disp(sprintf('.REC file read sucessfully'));
end

% generate the final matrix of images
dims = [stride nline nslice necho nphase ntype ndyn ngrad ndiffb];
image_data=zeros(dims,'single');
for I  = 1:nimg
    slice = parms.tags(I,1);
    slice_idx = find(slice_list == slice);
    phase = parms.tags(I,4);
    phase_idx = find(phase_list == phase);
    type = parms.tags(I,5);
    type_idx = find(types_list == type);
    echo = parms.tags(I,2);
    echo_idx = find(echo_list == echo);
    dyn = parms.tags(I,3);
    dyn_idx = find(dynamic_list == dyn);
    %seq = parms.tags(I,6);
    rec = parms.tags(I,7);
    if(strcmp(ver,'V4.1'))
        diffb = parms.tags(I,42);
        grad = parms.tags(I,43);
        diffb_idx = find(diffb_list == diffb);
        grad_idx = find(grad_list == grad);
        % Put anatomic view at end of first diffusion volumes
        if ((parms.tags(I,46) == 0)&&(parms.tags(I,47) == 0)&&(parms.tags(I,48) == 0)&&(parms.tags(I,42) == 1))
            diffb_idx = 1;
            grad_idx = find(grad_list == grad_list(end));
        end
    else
        diffb_idx = 1;
        grad_idx = 1;
    end
    start_image = 1+rec*nline*stride;
    end_image = start_image + stride*nline - 1;
    img = reshape(binary_1D(start_image:end_image),stride,nline);
    % rescale data to produce FP information (not SV, not DV)
    img = rescale_rec(img,parms.tags(I,:),ver,out_fmt);
    image_data(:,:,slice_idx,echo_idx,phase_idx,type_idx,dyn_idx,grad_idx,diffb_idx) = img;
end
return;

%==========================================================================

function img = rescale_rec(img,tag,ver,out_fmt)

% transforms SV data in REC files to SV, DV or FP data for output
switch( ver )
case { 'V3' }, ri_i = 8; rs_i = 9; ss_i = 10;
case { 'V4' }, ri_i = 12; rs_i = 13; ss_i = 14;
case { 'V4.1' }, ri_i = 12; rs_i = 13; ss_i = 14;
end;
RI = tag(ri_i);  % 'inter' --> 'RI'
RS = tag(rs_i);  % 'slope' --> 'RS'
SS = tag(ss_i);  % new var 'SS'
switch (out_fmt)
    case { 'FP' }, img = (RS*img + RI)./(RS*SS);
    case { 'DV' }, img = (RS*img + RI);
    case { 'SV' }, 
end

return;

%==========================================================================
function [template] = get_template_v3  % header information for V3 PAR files

template = { ...                                  
'.    Patient name                       :'    'char  scalar'    'patient';    ...   
'.    Examination name                   :'    'char  scalar'    'exam_name';   ... 
'.    Protocol name                      :'    'char  vector'    'protocol';   ... 
'.    Examination date/time              :'    'char  vector'    'exam_date';  ...
'.    Acquisition nr                     :'    'int   scalar'    'acq_nr';    ...
'.    Reconstruction nr                  :'    'int   scalar'    'recon_nr';  ...
'.    Scan Duration [sec]                :'    'float scalar'    'scan_dur';        ...
'.    Max. number of cardiac phases      :'    'int   scalar'    'max_card_phases'; ...
'.    Max. number of echoes              :'    'int   scalar'    'max_echoes'; ...
'.    Max. number of slices/locations    :'    'int   scalar'    'max_slices'; ... 
'.    Max. number of dynamics            :'    'int   scalar'    'max_dynamics'; ... 
'.    Max. number of mixes               :'    'int   scalar'    'max_mixes'; ... 
'.    Image pixel size [8 or 16 bits]    :'    'int   scalar'    'pixel_bits'; ... 
'.    Technique                          :'    'char  scalar'    'technique'; ...  
'.    Scan mode                          :'    'char  scalar'    'scan_mode'; ... 
'.    Scan resolution  (x, y)            :'    'int   vector'    'scan_resolution'; ... 
'.    Scan percentage                    :'    'int   scalar'    'scan_percentage'; ... 
'.    Recon resolution (x, y)            :'    'int   vector'    'recon_resolution'; ... 
'.    Number of averages                 :'    'int   scalar'    'num_averages'; ... 
'.    Repetition time [msec]             :'    'float scalar'    'repetition_time'; ...   
'.    FOV (ap,fh,rl) [mm]                :'    'float vector'    'fov'; ... 
'.    Slice thickness [mm]               :'    'float scalar'    'slice_thickness'; ...
'.    Slice gap [mm]                     :'    'float scalar'    'slice_gap'; ... 
'.    Water Fat shift [pixels]           :'    'float scalar'    'water_fat_shift'; ... 
'.    Angulation midslice(ap,fh,rl)[degr]:'    'float vector'    'angulation'; ...
'.    Off Centre midslice(ap,fh,rl) [mm] :'    'float vector'    'offcenter'; ... 
'.    Flow compensation <0=no 1=yes> ?   :'    'int   scalar'    'flowcomp'; ...
'.    Presaturation     <0=no 1=yes> ?   :'    'int   scalar'    'presaturation';... 
'.    Cardiac frequency                  :'    'int   scalar'    'card_frequency'; ...
'.    Min. RR interval                   :'    'int   scalar'    'min_rr_interval'; ...
'.    Max. RR interval                   :'    'int   scalar'    'max_rr_interval'; ...
'.    Phase encoding velocity [cm/sec]   :'    'float vector'    'venc'; ... 
'.    MTC               <0=no 1=yes> ?   :'    'int   scalar'    'mtc'; ...
'.    SPIR              <0=no 1=yes> ?   :'    'int   scalar'    'spir'; ...
'.    EPI factor        <0,1=no EPI>     :'    'int   scalar'    'epi_factor'; ...
'.    TURBO factor      <0=no turbo>     :'    'int   scalar'    'turbo_factor'; ...
'.    Dynamic scan      <0=no 1=yes> ?   :'    'int   scalar'    'dynamic_scan'; ...
'.    Diffusion         <0=no 1=yes> ?   :'    'int   scalar'    'diffusion'; ...
'.    Diffusion echo time [msec]         :'    'float scalar'    'diffusion_echo_time'; ...
'.    Inversion delay [msec]             :'    'float scalar'    'inversion_delay'; ...
};

return;

%==========================================================================
function [template] = get_template_v4    % header information for V4 PAR files

template = { ...                                  
'.    Patient name                       :' 'char  scalar' 'patient';    ...   
'.    Examination name                   :' 'char  vector' 'exam_name';   ... 
'.    Protocol name                      :' 'char  vector' 'protocol';   ... 
'.    Examination date/time              :' 'char  vector' 'exam_date';  ...
'.    Series Type                        :' 'char  vector' 'series_type';  ...
'.    Acquisition nr                     :' 'int   scalar' 'acq_nr';    ...
'.    Reconstruction nr                  :' 'int   scalar' 'recon_nr';  ...
'.    Scan Duration [sec]                :' 'float scalar' 'scan_dur';        ...
'.    Max. number of cardiac phases      :' 'int   scalar' 'max_card_phases'; ...
'.    Max. number of echoes              :' 'int   scalar' 'max_echoes'; ...
'.    Max. number of slices/locations    :' 'int   scalar' 'max_slices'; ... 
'.    Max. number of dynamics            :' 'int   scalar' 'max_dynamics'; ... 
'.    Max. number of mixes               :' 'int   scalar' 'max_mixes'; ... 
'.    Patient position                   :' 'char  vector' 'patient_position'; ... 
'.    Preparation direction              :' 'char  vector' 'preparation_dir'; ... 
'.    Technique                          :' 'char  scalar' 'technique'; ...  
'.    Scan resolution  (x, y)            :' 'int   vector' 'scan_resolution'; ... 
'.    Scan mode                          :' 'char  scalar' 'scan_mode'; ... 
'.    Repetition time [ms]               :' 'float scalar' 'repetition_time'; ...   
'.    FOV (ap,fh,rl) [mm]                :' 'float vector' 'fov'; ... 
'.    Water Fat shift [pixels]           :' 'float scalar' 'water_fat_shift'; ... 
'.    Angulation midslice(ap,fh,rl)[degr]:' 'float vector' 'angulation'; ...
'.    Off Centre midslice(ap,fh,rl) [mm] :' 'float vector' 'offcenter'; ... 
'.    Flow compensation <0=no 1=yes> ?   :' 'int   scalar' 'flowcomp'; ...
'.    Presaturation     <0=no 1=yes> ?   :' 'int   scalar' 'presaturation';... 
'.    Phase encoding velocity [cm/sec]   :' 'float vector' 'venc'; ... 
'.    MTC               <0=no 1=yes> ?   :' 'int   scalar' 'mtc'; ...
'.    SPIR              <0=no 1=yes> ?   :' 'int   scalar' 'spir'; ...
'.    EPI factor        <0,1=no EPI>     :' 'int   scalar' 'epi_factor'; ...
'.    Dynamic scan      <0=no 1=yes> ?   :' 'int   scalar' 'dynamic_scan'; ...
'.    Diffusion         <0=no 1=yes> ?   :' 'int   scalar' 'diffusion'; ...
'.    Diffusion echo time [msec]         :' 'float scalar' 'diffusion_echo_time'; ...
};

return;

%==========================================================================
function procpar = parsepp(pathname)
%----------------------------------------
% Maj Hedehus, Varian, Inc., Oct 2001.
% modified by J. Luci for more general use.
%----------------------------------------

fullname=pathname;
fid = fopen(fullname,'r');
if fid == -1
    error(sprintf('Can not open file %s',fullname));
end

while ~feof(fid)
    par  = fscanf(fid,'%s',1);
    type  = fscanf(fid,'%d',1);
    fgetl(fid); %skip rest of line
    nvals = fscanf(fid,'%d',1);

    switch type
        case 1  % float
            eval(['procpar.' par '= fscanf(fid,''%f'',nvals);']);
            fgetl(fid);
            fgetl(fid);
        case 2  % string
            fullstr=[];
            for ii = 1:nvals,
                str = uint8(fgetl(fid));
                if str(1) == 32, 
                    str=str(2:end);  %strip off leading space if it exists
                end
                str = str(2:size(str,2)-1); %strip off double quotes from the string
                fullstr = [fullstr str 10]; %10 is the ASCII code for LF
            end
            fullstr=char(fullstr);
            eval(['procpar.' par '= fullstr;']);
            fgetl(fid);
        case 3  % delay
            eval(['procpar.' par '= fscanf(fid,''%f'',nvals);']);
            fgetl(fid);
            fgetl(fid);
        case 4  % flag
            L = uint8(fgetl(fid));
            if L(1) == 32, 
                L=L(2:end); %strip off leading space if it exists
            end 
            L = L(2:size(L,2)-1);  %strip off double quotes from the string
            L=char(L);
            eval(['procpar.' par '= L;']);
            fgetl(fid);
        case 5  % frequency      
            eval(['procpar.' par '= fscanf(fid,''%f'',nvals);']);
            fgetl(fid);
            fgetl(fid);
        case 6  % pulse
            L = str2double(fgetl(fid));
            eval(['procpar.' par '= L;']);
            skip=fgetl(fid);
        case 7  % integer
            num = str2num(fgetl(fid));
            eval(['procpar.' par '= num;']);
            fgetl(fid);
    end
end
fclose(fid);

return;

function [data, dims, spc, origin, orient, parms] = readFID(fullFilename,transform)
% FID Reader
% Original By: Maj Hedehus, Varian, Inc., Sep 2001.
% Adapted By: Jeff Luci
% Adapted For vuOpenImage By: Kevin Wilson
% Date: May 7, 2007
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

fid = fopen(fullFilename,'r','ieee-be');
if fid == -1
    error('MATLAB:vuOpenImage:FileError',sprintf('Can not open file %s',fullname));
end

% Read datafileheader
nblocks   = fread(fid,1,'int32');
ntraces   = fread(fid,1,'int32');
np        = fread(fid,1,'int32');
ebytes    = fread(fid,1,'int32');
tbytes    = fread(fid,1,'int32');
bbytes    = fread(fid,1,'int32');
vers_id   = fread(fid,1,'int16');
status    = fread(fid,1,'int16');
nbheaders = fread(fid,1,'int32');
 
 
s_data    = bitget(status,1);
s_spec    = bitget(status,2);
s_32      = bitget(status,3);
s_float   = bitget(status,4);
s_complex = bitget(status,5);
s_hyper   = bitget(status,6);

% added to prevent dynamic reallocation of resources
RE=zeros(np/2, ntraces*nblocks);
IM=zeros(np/2, ntraces*nblocks);
 
 
inx = 1;
B = 1;
for b = 1:nblocks
    sprintf('read block %d\n',b);
    % Read a block header
    scale     = fread(fid,1,'int16');
    bstatus   = fread(fid,1,'int16');
    index     = fread(fid,1,'int16');
    mode      = fread(fid,1,'int16');
    ctcount   = fread(fid,1,'int32');
    lpval     = fread(fid,1,'float32');
    rpval     = fread(fid,1,'float32');
    lvl       = fread(fid,1,'float32');
    tlt       = fread(fid,1,'float32');
    
    for t = 1:ntraces
        % We have to read data every time in order to increment file pointer
        if s_float == 1
            data = fread(fid,np,'float32');
        elseif s_32 == 1
            data = fread(fid,np,'int32');
        else
            data = fread(fid,np,'int16');
        end
        
        RE(:,inx) = data(1:2:np);
        IM(:,inx) = data(2:2:np);
        inx = inx + 1;

    end %trace loop
    
end  % done reading one block

% Read procpar (if exists)
path = fileparts(fullFilename);
procparFile = fullfile(path,'procpar');
if (exist(procparFile,'file'))
    pp = parsepp(procparFile);
    parms = pp;
    
    % Dimension
    dims = pp.np/2;
    fov = [pp.lro];
    if (pp.ne>1)
        dims = [dims pp.ne];
    end 
    if (pp.nv>0)
        dims = [dims pp.nv];
        fov = [fov pp.lpe];
    end
    if (pp.nv2>0)
        dims = [dims pp.nv2];
        fov = [fov pp.lpe2];
    end
    if (pp.nv3>0)
        dims = [dims pp.nv3];
        fov = [fov pp.lpe3];
    end
    if (pp.ns>=1 && pp.nv2==0)
        dims = [dims pp.ns];
        % Note : thk is in mm (Convert to cm)
        fov = [fov pp.ns*pp.thk/10];
    end
    if (pp.ns>=1 && ~isempty(strfind(parms.seqfil,'csi')))
        dims = [dims pp.ns];
    end
    % Transpoose fov and dims
    fov(1:2) = [fov(2) fov(1)];
    dims(1:2) = [dims(2) dims(1)];

    % Spacing
    spc = fov./dims(1:3);

    if (~isempty(strfind(pp.seqfil,'ge3d')))
        origin = [pp.ppe -pp.pro pp.ppe2] - fov(1:3)./2;
        origin = origin + spc./2;
    else
        origin = [pp.ppe -pp.pro pp.pss0] - fov(1:3)./2;
        origin = origin + spc./2;
    end
    
    % Orientation
    orient = CalculateRotationMatrix(pp.phi,pp.theta,-pp.psi);
    orient = permute(orient,[2 1]);
else
    dims = size(RE);
    spc = [1 1 1];
    origin = [0 0 0];
    orient = eye(3);
    parms = 0;
end
multicycle = 0;
data = RE + IM*i;
if (pp.ns>=1 && pp.nv2==0)
    try
        data = reshape(data,[dims(2) dims(end) dims(1) dims(3:end-1)]);
    catch
        data = reshape(data,[dims(2) dims(end) dims(1) dims(3:end-1) numel(data)/prod(dims)]);
        dims = [dims size(data,ndims(data))];
        multicycle = 1;
    end
else
    try
        data = reshape(data,[dims(2) dims(1) dims(3:end)]);
    catch
        data = reshape(data,[dims(2) dims(1) dims(3:end) parms.acqcycles]);
        dims = [dims parms.acqcycles];
        multicycle = 1;
    end
end
if (multicycle == 0)
    if (pp.ne>1)
        if (pp.ns>=1 && pp.nv2==0)
            data = permute(data,[1 4:length(dims) 2 3]);
            dims = [dims(2) dims(3:end) dims(1)];
            spc = [spc(2) spc(3:end) spc(1)];
            origin = [origin(2) origin(3:end) origin(1)];
        else
            data = permute(data,[1 3:length(dims) 2]);
            dims = [dims(1) dims(3:end) dims(2)];
        end
    elseif (pp.ns>=1 && pp.nv2==0)
        data = permute(data,[1 3:length(dims) 2]);
    end
else
    if (pp.ne>1)
        if (pp.ns>=1 && pp.nv2==0)
            data = permute(data,[1 4:(length(dims)-1) 3 2 length(dims)]);
            dims = [dims(1) dims(4:end-1) dims(3) dims(2) dims(end)];
        else
            data = permute(data,[1 3:(length(dims)-1) 2 length(dims)]);
            dims = [dims(1) dims(3:end-1) dims(2) dims(end)];
        end
    elseif (pp.ns>=1 && pp.nv2==0)
        data = permute(data,[1 3:(length(dims)-1) 2 length(dims)]);
    end
end

if (strfind(parms.seqfil,'mgems'))
    % Flip every other image
    for k = 1:2:size(data,4)
        data(:,:,:,k) = flipdim(data(:,:,:,k),1);
    end
end

% Shift data if not CSI data
if (isempty(strfind(parms.seqfil,'csi')))
    if (ndims(data) < 4)
        data = peshift(data, procparFile);
    elseif (ndims(data) == 4)
        for j = 1:size(data,4)
            data(:,:,:,j) = peshift(data(:,:,:,j), procparFile);
        end
    end
else
    if (ndims(data)==3)
        % Reformat CSI data
        data = permute(data,[2 3 1]);
        data = flipdim(data,3);
        data = flipdim(data,1);
        spc = [spc(3) spc(1) spc(2)];
        dims = [dims(3) dims(1) dims(2)];
        origin = [pp.ppe -pp.ppe2 0] - spc.*dims./2;
        origin = origin + spc./2;
    elseif (ndims(data)==4)
        % Reformat CSI data
        data = permute(data,[2 3 4 1]);
        data = flipdim(data,4);
        data = flipdim(data,1);
        spc = [spc(3) spc(1) pp.ns*pp.thk/(10*dims(4)) spc(2)];
        dims = [dims(3) dims(1) dims(4) dims(2)];
        origin = [pp.ppe -pp.ppe2 pp.pss0 0] - spc.*dims./2;
        origin = origin + spc./2;
    end
        
end

% Transform to image space
if (isfield(parms,'seqfil')&&transform)
    if (strfind(parms.seqfil,'ge3d'))
        data = fftshift(abs(fftn(data)));
    elseif (strfind(parms.seqfil,'mgems'))
        data = fftshift(abs(fft2(data)),1);
        data = fftshift(data,2);
    elseif (strfind(parms.seqfil,'gems'))
        data = fftshift(abs(fft2(data)),1);
        data = fftshift(data,2);
    elseif (strfind(parms.seqfil,'sems'))
        data = fftshift(abs(fft2(data)),1);
        data = fftshift(data,2);  
    end
end


data = vuRowMajorColumnMajorSwitch(data);


fclose(fid);

return;

function [hdr,be] = read_hdr_raw(fname)
% Adapted in vuOpenImage
% By: Kevin Wilson
% Date: Sept 20, 2007
%
% Read a NIFTI-1 hdr file
% FORMAT [hdr,be] = read_hdr_raw(fname)
% fname - filename of image
% hdr   - a structure containing hdr info
% be    - whether big-endian or not
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: read_hdr_raw.m 253 2005-10-13 15:31:34Z guillaume $


hdr = [];
be  = [];
ok  = true;

% Name of header file
[pth,nam,ext] = fileparts(fname);
switch ext
case {'.img','.hdr'}
    hname = fullfile(pth,[nam '.hdr']);
case {'.nii'}
    hname = fullfile(pth,[nam '.nii']);
otherwise
    hname = fullfile(pth,[nam '.hdr']);
end;

% Open it if possible
fp  = fopen(hname,'r','native');
if fp==-1
    hdr = [];
    return;
end;

% Sort out endian issues of header
[unused,unused,mach] = fopen(fp);
if strcmp(mach,'ieee-be')
    be = true;
else
    be = false;
end;
fseek(fp,0,'bof');
fseek(fp,40,'bof');
nd = fread(fp,1,'int16')';
if isempty(nd),
    fclose(fp);
    hdr = [];
    return;
elseif nd<1 || nd>7 %#ok<BDSCI>
    be = ~be;
    fclose(fp);
    if be, mach = 'ieee-be';
    else   mach = 'ieee-le';
    end;
    fp = fopen(hname,'r',mach);
    if fp==-1
        hdr = [];
        return;
    end;
end;

% Is it NIFTI or not
fseek(fp,0,'bof');
fseek(fp,344,'bof');
mgc = deblank(char(fread(fp,4,'uint8')'));
switch mgc
case {'ni1','n+1'}
    org = niftistruc;
otherwise
    org = mayostruc;
end;
fseek(fp,0,'bof');
% Read the fields
for i=1:length(org)
    tmp = fread(fp,org(i).len,['*' org(i).dtype.prec])';
    if length(tmp) ~= org(i).len
disp([length(tmp) org(i).len]);
        tmp = org(i).def;
        ok  = false;
    end;
    tmp = feval(org(i).dtype.conv,tmp);
    hdr.(org(i).label) = tmp;
end;

fclose(fp);
if ~ok,
     fprintf('There was a problem reading the header of\n');
     fprintf('"%s".\n', fname);
     fprintf('It may be corrupted in some way.');
end;
return;

function o = niftistruc
% Create a data structure describing NIFTI headers
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: niftistruc.m 253 2005-10-13 15:31:34Z guillaume $


persistent org;
if ~isempty(org),
    o = org;
    return;
end;
t = struct('conv',{ @char , @uint8 , @int16 , @int32 , @single },...
           'prec',{'uint8', 'uint8', 'int16', 'int32', 'single'},...
           'size',{       1,      1,      2,      4,       4 });
c = t(1);
b = t(2);
s = t(3);
i = t(4);
f = t(5);

table = {...
    i, 1,'sizeof_hdr',348
    c,10,'data_type',[]
    c,18,'db_name',[]
    i, 1,'extents',[]
    s, 1,'session_error',[]
    c, 1,'regular','r'
    b, 1,'dim_info',[]
    s, 8,'dim',[3 0 0 0  1 1 1 1 1]
    f, 1,'intent_p1',0
    f, 1,'intent_p2',0
    f, 1,'intent_p3',0
    s, 1,'intent_code',0
    s, 1,'datatype',2
    s, 1,'bitpix',8
    s, 1,'slice_start',[]
    f, 8,'pixdim',[0 1 1 1]
    f, 1,'vox_offset',0
    f, 1,'scl_slope',1
    f, 1,'scl_inter',0
    s, 1,'slice_end',[]
    b, 1,'slice_code',[]
    b, 1,'xyzt_units',10
    f, 1,'cal_max',[]
    f, 1,'cal_min',[]
    f, 1,'slice_duration',[]
    f, 1,'toffset',[]
    i, 1,'glmax',[]
    i, 1,'glmin',[]
    c,80,'descrip','NIFTI-1 Image'
    c,24,'aux_file',''
    s, 1,'qform_code',0
    s, 1,'sform_code',0
    f, 1,'quatern_b',0
    f, 1,'quatern_c',0
    f, 1,'quatern_d',0
    f, 1,'qoffset_x',0
    f, 1,'qoffset_y',0
    f, 1,'qoffset_z',0
    f, 4,'srow_x',[1 0 0 0]
    f, 4,'srow_y',[0 1 0 0]
    f, 4,'srow_z',[0 0 1 0]
    c,16,'intent_name',''
    c, 4,'magic','ni1'};
org = struct('label',table(:,3),'dtype',table(:,1),'len',table(:,2),...
    'offset',0,'def',table(:,4));
os  = 0;
for j=1:length(org)
    os  = os + org(j).dtype.size*ceil(os/org(j).dtype.size);
    fun = org(j).dtype.conv;
    def = [org(j).def zeros(1,org(j).len-length(org(j).def))];
    org(j).def    = feval(fun,def);
    org(j).offset = os;
    os  = os + org(j).len*org(j).dtype.size;
end;
o = org;
return;

function o = mayostruc
% Create a data structure describing Analyze headers
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: mayostruc.m 253 2005-10-13 15:31:34Z guillaume $


persistent org;
if ~isempty(org),
    o = org;
    return;
end;
t = struct('conv',{ @char, @int16, @int32, @single },...
           'prec',{'uint8','int16','int32','single'},...
           'size',{      1,      2,      4,       4});
c = t(1);
s = t(2);
i = t(3);
f = t(4);
table = {...
    i, 1,'sizeof_hdr',348
    c,10,'data_type',[]
    c,18,'db_name',[]
    i, 1,'extents',[]
    s, 1,'session_error',[]
    c, 1,'regular','r'
    c, 1,'hkey_un0',[]
    s, 8,'dim',[3 1 1 1  1 1 1 1 1]
    c, 4,'vox_units',[]
    c, 8,'cal_units',[]
    s, 1,'unused1',[]
    s, 1,'datatype',[]
    s, 1,'bitpix',[]
    s, 1,'dim_un0',[]
    f, 8,'pixdim',[]
    f, 1,'vox_offset',0
    f, 1,'roi_scale',1
    f, 1,'funused1',0
    f, 1,'funused2',[]
    f, 1,'cal_max',[]
    f, 1,'cal_min',[]
    i, 1,'compressed',[]
    i, 1,'verified',[]
    i, 1,'glmax',[]
    i, 1,'glmin',[]
    c,80,'descrip','Analyze Image'
    c,24,'aux_file',''
    c, 1,'orient',[]
%    c,10,'originator',[]
    s, 5,'origin',[] % SPM version
    c,10,'generated',[]
    c,10,'scannum',[]
    c,10,'patient_id',[]
    c,10,'exp_date',[]
    c,10,'exp_time',[]
    c, 3,'hist_un0',[]
    i, 1,'views',[]
    i, 1,'vols_added',[]
    i, 1,'start_field',[]
    i, 1,'field_skip',[]
    i, 1,'omax',[]
    i, 1,'omin',[]
    i, 1,'smax',[]
    i, 1,'smin',[]};
org = struct('label',table(:,3),'dtype',table(:,1),'len',table(:,2),...
    'offset',0,'def',table(:,4));
os  = 0;
for j=1:length(org)
    os  = os + org(j).dtype.size*ceil(os/org(j).dtype.size);
    fun = org(j).dtype.conv;
    def = [org(j).def zeros(1,org(j).len-length(org(j).def))];
    org(j).def    = feval(fun,def);
    org(j).offset = os;
    os  = os + org(j).len*org(j).dtype.size;
end;
o = org;
return;

function M = Q2M(Q)
% Adapted in vuOpenImage
% By: Kevin Wilson
% Date: Sept 20, 2007
%
% Generate a rotation matrix from a quaternion xi+yj+zk+w,
% where Q = [x y z], and w = 1-x^2-y^2-z^2.
% See: http://skal.planet-d.net/demo/matrixfaq.htm
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: Q2M.m 253 2005-10-13 15:31:34Z guillaume $


Q = Q(1:3); % Assume rigid body
w = sqrt(1 - sum(Q.^2));
x = Q(1); y = Q(2); z = Q(3);
if w<1e-7,
    w = 1/sqrt(x*x+y*y+z*z);
    x = x*w;
    y = y*w;
    z = z*w;
    w = 0;
end;
xx = x*x; yy = y*y; zz = z*z; ww = w*w;
xy = x*y; xz = x*z; xw = x*w;
yz = y*z; yw = y*w; zw = z*w;
M = [...
(xx-yy-zz+ww)      2*(xy-zw)      2*(xz+yw) 0
    2*(xy+zw) (-xx+yy-zz+ww)      2*(yz-xw) 0
    2*(xz-yw)      2*(yz+xw) (-xx-yy+zz+ww) 0
           0              0              0  1];
return;

function data = readDAT(fullFilename)
% File reader for .dat files from Neuro Scan
% By: Kevin Wilson
% Date: October 12, 2007
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

fid = fopen(fullFilename,'rt');
if (fid == -1)
    disp('.DAT image open unsuccessful')
    data = -1;
    return;
end

% Get line until end of header
Parms = [];
line = fgetl(fid);

% Until Electrodes
while (isempty(findstr(line,'[Electrode Labels]')))
    % Find open and close brackets
    openB = findstr(line,'[');
    closeB = findstr(line,']');
    
    % Store parameters
    curParm = line(openB+1:closeB-1);
    
    % Get Value
    if (strcmp(curParm,'Channels')||strcmp(curParm,'Rate')||strcmp(curParm,'Points')||strcmp(curParm,'Xmin')|| ...
            strcmp(curParm,'Sweeps')||strcmp(curParm,'Accepted')||strcmp(curParm,'Rejected'))
        curValue = sscanf(line(closeB+1:end),'%f');
    else
        curValue = sscanf(line(closeB+1:end),'%s');
    end
    
    % Store
    Parms.(curParm) = curValue;
    
    % Next line
    line = fgetl(fid);
end

% Labels, XUnits, and YUnits
for ii = 1:3
    % Find open and close brackets
    openB = findstr(line,'[');
    closeB = findstr(line,']');
    
    % Store parameter name
    curParm = line(openB+1:closeB-1);
    curParm = curParm(~isspace(curParm));
    % Next line
    line = fgetl(fid);
    
    % Find open and close brackets
    openB = findstr(line,'[');
    closeB = findstr(line,']');
    
    % Check the lengths
    if (length(openB)~=length(closeB)||length(openB)~=Parms.Channels)
        error('MATLAB:vuOpenImage:DATRead','Number of Channels in header does not match actual number of channels');
    end
    
    % Get channel names
    for j = 1:Parms.Channels
        channel{j} = line(openB(j)+1:closeB(j)-1);
    end
    
    % Store Parm
    Parms.(curParm) = channel;
    
    % Next line
    line = fgetl(fid);
end

% Store containers
eTrialType = zeros(1,Parms.Sweeps);
eAccept = zeros(1,Parms.Sweeps);
eCorrect = zeros(1,Parms.Sweeps);
eRT = zeros(1,Parms.Sweeps);
eResp = zeros(1,Parms.Sweeps);
eData = zeros(Parms.Points,Parms.Channels,Parms.Sweeps);
count = 1;

% Read data
while (~feof(fid)||isempty(line))
    % Check
    if (isempty(findstr(line,'[Epoch Header]')))
        error('MATLAB:vuOpenImage:DATRead','Corrupt Epoch Header');
    else
        % Next line
        line = fgetl(fid);
        
        % Get Epoch Data Header
        while (isempty(findstr(line,'[Epoch Data]')))
            % Find open and close brackets
            openB = findstr(line,'[');
            closeB = findstr(line,']');

            if (~isempty(findstr(line,'[Trial Type]')))
                eTrialType(count) = sscanf(line(closeB+1:end),'%f');
            elseif (~isempty(findstr(line,'[Accept]')))
                eAccept(count) = sscanf(line(closeB+1:end),'%f');
            elseif (~isempty(findstr(line,'[Correct]')))
                eCorrect(count) = sscanf(line(closeB+1:end),'%f');
            elseif (~isempty(findstr(line,'[RT]')))
                eRT(count) = sscanf(line(closeB+1:end),'%f');
            elseif (~isempty(findstr(line,'[Resp]')))
                eResp(count) = sscanf(line(closeB+1:end),'%f');
            end

            % Next line
            line = fgetl(fid);
        end

        % Next line
        line = fgetl(fid);
        
        % Complex Data
        if (~isempty(findstr(line,'[Real Data]')))

            % Get Real data
            tempData = textscan(fid,'%f',(Parms.Points)*Parms.Channels);
            eData(:,:,count) = reshape(tempData{1},[(Parms.Points),Parms.Channels]);
            
            % Next line
            line = fgetl(fid);
            line = fgetl(fid);
            
            if (isempty(findstr(line,'[Imaginary Data]')))
                error('MATLAB:vuOpenImage:DATRead','Corrupt Epoch Data : No Imaginary Data Header');
            end

            % Get Imaginary data
            xxx = textscan(fid,'%f',(Parms.Points)*Parms.Channels);
            eData(:,:,count) = eData(:,:,count) + i*reshape(xxx{1},[(Parms.Points),Parms.Channels]);
            
        % Real Only 
        else
            % Get Data
            eData(1,:,count) = str2num(line);
            tempData = textscan(fid,'%f',(Parms.Points-1)*Parms.Channels);
            eData(2:end,:,count) = reshape(tempData{1},[(Parms.Points-1),Parms.Channels]);
        end
        
        % Next line
        line = fgetl(fid);
        line = fgetl(fid);
        count = count + 1;
    end
end

% Put together Data
data.Data = eData;
data.TrialType = eTrialType;
data.Accept = eAccept;
data.Correct = eCorrect;
data.RT = eRT;
data.Resp = eResp;
data.Parms = Parms;

fclose(fid);

return;

function xform = CalculateRotationMatrix(phi,theta,psi)
% Function to calculate a rotation matrix from euler angles
B = eye(3);
B(1:2,1:2) = [cosd(psi(1)) sind(psi(1));-sind(psi(1)) cosd(psi(1))];

C = eye(3);
C(2:3,2:3) = [cosd(theta(1)) sind(theta(1));-sind(theta(1)) cosd(theta(1))];

D = eye(3);
D(1:2,1:2) = [cosd(phi(1)) sind(phi(1));-sind(phi(1)) cosd(phi(1))];

xform = B*C*D;

function k=peshift(k, procpar)
%PESHIFT  Shifts the FOV of an image in the phase encode direction.
%
%   k=peshift(k, procpar); returns the linearly phase-ramped
%   k-space data that, when inverse Fourier transformed, will yield
%   the image shifted in the phase encode direction.
%
%   k is the k-space data for a stack of images.
%   procpar is the full path to the Varian procpar file.  This file
%      will be read by the embedded function parsepp, and the parameters
%      ppe and lpe will be used to calculate the phase encode increment.
%      If ppe is zero, the user will be prompted for the desired value.
%


%   VERSION 1.0, Released 16 April, 2007.  J. Luci
%   VERSION 1.1, Released 17 April, 2007.  J. Luci (multi-page data
%      supported)


%Read in the parameters from the procpar file
par=parsepp(procpar);

%Calculate the linear phase increment to affect the desired off-center FOV
phiIncrement = 2*pi*par.ppe*-1/par.lpe;
phiIncrement2 = 2*pi*par.ppe2*-1/par.lpe2;

if (par.ppe~=0)
    %construct the phase ramp for ppe
    shift=ones(size(k,1), size(k,2));
    for ii = 1:size(k,2),
        phi=ii*phiIncrement;
        shift(:,ii)=exp(i*phi)*shift(:,ii);
    end
    %Apply the phase correction
    for ii = 1:size(k,3)
        k(:,:,ii) = shift.*k(:,:,ii);
    end
end

if (par.ppe2~=0&&size(k,3)>1)
    %construct the phase ramp for ppe2
    shift2=ones(size(k,2), size(k,3));
    for ii = 1:size(k,3),
        phi2=ii*phiIncrement2;
        shift2(:,ii)=exp(i*phi2)*shift2(:,ii);
    end

    for ii = 1:size(k,1)
        k(ii,:,:) = shift2.*squeeze(k(ii,:,:));
    end
end