function vuWriteImage(X,filename,varargin)
% vuImageOpen allows users to select and open image files
%
%   SYNTAX:
%       vuWriteImage(X);
%       vuWriteImage(X,filename,options)
%
%   PARAMETERS (REQUIRED):
%       X : Input Image
%
%   OPTIONS & DEFAUTLS:
%       filename : dialog opens
%       datatype (Nifti only) : Percision of data to be written out
%               'int8','uint8','int16','uint16','int32','uint32','float32'
%       map = gray(256) (others include: hot, cool, bone, spring, summer,
%                   autumn, winter, hsv, jet, copper, pink)
%       Note: map only for Regular images.
%
%   OUTPUT:
%       X is the image to be saved
%
%   EXAMPLE:
%       im = imread('cameraman.tif');
%       vuWriteImage(im,'Jcameraman.jpg');
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 1
  warning('MATLAB:vuWriteImage:NotEnoughInputs', 'Not enough input arguments.');
end

% Get and parse options
p = inputParser;
p.addParamValue('map',gray(256),@(x) isa(x,'double'));
p.addParamValue('datatype','float32', @ischar);
p.FunctionName='vuWriteImage';
p.parse(varargin{:});

if nargin < 2
    [filename, pathname, filterindex] = uiputfile({'*.par *.PAR *.rec *.REC *.mha;*.MHA;*.nii;*.NII;*.nii.gz;*.NII.GZ','Medical Image Files'; ...
    '*.mat;*.MAT;*.dat;*.DAT','Matlab Data';'*.bmp;*.BMP;*.cur;*.CUR;*.gif;*.GIF;*.hdf;*.HDF;*.ico;*.ICP;*.jpg;*.JPG;*.jpeg;*.JPEG;*.pbm;*.PBM;*.pcx;*.PCX;*.pgm;*.PGM;*.png;*.PNG;*.pnm;*.PNM;*.ppm;*.PPM;*.ras;*.RAS;*.tif;*.TIF;*.tiff;*.TIFF;*.xwd;*.XWD', 'Image Files'; ...
    '*.*', 'All Files (*.*)'},'Save file as ...');

    % Check if file was selected
    if (~filterindex)
        return
    end
    
    % full file name
    fullFilename = fullfile(pathname,filename);
    
else
    pathname = [];
    fullFilename = filename;
end

% Get extension
len = length(filename);
ext = filename(len-3:len);

% Find file type
isMed = regexpi('.par .rec .mha .nii i.gz .hdr',ext);
isMatlab = regexpi('.mat .dat',ext);
isImage = regexpi('.bmp .cur .gif .hdf .ico .jpg .jpe .pbm .pcx .pgm .png .pnm .ppm .ras .tif .xwd',ext);


if (isMed)
    disp('Medical Image Write ...')
    if isstruct(X)
       if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
        disp('X is a valid MetaImage struct.')
       else
          error('MATLAB:vuWriteImage:InvalidStruct', 'The input image structure is not valid.');
       end
    else
        X = vuGenerateMetaImage(X);
    end
    % Row-major switch
    X = vuRowMajorColumnMajorSwitch(X);
    if (regexpi('.par .rec',ext))
        disp('PAR/REC write')
        writeParRec(X,pathname,filename);
    elseif (regexpi('.nii',ext))
        disp('NII write')
        writeNifti(X,fullFilename,p);
    elseif (~isempty(regexpi('.hdr',ext)))
        if (ndims(X.Data) ==2)
            MEX2DITKImageFileWriter(fullFilename,X)
        elseif (ndims(X.Data) == 3)
            if (~isfield(X,'Orientation'))
                X.Orientation = eye(3);
            end
            MEX3DITKImageFileWriter(fullFilename,X)
        else
            error('MATLAB:vuWriteImage:UnknownDimensions','Unsupported image dimensions.');
        end
    else
        if (ndims(X.Data) ==2)
            MEX2DITKImageFileWriter(fullFilename,X)
        elseif (ndims(X.Data) == 3)
            if (~isfield(X,'Orientation'))
                X.Orientation = eye(3);
            end
            MEX3DITKImageFileWriter(fullFilename,X)
        elseif (ndims(X.Data) == 4)
            if (~isfield(X,'Orientation'))
                X.Orientation = eye(3);
            end
            MEX4DITKImageFileWriter(fullFilename,X)
        else
            error('MATLAB:vuWriteImage:UnknownDimensions','Unsupported image dimensions.');
        end
    end
elseif (isMatlab)
    disp('Matlab Data Write (Column Major) ...')
    if isstruct(X)
        if(isfield(X,'Data'))
            tmp = X.Data;
            save fullFilename tmp
        else
            error('MATLAB:vuWriteImage:InvalidStucture','Invalid Image Structure.');
        end
    else
        save fullFilename X
    end
elseif (isImage)
    disp('Image Data Write ...')
    if isstruct(X)
        if(isfield(X,'Data'))
            imwrite(double(X.Data),p.Results.map,fullFilename,ext(2:4));
        else
            error('MATLAB:vuWriteImage:InvalidStucture','Invalid Image Structure.');
        end
    else
        imwrite(double(X),p.Results.map,fullFilename,ext(2:4));
    end
else
    disp('Binary Data Write (Row-Major) ...')
    if isstruct(X)
        if(isfield(X,'Data'))
            if(ndims(X.Data)==2)
                fwrite(fopen(fullFilename,'w'),permute(X.Data,[2 1]),'int16');
            elseif (ndims(X.Data)==3)
                fwrite(fopen(fullFilename,'w'),permute(X.Data,[2 1 3]),'int16');
            elseif (ndims(X.Data)==4)
                fwrite(fopen(fullFilename,'w'),permute(X.Data,[2 1 3 4]),'int16');
            else
                error('MATLAB:vuWriteImage:UnknownDimensions','Unknown Dimensions.');
            end
        else
            error('MATLAB:vuWriteImage:InvalidImageFormat','Invalid Image Format.');
        end
    else
        if(ndims(X)==2)
            fwrite(fopen(fullFilename,'w'),permute(X,[2 1]),'int16');
        elseif (ndims(X)==3)
            fwrite(fopen(fullFilename,'w'),permute(X,[2 1 3]),'int16');
        elseif (ndims(X)==4)
            fwrite(fopen(fullFilename,'w'),permute(X,[2 1 3 4]),'int16');
        else
            error('MATLAB:vuWriteImage:UnknownDimensions','Unknown Dimensions.');
        end
    end
end

function writeParRec(X,pathname,filename)
% By: Kevin Wilson
% Date: April 26, 2007
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if (~isstruct(X))
    warning ('MATLAB:vuWriteImage:Default','Generating default spacing and origin');
    X = vuGenerateMetaImage(X);
end

try
    par_file = X.Parms.filename;
catch
    par_file = [];
    warning ('MATLAB:vuWriteImage:ParFileNoDefined','The Original Par file was not defined.  No Par file generated');
end

if(~exist(par_file,'file'))
    copyPar = 0;
    warning('MATLAB:vuWriteImage:ParFileNoFound','The Original Par file was not found.  No Par file generated');
else
    copyPar = 1;
end

% Default
write_type = 'short';

% Block to determine whether to copy PAR file
try
    % Change spacing and origin back to mms
    X.Spc = X.Spc;
    X.Origin = X.Origin;
    
    % Version 3
    if (strcmp(X.Parms.version,'V3'))
        pixels = X.Parms.pixel_bits;
        sliceOri = X.Parms.tags(1,20);
    % Version 4 
    elseif (strcmp(X.Parms.version,'V4')||strcmp(X.Parms.version,'V4.1'))
        pixels = X.Parms.tags(1,8);
        sliceOri = X.Parms.tags(1,26);
    else
        copyPar = 0;
        warning('MATLAB:vuWriteImage:UnsupportedVersion','The Par File Version is Unsupported. No Par file Generated.');
    end

    switch (pixels)
        case { 8 }, write_type = 'int8';
        case { 16 }, write_type = 'short';
        otherwise, write_type = 'uchar';
    end
    
    % Get variables to compare whether im dimensions, spacing, or origin
    % has changed
    switch (sliceOri)
        case {1}
            dims = [X.Parms.fov(3) X.Parms.fov(1) X.Parms.fov(2)]./X.Spc;
            origin = vuCalculateImageCenter(X);
            origin = [origin(2) origin(3) origin(1)];
        case {2}
            dims = [X.Parms.fov(1) X.Parms.fov(2) X.Parms.fov(3)]./X.Spc;
            origin = vuCalculateImageCenter(X);
            origin = origin.*[1 -1 -1];
        case {3}
            dims = [X.Parms.fov(3) X.Parms.fov(2) X.Parms.fov(1)]./X.Spc;
            origin = vuCalculateImageCenter(X);
            origin = [origin(3) origin(2) origin(1)];
            origin = origin.*[1 -1 1];
        otherwise
            copyPar = 0;
            warning('MATLAB:vuWriteImage:UnknownOrientation','The Original Par File orientation is not reconized.  No Par file generated');
    end
    dims = round(dims);
    if (max(dims~=X.Dims(1:length(dims))))
        copyPar = 0;
        warning('MATLAB:vuWriteImage:ChangedDimensions','The dimensions/spacing of the image has changed for the original PAR file.  No Par file generated');
    % Allow some room for round off error
    elseif (max((origin - X.Parms.offcenter)>0.001))
        copyPar = 0;
        warning('MATLAB:vuWriteImage:ChangedOrigin','The origin of the image has changed for the original PAR file.  No Par file generated');
    end

catch
    copyPar = 0;
    warning('MATLAB:vuWriteImage:IncorrectTags','The image does not has sufficient information for PAR tags.  No PAR file generated');
end

% Copy Par
if (copyPar)
    try
        newPar = fullfile(pathname,filename);
        newPar = strcat(newPar(1:length(newPar)-3),'PAR');
        s=copyfile(par_file,newPar,'f');
        if (s==0)
            warning('MATLAB:vuWriteImage:ParFileError','The original PAR file could not be copied.  No PAR file generated.');
        end
    catch
        warning('MATLAB:vuWriteImage:ParFileError','The original PAR file could not be copied.  No PAR file generated.');
    end
end

% Write Rec
try
    newRec = fullfile(pathname,filename);
    newRec = strcat(newRec(1:length(newRec)-3),'REC');
    fid = fopen(newRec,'w','l');
    
    parms = X.Parms;
    types_list = unique(parms.tags(:,5)'); % to handle multiple types
    slice_list = unique(parms.tags(:,1)');
    phase_list = unique(parms.tags(:,4)');
    echo_list = unique(parms.tags(:,2)');
    dynamic_list = unique(parms.tags(:,3)');
    if(strcmp(parms.version,'V4.1'))
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
    if(strcmp(parms.version,'V4.1'))
        ndiffb = size(diffb_list,2);
        ngrad = size(grad_list,2);
    else
        ndiffb = 1;
        ngrad = 1;
    end
    ntype = size(types_list,2);
    
    if ( isfield(parms,'recon_resolution') )
        nline = parms.recon_resolution(1);
        stride = parms.recon_resolution(2);
    else
        nline = parms.tags(1,10);
        stride = parms.tags(1,11);
    end

    switch(parms.version)
    case {'V3'}, pixel_bits = parms.pixel_bits;
    case {'V4'}, pixel_bits = parms.tags(1,8);  % assume same for all imgs
    case {'V4.1'}, pixel_bits = parms.tags(1,8);  % assume same for all imgs
    end

    switch (pixel_bits)
        case { 8 }, write_type = 'int8';
        case { 16 },write_type = 'short';
        otherwise, write_type = 'uchar';
    end

    X.Data = reshape(X.Data,[stride nline nslice necho nphase ntype ndyn ngrad ndiffb]);
    
    switch( parms.version )
        case { 'V3' }, ri_i = 8; rs_i = 9; ss_i = 10;
        case { 'V4' }, ri_i = 12; rs_i = 13; ss_i = 14;
        case { 'V4.1' }, ri_i = 12; rs_i = 13; ss_i = 14;
    end;
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
        if(strcmp(parms.version,'V4.1'))
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
        img = X.Data(:,:,slice_idx,echo_idx,phase_idx,type_idx,dyn_idx,grad_idx,diffb_idx);
        % rescale data to produce FP information (not SV, not DV)
        RI = parms.tags(I,ri_i);  % 'inter' --> 'RI'
        RS = parms.tags(I,rs_i);  % 'slope' --> 'RS'
        SS = parms.tags(I,ss_i);  % new var 'SS'
        img = (img*(RS*SS)-RI)./RS;
        fwrite(fid,img,write_type);
    end
    
catch
    warning('MATLAB:vuWriteImage:RecFileScale','Unsufficient Data to rescale REC Data.  Given values written');
    
    try
        newRec = fullfile(pathname,filename);
        newRec = strcat(newRec(1:length(newRec)-3),'REC');
        fid = fopen(newRec,'w','l');
        fwrite(fid,X.Data,write_type);
    catch
        error('MATLAB:vuWriteImage:RecFileError','The Rec file could not be written.  No Rec file generated.');
    end
    
end


return;

function writeNifti(IM,fullFilename,p)
% Adapted in vuWriteImage
% By: Jeff Luci and Kevin Wilson
% Date: June 16, 2008

% Create a niftiHeader
% Check if one already exists by Checking for some known nifti parms
if (isfield(IM,'Parms')&&isfield(IM.Parms,'srow_x')&&isfield(IM.Parms,'aux_file')&&isfield(IM.Parms,'scl_slope'))
    % Assume niftiheader
    niftiHeader = IM.Parms;
else
    % Create one
    niftiHeader.sizeof_hdr = 348; %sizeof_hdr
    niftiHeader.data_type = [1 2 3 4 5 6 7 8 9 0]; %UNUSED data_type
    niftiHeader.db_name = [1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8]; %UNUSED db_name
    niftiHeader.extents = 0; %UNUSED extents
    niftiHeader.session_error = 0; %UNUSED session_error
    niftiHeader.regular = 0; %UNUSED regular
    dim_info=uint8(0);
    dim_info=bitset(dim_info,1);
    dim_info=bitset(dim_info,4);
    dim_info=bitset(dim_info,5);
    dim_info=bitset(dim_info,6);
    niftiHeader.dim_info = dim_info; %dim_info
    dim = ones(1,8);
    dim(1:length(IM.Dims)+1) = [ndims(IM.Data) IM.Dims];
    niftiHeader.dim = dim; %dim
    niftiHeader.intent_p1 = 0; %intent parameter 1 ***
    niftiHeader.intent_p2 = 0; %intent parameter 2 ***
    niftiHeader.intent_p3 = 0; %intent parameter 3 ***
    niftiHeader.intent_code = 0; %NIFIT intent code ***
    if (strcmp(p.Results.datatype,'int8'))
        datatypeN = 256;
    elseif (strcmp(p.Results.datatype,'uint16'))
        datatypeN = 512;
    elseif (strcmp(p.Results.datatype,'uint32'))
        datatypeN = 768;
    elseif (strcmp(p.Results.datatype,'int64'))
        datatypeN = 1024;
    elseif (strcmp(p.Results.datatype,'uint64'))
        datatypeN = 1280;
    elseif (strcmp(p.Results.datatype,'uint8'))
        datatypeN = 2;
    elseif (strcmp(p.Results.datatype,'int16'))
        datatypeN = 4;
    elseif (strcmp(p.Results.datatype,'int32'))
        datatypeN = 8;
    elseif (strcmp(p.Results.datatype,'float32'))
        datatypeN = 16;
    else
        error('Data type not supported');
    end 
    niftiHeader.datatype = datatypeN; %datatype
    niftiHeader.bitpix = 32; %bitpix
    niftiHeader.slice_start = 1; %first slice number
    pix_dim = ones(1,8);
    if (~isfield(IM,'Orientation'))
        IM.Orientation = eye(3);
    end
    if (det(IM.Orientation)<1)
        pix_dim(1) = -1;
        IM.Orientation = IM.Orientation*[1 0 0;0 1 0;0 0 -1];
    end
    pix_dim(2:length(IM.Spc)+1) = IM.Spc;
    niftiHeader.pix_dim = pix_dim; %grid spacings ***
    niftiHeader.vox_offset = 348; %vox_offset
    niftiHeader.scl_slope = 1; %scl_slope
    niftiHeader.scl_inter = 0; %scl_inter
    try
        numSlice = IM.Dims(3);
    catch
        numSlice = 1;
    end
    niftiHeader.slice_end = numSlice; %last slice number
    niftiHeader.slice_code = 1; %slice timing code *** WHY CHAR?
    xyzt_units=uint8(0);
    xyzt_units=bitset(xyzt_units,2);
    xyzt_units=bitset(xyzt_units,4);
    niftiHeader.xyzt_units = xyzt_units; %xyzt_units
    niftiHeader.cal_max = max(IM.Data(:)); %display max
    niftiHeader.cal_min = min(IM.Data(:)); %display min
    if (isfield(IM,'Parms')&&isfield(IM.Parms,'repetition_time'))
        niftiHeader.slice_duration = IM.Parms.repetition_time; %slice time
    else
        niftiHeader.slice_duration = 1; %slice time
    end
    niftiHeader.toffset = 0; %toffset
    niftiHeader.glmax = 0; %unused
    niftiHeader.glmin = 0; %unused
    descrip = 'Written using vuTools';
    for ii = size(descrip,2)+1:80
        descrip = [descrip ' '];
    end
    niftiHeader.descrip = descrip; %description, use Varian seqfil
    descrip = 'niftiFilename';
    for ii = size(descrip,2)+1:24
        descrip = [descrip ' '];
    end
    niftiHeader.aux_file = descrip; %alternate filename
    niftiHeader.qform_code = 1; %qform code ***
    niftiHeader.sform_code = 0; %sform code ***
    a = 1 + IM.Orientation(1,1) + IM.Orientation(2,2) + IM.Orientation(3,3);
    if (a>0.5)
        a = 0.5*sqrt(a);
        b = 0.25*(IM.Orientation(3,2)-IM.Orientation(2,3)) / a;
        c = 0.25*(IM.Orientation(1,3)-IM.Orientation(3,1)) / a;
        d = 0.25*(IM.Orientation(2,1)-IM.Orientation(1,2)) / a;
    else
        xd = 1 + IM.Orientation(1,1) - (IM.Orientation(2,2)+IM.Orientation(3,3));
        yd = 1 + IM.Orientation(2,2) - (IM.Orientation(1,1)+IM.Orientation(3,3));
        zd = 1 + IM.Orientation(3,3) - (IM.Orientation(1,1)+IM.Orientation(2,2));
        if xd > 1
            b = 0.5 * sqrt(xd);
            c = 0.25 * (IM.Orientation(1,2)+IM.Orientation(2,1)) / b;
            d = 0.25 * (IM.Orientation(1,3)+IM.Orientation(3,1)) / b;
            a = 0.25 * (IM.Orientation(3,2)-IM.Orientation(2,3)) / b;
        elseif yd > 1
            c = 0.5 * sqrt(yd);
            b = 0.25 * (IM.Orientation(1,2)+IM.Orientation(2,1)) / c;
            d = 0.25 * (IM.Orientation(2,3)+IM.Orientation(3,2)) / c;
            a = 0.25 * (IM.Orientation(1,3)-IM.Orientation(3,1)) / c;
        else
            d = 0.5 * sqrt(zd);
            b = 0.25 * (IM.Orientation(1,3)+IM.Orientation(3,1)) / d;
            c = 0.25 * (IM.Orientation(2,3)+IM.Orientation(3,2)) / d;
            a = 0.25 * (IM.Orientation(2,1)-IM.Orientation(1,2)) / d;
        end
        if a < 0
            b = -b;
            c = -c;
            d = -d;
            a = -a;
        end
    end
    niftiHeader.quartern_b = b; %Quarternion b param ***
    niftiHeader.quartern_c = c; %Quarternion c param ***
    niftiHeader.quartern_d = d; %Quarternion d param ***
    niftiHeader.qoffset_x = IM.Origin(1); %Quarternion x param ***
    niftiHeader.qoffset_y = IM.Origin(2); %Quarternion y param ***
    try
        niftiHeader.qoffset_z = IM.Origin(3); %Quarternion z param ***
    catch
        niftiHeader.qoffset_z = 0; %Quarternion z param ***
    end
    srow_x = [IM.Orientation(1,1:3) IM.Origin(1)];
    srow_y = [IM.Orientation(2,1:3) IM.Origin(2)];
    srow_z = [IM.Orientation(3,1:3) IM.Origin(3)];
    niftiHeader.srow_x = srow_x; %1st row affine transform
    niftiHeader.srow_y = srow_y; %2nd row affine transform
    niftiHeader.srow_z = srow_z; %3rd row affine transform
    descrip = 'Nifti File';
    for ii = size(descrip,2)+1:16
        descrip = [descrip ' '];
    end
    niftiHeader.intent_name = descrip; %name of data
    niftiHeader.magic = [uint8('n+1') 0]; %Magic string
end


% Added to write out the data
% Get data type
if isfield(niftiHeader,'datatype')
    if (niftiHeader.datatype==256)
        datatype = 'int8';
    elseif (niftiHeader.datatype==512)
        datatype = 'uint16';
    elseif (niftiHeader.datatype==768)
        datatype = 'uint32';
    elseif (niftiHeader.datatype==1024)
        datatype = 'int64';
    elseif (niftiHeader.datatype==1280)
        datatype = 'uint64';
    elseif (niftiHeader.datatype==2)
        datatype = 'uint8';
    elseif (niftiHeader.datatype==4)
        datatype = 'int16';
    elseif (niftiHeader.datatype==8)
        datatype = 'int32';
    elseif (niftiHeader.datatype==16)
        datatype = 'float32';
    else
        error('Data type not supported');
    end          
else
    datatype = 'uint8';
end

% Rescale
if isfield(niftiHeader,'scl_slope')||isfield(niftiHeader,'scl_inter')
    slope = 1;
    inter = 0;
    if isfield(niftiHeader,'scl_slope')&&niftiHeader.scl_slope~=0
        slope = niftiHeader.scl_slope;
    end
    if isfield(niftiHeader,'scl_inter')
        inter = niftiHeader.scl_inter;
    end
    IM.Data = (IM.Data - inter)/slope;
end

% Write the Nifti File
fid=fopen(fullFilename, 'wb');
fwrite(fid, niftiHeader.sizeof_hdr, 'int'); %sizeof_hdr
fwrite(fid, niftiHeader.data_type, 'char'); %UNUSED data_type
fwrite(fid, niftiHeader.db_name, 'char'); %UNUSED db_name
fwrite(fid, niftiHeader.extents, 'int'); %UNUSED extents
fwrite(fid, niftiHeader.session_error, 'short'); %UNUSED session_error
fwrite(fid, niftiHeader.regular, 'char'); %UNUSED regular
fwrite(fid, niftiHeader.dim_info, 'char'); %dim_info
fwrite(fid, niftiHeader.dim, 'short'); %dim
fwrite(fid, niftiHeader.intent_p1, 'float'); %intent parameter 1 ***
fwrite(fid, niftiHeader.intent_p2, 'float'); %intent parameter 2 ***
fwrite(fid, niftiHeader.intent_p3, 'float'); %intent parameter 3 ***
fwrite(fid, niftiHeader.intent_code, 'short'); %NIFIT intent code ***
fwrite(fid, niftiHeader.datatype, 'short'); %datatype
fwrite(fid, niftiHeader.bitpix, 'short'); %bitpix
fwrite(fid, niftiHeader.slice_start, 'short'); %first slice number
fwrite(fid, niftiHeader.pix_dim, 'float'); %grid spacings ***
fwrite(fid, niftiHeader.vox_offset, 'float'); %vox_offset
fwrite(fid, niftiHeader.scl_slope, 'float'); %scl_slope
fwrite(fid, niftiHeader.scl_inter, 'float'); %scl_inter
fwrite(fid, niftiHeader.slice_end, 'short'); %last slice number
fwrite(fid, niftiHeader.slice_code, 'char'); %slice timing code *** WHY CHAR?
fwrite(fid, niftiHeader.xyzt_units, 'char'); %xyzt_units
fwrite(fid, niftiHeader.cal_max, 'float'); %display max
fwrite(fid, niftiHeader.cal_min, 'float'); %display min
fwrite(fid, niftiHeader.slice_duration, 'float'); %slice time
fwrite(fid, niftiHeader.toffset, 'float'); %UNUSED slice_toffset
fwrite(fid, niftiHeader.glmax, 'int'); %unused
fwrite(fid, niftiHeader.glmin, 'int'); %unused
fwrite(fid, niftiHeader.descrip, 'char'); %description
fwrite(fid, niftiHeader.aux_file, 'char'); %alternate filename
fwrite(fid, niftiHeader.qform_code, 'short'); %qform code ***
fwrite(fid, niftiHeader.sform_code, 'short'); %sform code ***
fwrite(fid, niftiHeader.quartern_b, 'float'); %Quarternion b param 
fwrite(fid, niftiHeader.quartern_c, 'float'); %Quarternion c param 
fwrite(fid, niftiHeader.quartern_d, 'float'); %Quarternion d param 
fwrite(fid, niftiHeader.qoffset_x, 'float'); %Quarternion x param ***
fwrite(fid, niftiHeader.qoffset_y, 'float'); %Quarternion y param ***
fwrite(fid, niftiHeader.qoffset_z, 'float'); %Quarternion z param ***
fwrite(fid, niftiHeader.srow_x, 'float'); %1st row affine transform
fwrite(fid, niftiHeader.srow_y, 'float'); %2nd row affine transform
fwrite(fid, niftiHeader.srow_z, 'float'); %3rd row affine transform
fwrite(fid, niftiHeader.intent_name, 'char'); %name of data
fwrite(fid, niftiHeader.magic, 'uint8'); %Magic string
fwrite(fid, '    ', 'char'); %Blank characters
fwrite(fid, IM.Data, datatype);
fclose(fid);

return;
