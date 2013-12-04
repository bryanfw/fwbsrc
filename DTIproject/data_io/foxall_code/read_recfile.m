%
% READ RECFILE
%
% Read the binary data from the .REC file of a Philips  
% .PAR/.REC pair. Format the data into a multidimensional
% array that contains an indexed set of images.
%
% FORMAT:
%
% [image_data,TNI,NLINE,stride] = read_recfile(recfile_name,parfile);
%
% INPUT:
%
% recfile_name     string        File path name of the .REC file
% 
% parfile          structure     MATLAB structure read from the .PAR
%                                file.
%
% OUTPUT:
%
% image_data       int[]         Multidimensional integer array that
%                                contains data from the .REC file.
%                                The array has three indices:
%
%                                1. Image index running 1...TNI
%                                2. Phase encode line index 1...NLINE
%                                3. Read out pixel index 1...stride
%
% TNI              int           Total number of images
%
% NLINE            int           Total number of phase encode lines
%
% stride           int           Total number of read pixels per line.
%
% $Id: read_recfile.m,v 1.2 2004/08/24 17:07:02 dfoxall Exp $
%

function [image_data,TNI,NLINE,stride] = read_recfile(recfile_name,parfile);

%-------------------------------------------------------------------
% INTERPRET THE .PAR FILE
%
% The information needed to interpret the format of the .REC file
% correctly is stored in the MATLAB structure parfile passed as
% an argument to this function. It expects the parfile fields to have
% the standard names defined in the V3.m and V4.m M-FILES. If
% user selected field names have been applied when reading the
% .PAR file expect trouble. 
%
%-------------------------------------------------------------------
V4_PIXEL_BIT_INDEX   = 8;
V4_RECON_RES_INDEX_1 = 10;
V4_RECON_RES_INDEX_2 = 11;
pixel_bits           = 0;

if ( isfield(parfile,'version') )
   switch( char(parfile.version) )
   case { 'V3' }
        pixel_bits = parfile.pixel_bits;
   case { 'V4' }
        pixel_bits = parfile.scan_tags(1,V4_PIXEL_BIT_INDEX);
   otherwise
        pixel_bits = -1;
   end
end

scan_tag_size  = size(parfile.scan_tags);
TNI            = scan_tag_size(1);

if ( isfield(parfile,'recon_resolution') )
    NLINE      = parfile.recon_resolution(1);
    stride     = parfile.recon_resolution(2);
else
    NLINE      = parfile.scan_tags(1,V4_RECON_RES_INDEX_1);
    stride     = parfile.scan_tags(1,V4_RECON_RES_INDEX_2);
end

fprintf(1,'\n');
fprintf(1,'total number of (%d X %d) images expected is:  %d\n', NLINE, stride, TNI);


switch (pixel_bits)
    case { 8 }
         read_type = 'int8';
    case { 16 }
         read_type = 'int16';
    otherwise
         read_type = 'uchar';
end

fprintf(1,'data pixel size  : %d\n',pixel_bits);
fprintf(1,'reading data type: %s\n',read_type);

%-------------------------------------------------------------------
% READ THE .REC FILE
%-------------------------------------------------------------------
fid                         = fopen(recfile_name,'r');
[binary_1D,read_size_act]   = fread(fid,inf,read_type);
fclose(fid);

fprintf(1,'Expecting to read %d integers.  Found %d integers \n',TNI*NLINE*stride,read_size_act);

if (read_size_act ~= TNI*NLINE*stride)

    if (read_size_act > TNI*NLINE*stride)
        error('.REC file has more data than expected from .PAR file')
    else
        error('.REC file has less data than expected from .PAR file')
    end

else
    fprintf(1,'.REC file read sucessfully \n');
end

fprintf(1,'\n');

%-------------------------------------------------------------------
% FORMAT THE BINARY DATA FOR OUTPUT
%-------------------------------------------------------------------
for I  = 1:TNI
    start_image       = (I - 1)*NLINE*stride;
for L  = 1:NLINE
    start_line        = (L - 1)*stride;
    start_pix         = start_image + start_line + 1;
    end_pix           = (start_pix - 1) + stride;
    image_data(I,L,:) = binary_1D(start_pix:end_pix)';
end
end

