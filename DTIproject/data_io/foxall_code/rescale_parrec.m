%
% RESCALE_PARREC
%
% Rescale image data read from a Philips .PAR/.REC file pair.
%
% The data stored in the .REC file is in a compressed integer 
% format derived by applying a linear transformation to the
% orignal floating point MR data as it was stored off into 
% scanner the database. The integer form is a convenient
% display format, but not useful for image processing 
% operations. The rescale process restores the data to 
% a floating point form suitable for onward processing.
%
% FORMAT:
%
% rescaled_data = rescale_parrec(image_data,parfile)
%
% INPUT:
%
% image_data       int[]       Multidimensional image data 
%                              array read from the .REC file.
%
% parfile          structure   .PAR file information held in
%                               MATLAB structure format.
%  
% OUTPUT:  
%  
% rescaled_data    float[]     Converted image series as a  
%                              multidimensional MATLAB array.  
%  
%  
% $Id: rescale_parrec.m,v 1.3 2004/08/24 17:52:11 dfoxall Exp $  
%  

function rescaled_data = rescale_parrec(image_data,parfile)

SLOPE_INDEX     = 0;
INTER_INDEX     = 0;

if ( isfield(parfile,'version') )
    switch( char(parfile.version) )
    case { 'V3' }
         SLOPE_INDEX  = 9;
         INTER_INDEX  = 8;
    case { 'V4' }
         SLOPE_INDEX  = 13;
         INTER_INDEX  = 12;
    otherwise
         SLOPE_INDEX  = -1;
         INTER_INDEX  = -1;
    end
end

size_image_data = size(image_data);
TNI             = size_image_data(1);
NLINE           = size_image_data(2);
stride          = size_image_data(3);

for I = 1:TNI
    slope                = parfile.scan_tags(I,SLOPE_INDEX);
    inter                = parfile.scan_tags(I,INTER_INDEX);
%   fprintf(1,'rescaling image %d with slope = %f and intercept = %f \n',I,slope,inter);
    rescaled_data(I,:,:) = slope*image_data(I,:,:) + inter;
end
