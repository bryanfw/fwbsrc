%% CLONEDICOM     Clone a Philips DICOM file
%
% [SUCCESS] = CLONEDICOM(OUTPUTNAME,DATA,INFO)
%
%   OUTPUTNAME is the output DICOM filename
%
%   DATA is an N-dimensional array holding the image data.
%
%   INFO is a structure containing details from DICOM header.  The name of
%   the DICOM file to clone should already be in the INFO structure.
%
%   DATA and INFO should have the format as that returned by LOADDICOM.
%
%  See also: LOADDICOM
%
%  Dependencies: none
%

%% Revision History
% * 2008.05.22    initial version - brianwelch

%% Function definition
function success = cloneDICOM(outputname,data,info)

% Assume successful unless proven to be otherwise
success = 1;

%% Grab name of DICOM file to clone
clonedname = info.filename;

%% Open cloned file and read header
fidcloned = fopen(clonedname, 'r');
header = fread(fidcloned,info.pixel_data_offset_bytes,'uint8');
fclose(fidcloned);

%% Open output file 
fidoutput = fopen(outputname,'wb');

%% Write DICOM header
write_count = fwrite(fidoutput,header,'uint8');

%% Determine order for writing
% sort linearized info.table_row_index_array and keep sorted indices
[dummy, write_order_indices] = sort(info.table_row_index_array(:)); 

%% Write DICOM image data

% temporarily reshape data to a continuous stack of images
size_data = size(data);
data = reshape(data,[size_data(1) size_data(2) info.n_data_imgs]);

% loop through all image dimensions
for k=1:length(write_order_indices),
    
    % find the table_row_index that is associated with this image
    table_row_index = info.table_row_index_array( write_order_indices(k) );
   
    if(table_row_index>0),        
        % rescale intercept
        ri = info.imgdef.rescale_intercept.vals(table_row_index);

        % rescale slope
        rs = info.imgdef.rescale_slope.vals(table_row_index);

        % scale slope
        ss = info.imgdef.scale_slope.vals(table_row_index);

        % convert data for writing
        tmpimg = data(:,:, write_order_indices(k) );
        
        % Apply image scaling by using info.table_index & info.table
        switch info.loadopts.scale
            case 'FP',
                %data(:,:,k) = (data(:,:,k) * rs + ri)/(rs * ss);
                tmpimg = ( tmpimg * (rs * ss) - ri ) / rs;
                if(k==12+24*2),
                    disp('hello');
                end
            case 'DV',
                %data(:,:,k) = data(:,:,k) * rs + ri;
                tmpimg = ( tmpimg - ri ) / rs;
            case 'PV',
                % do nothing
                % values are already the pixel value stored in the REC file
                % tmpimg = tmpimg;
            otherwise,
                if(k==1),
                    warning( sprintf('Unkown scale type option : ''%s''.  Will assumme floating point (''FP'') instead',info.loadopts.scale) );
                end
                tmpimg = ( tmpimg * (rs * ss) - ri ) / rs;
        end
        
        % transpose image
        tmpimg = tmpimg.';       
        if info.pixel_bits==16,
            count = fwrite(fidoutput,tmpimg(:),'int16','ieee-le');
        elseif info.pixel_bits==8,
            count = fwrite(fidoutput,tmpimg(:),'int8','ieee-le');
        else
            error( sprintf('Unsupported image pixel bits : %d', info.pixel_bits) );
        end 
        
    end 
end
fclose(fidoutput);

