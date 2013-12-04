function dicomheader = dicominfofast(filename,parrec_version)

%% Open the file.
fid = fopen(filename, 'r');
if (fid < 0)
  error( sprintf('Cannot open file "%s" for reading.', filename) );
end

%% Initialize DICOMHEADER structure
% Serves to fix the display order
dicomheader.parrec_version = [];
dicomheader.series = [];
dicomheader.stack = [];
dicomheader.frame = [];
dicomheader.stack_cnt = [];
dicomheader.frame_cnt = [];
dicomheader.tags = [];
dicomheader.pixel_data_offset_bytes = [];
dicomheader.pixel_data_length_bytes = [];

%% Initialize DICOM file metadata header entries
dicomheader.series.TransferSyntaxUID = [];
dicomheader.series.MediaStorageSOPClassUID = [];

%% Load PHILIPS DICOM details
philips = dictPARREC;
pardef_fieldnames = fieldnames(philips.pardef);
imgdef_fieldnames = fieldnames(philips.imgdef);

%% set parrec_version
if nargin==1,
    dicomheader.parrec_version = 0;
    for k=1:length(pardef_fieldnames),
        if ( (philips.pardef.(pardef_fieldnames{k}).par_min_version < Inf) & (philips.pardef.(pardef_fieldnames{k}).par_min_version > dicomheader.parrec_version) ),
            dicomheader.parrec_version = philips.pardef.(pardef_fieldnames{k}).par_min_version;
        end
    end
else
    dicomheader.parrec_version = parrec_version;
end

%% Create series, stack and frame level structures
dicomheader.tags.series.tag_ids = sparse(65536,65536);
dicomheader.tags.stack.tag_ids  = sparse(65536,65536);
dicomheader.tags.frame.tag_ids  = sparse(65536,65536);

dicomheader.tags.series.tag_cols = sparse(65536,65536);
dicomheader.tags.stack.tag_cols  = sparse(65536,65536);
dicomheader.tags.frame.tag_cols  = sparse(65536,65536);

dicomheader.tags.series.tag_cnt = 0;
dicomheader.tags.stack.tag_cnt = 0;
dicomheader.tags.frame.tag_cnt = 0;

dicomheader.tags.series.tag_names = {};
dicomheader.tags.stack.tag_names = {};
dicomheader.tags.frame.tag_names = {};

% parameter definitions
for k=1:length(pardef_fieldnames),

    if ( (philips.pardef.(pardef_fieldnames{k}).par_min_version <= dicomheader.parrec_version) & (philips.pardef.(pardef_fieldnames{k}).par_max_version >= dicomheader.parrec_version) ),
        
        group = philips.pardef.(pardef_fieldnames{k}).group;
        dicom_sq = philips.pardef.(pardef_fieldnames{k}).dicom_sq;
        for n = 1:length(philips.pardef.(pardef_fieldnames{k}).element),
            element = philips.pardef.(pardef_fieldnames{k}).element(n);
            dicomheader.(dicom_sq).(pardef_fieldnames{k}) = [];
            dicomheader.tags.(dicom_sq).tag_cnt = dicomheader.tags.(dicom_sq).tag_cnt + 1;
            dicomheader.tags.(dicom_sq).tag_ids(group+1,element+1) = dicomheader.tags.(dicom_sq).tag_cnt;
            dicomheader.tags.(dicom_sq).tag_cols(group+1,element+1) = n;
            dicomheader.tags.(dicom_sq).tag_names{ dicomheader.tags.(dicom_sq).tag_cnt } = pardef_fieldnames{k};
        end

    end

end

% image definitions
for k=1:length(imgdef_fieldnames),
    
    if ( (philips.imgdef.(imgdef_fieldnames{k}).par_min_version <= dicomheader.parrec_version) & (philips.imgdef.(imgdef_fieldnames{k}).par_max_version >= dicomheader.parrec_version) ),
    
        group = philips.imgdef.(imgdef_fieldnames{k}).group;
        dicom_sq = philips.imgdef.(imgdef_fieldnames{k}).dicom_sq;
        for n = 1:length(philips.imgdef.(imgdef_fieldnames{k}).element),
            element = philips.imgdef.(imgdef_fieldnames{k}).element(n);
            dicomheader.(dicom_sq).(imgdef_fieldnames{k}) = [];
            dicomheader.tags.(dicom_sq).tag_cnt = dicomheader.tags.(dicom_sq).tag_cnt + 1;
            dicomheader.tags.(dicom_sq).tag_ids(group+1,element+1) = dicomheader.tags.(dicom_sq).tag_cnt;
            dicomheader.tags.(dicom_sq).tag_cols(group+1,element+1) = n;
            dicomheader.tags.(dicom_sq).tag_names{ dicomheader.tags.(dicom_sq).tag_cnt } = imgdef_fieldnames{k};
        end
    
    end
    
end

%% Get the possible DICOM header and inspect it for DICOM-like data.
header = fread(fid, 132, 'uint8');
if ~( ( numel(header)==132 ) && ( sum(header(1:128))==0 ) && isequal(char(header(129:132))','DICM') )
    % It's not a proper DICOM file.
    data_details = [];
    fclose(fid);
    return;
end

%% Read file meta data
file_endian = 'ieee-le';
group   = fread(fid, 1, 'uint16', file_endian);
element = fread(fid, 1, 'uint16', file_endian);
VR      = fread(fid, 2, 'uint8=>char', file_endian);
VL      = fread(fid, 1, 'uint16', file_endian);
GL0002  = fread(fid, 1, 'uint32', file_endian); 

%% Precalculate group and element values
g0002 = hex2dec('0002'); e0010 = hex2dec('0010'); % TransferSyntaxUID
g0002 = hex2dec('0002'); e0002 = hex2dec('0002'); % MediaStorageSOPClassUID
g7FE0 = hex2dec('7FE0'); e0010 = hex2dec('0010');
gFFFE = hex2dec('FFFE'); eE000 = hex2dec('E000');
gFFFE = hex2dec('FFFE'); eE00D = hex2dec('E00D');
gFFFE = hex2dec('FFFE'); eE0DD = hex2dec('E0DD');

% Stack SQ, begin stack section, section_current = section_stack
g2001 = hex2dec('2001'); e105F = hex2dec('105F'); 
% "item begin tag" at SQ_depth = 1 increments
% "end of SQ tag" at SQ_depth = 1 sets section_current = section_series

% PerFrameFunctionalGroupsSequence, begins frame section, section_current = section_frame
g5200 = hex2dec('5200'); e9230 = hex2dec('9230');
% "item begin tag" at SQ_depth = 1 increments
% "end of SQ tag" at SQ_depth = 1 sets section_current = section_series

%% VR names groups
VR_names_with_int32_VL = { ('OB')' ('OW')' ('SQ')' ('UN')' };
VR_names_char    = {('AE')' ('AS')' ('CS')' ('DA')' ('DS')' ('DT')' ('IS')' ('LO')' ('LT')' ('PN')' ('SH')' ('ST')' ('TM')' ('UI')' ('UT')' };
VR_names_uint8   = { ('OB')' };
VR_names_uint16  = { ('AT')' ('OW')' ('US')' };
VR_names_int16   = { ('SS')' };
VR_names_uint32  = { ('UL')' };            
VR_names_int32   = { ('SL')' };
VR_names_float32 = { ('FL')' ('OF')' };
VR_names_float64 = { ('FD')' };
VR_names_unknown = { ('UN')' ('SQ')' };

%% Read file meta data
read_count = 132 + 2 + 2 + 2 + 2 + 4;
while read_count < (GL0002 + 144),
    group   = fread(fid, 1, 'uint16', file_endian);
    element = fread(fid, 1, 'uint16', file_endian);
    VR      = fread(fid, 2, 'uint8=>char', file_endian);
    read_count = read_count + 6;
    
    switch VR,
        case VR_names_with_int32_VL,
            dummy = fread(fid, 2, 'uint8', file_endian);
            VL = fread(fid, 1, 'int32', file_endian);
            read_count = read_count + 6;    
        otherwise,
            VL = fread(fid, 1, 'int16', file_endian);
            read_count = read_count + 2;
    end
    value = fread(fid, VL, 'uint8', file_endian);
    read_count = read_count + VL;
    
    if group==g0002 && element==e0010,
        dicomheader.series.TransferSyntaxUID = char(value(1:(VL-1))');
    end
    
    if group==g0002 && element==e0002,
        dicomheader.series.MediaStorageSOPClassUID = char(value(1:(VL-1))');
    end
end

%% Set file endian
switch dicomheader.series.TransferSyntaxUID,
    case '1.2.840.10008.1.2.1',
        file_endian = 'ieee-le';
    case '1.2.840.10008.1.2.2',
        file_endian = 'ieee-be';
    otherwise,
        fclose(fid);
        error( sprintf('Unsupported transfer syntax : %s', transfer_syntax) );
end

%%  Section identifiers
section_series = 1;
section_stack  = 2;
section_frame  = 3;
section_current = section_series;

%% Section counters
dicomheader.stack_cnt = 0;
dicomheader.frame_cnt = 0;
SQ_depth = 0;

hwaitbar = waitbar(0,'Parsing DICOM header ...');

%% Iterative read buffer
while true,
    
    [group_element, count] = fread(fid, 2, 'uint16', file_endian); 
    group   = group_element(1)
    element = group_element(2)
    VR      = fread(fid, 2, 'uint8=>char', file_endian);
    read_count = read_count + 6

    if VR(1)=='S' && VR(2)=='Q',
        SQ_depth = SQ_depth + 1;
    end
    
    if count<2,
        
        fclose(fid);
        error('reached end of the file');
  
    elseif group==g7FE0 && element==e0010,
        
        dummy = fread(fid, 2, 'uint8', file_endian);
        VL = fread(fid, 1, 'int32', file_endian);
        read_count = read_count + 6;
        
        dicomheader.pixel_data_offset_bytes = read_count;        
        dicomheader.pixel_data_length_bytes = VL;
        break;
 
    elseif group==gFFFE,
        
        switch element,
            case eE000, % item begin
                VL = fread(fid, 1, 'int16', file_endian);
                read_count = read_count + 2;      

                if SQ_depth==1,
                    if section_current==section_frame,
                        dicomheader.frame_cnt = dicomheader.frame_cnt + 1;
                        % pre-allocate frames
                        if dicomheader.frame_cnt==2,
                            nframes = str2num(dicomheader.series.NumberOfFrames{1});
                            dicomheader.frame(2:nframes) = dicomheader.frame(1);
                        end
                        
                        if mod(dicomheader.frame_cnt,25)==0,
                            waitbar(dicomheader.frame_cnt/nframes, hwaitbar, sprintf('Parsing DICOM header ... frames %5d of %5d', dicomheader.frame_cnt, nframes) );
                        end
                        
                    elseif section_current==section_stack,
                        dicomheader.stack_cnt = dicomheader.stack_cnt + 1;
                    end
                end
                
            case eE00D, % item end
                VL = fread(fid, 1, 'int16', file_endian);
                read_count = read_count + 2;
                
            case eE0DD, % end of SQ              
                VL = fread(fid, 1, 'int16', file_endian);
                read_count = read_count + 2;
                
                SQ_depth = SQ_depth - 1;
                
                if SQ_depth==0,
                    section_current = section_series;
                end
        end
        
    else
    
        if group==g2001 && element==e105F, % Stack
            section_current = section_stack;
        end
                
        if group==g5200 && element==e9230, % PerFrameFunctionalGroupsSequence
            section_current = section_frame;
        end
        
        switch VR,
            case VR_names_with_int32_VL,
                dummy = fread(fid, 2, 'uint8', file_endian);
                VL = fread(fid, 1, 'int32', file_endian);
                read_count = read_count + 6;    
            otherwise,
                VL = fread(fid, 1, 'int16', file_endian);
                read_count = read_count + 2;
        end

        if VL>0,
            
            switch VR,
                case VR_names_char,
                    value = strtrim( fread(fid, VL, 'uint8=>char', file_endian)' );
                case VR_names_uint8,
                    value = fread(fid, VL, 'uint8', file_endian)';
                case VR_names_uint16,
                    value = fread(fid, VL/2, 'uint16', file_endian)';
                case VR_names_int16,
                    value = fread(fid, VL/2, 'int16', file_endian)';
                case VR_names_uint32,
                    value = fread(fid, VL/4, 'uint32', file_endian)';                   
                case VR_names_int32,
                    value = fread(fid, VL/4, 'int32', file_endian)';
                case VR_names_float32
                    value = fread(fid, VL/4, 'float32', file_endian)';
                case VR_names_float64
                    value = fread(fid, VL/8, 'float64', file_endian)';
                otherwise % VR_names_unknown
                    value = fread(fid, VL, 'uint8', file_endian)';                   
            end
            read_count = read_count + VL;
            
            switch section_current,
                case section_series,
                    tag_id  = dicomheader.tags.series.tag_ids(group+1,element+1);
                    tag_col = dicomheader.tags.series.tag_cols(group+1,element+1);
                    if tag_id>0,
                        tag_name = dicomheader.tags.series.tag_names{tag_id};
                        dicomheader.series.(tag_name){tag_col} = value;
                    end                  
                case section_stack,
                    tag_id  = dicomheader.tags.stack.tag_ids(group+1,element+1);
                    tag_col = dicomheader.tags.stack.tag_cols(group+1,element+1);
                    if tag_id>0,
                        tag_name = dicomheader.tags.stack.tag_names{tag_id};
                        dicomheader.stack(dicomheader.stack_cnt).(tag_name){tag_col} = value;
                    end
                case section_frame,
                    tag_id  = dicomheader.tags.frame.tag_ids(group+1,element+1);
                    tag_col = dicomheader.tags.frame.tag_cols(group+1,element+1);
                    if tag_id>0,
                        tag_name = dicomheader.tags.frame.tag_names{tag_id};
                        dicomheader.frame(dicomheader.frame_cnt).(tag_name){tag_col} = value;
                    end    
            end
                        
        end
        
    end
end
close(hwaitbar);
fclose(fid);




end