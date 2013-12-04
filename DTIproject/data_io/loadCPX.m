%% LOADCPX     Load a Philips CPX file
%
% [DATA,INFO] = LOADCPX(FILENAME)
%
%   FILENAME is a string containing a file prefix or name of the CPX
%   image file, e.g. MYDATA or MYDATA.CPX
%
%   DATA is an N-dimensional array holding the complex image-space data.
%
%   INFO is a structure containing data details from the CPX file
%
% [DATA,INFO] = LOADCPX([])
%
%   When the passed FILENAME is not provided or is an empty array or empty 
%   string.  The user chooses a file using UIGETFILE.
%
% [DATA,INFO] = LOADCPX(FILENAME,'OptionName1',OptionValue1,...)
%
%   Options can be passed to LOADCPX to control the range/pattern of
%   loaded data, verbose output, etc.  The list below shows the avialable 
%   options.  Names are case-sensitive
%
%       OptionName          OptionValue       Description    
%       ----------          -----------     ---------------  
%       'x'                 numeric         image rows      
%       'y'                 numeric         image columns   
%       'z'                 numeric         image slices
%       'ec'                numeric         echoes          
%       'dyn'               numeric         dynamics        
%       'ph'                numeric         cardiac phases  
%       'row'               numeric         rows  
%       'mix'               numeric         mixes 
%       'coil'              numeric         coils
%       'verbose'           logical         [ true |{false}]
%       'savememory'        logical         [{true}| false ]
%
%       When 'savememory' is true, SINGLE precision is used instead of DOUBLE
%
%   Example:
%       myfile = 'MYDATA.CPX';
%       [data,info] = loadCPX(myfile,'coil',[1 5],'verbose',true);
%
% [DATA,INFO] = LOADCPX(FILENAME,LOADOPTS)
%
%   LOADOPTS is a structure with fieldnames equal to any of the possible
%   OptionNames.
%
%   Example:
%       loadopts.coil = [1 5];
%       loadopts.verbose = true;
%       [data,info] = loadCPX(myfile,loadopts);
%
%   For any dimension, values may be repeated and appear in any order.
%   Values that do not intersect with the available values for that
%   dimension will be ignored.  If the intersection of the user-defined
%   dimension values and the available dimension range has length zero, an
%   error is generated.  The order of the user-defined pattern is preserved.
%
%   Example:
%       % load a specific pattern of locations (-1 will be ignored)
%       loadopts.loc = [1 1 2 1 1 2 -1];
%       [data,info] = loadCPX(myfile,loadopts);
%
% INFO = LOADCPX(FILENAME)
%
%   If only one return argument is provided, the INFO structure will be
%   returned.  DATA will not be loaded (fast execution).
%
% INFO structure contents
%
%   The INFO structure contains all the information from the CPX image
%   512 byte (128 32-bit signed integers or floats) preambles appearing before each 
%   complex image.  The list below describes some the additional fields 
%   found within INFO
%
%   Preamble  Dimension CPX Field                
%   Location  Name      Name                                Description
%   --------  --------- ----------------------------------- ------------------------------------------------------------------
%    1        mix       mix                                 mixed sequence number
%    2        z         location                            (?) (ignored)
%    3                  perp_dir_coord                      slice number
%    4                  kz_profile_nr                       (used for 3D CSI?) (ignored)
%    5        ec        echo                                echo number
%    6        ph        card_phase                          cardiac phase number        
%    7        dyn       dyn_scan                            dynamic scan number
%    8        row       row                                 row number
% 
%    9                  COMPLEX_MATRIX_EXISTENCE            (0 in END_OF_DATA_MARK)
%   10                  COMPLEX_MATRIX_BYTE_OFFSET_32_BITS
%   11-12     x y       COMPLEX_MATRIX_RESOLUTIONS
%   13                  COMPLEX_MATRIX_DATA_SIZE            (in units of disc blocks)  
%   14                  DATA_COMPRESSION_FACTOR             (scaling values for inverse data compression)
%   
%   15                  SEQUENCE_NR                         sequence number (only needed for compatibilty with older releases) (ignored)
%
%   16                  COMPLEX_SCALING_FACTOR              (float)
%   17                  COMPLEX_SCALING_MINIMUM             (float)
%   18                  EXTRA_ATTRIBUTE_VALUE
%
%   19                  RECEIVE_CHANNEL_NR                  (nothing to do with channel numbers of received data)
%
%   20                  UNKNOWN_ATTR_NR
%   21                  CPX_FILE_LAYOUT
%   22        coil      EXTERN_CHANNEL_NR                   (absolute channel index, i.e. coil)
%   23-32?              AD_HOC_INT_ARR                      (length 10?)
%   33-42?              AD_HOC_FLOAT_ARR                    (float) (length 10?)
%   43-44               COMPLEX_MATRIX_BYTE_OFFSET_64_BITS  (allows CPX files > 2GB, CPX layout 2 and higher)
%   45-128              BLOCK_FILLER
%
%
%  References:
%  XJS-154-1528 "RS: Reconstruction" by G. van Ensbergen (rcrsrec.pdf, p. 150) 
%  XJS-154-1529 "IS: Reconstruction" by Adri Duijndam (rcisrec.pdf, pp. 434-440)

%   The INFO.ROW_INDEX_ARRAY is a special array that is the same 
%   size as the DATA array (minus the first two dimensions used to store 
%   X, Y).  A given index for a given CPX 2D image  in the DATA 
%   array will return the label index number describing the details of that 
%   raw data vector in the INFO.LABELS array when that same index is used 
%   with the INFO.ROW_INDEX_ARRAY.  This provides a quick way to 
%   recall the details of any individual CPX 2D image contained within DATA.  
%   If the INFO.ROW_INDEX_ARRAY holds a ZERO for a given index, there 
%   was no label from the LAB file that matched the dimension location in DATA.  
%
%  See also: LOADPARREC, LOADXMLREC, LOADLABRAW
%

%% Revision History
% * 2010.11.03    initial version - welcheb
% * 2012.07.18    major update using actual documentation - welcheb
% * 2013.09.05    changed slice number to come from perp_dir_coord instead of location - welcheb

%% Function definition
function [data,info] = loadCPX(filename,varargin)

%% Start execution time clock and initialize DATA and INFO to empty arrays
tic;
data=[];
info=[];

%% Initialize INFO structure
% Serves to fix the display order
info.filename = [];
info.loadopts = [];
info.dims = [];
info.preamble_table_index_array = [];
info.preamble_table = [];
info.preambles = [];
info.fseek_offsets = [];
info.nPreambles = [];
info.nLoadedImages = [];
info.datasize = [];

%% Allow user to select a file if input FILENAME is not provided or is empty
if nargin<1 || isempty(filename),
    [fn, pn] = uigetfile({'*.CPX'},'Select a CPX file');
    if fn~=0,
        filename = sprintf('%s%s',pn,fn);
    else
        disp('LOADCPX cancelled');
        return;
    end
end

%% Parse the filename.
% It may be the LAB filename, RAW filename or just the filename prefix
% Instead of REGEXP, use REGEXPI which igores case
toks = regexpi(filename,'^(.*?)(\.cpx)?$','tokens');
prefix = toks{1}{1};
cpxname = sprintf('%s.CPX',prefix);
info.filename = filename;

%% Pre-allocate preamble_table
preallocate_count = 20e3;
info.preamble_table = zeros(preallocate_count,127);

%% Open CPX for reading
cpxfid = fopen(cpxname,'rb','ieee-le');

%% Read all CPX preambles
preamble_length = 128; % 128 4-byte integers (512 bytes)
preamble_cnt = 0;
read_flag = 1;
while(read_flag)
    [preamble_001_015, readsize_001_015] = fread(cpxfid,[1 15],'int32');
    [preamble_016_017, readsize_016_017] = fread(cpxfid,[1  2],'float32');
    [preamble_018_032, readsize_018_032] = fread(cpxfid,[1 15],'int32');
    [preamble_033_042, readsize_033_042] = fread(cpxfid,[1 10],'float32');
    [preamble_043_044, readsize_043_044] = fread(cpxfid,[1  1],'int64');
    [preamble_045_128, readsize_045_128] = fread(cpxfid,[1 84],'int32');
    
    CPX_FILE_LAYOUT = preamble_018_032(21-18+1);
    if CPX_FILE_LAYOUT~=2,
        error('CPX_FILE_LAYOUT = %d (only 2 is supported)', CPX_FILE_LAYOUT);
    end
    
    if ( (readsize_045_128==84) && (preamble_001_015(9)==1) ),
        preamble_cnt = preamble_cnt + 1;
        info.preamble_table(preamble_cnt,:) = [preamble_001_015 preamble_016_017 preamble_018_032 preamble_033_042 preamble_043_044 preamble_045_128];    
        x = info.preamble_table(preamble_cnt,11);
        y = info.preamble_table(preamble_cnt,12);
        data_compression_factor = info.preamble_table(preamble_cnt,14);
        bytes_to_next_preamble = x * y * 4 * 2 / data_compression_factor;
        fseek(cpxfid, bytes_to_next_preamble, 'cof');
    else
        read_flag = 0;
    end
    
end
fclose(cpxfid);
info.nPreambles = preamble_cnt;

%% Delete unused preallocated columns
info.preamble_table((preamble_cnt+1):end,:) = [];

%% Parse preamble table
info.preambles.mix.vals  = info.preamble_table(:,1)';
%info.preambles.z.vals    = info.preamble_table(:,2)'; % location
info.preambles.z.vals    = info.preamble_table(:,3)'; % perp_dir_coord
info.preambles.ec.vals   = info.preamble_table(:,5)';
info.preambles.ph.vals   = info.preamble_table(:,6)';
info.preambles.dyn.vals  = info.preamble_table(:,7)';
info.preambles.row.vals  = info.preamble_table(:,8)';
info.preambles.x.vals = info.preamble_table(:,11)';
info.preambles.y.vals = info.preamble_table(:,12)';
info.preambles.coil.vals   = info.preamble_table(:,22)';

info.preambles.data_compression_factor.vals = info.preamble_table(:,14)';
info.preambles.byte_offset.vals = info.preamble_table(:,43)';

%% Find unique values of each label field
info.preamble_fieldnames = fieldnames(info.preambles);
for k=1:length(info.preamble_fieldnames),
    info.preambles.(info.preamble_fieldnames{k}).uniq = unique( info.preambles.(info.preamble_fieldnames{k}).vals ); 
end

%% Dimension names
dimnames  = {'x','y','z','ec','dyn','ph','mix','row','coil'};
dimfields = {'x','y','z','ec','dyn','ph','mix','row','coil'};

%% Calculate dimensions of normal data
info.dims.nX             = max(info.preambles.x.uniq);
info.dims.nY             = max(info.preambles.y.uniq);
info.dims.nZ             = length(info.preambles.z.uniq);
info.dims.nEc            = length(info.preambles.ec.uniq);
info.dims.nDyn           = length(info.preambles.dyn.uniq);
info.dims.nPh            = length(info.preambles.ph.uniq);
info.dims.nMix           = length(info.preambles.mix.uniq);
info.dims.nRow           = length(info.preambles.row.uniq);
info.dims.nCoil          = length(info.preambles.coil.uniq);

%% With known possible dimension names, the load options can now be parsed
p = inputParser;
p.StructExpand = true;
p.CaseSensitive = true;
p.KeepUnmatched = false; % throw an error for unmatched inputs
p.addRequired('filename', @ischar);
for k=1:length(dimnames),
    p.addParamValue(dimnames{k}, [], @isnumeric);
end
p.addParamValue('verbose', false, @islogical);
p.addParamValue('savememory', true, @islogical);
p.parse(filename, varargin{:});

%% Return loadopts structure inside INFO structure
% remove filename field - it is passed as the first required argument
info.loadopts = rmfield(p.Results,'filename');

%% Find the unique set of values for each dimension name
info.dims.x    = 1:info.dims.nX;
info.dims.y    = 1:info.dims.nY;
info.dims.z    = (1:info.dims.nZ)-1;
info.dims.ec   = info.preambles.ec.uniq;
info.dims.dyn  = info.preambles.dyn.uniq;
info.dims.ph   = info.preambles.ph.uniq;
info.dims.mix  = info.preambles.mix.uniq;
info.dims.row  = info.preambles.row.uniq;
info.dims.coil = info.preambles.coil.uniq;

%% Find intersection of available dimensions with LOADOPTS dimensions
for k=1:length(dimnames),
    if ~isempty(info.loadopts.(dimnames{k})),
        info.dims.(dimnames{k}) = intersect_a_with_b(info.loadopts.(dimnames{k}),info.dims.(dimnames{k}));
    end
end

%% Calculate data size
datasize = []; 
for k=1:length(dimnames),
    datasize = [datasize length(info.dims.(dimnames{k}))];
end
info.datasize = datasize;

% throw error if any dimension size is zero
if any(info.datasize==0),
    zero_length_str = sprintf(' ''%s'' ', dimnames{find(info.datasize==0)});
    error('size of selected data to load has zero length along dimension(s): %s', zero_length_str);
end

%% Skip data loading if only one output argument is provided, return INFO
if nargout==1,
    info.preamble_index_array = [1:size(info.preambles,1)];
    data=info;
    return;
end

%% Create array to hold label row numbers for loaded data
% skip the coil and kx dimensions
info.preamble_row_index_array = zeros(datasize(3:end));

%% Pre-allocate DATA array
if info.loadopts.savememory==true,
    data = complex(zeros(info.datasize,'single'),zeros(info.datasize,'single'));
else
    data = complex(zeros(info.datasize),zeros(info.datasize));
end

%% Read RAW data for selected dimension ranges
cpxfid = fopen(cpxname,'r','ieee-le');
if cpxfid<0,
    error(sprintf('cannot open CPX file: %s', cpxname));
end

info.nLoadedImages=0;

cpx_data_fread_size = double(info.dims.nX * info.dims.nY * 2);
cpxdata_2d = complex(zeros(info.dims.nX,info.dims.nY),zeros(info.dims.nX,info.dims.nY));

for n=1:info.nPreambles,
    
    load_flag=1;
    dim_assign_indices_full_array = [];
           
    for k=3:length(dimfields), % 'x' and 'y' are the first 3 dimfields
        
        dimval = info.preambles.(dimfields{k}).vals(n);
        
        % it is allowed that the dimval appears more than once 
        % in the requested dimension ranges to be loaded
        dim_assign_indices = find(dimval==info.dims.(dimnames{k}));
        
        if isempty(dim_assign_indices),
            load_flag=0;
            break;
        else
           
            if k>3,
                
                dim_assign_indices_full_array_new = zeros( size(dim_assign_indices_full_array,1)*length(dim_assign_indices), size(dim_assign_indices_full_array,2)+1);
                
                mod_base_a = size(dim_assign_indices_full_array,1);
                mod_base_b = length(dim_assign_indices);
                
                for d=1:size(dim_assign_indices_full_array_new,1),
                    dim_assign_indices_full_array_new(d,:) = [dim_assign_indices_full_array(mod(d,mod_base_a)+1,:) dim_assign_indices(mod(d,mod_base_b)+1)];
                end
                
            else
                dim_assign_indices_full_array_new = dim_assign_indices(:);
            end
            
            dim_assign_indices_full_array = dim_assign_indices_full_array_new;
            
        end
    end
    
    if load_flag,
        
        info.nLoadedImages = info.nLoadedImages+1;
        
        byte_offset = double( info.preambles.byte_offset.vals(n) );
        data_compression_factor = double( info.preambles.data_compression_factor.vals(n) );
        
        status = fseek(cpxfid, byte_offset, 'bof');
        
        switch data_compression_factor,
            case 1,
                read_type = 'float32';
            case 2,
                read_type = 'int16';
            case 4,
                read_type = 'int8';
        end
        
        cpxdata_1d = double(fread(cpxfid, cpx_data_fread_size, read_type));
        
        cpxdata_2d = reshape( complex(cpxdata_1d(1:2:end),cpxdata_1d(2:2:end)),[info.dims.nX info.dims.nY]);
        
        % insert cpxdata_2d into proper locations of the data array
        for d=1:size(dim_assign_indices_full_array,1),
            
            dim_assign_str = sprintf(',%d', dim_assign_indices_full_array(d,:) );
            
            % delete initial comma
            dim_assign_str(1) = [];
            
            % assign index to table_index table
            eval( sprintf('info.preamble_row_index_array(%s)=%d;', dim_assign_str, n) );
        
            % assign read image to correct location in data array
            eval( sprintf('data(:,:,%s)=cpxdata_2d;', dim_assign_str) );
                    
        end
    else
        %if info.loadopts.verbose==true,
        %    disp( sprintf('%d not loaded', n) );
        %end
    end
    
end

%% Close CPX file
fclose(cpxfid);

%% If VERBOSE, display execution information
if info.loadopts.verbose==true,
    disp( sprintf('Loaded %d of %d available CPX images', info.nLoadedImages, info.nPreambles) );
    tmpstr = '';
    for k=1:length(dimnames),
        tmpstr = sprintf('%s, # %s: %d', tmpstr, dimnames{k}, length(info.dims.(dimnames{k})) );
    end
    disp( sprintf('Total execution time = %.3f seconds', toc) );
end

%% Find intersection of vector a with vector b without sorting 
function c = intersect_a_with_b(a,b)
c = a;
% work backwards in order to use [] assignment
for k=length(a):-1:1,
    if isempty(find(a(k)==b)),
        c(k)=[]; 
    end
end

% force c to be a row vector for easier display
c = c(:).';


