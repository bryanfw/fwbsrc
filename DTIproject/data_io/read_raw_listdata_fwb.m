function [data_bundle info]=read_raw_listdata_fwb(list_fname)
% [kspace_data_struct info_struct]=READ_RAW_LISTDATA_FWB(dotlist_filename)
% reads raw_xxx.list/.data or cpx_xxx.list/.data files
% 
% cpx_xxx.list/.data notes:
%   The output is a sub-set of the outputs listed below (outputs marked by ^^
%   are only for raw_xxx.list/.data)
%   
%
% raw_xxx.list/.data notes:
%    it also applies FRX and PHX phase correction
%    designed & tested on EPI, may also work for TSE / GRaSE
%    Does not FFT anything or partial-fourier ("half-scan") correct
%    Data is provided according to its type (NOI,FRX,PHX,NAV,STD)
%    Data is arranged in k-space according to the dimensions given in the .list
%    file header. Therefore, there will be some all zero-lines, due to :
%        half-scan 
%        sense
%        partial-echo 
%        because the collected voxel dimensions didn't fit perfectly into recon
%           voxel dimensions (fixed via zero-padding & k-space interpolation)
%
% Input: 
%   -.list filename
%       If you don't give it an input argument, it will bring up a GUI
%
% Output: 
%   -data_bundle (1st output) 
%       all the data in the .data file parse and organized to the 
%       best of my ability. also a couple of masks that you *may* find helpful 
%       (don't hope for much there)
%
%       For each type of data line present (STD,PHX^^,NOI^^,etc^^) you get a
%       struct.The fields for each struct can differ, but they are all listed
%       here:
%         +   .info - a struct of arrays containing all the information from the
%               .list textfile lines - unmodified
%         +   .info_adjusted - a similar struct, but modified so all indexing 
%               now uses positive numbers (MATLAB compatible)
%         +   .data - the parsed and organized data, 11-D or so
%         +   .data_corrected^^ - all FRX/PHX corrections applied and reduced
%              readout dimension to just FOV
%         +   .dim_labels - cell array of strings telling the user what each
%               dimension in .data is for
%
%    -info (2nd output)
%       This struct that contains unmodified info parsed from the .list file's
%       header
%       It may also contain other useful information that was calculated in here 
%       (calculated fields suffixed are with __calc)
%
%
%
%          
% Frederick W Bryan 2013
% 
% Put together with parts and tricks from:
% Brian Welch        readListData.m
% Hui Wang           read_datafile_complex.m
%                    read_listfile_complex.m
%                    read_listdata_complex.m
% Amol Pednekar      ReadLIST_withPhaseCorrection.m
% These are all (distantly) descended from David Foxall's scripts
%
%
% Revision History
% 2013-10-30 fwb - created
% 2013-11-04 fwb - continual changes
% 2013-11-06 fwb - made compatible to read both .list and .data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GET/LOOK FOR & PARSE .LIST FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check outputs/inputs
if nargout<2
    error('read_raw_listdata:nargout_too_small','Not enough output arguments');
end

% call gui to select list file 
if nargin<1 
    [fname,pname] = uigetfile('*.LIST','Select *.LIST file');
    list_fname=[pname fname];
else
    [pname fname ext] =fileparts(list_fname);
    fname = [fname ext];
end

if ~strcmpi(list_fname(end-4:end),'.LIST'); % strcmpi is case insensitive
    error('read_raw_listdat:not_list_file',...
        'File not selected or not a .LIST file.');
end

% look for .data file in same directory
files = dir([list_fname(1:end-5) '.*']); files = {files.name};

matching_files = strcmpi([fname(1:end-5) '.DATA'],files);
if sum(matching_files) == 1
    data_fname = [pname filesep files{matching_files}];
else
    error('read_raw_listdata:datamissing','.data file not found');
end
  
tmp_fname = [list_fname(1:end-5) '_tmp.mat'];
if exist(tmp_fname,'file'); % if this file exists, skip a bunch of stuff
    load(tmp_fname)
    fprintf('Saved time by loading temporary file: %s\n',tmp_fname);
else

    % read the list file
    fid = fopen(list_fname,'r');
    if fid==-1,
        error('read_raw_listdata:fileerror',...
            'cannot open %s for reading', list_fname );
    end

    % parse the list file until the data vector lines
    % this is so that we can initialize the matrix sizes and determine
    % the structure of our data struct
    stop_condition = 0; 
    ind = 0;
    while stop_condition == 0;
        line = fgetl(fid); 
        ind = ind+1;

        if line == -1;
            error('got to end of file before seeing data lines'); 
        end

        switch line(1)
            case '#' % comment line - do nothing
                % use a big threshold on that length so that it goes faster by
                % skipping the strcmp on shorter lines
                if length(line)>100 && strcmp(line(1:5), '# typ');
                    header_line = line;
                end
            case '.' % info line
                fields = textscan(line, '%*s %*s %*s %*s %[^:] %*s %f%f');
                name = strtrim(fields{1}{:});
                val = fields{2};
                val2 = fields{3}; % only non-empty for a few field types
                info.(cleanFieldname(name)) = [val val2];
            case ' ' % data line - make different data structures for each type
                stop_condition = 1;
            otherwise
                disp('hmm.... should never have gone here');
                error('Call fred 803.983.9352');
        end

    end

    % read body of text once again, close file
    frewind(fid); % rewind file
    list_text = fread(fid,inf,'uint8=>char')'; 
    fclose(fid);
    
    % define & clean attr names (will be used as var names and fieldnames)
    attr = textscan(header_line,'%s','Delimiter',' ','MultipleDelimsAsOne',1);
    attr = attr{1}(2:end); % get rid of '#'
    for k=1:length(attr)
        attr{k}=cleanFieldname(attr{k});
    end
    
    % define a regexp pattern
    pattern = '(?<typ>\w+)'; % note that the first attr is right here!
    for k=2:length(attr),
        pattern = sprintf('%s\\s+(?<%s>-?\\d+)',pattern,attr{k});
    end

    % use regexp to match up the data lines (which match the defined pattern) and
    % put them into a struct
    % this is kind of a trick - instead of returning the indicies, the 'names'
    % argument tells it to return the actual string 
    lines_for_data = regexp(list_text, pattern, 'names');
    print_status_bar(1,length(fields)+2,20,'Processing .list file ');
    print_status_bar(2,length(fields)+2,20,'Processing .list file ');
    
    % free that significant (~70 MB) chunk of mem
    clear list_text 
    
    % drop n_a_ field (save memory and time when str2double comes in)
    lines_for_data = rmfield(lines_for_data,'n_a_');
    
    % the regexp command is quick, but left us with strings in numeric fields.
    % all of our numbers, which were matched as numbers in the
    % regexp come back as strings - which we need to put back as numbers
    % to fix this im going to change from the current struct array to a size==1
    % struct with longs arrays 
    % THIS TAKES MOST OF THE TIME - FIXME
    fields = fieldnames(lines_for_data);
    for ii = 1:length(fields);
        if strcmp(fields(ii),'typ');
            data_info_struct.(fields{ii}) = {lines_for_data.(fields{ii})};
        else % everything else is numeric
            % using tmp makes this process a little faster than the one-liner
            % version
            tmp = {lines_for_data.(fields{ii})};
            data_info_struct.(fields{ii}) = ...
                str2double(tmp);
        end
        print_status_bar(ii+2,length(fields)+2,20,'Processing .list file ');
    end

    % clear away that old struct (~550MB) new one is ~65MB
    clear lines_for_data tmp

    % get the different types of data
    % % assumes that the order is NOI>FRX>PHX>STD
    indNOI = find(strcmp('NOI',data_info_struct.typ));
    indFRX = find(strcmp('FRX',data_info_struct.typ));
    indPHX = find(strcmp('PHX',data_info_struct.typ));
    indSTD = find(strcmp('STD',data_info_struct.typ));
    indNAV = find(strcmp('NAV',data_info_struct.typ));

    % make substructs
    for ff = 1:length(fields)
        NOIstr.info.(fields{ff}) = data_info_struct.(fields{ff})(indNOI);
        FRXstr.info.(fields{ff}) = data_info_struct.(fields{ff})(indFRX);
        PHXstr.info.(fields{ff}) = data_info_struct.(fields{ff})(indPHX);
        STDstr.info.(fields{ff}) = data_info_struct.(fields{ff})(indSTD);
        NAVstr.info.(fields{ff}) = data_info_struct.(fields{ff})(indNAV);
    end


    % Open and read binary complex data
    fid = fopen(data_fname,'r');
    if fid~=-1
        [raw_data,raw_data_size] = fread(fid,inf,'float32');
    else
        error('read_raw_listdata:fileerror',...
            'cannot open %s for reading', list_fname );
    end
    fclose(fid);

    % Validate file size
    % Calculate the expected size of the file from the last row's offset and sz
    %   numbers
    % (last_offset + last_size)/FLOAT32_BYTES should = raw_data_size
    FLOAT32_BYTES = 4; % this is a constant and just here for clarity
    exp_size = (data_info_struct.offset(end) + ...
        data_info_struct.size(end))/FLOAT32_BYTES;
    if raw_data_size ~= exp_size
        error('read_raw_listdata:datafilesize',...
            '.data file size: expected %i bytes, read %i bytes',...
            exp_size,raw_data_size);
    end

    % done with a data_info_struct now (~60MB)
    clear data_info_struct

    save(tmp_fname);
    fprintf('Saved temporary file: %s\n', tmp_fname);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOW ACTUALLY DEAL WITH THE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General 
% these should always be there
X_resolution = info.X_resolution;
Y_resolution = info.Y_resolution;

% this stuff varies depending on whether we are reading raw_xxx.list/data or
% cpx_xxx.list/data
if ~exist('info.X_direction_SENSE_factor','var');
    X_direction_SENSE_factor=1;
    Y_direction_SENSE_factor=1;
else
    X_direction_SENSE_factor = info.X_direction_SENSE_factor;
    Y_direction_SENSE_factor = info.Y_direction_SENSE_factor;
end

if ~isfield(info,'kx_oversample_factor');
    kx_oversample_factor = 1;
    ky_oversample_factor = 1;
    kx_min = 0;
    kx_max = X_resolution-1;
    ky_min = 1;
    ky_max = Y_resolution;
    X_range = [0 X_resolution-1];
    iscpxfile = true; 
    fprintf('\nFile is detected to be cpx_XXX.list/.data file.\n');
else
    kx_oversample_factor = info.kx_oversample_factor;
    ky_oversample_factor = info.ky_oversample_factor;
    kx_min = info.kx_range(1);
    kx_max = info.kx_range(2);
    ky_min = info.ky_range(1);
    ky_max = info.ky_range(2);
    X_range = info.X_range;
    iscpxfile = false; 
    fprintf('\nFile is detected to be raw_XXX.list/.data file.\n');
end

% For each datatype the fields included are based on that chart at the end of
% this script, which shows which fields could ever be populated for each
% datatype
% I don't populate any that wont ever be populated
% I also don't populate all of the ones that could be poluated (possible FIXME)

%%%%%%%%%%%%%%% deal with NOI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    n.a.  aver  sign  rf    grad  
% --- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
% NOI                           *     *                                                      
if ~isempty(indNOI)
    % make a second version that is good for indexing
    ky_min = 0; % for NOI lines - this is definately non-trivial for STD lines
    tstr = fix_str_for_indexing(NOIstr.info,ky_min); % temporary struct
    NOIstr.info_adjusted = tstr;
    % preallocate data array
    loca = unique(tstr.loca); num_loca = length(loca);
    chan = unique(tstr.chan); num_chan = length(chan);
    num_kx = NOIstr.info.size(1)/FLOAT32_BYTES/2; 
    NOIstr.data = complex(nan(num_kx, num_loca, num_chan));
    NOIstr.dim_labels = {'kx','loca','chan'};
    for ii=1:length(indNOI)
        beg = NOIstr.info.offset(ii)/FLOAT32_BYTES + 1;
        stp = (NOIstr.info.offset(ii) + NOIstr.info.size(ii))/FLOAT32_BYTES;
        tmp = raw_data(beg:stp);
        tmp = tmp(1:2:end) + 1i*tmp(2:2:end); % real and imag are interleaved 
        NOIstr.data(:,tstr.loca(ii), tstr.chan(ii)) = tmp;
        print_status_bar(ii,length(indNOI),40,'Processing NOI data ');
    end
    clear tstr    
end

%%%%%%%%%%%%%%% deal with FRX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    n.a.  aver  sign  rf    grad  
% --- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
% FRX   *                 *     *     *                                         *            
if ~isempty(indFRX)
    ky_min = 0; % also zero in the FRX case
    tstr = fix_str_for_indexing(FRXstr.info,ky_min); % temporary struct
    FRXstr.info_adjusted = tstr;
    mix = unique(tstr.mix);   num_mix = length(mix);
    echo = unique(tstr.echo);   num_echo = length(echo);
    loca = unique(tstr.loca); num_loca = length(loca);
    chan = unique(tstr.chan); num_chan = length(chan);
    sign = unique(tstr.sign); num_sign = length(sign);
    num_kx = tstr.size(1)/FLOAT32_BYTES/2; % this may not be the be the best way 
    FRXstr.data = complex(nan(num_kx,num_mix, num_echo,num_loca,num_chan,num_sign));
    FRXstr.dim_labels = {'kx','mix','echo','loca','chan','sign'};
    for ii=1:length(indFRX)
        beg = FRXstr.info.offset(ii)/FLOAT32_BYTES + 1;
        stp = (FRXstr.info.offset(ii) + FRXstr.info.size(ii))/FLOAT32_BYTES;
        tmp = raw_data(beg:stp);
        tmp = tmp(1:2:end) + 1i*tmp(2:2:end); % real and imag are interleaved 
        FRXstr.data(:,tstr.mix(ii), tstr.echo(ii),tstr.loca(ii), tstr.chan(ii), ....
            tstr.sign(ii)) = tmp;
        print_status_bar(ii,length(indFRX),40,'Processing FRX data ');
    end
    clear tstr
end
%%%%%%%%%%%%%%% deal with PHX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    n.a.  aver  sign  rf    grad  
% --- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
% PHX   *                 *     *     *                                         *     *     *
if ~isempty(indPHX)
    ky_min = info.ky_range(1); % NOT TRIVIAL (FIXME - is this the general case?)
    tstr = fix_str_for_indexing(PHXstr.info,ky_min); % temporary struct
    PHXstr.info_adjusted = tstr;
    num_kx = tstr.size(1)/FLOAT32_BYTES/2; %FIXME; can be buggy is partial echo aq'd
    mix = unique(tstr.mix);   num_mix = length(mix);
    echo = unique(tstr.echo);   num_echo = length(echo);
    loca = unique(tstr.loca); num_loca = length(loca);
    chan = unique(tstr.chan); num_chan = length(chan);
    sign = unique(tstr.sign); num_sign = length(sign);
    grad = unique(tstr.grad); num_grad = length(grad);
    PHXstr.data = complex(nan(num_kx,num_mix,num_echo,num_loca,num_chan,...
        num_sign,num_grad));
    PHXstr.dim_labels = {'kx','mix','echo','loca','chan','sign','grad'};
    for ii=1:length(indPHX)
        beg = PHXstr.info.offset(ii)/FLOAT32_BYTES + 1;
        stp = (PHXstr.info.offset(ii) + PHXstr.info.size(ii))/FLOAT32_BYTES;
        tmp = raw_data(beg:stp);
        tmp = tmp(1:2:end) + 1i*tmp(2:2:end); % real and imag are interleaved 
        PHXstr.data(:,tstr.mix(ii),tstr.echo(ii), tstr.loca(ii), tstr.chan(ii),...
            tstr.sign(ii), tstr.grad(ii)) = tmp;
        print_status_bar(ii,length(indPHX),40,'Processing PHX data ');
    end
    clear tstr
end



%%%%%%%%%%%%%%% deal with STD (the real data) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    n.a.  aver  sign  rf    grad  
% --- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
% STD   *     *     *     *     *     *     *     *     *     *           *     *     *     *
% make a second version that is good for indexing and preallocation
% read the bit about Fourier Transofrmation in the long comment block below

% in the raw_xxx case, this variable is appropriate, in the cpx_xxx case, they
% clearly are not fourier transform lengths, so the name is a little silly. in
% either case, they exist just for size checking and pre-allocation
fftx_size = X_resolution *  kx_oversample_factor / X_direction_SENSE_factor;
ffty_size = Y_resolution * ky_oversample_factor / Y_direction_SENSE_factor;

if ~iscpxfile
    info.fftx_size__calc = fftx_size;
    info.ffty_size__calc = ffty_size;
end

% make sure the min/max kx/ky fits in that fftx_size
if ~( (kx_max - kx_min +1) <= fftx_size)
    error('read_raw_listdata:bad_fftx_math',...
        'kx data indicies larger than expected. check code.');
end
if ~( (ky_max - ky_min +1) <= ffty_size)
    error('read_raw_listdata:bad_ffty_math',...
        'ky data indicies larger than expected. check code.');
end

if iscpxfile
    kx_offset=0;
    ky_offset=0;
else
    kx_offset = -floor(fftx_size/2); % subtract this the kx values
    ky_offset = -floor(ffty_size/2); % used in the fix_str_for_indexing to do an
                                     % an offset correction
    info.kx_offset__calc = kx_offset;
    info.ky_offset__calc = ky_offset;
end


% calculate the indicies that are needed to extract the smaller FOV from teh
% oversampled read-out direction
x_lims = [X_range(1)-kx_offset+1, X_range(2)-kx_offset+1];

tstr = fix_str_for_indexing(STDstr.info,ky_offset); % temporary struct
STDstr.info_adjusted = tstr;

if iscpxfile
    % cheat on names here and use "kz" instead of rewriting other parts
    tstr.aver = ones(size(tstr.mix));
    tstr.kz = tstr.z;
    tstr.ky = tstr.y;
    STDstr.dim_labels = {'x','y','z','mix','dyn','card','echo','loca','chan',...
        'extr_combined'};
else % only present in raw_xxx
    STDstr.dim_labels = {'kx','ky','kz','mix','dyn','card','echo','loca','chan',...
        'extr_combined','aver'};
end
% these fields are in both 
kz = unique(tstr.kz); num_kz = length(kz); 
aver = unique(tstr.aver); num_aver = length(aver); 
mix = unique(tstr.mix);   num_mix = length(mix);
dyn = unique(tstr.dyn);   num_dyn = length(dyn);
card = unique(tstr.card);   num_card = length(card);
echo = unique(tstr.echo);   num_echo = length(echo);
loca = unique(tstr.loca); num_loca = length(loca);
chan = unique(tstr.chan); num_chan = length(chan);
extr = unique(tstr.extr_combined); num_extr = length(extr);
STDstr.info.extr_combined = tstr.extr_combined - 1; % -1 to put in 0 reference
% pre-allocate
% sign, grad, are not used and
% extr1/2 are the same dimension (making extr_combined)
% .data is uncorrected
STDstr.data = complex(nan(fftx_size, ffty_size, num_kz,...
    num_mix,num_dyn,num_card,num_echo,num_loca,num_chan, num_extr,num_aver));

if ~iscpxfile
    % if we are reading raw_xxx.list/.data then,
    % .data_corrected has FRX/PHX corrections and is the final size in readout
    STDstr.data_corrected = complex(nan(info.X_resolution,ffty_size, num_kz,...
        num_mix,num_dyn,num_card,num_echo,num_loca,num_chan, num_extr,num_aver));
end
    
for ii=1:length(indSTD)
    beg = STDstr.info.offset(ii)/FLOAT32_BYTES + 1;
    stp = (STDstr.info.offset(ii) + STDstr.info.size(ii))/FLOAT32_BYTES;
    tmp = raw_data(beg:stp);
    tmp = tmp(1:2:end) + 1i*tmp(2:2:end); % real and imag are interleaved 
    STDstr.data((kx_min:kx_max)-kx_offset+1, tstr.ky(ii),tstr.kz(ii),...
        tstr.mix(ii),tstr.dyn(ii),tstr.card(ii),tstr.echo(ii),tstr.loca(ii),...
        tstr.chan(ii),tstr.extr_combined(ii),tstr.aver(ii)) = tmp;
    
    if ~iscpxfile
        % create the corrected data

        % fft the tmp_line (readout direction fft)
        tmp = ft(tmp);
        % reduce the FOV to just the expected
        tmp = tmp(x_lims(1):x_lims(2));

        % find the right FRX correction line (line will get accessed many time)
        FRXline = FRXstr.data(:,tstr.mix(ii), tstr.echo(ii), ...
            tstr.loca(ii), tstr.chan(ii), tstr.sign(ii) );
        % find the right PHX correction line (line will get accessed many time)
        PHXline = PHXstr.data(:,tstr.mix(ii), tstr.echo(ii), ...
            tstr.loca(ii), tstr.chan(ii), tstr.sign(ii),:);  

        % apply corrections
        tmp = tmp.*FRXline(:).*PHXline(:);
        tmp = ift(tmp); % ift to k-space to avoid the half-way-fft'd wierdness

        % put corrected data in corrected struct field
        STDstr.data_corrected(:, tstr.ky(ii),tstr.kz(ii),tstr.mix(ii),...
            tstr.dyn(ii),tstr.card(ii),tstr.echo(ii),tstr.loca(ii),...
            tstr.chan(ii),tstr.extr_combined(ii),tstr.aver(ii)) = tmp;
    end
    
    print_status_bar(ii,length(indSTD),40,'Processing STD data '); 
end

% currently in-slice data that was uncollected due to half-scan or partial echo
% is due to the initialization NaN, when it really should be 0. That's fixed
% here.
STDstr.data = convert_in_slice_nans_to_zeros(STDstr.data);
if ~iscpxfile
    STDstr.data_corrected = convert_in_slice_nans_to_zeros(STDstr.data_corrected);
end

if any(isnan(STDstr.data(:)))
    fprintf('\n');% newline
    disp(...
        'There are some NaNs left in STDstr.data, indicating empty dimensions.');
    disp('This can happen if one dynamic has two averages while another has');
    disp('only one, or something like that so that the data doesn''t fit');
    disp('perfectly into the N-dimensional hyper-rectangular matrix');
    fprintf('\n');% newline
    
end

%%%%%%%%%%%%%%% deal with NAV data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    n.a.  aver  sign  rf    grad  
% --- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
% NAV   *     *     *     *     *     *     *     *     *     *           *     *     *     *
if ~isempty(indNAV)
    kx_offset = -floor(fftx_size/2) ; % subtract this the kx values
    ky_offset = min(NAVstr.info.ky); % used in the fix_str_for_indexing to do an 
                                      % an offset correction
    tstr = fix_str_for_indexing(NAVstr.info,ky_offset); % temporary struct
    NAVstr.info_adjusted = tstr;
    kz = unique(tstr.kz); num_kz = length(kz); 
    ky = unique(tstr.ky); num_ky = length(ky);
    mix = unique(tstr.mix);   num_mix = length(mix);
    dyn = unique(tstr.dyn);   num_dyn = length(dyn);
    card = unique(tstr.card);   num_card = length(card);
    echo = unique(tstr.echo);   num_echo = length(echo);
    loca = unique(tstr.loca); num_loca = length(loca);
    chan = unique(tstr.chan); num_chan = length(chan);
    extr = unique(tstr.extr_combined); num_extr = length(extr);
    aver = unique(tstr.aver); num_aver = length(aver);
    sign = unique(tstr.aver); num_sign = length(sign);
    NAVstr.info.extr_combined = tstr.extr_combined - 1;
    NAVstr.dim_labels = {'kx','ky(modified)','kz','mix','dyn','card','echo',...
        'loca','chan','extr_combined','aver','sign'};
    NAVstr.data = complex(nan(fftx_size, num_ky, num_kz,num_mix,num_dyn,num_card,...
        num_echo,num_loca,num_chan, num_extr,num_aver,num_sign));
    % grad is not used
    % extr1/2 are the same dimension (extr_combined)
    for ii=1:length(indNAV)
        beg = NAVstr.info.offset(ii)/FLOAT32_BYTES + 1;
        stp = (NAVstr.info.offset(ii) + NAVstr.info.size(ii))/FLOAT32_BYTES;
        tmp = raw_data(beg:stp);
        tmp = tmp(1:2:end) + 1i*tmp(2:2:end); % real and imag are interleaved 
        NAVstr.data((kx_min:kx_max)-kx_offset+1, tstr.ky(ii),tstr.kz(ii),tstr.mix(ii),...
            tstr.dyn(ii),tstr.card(ii),tstr.echo(ii),tstr.loca(ii),tstr.chan(ii),...
            tstr.extr_combined(ii),tstr.aver(ii),tstr.aver(ii)) = tmp;
        print_status_bar(ii,length(indNAV),40,'Processing NAV data ');
    end
    NAVstr.data = convert_in_slice_nans_to_zeros(NAVstr.data);
end

if ~iscpxfile
    %%%%%%%%%%%%%%% make shot mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assumes that the shots collected the same space in each slice
    % this is useful so that we can figure out how to take apart the k-space  
    % data in reconstruction algorithms that do something on a per-shot basis
    if (ky_max - ky_min +1) ~= info.number_of_gradient_echoes;
        fprintf('Number of ky ~= to EPI factor. Must be multi-shot.\n');
        fprintf('Creating shot mask\n');

        % get unique ky-gradient number pairs
        tmp = unique([tstr.ky(:), tstr.grad(:)],'rows');
        % this will be as long as the number of ky lines (assumes that ky lines
        % are collected on the same gradient "pass" in every slice/echo/dyn/etc)
        if size(tmp,1) ~= length(unique(tstr.ky))
            error('read_raw_listdata:bad_shot_mask_logic',...
                'Logic is faulty/not robust enough for multishot dti.');
        end    
        shot_mask = zeros(ffty_size,1);
        shot_mask(tmp(:,1)) = tmp(:,2); % one col is ky ind,  other is grad ind
        shot_mask = repmat(shot_mask',[fftx_size,1]);

        % two ways to figure out num_shots
        num_shots =...
            round((info.ky_range(2)-info.ky_range(1)+1)/...
            info.number_of_gradient_echoes);
        num_shots2 = sum(tmp(:,2)==1);
        if num_shots ~= num_shots2
            error('read_raw_listdata:bad_shot_mask_logic', ...
                'Logic is faulty for calculating number of shots');
        end
        shot_mask_by_shot = repmat(1:num_shots,[1 ceil(ffty_size/num_shots)]); 
        start = find(shot_mask_by_shot>0,1); % may have uncollected lines at beginning
        % zero any places that were uncollected
        % (get this info from first line of shot mask)
        shot_mask_by_shot = shift_array(shot_mask_by_shot,start).*(shot_mask(1,:)>0);

        % make it a full matrix
        shot_mask_by_shot = repmat(shot_mask_by_shot,[fftx_size,1]);

    end

    %%%%%%%%%%%%%%% make sign mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assumes that this is constant for all slices acq'd
    tmp = unique([tstr.ky(:),STDstr.info.sign(:)],'rows');
    if size(tmp,1) == length(unique(tstr.ky))
        sign_mask = zeros(ffty_size,1);
        sign_mask(tmp(:,1)) = tmp(:,2); 
        sign_mask = repmat(sign_mask',[fftx_size,1]);
    elseif size(tmp,1) > length(unique(tstr.ky))

        % lets see if we can figure out why it's not consistent in all cases
        fields = {'mix','dyn','card','echo','extr_combined','aver'};

        % ditch singleton dimensions for clarity
        keep = zeros(length(fields),1);
        for ii=1:length(fields);
            if length(unique(tstr.(fields{ii})))>1
                keep(ii)=1;
            end
        end
        fields = fields(logical(keep)); % these are the non-singleton dims

        tmp = [tstr.ky(:) STDstr.info.sign(:)];
        for ii = 1:length(fields)
            tmp = [tmp tstr.(fields{ii})(:)];
        end
        tmp = unique(tmp,'rows');
        sign_mask.ky_line_sign_combos = tmp;
        sign_mask.fields_shown = [{'ky','sign'},fields];

        fprintf('\n');
        disp('Some lines are covered with both positive and negative gradients.');
        disp('Current logic for creation of sign mask isn''t ready for that. ');
        disp('A fully-usuable Sign mask is not being created');
        disp('However, the information to make your own sign mask is provided');
        disp('in the data bundle. This shows you the sign for each k-space line');
        disp('for every unique case. The data dimensions considered for these');
        fprintf('combinations is given as .fields_shown\n\n');
    else
        error('read_raw_listdata:bad_sign_mask_logic',...
            'Logic is faulty/not robust for creation of sign mask.');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CREAT OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% core stuff
if ~iscpxfile
    data_bundle.NOIstr = NOIstr;
    data_bundle.FRXstr = FRXstr;
    data_bundle.PHXstr = PHXstr;
    data_bundle.NAVstr = NAVstr;
end
% both cpx/raw .list/.data
data_bundle.STDstr = STDstr;

% other stuff
if exist('shot_mask','var');
    data_bundle.shot_mask = shot_mask;
    data_bundle.shot_mask_by_shot = shot_mask_by_shot;
end
if exist('sign_mask','var');
    data_bundle.sign_mask = sign_mask;
end

tmp = whos('data_bundle'); 
fprintf('\nOutput data is %.2f MB\n\n',tmp.bytes/2^20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = convert_in_slice_nans_to_zeros(in)
% assumes first two dimensions are slice dimensions
dsz = size(in);
reshaped = reshape(in,[dsz(1)*dsz(2) prod(dsz(3:end))]);
for ii = 1:size(reshaped,2);
    nan_mat = isnan(reshaped(:,ii));
    % has both nan and non-nan elements
    if sum(nan_mat(:)) < numel(nan_mat) && sum(nan_mat(:)) > 0 
        reshaped(nan_mat,ii) = 0;
    end
end
out = reshape(reshaped,dsz);

function out = shift_array(in, shift)
% shift is an integer. A positive number shifts the array right 
% shift of zero does nothing
% new spots are zero padded
% for circshift use circshift
sz = size(in);
in = in(:); % vectorize
if abs(shift)>=length(in)
    out = zeros(size(in));
    return;
end
if shift>=0
    ind = 1:length(in)-shift;
    out = reshape([zeros(shift,1);in(ind)],sz);
else
    ind = abs(shift)+1:length(in);
    out = reshape([in(ind); zeros(abs(shift),1)],sz);
end
    
function clean_str = fix_str_for_indexing(improper_str,ky_offset)
% this function creates a struct with vectors that can be used to index into a
% valid data matrix, fixing things like zero/negative-indexing into arrays

fields = fieldnames(improper_str);

for ff = 1:length(fields)
    switch fields{ff}
        case {'typ','size'} % keep but don't touch these
            clean_str.(fields{ff}) = improper_str.(fields{ff});
        case 'sign' % make this an array 1:1 -or- 1:2 (if differing signs presnet)
            [~,~,IC] = unique(improper_str.sign);
            clean_str.sign = IC;
        case 'chan' % these are kind of meaningless numbers, so replace them
            [~,~,IC] = unique(improper_str.chan);
            clean_str.chan = IC;
        case 'ky'; % ky is the only one with negative values
            clean_str.ky = improper_str.ky - ky_offset + 1;
        otherwise
            clean_str.(fields{ff}) = improper_str.(fields{ff}) + 1;
    end
end

% create a new field that combines extr1 and extr2
[~,~,IC] = ...
    unique([improper_str.extr1(:),improper_str.extr2(:)],'rows');
clean_str.extr_combined = IC;

function [s] = cleanFieldname(s)

illegal_chars = {...
    '+','-','*','.',...
    '^','\','/','.',...
    '=','~','<','>',...
    '&','|',':',';',...
    '(',')','{','}',...
    '[',']','{','}',...
    '''','%',' ','!', ...
    '@','#','$','`',...
    '?',',','"',...
    };

general_replacement_char = '_';
firstchar_replacement_char = 'x'; % cannot be an underscore

for k=1:length(illegal_chars),
    s = strrep(s,illegal_chars{k},general_replacement_char);
end

% first character cannot be a number
firstchar_code = double(s(1));
if ( (firstchar_code>=double('0')) && (firstchar_code<=double('9')) )
    s(1) = firstchar_replacement_char;
end

% first character cannot be an underscore
if(s(1)=='_'),
    s(1) = firstchar_replacement_char;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INTERESTING INFORMATION* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *copied from inside a .list file
%
% Fourier transformations:
% ------------------------
% The Fourier transformation lengths for STD, NAV and REJ data of a SENSE scan are (by definition) equal to
% reconstruction matrix length * k-space oversample sample factor / SENSE factor
% in each applicable direction.
%
% The centre point of the k-space must be positioned in the middle of the interval with this length,
% (or just right from the middle if the Fourier transformation length is even), according to the
% specified k-space coordinate ranges. Apply a cyclic wrap-around to data that extends outside the
% interval, and zerofill missing data before Fourier transformation.
% After Fourier transformation, the imaging space coordinate ranges specify the position, length 
% and direction of the relevant part of the frequency spectrum. Apply a cyclic wrap-around if
% needed. If start value > end value, read the spectrum part between these positions in reverse order.
%
% Corrections:
% ------------
% Depending on the type of scan, one or more of the following corrections may have to be applied:
% - frequency spectrum correction (using FRX data)
% - EPI/GraSE phase correction    (using PHX data)
% - rf echo amplitude correction  (using PHX data, only for GraSE scans)
%
% These corrections can be applied to STD, NAV and REJ data vectors after forward Fourier transformation
% in the measurement direction, by multiplying the selected part of the frequency spectrum of
% these data vectors with the PHX and FRX correction data vector that has corresponding identifying
% complex data vector attributes.
% 
% Remarks:
% - In GraSE scans you may get a PHX correction data vector for each rf echo. These correction vectors
%   differ only a scale factor from each other, and can be used for rf echo amplitude correction.
%   If you don't want this correction, use the correction data vectors with rf echo = 0 only.
% - Depending on the type of scan, several of the FRX correction data vectors may be identical.
%
% Identifying attributes of complex data vectors:
% -----------------------------------------------
% The next table specifies the identifying attributes for each type of complex data vector:
%
% typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    n.a.  aver  sign  rf    grad  
% --- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
%
% STD   *     *     *     *     *     *     *     *     *     *           *     *     *     *
% REJ   *     *     *     *     *     *     *     *     *     *           *     *     *     *
% PHX   *                 *     *     *                                         *     *     *
% FRX   *                 *     *     *                                         *            
% NOI                           *     *                                                      
% NAV   *     *     *     *     *     *     *     *     *     *           *     *     *     *


% BONUS/UNUSED/SAVE_FOR_LATER CODE



% 
% 
%     attr = {'mix','dyn','card','echo','loca','chan','extr1','extr2','ky','kz',...
%         'n.a.','aver','sign','rf','grad','enc','rtop','rr','size','offset'};
%     for k=1:length(attr),
%         attr{k} = cleanFieldname( attr{k} );
%     end



