function [lstparms, DVattribs] = read_listfile_either_fwb(file_name)
%--------------------------------------------------------------------------
% read all the lines (general info, comments, data vector) of the LIST file
% for complex data
% This file is slightly different from LIST file for raw data
% Input: 
%             *.list
% Outputs:
%             lstparms: a 1*1 structure containing all the general information
%             and data vector, e.g., lstparms.label(NV,:),
%             lstparms.attrib(NV,:), lstparms.(fieldname)
%             DVattribs: a 1*1 cell containing 22*1 cell array
%                        '#', 'tpy', 'mix', 'dyn'....'offset'
%                         # typ mix   dyn   card  echo  loca  chan  extr1 
%                         extr2 y     z     n.a.  n.a.  n.a.  n.a.  n.a.  
%                         n.a.  n.a.  n.a.  size  offset
% Example:
% [lstparms, DVattribs] = read_listfile_complex('cpx_040(1028).list');
%==========================================================================
% Author        Hui Wang, Ph.D student
%               University of Louisville, KY
%               12/20/2010
%==========================================================================
fprintf('Reading .LIST file ...');

% stuff to determine whether to deal w/ cpx or raw data (changes template below)
ftype = '';
if strcmp(file_name(1:3),'cpx');
    ftype = 'cpx'; fprintf(' cpx');
elseif strcmp(file_name(1:3),'raw');
    ftype = 'raw'; fprintf(' raw');
end

while ~any(strcmp(ftype,{'cpx','raw'})) && ~isequal(ftype, 0);
    fprintf('Could not determine wheter file is cpx_XXX.list or raw_XXX.list \n style from file name.\n');
    ftype = input('Which is it [''raw'' or ''cpx'']? ');
end
if ftype == 0;
    disp('0 input, exiting');
    return;
end

nlines = 0; NG=0; NC=0; NV=0;
fid = fopen(file_name,'r');
if (fid < 1), error(['.LIST file ', file_name, ' not found.']); end;
% Scan the .list file line by line, lines (a 1*nlines cell array) stores
% every line. geninfo (a 1*NG cell array) stores all the general information
% comment (a 1*NC cell array) stores all comments. Each cell in those cell
% array is one line from the .list file. lstparms.lable (a char array) stores
% the labels for the data vector, lstparms.attrib (a double array) stores
% the attributes of the data vector. 
while ~feof(fid)
    curline = fgetl(fid);
    if ~isempty(curline)
        nlines = nlines + 1;  % Recording number of lines
        lines(nlines) = cellstr(curline); % Allows variable line size
        if (strmatch('# typ mix',curline))
            DVattribs = textscan(curline,'%s'); % Data vector attributes
        end
        % Geninfo, Comments and data vectors
        curline = char(lines(nlines)); firstc = curline(1);
        if (size(curline,2) > 0)
            switch (firstc)
                case '.'
                    NG = NG + 1; geninfo(NG) = lines(nlines); % Number of general information.
                case '#'
                    NC = NC + 1; comment(NC) = lines(nlines); % Number of comments.
                otherwise
                    NV = NV + 1;                              % Number of data vector.
                    [token, remainder] = strtok(curline);
                    lstparms.lable(NV,:) = token; lstparms.attrib(NV,:) = str2num(remainder);
            end
        end
    end
end
fclose(fid);
if (NG < 1 || NC < 1), error('.LIST file has invalid format'); end;


% Get LIST information template
% this is the only difference between the two (cpx_XXX.list and regualaw raw_XXX.list)
if strcmp(ftype,'cpx')
    info_template = get_info_template_complex;
else % must be raw
    info_template = get_info_template;
end
% parse the data information, i.e., extract the general information stored
% in geninfo and write them into lstparms.(fieldname). The fieldname is
% defined in get_info_template_complex by the programmer
for S=1:size(info_template,1)
    line_key = char(info_template(S,1));
    value_type = char(info_template(S,2));
    field_name = char(info_template(S,3));
    L = strmatch(line_key,geninfo);
    if ~isempty(L)
        curline = char(geninfo(L));
        value_start = 1 + strfind(curline,':');
        value_end = size(curline,2);
    else
        value_type = ':-( VALUE NOT FOUND )-:';
    end
    switch (value_type)
        case {'float scalar' 'int   scalar' 'float vector' 'int   vector'}
            lstparms.(field_name) = str2num(curline(value_start:value_end));
        case {'char  scalar' 'char  vector'}
            lstparms.(field_name) = strtrim(curline(value_start:value_end));
        otherwise
            lstparms.(field_name) = '';
    end
end
fprintf(' complete.\n');


function [info_template] = get_info_template_complex    
% General information template from exported LIST file for complex data
% Purpose: for parsing the general information
% Modified from Amol's codes, hui
% Compared to the template for raw LIST file, it doesn't have the
% following: k-space coordinate ranges, k-space oversample factors
% imaging space coordinate ranges.
% Output: info_template: a 16*3 cell array. 
info_template = { ...
    '.    0    0    0  number_of_mixes                    :' 'int   scalar' 'number_of_mixes'; ...
    '.    0    0    0  number_of_encoding_dimensions      :' 'int   scalar' 'number_of_encoding_dimensions'; ...
    '.    0    0    0  number_of_dynamic_scans            :' 'int   scalar' 'number_of_dynamic_scans'; ...
    '.    0    0    0  number_of_cardiac_phases           :' 'int   scalar' 'number_of_cardiac_phases'; ...
    '.    0    0    0  number_of_echoes                   :' 'int   scalar' 'number_of_echoes'; ...
    '.    0    0    0  number_of_locations                :' 'int   scalar' 'number_of_locations'; ...
    '.    0    0    0  number_of_extra_attribute_1_values :' 'int   scalar' 'number_of_extra_attribute_1_values'; ...
    '.    0    0    0  number_of_extra_attribute_2_values :' 'int   scalar' 'number_of_extra_attribute_2_values'; ...
    '.    0    0    0  number_of_signal_averages          :' 'int   scalar' 'number_of_signal_averages'; ...
    '.    0    0    0  number of coil channels            :' 'int   scalar' 'number_of_coil_channels'; ...
    '.    0    0    0  X-resolution                       :' 'int   scalar' 'X_resolution'; ...
    '.    0    0    0  Y-resolution                       :' 'int   scalar' 'Y_resolution'; ...
    '.    0    0    0  Z-resolution                       :' 'int   scalar' 'Z_resolution'; ...
    '.    0    0    0  X-direction SENSE factor           :' 'float scalar' 'X_direction_SENSE_factor'; ...
    '.    0    0    0  Y-direction SENSE factor           :' 'float scalar' 'Y_direction_SENSE_factor'; ...
    '.    0    0    0  Z-direction SENSE factor           :' 'float scalar' 'Z_direction_SENSE_factor'; ...
    };
return;


function [info_template] = get_info_template    % Information from exported LIST file
info_template = { ...
    '.    0    0    0  number_of_mixes                    :' 'int   scalar' 'number_of_mixes'; ...
    '.    0    0    0  number_of_encoding_dimensions      :' 'int   scalar' 'number_of_encoding_dimensions'; ...
    '.    0    0    0  number_of_dynamic_scans            :' 'int   scalar' 'number_of_dynamic_scans'; ...
    '.    0    0    0  number_of_cardiac_phases           :' 'int   scalar' 'number_of_cardiac_phases'; ...
    '.    0    0    0  number_of_echoes                   :' 'int   scalar' 'number_of_echoes'; ...
    '.    0    0    0  number_of_locations                :' 'int   scalar' 'number_of_locations'; ...
    '.    0    0    0  number_of_extra_attribute_1_values :' 'int   scalar' 'number_of_extra_attribute_1_values'; ...
    '.    0    0    0  number_of_extra_attribute_2_values :' 'int   scalar' 'number_of_extra_attribute_2_values'; ...
    '.    0    0    0  number_of_signal_averages          :' 'int   scalar' 'number_of_signal_averages'; ...
    '.    0    0    0  number of coil channels            :' 'int   scalar' 'number_of_coil_channels'; ...
    '.    0    0    0  kx_range                           :' 'int   vector' 'kx_range'; ...
    '.    0    0    0  ky_range                           :' 'int   vector' 'ky_range'; ...
    '.    0    0    0  kz_range                           :' 'int   vector' 'kz_range'; ...
    '.    0    0    0  kx_oversample_factor               :' 'float scalar' 'kx_oversample_factor'; ...
    '.    0    0    0  ky_oversample_factor               :' 'float scalar' 'ky_oversample_factor'; ...
    '.    0    0    0  kz_oversample_factor               :' 'float scalar' 'kz_oversample_factor'; ...
    '.    0    0    0  X-resolution                       :' 'int   scalar' 'X_resolution'; ...
    '.    0    0    0  Y-resolution                       :' 'int   scalar' 'Y_resolution'; ...
    '.    0    0    0  Z-resolution                       :' 'int   scalar' 'Z_resolution'; ...
    '.    0    0    0  X-direction SENSE factor           :' 'float scalar' 'X_direction_SENSE_factor'; ...
    '.    0    0    0  Y-direction SENSE factor           :' 'float scalar' 'Y_direction_SENSE_factor'; ...
    '.    0    0    0  Z-direction SENSE factor           :' 'float scalar' 'Z_direction_SENSE_factor'; ...
    '.    0    0    0  X_range                            :' 'int   vector' 'X_range'; ...
    '.    0    0    0  Y_range                            :' 'int   vector' 'Y_range'; ...
    '.    0    0    0  Z_range                            :' 'int   vector' 'Z_range'; ...
    '.    0    0    0  0th_order_phase_error_X            :' 'float scalar' 'zero_order_phase_error_X'; ...
    '.    0    0    0  1st_order_phase_error_X            :' 'float scalar' 'first_order_phase_error_X'; ...
    };
return;
