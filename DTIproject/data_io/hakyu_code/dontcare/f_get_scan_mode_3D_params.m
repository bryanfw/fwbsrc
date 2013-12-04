
function scan_mode_3D_related_params = f_get_scan_mode_3D_params(fname)
%[f_get_scan_mode_3D_params] 3-D scan mode related parameters, mostly kz
%directional parameters.
%
% USAGE:
%   scan_mode_3D_related_params = f_get_scan_mode_3D_params(fname)
%
% INPUT:
%   fname: Raw data name, e.g. 'raw_006'
%
% OUTPUT:
%   kz_range
%   kz_oversample_factor
%   Z-resolution
%   Z-direction SENSE factor
%
% NOTE:
%   This routine is reading kz directional parameters when scan mode of 3-D
%   is used.
%
%
% Last modified
% 2012.02.01.
%   Generated to get 3-D scan mode related kz directional parameters.
% 2012.07.23.
%   Add ky and kx oversample factor.
%
% Ha-Kyu



%% Preliminary
if ~isempty(regexp(fname,'(\\)', 'once')) && ...
        isempty(regexp(fname,'(\.list)', 'once')) % fullpath w/o .list
    %ind = regexp(fname,'(\\.)'); 
    f_s = [fname '.list'];
elseif ~isempty(regexp(fname,'(\\)', 'once')) && ...
        ~isempty(regexp(fname,'(\.list)', 'once')) % fullpath w/ .list
    %ind = regexp(fname,'(\\+\w+\.list)$'); 
    f_s = fname;
elseif isempty(regexp(fname,'(\\)', 'once')) && ...
        isempty(regexp(fname,'(\.list)', 'once')) % only raw_###
    f_s = [fname '.list'];
else % raw_###.list
    f_s = fname;
end


%% Main
scan_mode_3D_related_params = []; % struct output
%search_params = {'kz_range','kz_oversample_factor','Z-resolution','Z-direction SENSE factor'};

tic
fid=fopen(f_s);
while 1    
    tline = fgetl(fid);
    
    if ~ischar(tline), break, end
    
    % Find kz_range
    idx = strfind(tline,'kz_range');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        kz_range = str2num(s);
    end
    
    % Find kx_oversample_factor
    idx = strfind(tline,'kx_oversample_factor');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        kx_oversample_factor = str2num(s);
    end
    
    % Find ky_oversample_factor
    idx = strfind(tline,'ky_oversample_factor');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        ky_oversample_factor = str2num(s);
    end
    
    % Find kz_oversample_factor
    idx = strfind(tline,'kz_oversample_factor');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        kz_oversample_factor = str2num(s);
    end
    
    % Find Z-resolution
    idx = strfind(tline,'Z-resolution');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        Z_resolution = str2num(s);
    end
    
    % Find Z-direction SENSE factor
    idx = strfind(tline,'Z-direction SENSE factor');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        Z_direction_SENSE_factor = str2num(s);
    end
    
    % Finish if all search_params were found.
    if exist('kz_range','var') && exist('kz_oversample_factor','var') && ...
            exist('Z_resolution','var') && exist('Z_direction_SENSE_factor','var')
        break
    end
    
end
fclose(fid);

if ~isempty(kz_range)
    scan_mode_3D_related_params.kz_range = kz_range;
    clear  kz_range
end
if ~isempty(kz_oversample_factor)
    scan_mode_3D_related_params.kz_oversample_factor = kz_oversample_factor;
    clear  kz_oversample_factor
end
if ~isempty(ky_oversample_factor)
    scan_mode_3D_related_params.ky_oversample_factor = ky_oversample_factor;
    clear  ky_oversample_factor
end
if ~isempty(kx_oversample_factor)
    scan_mode_3D_related_params.kx_oversample_factor = kx_oversample_factor;
    clear  kx_oversample_factor
end
if ~isempty(Z_resolution)
    scan_mode_3D_related_params.Z_resolution = Z_resolution;
    clear  kz_range
end
if ~isempty(Z_direction_SENSE_factor)
    scan_mode_3D_related_params.Z_direction_SENSE_factor = Z_direction_SENSE_factor;
    clear  kz_range
end

t = toc;
fprintf('kz directional scan parameters are read in %f sec\n',t)
fprintf('\n')















