
function params = f_get_params_from_listflie(fname)
%[f_get_params_from_listflie] reads scan parameters from .LIST file.
%
% USAGE:
%   params = f_get_params_from_listflie(fname)
%
% INPUT:
%   fname: Raw data name, e.g. 'raw_006'
%
%
% Last modified
% 2012.07.23.
%   Read scan parameters from .LIST file. This is generated from
%   [f_get_scan_mode_3D_params.m].
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
params = []; % struct output
%search_params = {'kz_range','kz_oversample_factor','Z-resolution','Z-direction SENSE factor'};

tic
fid=fopen(f_s);
while 1    
    tline = fgetl(fid);
    
    if ~ischar(tline), break, end
    
    %------------------------------------------------ k-space range
    % Find kz_range
    idx = strfind(tline,'kz_range');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        kz_range = str2num(s);
    end
    
    % Find ky_range
    idx = strfind(tline,'ky_range');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        ky_range = str2num(s);
    end
    
    % Find kx_range
    idx = strfind(tline,'kx_range');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        kx_range = str2num(s);
    end
    
    %------------------------------------------------ k-space ovs_fac
    
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
    
    %------------------------------------------------ resolution
    
    % Find Z-resolution
    idx = strfind(tline,'Z-resolution');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        Z_resolution = str2num(s);
    end
    
    % Find Y-resolution
    idx = strfind(tline,'Y-resolution');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        Y_resolution = str2num(s);
    end
    
    % Find X-resolution
    idx = strfind(tline,'X-resolution');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        X_resolution = str2num(s);
    end
    
    %------------------------------------------------ SENSE factor
    
    % Find Z-direction SENSE factor
    idx = strfind(tline,'Z-direction SENSE factor');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        Z_direction_SENSE_factor = str2num(s);
    end
    
    % Find Y-direction SENSE factor
    idx = strfind(tline,'Y-direction SENSE factor');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        Y_direction_SENSE_factor = str2num(s);
    end
    
    % Find X-direction SENSE factor
    idx = strfind(tline,'X-direction SENSE factor');
    if ~isempty(idx)
        idx = strfind(tline,':');
        s = tline(idx+1:end);
        X_direction_SENSE_factor = str2num(s);
    end
    
    %------------------------------------------------ 
    
    % Finish if all search_params were found.
    if exist('kx_range','var') && exist('kx_oversample_factor','var') && ...
            exist('X_resolution','var') && exist('X_direction_SENSE_factor','var')
        break
    end
    
end
fclose(fid);

if exist('kz_range','var')
    params.kz_range = kz_range;
    clear  kz_range
end
if exist('ky_range','var')
    params.ky_range = ky_range;
    clear  ky_range
end
if exist('kx_range','var')
    params.kx_range = kx_range;
    clear  kx_range
end
if exist('kz_oversample_factor','var')
    params.kz_oversample_factor = kz_oversample_factor;
    clear  kz_oversample_factor
end
if exist('ky_oversample_factor','var')
    params.ky_oversample_factor = ky_oversample_factor;
    clear  ky_oversample_factor
end
if exist('kx_oversample_factor','var')
    params.kx_oversample_factor = kx_oversample_factor;
    clear  kx_oversample_factor
end
if exist('Z_resolution','var')
    params.Z_resolution = Z_resolution;
    clear  Z_resolution
end
if exist('Y_resolution','var')
    params.Y_resolution = Y_resolution;
    clear  Y_resolution
end
if exist('X_resolution','var')
    params.X_resolution = X_resolution;
    clear  X_resolution
end
if exist('Z_direction_SENSE_factor','var')
    params.Z_direction_SENSE_factor = Z_direction_SENSE_factor;
    clear  Z_direction_SENSE_factor
end
if exist('Y_direction_SENSE_factor','var')
    params.Y_direction_SENSE_factor = Y_direction_SENSE_factor;
    clear  Y_direction_SENSE_factor
end
if exist('X_direction_SENSE_factor','var')
    params.X_direction_SENSE_factor = X_direction_SENSE_factor;
    clear  X_direction_SENSE_factor
end

t = toc;
fprintf('scan parameters are read from [%s] in %f sec\n',f_s,t)
fprintf('\n')















