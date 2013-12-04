
function start_sign = f_get_grad_start_sign(fname)
%[f_get_grad_start_sign] get readout gradient starting sign for STD and PHC
%or PHX.
%
% USAGE:
%   start_sign = f_get_grad_start_sign(fname)
%
% INPUT:
%   fname: Raw data name, e.g. 'raw_006'
%
% OUTPUT:
%   start_sign: A struct with fields as,
%       STD_start_sign, if exist
%       PHC_start_sign, if exist
%       PHX_start_sign, if exist
%
% NOTE:
%   This routine is much faster than using regexp used in [readListData.m].
%
%
% Last modified
% 2010.08.06.
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
tic
flag_readSTD = 1;
flag_readPHC = 1;
flag_readPHX = 1;

STD_start_sign = [];
PHC_start_sign = [];
PHX_start_sign = [];
start_sign = [];

fid=fopen(f_s);
while 1    
    tline = fgetl(fid);
    
    if ~ischar(tline), break, end
    
    if length(tline) > 5 && strcmpi(tline(3:5),'STD') && flag_readSTD && ~strcmpi(tline(1),'#')
        numline = str2num((tline(6:end)));
        if ~isempty(numline)
            STD_start_sign = numline(13);
            flag_readSTD = 0;
        end
    end
    
    if length(tline) > 5 && strcmpi(tline(3:5),'PHC') && flag_readPHC && ~strcmpi(tline(1),'#')
        numline = str2num((tline(6:end)));        
        if ~isempty(numline)
            PHC_start_sign = numline(13);
            flag_readPHC = 0;
        end
    end
    
    if length(tline) > 5 && strcmpi(tline(3:5),'PHX') && flag_readPHX && ~strcmpi(tline(1),'#')
        numline = str2num((tline(6:end)));
        if ~isempty(numline)
            PHX_start_sign = numline(13);
            flag_readPHX = 0;
        end
    end
    
    if flag_readSTD==0 && (flag_readPHC==0 || flag_readPHX==0)
        break;
    end
    
end
fclose(fid);

if ~isempty(STD_start_sign)
    start_sign.STD_start_sign = STD_start_sign;
    clear  STD_start_sign
end
if ~isempty(PHC_start_sign)
    start_sign.PHC_start_sign = PHC_start_sign;
    clear  PHC_start_sign
end
if ~isempty(PHX_start_sign)
    start_sign.PHX_start_sign = PHX_start_sign;
    clear  PHX_start_sign
end

t = toc;
fprintf('Gradient starting sign is read in %f sec\n',t)
fprintf('\n')















