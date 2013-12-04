
%[s_resize_SMAP] resizes SMAP by resize, zeropad and/or cut.
%
%
% Last modified
% 2010.08.09
% 2010.08.18.
%   Pack memory.
% 2010.09.22.
%   Add REFparams to f_resize_SMAP.
% 2010.09.24.
%   Save I_smap_4d to saveDataDir_s from sharedDir_s.
% 2010.12.09.
%   Load 'I_smap__(REFparams.filename)' to distinguish REF files from
%   different slice orientations.
% 2012.05.23.
%   Take care of volume coil case (REFparams.filename is empty).
%
% Ha-Kyu



%% Prepare SMAP

cd(saveDataDir_s)
if ~(exist('I_smap_4d.mat','file')==2) && ~isempty(REFparams.filename)
    
    % Load SMAP.
    cd(sharedDir_s)
    %load  I_smap
    eval(sprintf('load  I_smap__%s',REFparams.filename))
    
    % Process SMAP.
    [I_smap_4d,I_mask_S_3d] = f_resize_SMAP(I_smap, DWIparams, REFparams, RECONparams);
    [sny,snx,snz,snc] = size(I_smap_4d);
    %clear  I_smap
    eval(sprintf('clear  I_smap__%s',REFparams.filename))
    
    % Report.
    fprintf('\n')
    fprintf('    I_smap_4d,     [sny,snx,snz,snc]  = [%d,%d,%d,%d]\n',sny,snx,snz,snc)
    fprintf('    I_mask_S_3d,   [sny,snx,snz]      = [%d,%d,%d]\n',sny,snx,snz)
    
    % Save SMAP.
    cd(saveDataDir_s)
    eval(sprintf('save  I_smap_4d  I_smap_4d'))
    eval(sprintf('save  I_mask_S_3d  I_mask_S_3d'))
    
    % Clear data after processing.
    clear  I_smap_4d  I_mask_S_3d
end

% Pack memory.
pack

fprintf('\n\n')





