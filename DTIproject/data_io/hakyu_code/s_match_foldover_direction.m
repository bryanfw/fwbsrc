
%[s_match_foldover_direction] match foldover direction between REF and IMG
%data.
%
%
% Last modified
% 2012.07.14.
%
%
% Ha-Kyu


%% Match foldover direction.
% This is to match the object position of SMAP to foldover and fat shift
% direction of DWI.

if ~isempty(REFparams.filename)
    if strcmpi(DWIparams.slice_orientation,'transverse') && ...
            strcmpi(REFparams.slice_orientation,'transverse')
        
        % Reserve for image echo.
        I_smap_echo_4d = I_smap_img_4d;
        I_mask_S_echo_3d = I_mask_S_img_3d;
        
        s_match_transverse_REF_DWI
        
        % Reserve for image echo.
        I_smap_img_4d = I_smap_echo_4d;
        I_mask_S_img_3d = I_mask_S_echo_3d;
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        
    elseif strcmpi(DWIparams.slice_orientation,'sagittal') && ...
            strcmpi(REFparams.slice_orientation,'sagittal')
        
        % Reserve for image echo.
        I_smap_echo_4d = I_smap_img_4d;
        I_mask_S_echo_3d = I_mask_S_img_3d;
        
        s_match_sagittal_REF_DWI
        
        % Reserve for image echo.
        I_smap_img_4d = I_smap_echo_4d;
        I_mask_S_img_3d = I_mask_S_echo_3d;
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        
    elseif strcmpi(DWIparams.slice_orientation,'coronal') && ...
            strcmpi(REFparams.slice_orientation,'coronal')
        
        % Reserve for image echo.
        I_smap_echo_4d = I_smap_img_4d;
        I_mask_S_echo_3d = I_mask_S_img_3d;
        
        s_match_coronal_REF_DWI
        
        % Reserve for image echo.
        I_smap_img_4d = I_smap_echo_4d;
        I_mask_S_img_3d = I_mask_S_echo_3d;
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        
    else
        error('s_recon_IMG:main','Unknown REF and DWI slice_orientation')
    end
end