
%[s_match_coronal_REF_DWI] matches object orientation between
%REF and DWI based on fold-over and fat-shift direction when both REF and
%DWI are acquired at coronal orientation.
%
%
% Last modified
% 2010.11.03.



%% Match object position
[sny,snx,snz,snc] = size(I_smap_echo_4d);

if strcmpi(REFparams.fold_over_dir,'RL') && ...
        strcmpi(DWIparams.fold_over_dir,'RL')
    
    if strcmpi(REFparams.fat_shift_dir,'F')
        I_smap_echo_temp_4d = zeros(sny,snx,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(sny,snx,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'L')
                    I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
                        fliplr(flipud(I_smap_echo_4d(:,:,ind_slice,ind_coil)));
                    if ind_coil==1
                        I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
                            fliplr(flipud(I_mask_S_echo_3d(:,:,ind_slice)));
                    end
                elseif  strcmpi(DWIparams.fat_shift_dir,'R')
                    % not set up yet
                else
                    error('s_match_coronal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
        
    elseif strcmpi(REFparams.fat_shift_dir,'H')
        
        I_smap_echo_temp_4d = zeros(sny,snx,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(sny,snx,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'L')
                    % not set up yet
                elseif  strcmpi(DWIparams.fat_shift_dir,'R')
                    % not set up yet
                else
                    error('s_match_coronal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
    else
        error('s_match_coronal_REF_DWI:main','Unknown REF fat_shift_dir')
    end
    
elseif strcmpi(REFparams.fold_over_dir,'RL') && ...
        strcmpi(DWIparams.fold_over_dir,'FH')
    
    if strcmpi(REFparams.fat_shift_dir,'F')
        I_smap_echo_temp_4d = zeros(snx,sny,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(snx,sny,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'F')
                    % not set up yet
                elseif  strcmpi(DWIparams.fat_shift_dir,'H')
                    % not set up yet
                else
                    error('s_match_coronal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
        
    elseif strcmpi(REFparams.fat_shift_dir,'H')
        
        I_smap_echo_temp_4d = zeros(snx,sny,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(snx,sny,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'F')
                    % not set up yet
                elseif  strcmpi(DWIparams.fat_shift_dir,'H')
                    % not set up yet
                else
                    error('s_match_coronal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
    else
        error('s_match_coronal_REF_DWI:main','Unknown REF fat_shift_dir')
    end
    
elseif strcmpi(REFparams.fold_over_dir,'FH') && ...
        strcmpi(DWIparams.fold_over_dir,'RL')
    
    if strcmpi(REFparams.fat_shift_dir,'L')
        I_smap_echo_temp_4d = zeros(snx,sny,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(snx,sny,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'L')
                    % not set up yet
                elseif  strcmpi(DWIparams.fat_shift_dir,'R')
                    % not set up yet
                else
                    error('s_match_coronal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
        
    elseif strcmpi(REFparams.fat_shift_dir,'R')
        
        I_smap_echo_temp_4d = zeros(snx,sny,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(snx,sny,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'L')
                    % not set up yet
                elseif  strcmpi(DWIparams.fat_shift_dir,'R')
                    % not set up yet
                else
                    error('s_match_coronal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
    else
        error('s_match_coronal_REF_DWI:main','Unknown REF fat_shift_dir')
    end
    
elseif strcmpi(REFparams.fold_over_dir,'FH') && ...
        strcmpi(DWIparams.fold_over_dir,'FH')
    
    if strcmpi(REFparams.fat_shift_dir,'L')
        I_smap_echo_temp_4d = zeros(sny,snx,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(sny,snx,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'F')
                    % not set up yet
                elseif  strcmpi(DWIparams.fat_shift_dir,'H')
                    % not set up yet
                else
                    error('s_match_coronal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
        
    elseif strcmpi(REFparams.fat_shift_dir,'R')
        
        I_smap_echo_temp_4d = zeros(sny,snx,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(sny,snx,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'F')
                    % not set up yet
                elseif  strcmpi(DWIparams.fat_shift_dir,'H')
                    % not set up yet
                else
                    error('s_match_coronal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
    else
        error('s_match_coronal_REF_DWI:main','Unknown REF fat_shift_dir')
    end
else
    error('s_match_coronal_REF_DWI:main','Unknown REF and DWI fold_over_dir match')
end






