
%[s_match_sagittal_REF_DWI] matches object orientation between
%REF and DWI based on fold-over and fat-shift direction when both REF and
%DWI are acquired at sagittal orientation.
%
%
% Last modified
% 2010.09.26.
% 2010.10.06.
%   Use fliplr(flipud()) for REF (AP-F) and DWI (AP-P). This is for recon
%   of [scan20101005_r2111__Yankeelov_2351_3T] data.
% 2010.12.22.
%   REF(AP-F)-DWI(AP-P) is verified.



%% Match object position
[sny,snx,snz,snc] = size(I_smap_echo_4d);

if strcmpi(REFparams.fold_over_dir,'AP') && ...
        strcmpi(DWIparams.fold_over_dir,'AP')
    
    if strcmpi(REFparams.fat_shift_dir,'F')
        I_smap_echo_temp_4d = zeros(sny,snx,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(sny,snx,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'P')
                    % 2010.12.22.
                    I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
                        fliplr(flipud(I_smap_echo_4d(:,:,ind_slice,ind_coil))); % 7T and 3T (scan20101005_r2111__Yankeelov_2351_3T)
                    %I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
                    %    fliplr((I_smap_echo_4d(:,:,ind_slice,ind_coil))); % scan20100916_r2111__Lori_Phantom_3T                    
                    if ind_coil==1
                        %I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
                        %    fliplr((I_mask_S_echo_3d(:,:,ind_slice)));
                        I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
                            fliplr(flipud(I_mask_S_echo_3d(:,:,ind_slice)));
                    end
                elseif  strcmpi(DWIparams.fat_shift_dir,'A')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         flipud(fliplr(I_smap_echo_4d(:,:,ind_slice,ind_coil)));
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             flipud(fliplr(I_mask_S_echo_3d(:,:,ind_slice)));
%                     end
                else
                    error('s_match_sagittal_REF_DWI:main','Unknown DWI fat_shift_dir')
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
                if strcmpi(DWIparams.fat_shift_dir,'P')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         flipud(I_smap_echo_4d(:,:,ind_slice,ind_coil));
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             flipud(I_mask_S_echo_3d(:,:,ind_slice));
%                     end
                elseif  strcmpi(DWIparams.fat_shift_dir,'A')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         fliplr(I_smap_echo_4d(:,:,ind_slice,ind_coil));
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             fliplr(I_mask_S_echo_3d(:,:,ind_slice));
%                     end
                else
                    error('s_match_sagittal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
    else
        error('s_match_sagittal_REF_DWI:main','Unknown REF fat_shift_dir')
    end
    
elseif strcmpi(REFparams.fold_over_dir,'AP') && ...
        strcmpi(DWIparams.fold_over_dir,'RL')
    
    if strcmpi(REFparams.fat_shift_dir,'L')
        I_smap_echo_temp_4d = zeros(snx,sny,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(snx,sny,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'L')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         permute(flipud(I_smap_echo_4d(:,:,ind_slice,ind_coil)),[2,1]);
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             permute(flipud(I_mask_S_echo_3d(:,:,ind_slice)),[2,1]);
%                     end
                elseif  strcmpi(DWIparams.fat_shift_dir,'R')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         permute(fliplr(I_smap_echo_4d(:,:,ind_slice,ind_coil)),[2,1]);
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             permute(fliplr(I_mask_S_echo_3d(:,:,ind_slice)),[2,1]);
%                     end
                else
                    error('s_match_sagittal_REF_DWI:main','Unknown DWI fat_shift_dir')
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
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         permute(I_smap_echo_4d(:,:,ind_slice,ind_coil),[2,1]);
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             permute(I_mask_S_echo_3d(:,:,ind_slice),[2,1]);
%                     end
                elseif  strcmpi(DWIparams.fat_shift_dir,'R')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         permute(fliplr(flipud(I_smap_echo_4d(:,:,ind_slice,ind_coil))),[2,1]);
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             permute(fliplr(flipud(I_mask_S_echo_3d(:,:,ind_slice))),[2,1]);
%                     end
                else
                    error('s_match_sagittal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
    else
        error('s_match_sagittal_REF_DWI:main','Unknown REF fat_shift_dir')
    end
    
elseif strcmpi(REFparams.fold_over_dir,'RL') && ...
        strcmpi(DWIparams.fold_over_dir,'AP')
    
    if strcmpi(REFparams.fat_shift_dir,'P')
        I_smap_echo_temp_4d = zeros(snx,sny,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(snx,sny,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'P')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         permute(flipud(I_smap_echo_4d(:,:,ind_slice,ind_coil)),[2,1]);
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             permute(flipud(I_mask_S_echo_3d(:,:,ind_slice)),[2,1]);
%                     end
                elseif  strcmpi(DWIparams.fat_shift_dir,'A')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         permute(fliplr(I_smap_echo_4d(:,:,ind_slice,ind_coil)),[2,1]);
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             permute(fliplr(I_mask_S_echo_3d(:,:,ind_slice)),[2,1]);
%                     end
                else
                    error('s_match_sagittal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
        
    elseif strcmpi(REFparams.fat_shift_dir,'A')
        
        I_smap_echo_temp_4d = zeros(snx,sny,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(snx,sny,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'P')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         permute(flipud(fliplr(I_smap_echo_4d(:,:,ind_slice,ind_coil))),[2,1]);
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             permute(flipud(fliplr(I_mask_S_echo_3d(:,:,ind_slice))),[2,1]);
%                     end
                elseif  strcmpi(DWIparams.fat_shift_dir,'A')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         permute(I_smap_echo_4d(:,:,ind_slice,ind_coil),[2,1]);
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             permute(I_mask_S_echo_3d(:,:,ind_slice),[2,1]);
%                     end
                else
                    error('s_match_sagittal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
    else
        error('s_match_sagittal_REF_DWI:main','Unknown REF fat_shift_dir')
    end
    
elseif strcmpi(REFparams.fold_over_dir,'RL') && ...
        strcmpi(DWIparams.fold_over_dir,'RL')
    
    if strcmpi(REFparams.fat_shift_dir,'P')
        I_smap_echo_temp_4d = zeros(sny,snx,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(sny,snx,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'L')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         fliplr(flipud(I_smap_echo_4d(:,:,ind_slice,ind_coil)));
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             fliplr(flipud(I_mask_S_echo_3d(:,:,ind_slice)));
%                     end
                elseif  strcmpi(DWIparams.fat_shift_dir,'R')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         I_smap_echo_4d(:,:,ind_slice,ind_coil);
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             I_mask_S_echo_3d(:,:,ind_slice);
%                     end
                else
                    error('s_match_sagittal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
        
    elseif strcmpi(REFparams.fat_shift_dir,'A')
        
        I_smap_echo_temp_4d = zeros(sny,snx,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(sny,snx,snz,'single');
        for ind_slice = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'L')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         flipud(I_smap_echo_4d(:,:,ind_slice,ind_coil));
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             flipud(I_mask_S_echo_3d(:,:,ind_slice));
%                     end
                elseif  strcmpi(DWIparams.fat_shift_dir,'R')
%                     I_smap_echo_temp_4d(:,:,ind_slice,ind_coil) = ...
%                         fliplr(I_smap_echo_4d(:,:,ind_slice,ind_coil));
%                     if ind_coil==1
%                         I_mask_S_echo_temp_3d(:,:,ind_slice) = ...
%                             fliplr(I_mask_S_echo_3d(:,:,ind_slice));
%                     end
                else
                    error('s_match_sagittal_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
    else
        error('s_match_sagittal_REF_DWI:main','Unknown REF fat_shift_dir')
    end
else
    error('s_match_sagittal_REF_DWI:main','Unknown REF and DWI fold_over_dir match')
end






