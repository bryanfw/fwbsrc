
%[s_match_transverse_REF_DWI] matches object orientation between
%REF and DWI based on fold-over and fat-shift direction when both REF and
%DWI are acquired at transverse orientation.
%
%
% Last modified
%   2010.09.23.
%   2012.07.13.
%
% Ha-Kyu



%% Match object position
[sny,snx,snz,snc] = size(I_smap_echo_4d);
[sny1,snx1,snz1,snc1] = size(I_mask_S_echo_3d);

if strcmpi(REFparams.fold_over_dir,'AP') && ...
        strcmpi(DWIparams.fold_over_dir,'AP')
    
    if strcmpi(REFparams.fat_shift_dir,'L')
        I_smap_echo_temp_4d = zeros(sny,snx,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(sny1,snx1,snz1,'single');
        for ind_sl = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'P')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        I_smap_echo_4d(:,:,ind_sl,ind_coil);
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        I_mask_S_echo_3d(:,:,ind_sl);
                    %end
                elseif  strcmpi(DWIparams.fat_shift_dir,'A')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        flipud(fliplr(I_smap_echo_4d(:,:,ind_sl,ind_coil)));
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        flipud(fliplr(I_mask_S_echo_3d(:,:,ind_sl)));
                    %end
                else
                    error('s_match_transverse_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        for ind_sl = 1:snz1
            if strcmpi(DWIparams.fat_shift_dir,'P')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            I_mask_S_echo_3d(:,:,ind_sl);
            elseif  strcmpi(DWIparams.fat_shift_dir,'A')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            flipud(fliplr(I_mask_S_echo_3d(:,:,ind_sl)));
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
        
    elseif strcmpi(REFparams.fat_shift_dir,'R')
        
        I_smap_echo_temp_4d = zeros(sny,snx,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(sny1,snx1,snz1,'single');
        for ind_sl = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'P')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        flipud(I_smap_echo_4d(:,:,ind_sl,ind_coil));
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        flipud(I_mask_S_echo_3d(:,:,ind_sl));
                    %end
                elseif  strcmpi(DWIparams.fat_shift_dir,'A')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        fliplr(I_smap_echo_4d(:,:,ind_sl,ind_coil));
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        fliplr(I_mask_S_echo_3d(:,:,ind_sl));
                    %end
                else
                    error('s_match_transverse_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        for ind_sl = 1:snz1
            if strcmpi(DWIparams.fat_shift_dir,'P')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            flipud(I_mask_S_echo_3d(:,:,ind_sl));
            elseif  strcmpi(DWIparams.fat_shift_dir,'A')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            fliplr(fliplr(I_mask_S_echo_3d(:,:,ind_sl)));
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
    else
        error('s_match_transverse_REF_DWI:main','Unknown REF fat_shift_dir')
    end
    
elseif strcmpi(REFparams.fold_over_dir,'AP') && ...
        strcmpi(DWIparams.fold_over_dir,'RL')
    
    if strcmpi(REFparams.fat_shift_dir,'L')
        I_smap_echo_temp_4d = zeros(snx,sny,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(snx1,sny1,snz1,'single');
        for ind_sl = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'L')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        permute(flipud(I_smap_echo_4d(:,:,ind_sl,ind_coil)),[2,1]);
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        permute(flipud(I_mask_S_echo_3d(:,:,ind_sl)),[2,1]);
                    %end
                elseif  strcmpi(DWIparams.fat_shift_dir,'R')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        permute(fliplr(I_smap_echo_4d(:,:,ind_sl,ind_coil)),[2,1]);
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        permute(fliplr(I_mask_S_echo_3d(:,:,ind_sl)),[2,1]);
                    %end
                else
                    error('s_match_transverse_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        for ind_sl = 1:snz1
            if strcmpi(DWIparams.fat_shift_dir,'L')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            permute(flipud(I_mask_S_echo_3d(:,:,ind_sl)),[2,1]);
            elseif  strcmpi(DWIparams.fat_shift_dir,'R')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            permute(fliplr(I_mask_S_echo_3d(:,:,ind_sl)),[2,1]);
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
        
    elseif strcmpi(REFparams.fat_shift_dir,'R')
        
        I_smap_echo_temp_4d = zeros(snx,sny,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(snx1,sny1,snz1,'single');
        for ind_sl = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'L')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        permute(I_smap_echo_4d(:,:,ind_sl,ind_coil),[2,1]);
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        permute(I_mask_S_echo_3d(:,:,ind_sl),[2,1]);
                    %end
                elseif  strcmpi(DWIparams.fat_shift_dir,'R')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        permute(fliplr(flipud(I_smap_echo_4d(:,:,ind_sl,ind_coil))),[2,1]);
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        permute(fliplr(flipud(I_mask_S_echo_3d(:,:,ind_sl))),[2,1]);
                    %end
                else
                    error('s_match_transverse_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        for ind_sl = 1:snz1
            if strcmpi(DWIparams.fat_shift_dir,'L')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            permute(I_mask_S_echo_3d(:,:,ind_sl),[2,1]);
            elseif  strcmpi(DWIparams.fat_shift_dir,'R')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            permute(fliplr(flipud(I_mask_S_echo_3d(:,:,ind_sl))),[2,1]);
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
    else
        error('s_match_transverse_REF_DWI:main','Unknown REF fat_shift_dir')
    end
    
elseif strcmpi(REFparams.fold_over_dir,'RL') && ...
        strcmpi(DWIparams.fold_over_dir,'AP')
    
    if strcmpi(REFparams.fat_shift_dir,'P')
        I_smap_echo_temp_4d = zeros(snx,sny,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(snx1,sny1,snz1,'single');
        for ind_sl = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'P')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        permute(flipud(I_smap_echo_4d(:,:,ind_sl,ind_coil)),[2,1]);
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        permute(flipud(I_mask_S_echo_3d(:,:,ind_sl)),[2,1]);
                    %end
                elseif  strcmpi(DWIparams.fat_shift_dir,'A')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        permute(fliplr(I_smap_echo_4d(:,:,ind_sl,ind_coil)),[2,1]);
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        permute(fliplr(I_mask_S_echo_3d(:,:,ind_sl)),[2,1]);
                    %end
                else
                    error('s_match_transverse_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        for ind_sl = 1:snz1
            if strcmpi(DWIparams.fat_shift_dir,'P')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            permute(flipud(I_mask_S_echo_3d(:,:,ind_sl)),[2,1]);
            elseif  strcmpi(DWIparams.fat_shift_dir,'A')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            permute(fliplr(I_mask_S_echo_3d(:,:,ind_sl)),[2,1]);
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
        
    elseif strcmpi(REFparams.fat_shift_dir,'A')
        
        I_smap_echo_temp_4d = zeros(snx,sny,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(snx1,sny1,snz1,'single');
        for ind_sl = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'P')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        permute(flipud(fliplr(I_smap_echo_4d(:,:,ind_sl,ind_coil))),[2,1]);
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        permute(flipud(fliplr(I_mask_S_echo_3d(:,:,ind_sl))),[2,1]);
                    %end
                elseif  strcmpi(DWIparams.fat_shift_dir,'A')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        permute(I_smap_echo_4d(:,:,ind_sl,ind_coil),[2,1]);
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        permute(I_mask_S_echo_3d(:,:,ind_sl),[2,1]);
                    %end
                else
                    error('s_match_transverse_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        for ind_sl = 1:snz1
            if strcmpi(DWIparams.fat_shift_dir,'P')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            permute(flipud(fliplr(I_mask_S_echo_3d(:,:,ind_sl))),[2,1]);
            elseif  strcmpi(DWIparams.fat_shift_dir,'A')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            permute(I_mask_S_echo_3d(:,:,ind_sl),[2,1]);
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
    else
        error('s_match_transverse_REF_DWI:main','Unknown REF fat_shift_dir')
    end
    
elseif strcmpi(REFparams.fold_over_dir,'RL') && ...
        strcmpi(DWIparams.fold_over_dir,'RL')
    
    if strcmpi(REFparams.fat_shift_dir,'P')
        I_smap_echo_temp_4d = zeros(sny,snx,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(sny1,snx1,snz1,'single');
        for ind_sl = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'L')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        fliplr(flipud(I_smap_echo_4d(:,:,ind_sl,ind_coil)));
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        fliplr(flipud(I_mask_S_echo_3d(:,:,ind_sl)));
                    %end
                elseif  strcmpi(DWIparams.fat_shift_dir,'R')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        I_smap_echo_4d(:,:,ind_sl,ind_coil);
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        I_mask_S_echo_3d(:,:,ind_sl);
                    %end
                else
                    error('s_match_transverse_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        for ind_sl = 1:snz1
            if strcmpi(DWIparams.fat_shift_dir,'L')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            fliplr(flipud(I_mask_S_echo_3d(:,:,ind_sl)));
            elseif  strcmpi(DWIparams.fat_shift_dir,'R')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            I_mask_S_echo_3d(:,:,ind_sl);
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
        
    elseif strcmpi(REFparams.fat_shift_dir,'A')
        
        I_smap_echo_temp_4d = zeros(sny,snx,snz,snc,'single');
        I_mask_S_echo_temp_3d = zeros(sny,snx,snz,'single');
        for ind_sl = 1:snz
            for ind_coil = 1:snc
                if strcmpi(DWIparams.fat_shift_dir,'L')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        flipud(I_smap_echo_4d(:,:,ind_sl,ind_coil));
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        flipud(I_mask_S_echo_3d(:,:,ind_sl));
                    %end
                elseif  strcmpi(DWIparams.fat_shift_dir,'R')
                    I_smap_echo_temp_4d(:,:,ind_sl,ind_coil) = ...
                        fliplr(I_smap_echo_4d(:,:,ind_sl,ind_coil));
                    %if ind_coil==1
                    %    I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                    %        fliplr(I_mask_S_echo_3d(:,:,ind_sl));
                    %end
                else
                    error('s_match_transverse_REF_DWI:main','Unknown DWI fat_shift_dir')
                end
            end
        end
        for ind_sl = 1:snz1
            if strcmpi(DWIparams.fat_shift_dir,'L')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            flipud(I_mask_S_echo_3d(:,:,ind_sl));
            elseif  strcmpi(DWIparams.fat_shift_dir,'R')
                I_mask_S_echo_temp_3d(:,:,ind_sl) = ...
                            fliplr(I_mask_S_echo_3d(:,:,ind_sl));
            end
        end
        clear  I_smap_echo_4d  I_mask_S_echo_3d
        I_smap_echo_4d = I_smap_echo_temp_4d;
        I_mask_S_echo_3d = I_mask_S_echo_temp_3d;
        clear  I_smap_echo_temp_4d  I_mask_S_echo_temp_3d
    else
        error('s_match_transverse_REF_DWI:main','Unknown REF fat_shift_dir')
    end
else
    error('s_match_transverse_REF_DWI:main','Unknown REF and DWI fold_over_dir match')
end






