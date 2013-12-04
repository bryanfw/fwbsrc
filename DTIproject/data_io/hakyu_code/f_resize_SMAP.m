
function [I_smap_4d,I_mask_3d] = f_resize_SMAP(I_smap, DWIparams, REFparams, RECONparams)
%[f_resize_SMAP] resizes coil sensitivity map to be ready for
%reconstruction.
%
% INPUT
%   I_smap          : Coil sensitivity map. A struct including I_smap.smap.
%   DWIparams       : Parameters for DWI scan.
%   REFparams       : Parameters for REF scan.
%   RECONparams     : Parameters for reconstruction.
%
% OUTPUT
%   I_smap_4d   :   Coil sensitivity map ready for reconstruction.
%   I_mask_3d   :   Mask of S. Generated from root-sum-of-square (NORM) of original S.
%
%
% Last modified
% 2010.07.12. Modified from [prepare_S_v2.m].
% 2010.08.09
%   This is generated from [prepare_S_v3.m].
%   Removed REFparams from input, because it is not used here.
% 2010.09.23.
%   Consider REFparams.fold_over_dir.
% 2010.11.06.
%   Add coronal object repositioning.



%% Check intput
if ~isstruct(I_smap)
    error('f_resize_SMAP:main','Input [I_smap] must be a STRUCT.')
end
if ~isfield(I_smap,'smap')
    error('f_resize_SMAP:main','Input [I_smap] must have a FIELD named ''smap''.')
end
[sny,snx,snz,snc] = size(I_smap.smap);



%% Prepare SMAP

% Take intermediate and final SMAP size.
s_N_res_m = RECONparams.s_N_res_m;
s_N_res_p = RECONparams.s_N_res_p;
s_N_final_m = RECONparams.s_N_final_m;
s_N_final_p = RECONparams.s_N_final_p;


% Reserve output depending on the foldover direction.
% Transverse.
if strcmpi(DWIparams.slice_orientation,'transverse') && ...
        strcmpi(REFparams.slice_orientation,'transverse')
    if (strcmpi(REFparams.fold_over_dir,'AP') && ...
            strcmpi(DWIparams.fold_over_dir,'RL')) || ...
            (strcmpi(REFparams.fold_over_dir,'RL') && ...
            strcmpi(DWIparams.fold_over_dir,'AP'))
        flag_fold_over_match = false;
        I_smap_4d = zeros(s_N_final_m,s_N_final_p,snz,snc,'single');
        I_mask_3d = zeros(s_N_final_m,s_N_final_p,snz);        
    elseif (strcmpi(REFparams.fold_over_dir,'AP') && ...
            strcmpi(DWIparams.fold_over_dir,'AP')) || ...
            (strcmpi(REFparams.fold_over_dir,'RL') && ...
            strcmpi(DWIparams.fold_over_dir,'RL'))
        flag_fold_over_match = true;
        I_smap_4d = zeros(s_N_final_p,s_N_final_m,snz,snc,'single');
        I_mask_3d = zeros(s_N_final_p,s_N_final_m,snz);
    else
        error('f_resize_SMAP:main','Unknown REF fold_over_dir and DWI fold_over_dir')
    end
end

% Sagittal.
if strcmpi(DWIparams.slice_orientation,'sagittal') && ...
        strcmpi(REFparams.slice_orientation,'sagittal')
    if (strcmpi(REFparams.fold_over_dir,'AP') && ...
            strcmpi(DWIparams.fold_over_dir,'AP')) || ...
            (strcmpi(REFparams.fold_over_dir,'FH') && ...
            strcmpi(DWIparams.fold_over_dir,'FH'))
        flag_fold_over_match = true;
        I_smap_4d = zeros(s_N_final_p,s_N_final_m,snz,snc,'single');
        I_mask_3d = zeros(s_N_final_p,s_N_final_m,snz);        
    elseif (strcmpi(REFparams.fold_over_dir,'AP') && ...
            strcmpi(DWIparams.fold_over_dir,'FH')) || ...
            (strcmpi(REFparams.fold_over_dir,'FH') && ...
            strcmpi(DWIparams.fold_over_dir,'AP'))
        %flag_fold_over_match = false;
        %I_smap_4d = zeros(s_N_final_m,s_N_final_p,snz,snc,'single');
        %I_mask_3d = zeros(s_N_final_m,s_N_final_p,snz);        
    else
        error('f_resize_SMAP:main','Unknown REF fold_over_dir and DWI fold_over_dir')
    end
end

% Coronal.
if strcmpi(DWIparams.slice_orientation,'coronal') && ...
        strcmpi(REFparams.slice_orientation,'coronal')
    if (strcmpi(REFparams.fold_over_dir,'RL') && ...
            strcmpi(DWIparams.fold_over_dir,'RL')) || ...
            (strcmpi(REFparams.fold_over_dir,'FH') && ...
            strcmpi(DWIparams.fold_over_dir,'FH'))
        flag_fold_over_match = true;
        I_smap_4d = zeros(s_N_final_p,s_N_final_m,snz,snc,'single');
        I_mask_3d = zeros(s_N_final_p,s_N_final_m,snz);        
    elseif (strcmpi(REFparams.fold_over_dir,'RL') && ...
            strcmpi(DWIparams.fold_over_dir,'FH')) || ...
            (strcmpi(REFparams.fold_over_dir,'FH') && ...
            strcmpi(DWIparams.fold_over_dir,'RL'))
        %flag_fold_over_match = false;
        %I_smap_4d = zeros(s_N_final_m,s_N_final_p,snz,snc,'single');
        %I_mask_3d = zeros(s_N_final_m,s_N_final_p,snz);        
    else
        error('f_resize_SMAP:main','Unknown REF fold_over_dir and DWI fold_over_dir')
    end
end



%% Resize and cut and/or zeropad SMAP
for ind_slice = 1:snz
    for ind_coil = 1:snc
        smap = I_smap.smap(:,:,ind_slice,ind_coil);
        
        if flag_fold_over_match==false
            
            %-----------
            % 1. Resize.
            smap_res = imresize(smap,[s_N_res_m,s_N_res_p]);
            
            %-------------------
            % 2. Cut or zeropad.            
            % Frequency-encoding direction.
            if RECONparams.s_N_res_m > RECONparams.s_N_final_m % cut
                start_vox = RECONparams.s_cut_vox_m+1;
                end_vox = s_N_res_m - RECONparams.s_cut_vox_m;
                smap_res_extr_m = smap_res(start_vox:end_vox,:);
            else % zeropad
                s_zeropad_vox_m = RECONparams.s_zeropad_vox_m;
                smap_res_extr_m = padarray(padarray(smap_res,[s_zeropad_vox_m,0],0,'pre'), ...
                    [s_zeropad_vox_m,0],0,'post');
            end
            
            % Phase-encoding direction.
            if RECONparams.s_N_res_p > RECONparams.s_N_final_p % cut
                start_vox = RECONparams.s_cut_vox_p+1;
                end_vox = s_N_res_p - RECONparams.s_cut_vox_p;
                smap_res_extr_mp = smap_res_extr_m(:,start_vox:end_vox);
            else % zeropad
                s_zeropad_vox_p = RECONparams.s_zeropad_vox_p;
                smap_res_extr_mp = padarray(padarray(smap_res_extr_m,[0,s_zeropad_vox_p],0,'pre'), ...
                    [0,s_zeropad_vox_p],0,'post');
            end
            
            % Check output.
            if any(size(smap_res_extr_mp)~=[RECONparams.s_N_final_m,RECONparams.s_N_final_p])
                error('f_resize_SMAP:main','Check the size of ''smap_res_extr_mp''.')
            end
            
            % Clear data.
            clear  smap_res_extr_m  smap_res
            
        elseif flag_fold_over_match==true
            
            %-----------
            % 1. Resize.
            smap_res = imresize(smap,[s_N_res_p,s_N_res_m]);
            
            %-------------------
            % 2. Cut or zeropad.            
            % Frequency-encoding direction.
            if RECONparams.s_N_res_m > RECONparams.s_N_final_m % cut
                start_vox = RECONparams.s_cut_vox_m+1;
                end_vox = s_N_res_m - RECONparams.s_cut_vox_m;
                smap_res_extr_m = smap_res(:,start_vox:end_vox);
            else % zeropad
                s_zeropad_vox_m = RECONparams.s_zeropad_vox_m;
                smap_res_extr_m = padarray(padarray(smap_res,[0,s_zeropad_vox_m],0,'pre'), ...
                    [0,s_zeropad_vox_m],0,'post');
            end
            
            % Phase-encoding direction.
            if RECONparams.s_N_res_p > RECONparams.s_N_final_p % cut
                start_vox = RECONparams.s_cut_vox_p+1;
                end_vox = s_N_res_p - RECONparams.s_cut_vox_p;
                smap_res_extr_mp = smap_res_extr_m(start_vox:end_vox,:);
            else % zeropad
                s_zeropad_vox_p = RECONparams.s_zeropad_vox_p;
                smap_res_extr_mp = padarray(padarray(smap_res_extr_m,[s_zeropad_vox_p,0],0,'pre'), ...
                    [s_zeropad_vox_p,0],0,'post');
            end
            
            % Check output.
            if any(size(smap_res_extr_mp)~=[RECONparams.s_N_final_p,RECONparams.s_N_final_m])
                error('f_resize_SMAP:main','Check the size of ''smap_res_extr_mp''.')
            end
            
            % Clear data.
            clear  smap_res_extr_m  smap_res
        end
        
        % Reserve output.
        I_smap_4d(:,:,ind_slice,ind_coil) = smap_res_extr_mp;
    end % for ind_coil
        
    % Generate mask.
    I_smap_ssq = abs(sqrt(sum(squeeze(I_smap_4d(:,:,ind_slice,:)) .* ...
        conj(squeeze(I_smap_4d(:,:,ind_slice,:))),3)));
    %thresh = graythresh(I_smap_ssq_m);
    thresh = 0;
    I_mask_m = I_smap_ssq > thresh;
    %I_mask_m = imfill(double(I_mask_m),'holes');
    I_mask_m = imfill(double(imclose(I_mask_m,strel('disk',5))),'holes');
    I_mask_3d(:,:,ind_slice) = I_mask_m > thresh; % logical
end % for ind_slice
clear  I_smap  I_mask_m  thresh  I_smap_ssq





%% END






