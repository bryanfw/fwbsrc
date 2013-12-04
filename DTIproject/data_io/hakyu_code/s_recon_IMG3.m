
%[s_recon_IMG3] reconstructs IMG data for b-value and slices with no phase
%correction. This is for comparison with recon with phase correction.
%
%
% Last modified
% 2010.02.11.
% 2010.02.17. saveReconDir_s -> saveDataDir_s to separate DW data and recon data.
% 2010.06.30. Decorrelate g_m.
% 2010.07.12. Consider the case when foldover is 'RL'.
%             Put '% Match foldover direction.' part.
% 2010.07.22. Consider the case when foldover is 'RL'.
%             Put '% Adjust reconstructed object position variable with foldover and fatshift.' part.
% 2010.07.26. Change psi to PSI to remove a conflict to Matlab function
%             'psi()'.
% 2010.07.28. Modified generating 'E_m' for each ind_avg.
% 2010.08.09
%   This is generated from [IMG_recon.m].
% 2010.08.18.
%   Pack memory.
% 2010.09.20.
%   Check slice_orientation.
% 2010.09.21.
%   Use K_alias_6d and K_nav_6d instead of K_alias_5d and K_nav_5d
%   considering multiple non-zero b-values.
% 2010.09.22.
%   Consider REFparams.fold_over_dir='RL'.
% 2010.09.24.
%   Multiply (-1) to odd line of k-space (undersampled) data. This is only
%   the case for 3T.
% 2010.09.28.
%   Take care old format of k-space data, K_alias_5d and K_nav_5d.
% 2010.09.30.
%   Erode mask only SMAPparams.vox_erode voxels.
% 2010.10.04.
%   Erode mask by round(SMAPparams.vox_erode/1.5) voxels as NAV recon does.
% 2010.11.03.
%   Adjust coronal reconstructed image object position.
% 2010.11.07.
%   Add silce by slice saving in '%% Save I_alias_sense_6d for all z' cell.
% 2010.12.07.
%   Multiply (-1) for 3T data. Adjust recon mask zeroing range.
% 2010.12.21.
%   Use I_alias_sense_sl%.2d_7d rather I_alias_sense_7d for reducing
%   memories being used.
% 2011.03.03.
%   Apply  'K_alias_shot_m(1:2:end,:) = K_alias_shot_m(1:2:end,:)*(-1);'
%   for 7T data when ZOOM is used.
% 2011.04.16.
%   Recon selectively only for b=0 images to check recon images, because
%   this is with no phase correction.
% 2012.01.23
%   Use improved recon algorithm as used in [s_recon_DWI2.m].
% 2012.02.01
%   Modify for the case when scan mode = 3D.
% 2012.05.17
%   To solve object position difference between SMAP and IMG data, first
%   recon central slice using original SMAP. Then compare image centroid
%   between the reconstructed IMG image and sum-of-square SMAP image. Then
%   shift IMG data and proceed recon for whole IMG image slices. The
%   difference in object position was found and verified using the 16
%   channel spine coil and phantom study
%   [s_batch__scan20120515_r4140__Jeong_999999.m]. This will be handled in
%   DATAFLAGparams.matchObjPosSMAP_IMG.
% 2012.05.23.
%   Take care of volume coil case (REFparams.filename is empty).
% 2012.06.16.
%   Take care of scan_mode == 3D.
% 2012.07.05.
%   Keep Head-to-Foot slice order of SMAP by setting
%   SMAPparams.flip_slice_order == false for 3D case. Then regenerate
%   ind_slice1 accordingly.
%
% Ha-Kyu



%% Recon IMG data for all DW and SLICES
warning off

% Load smap data.
if ~isempty(REFparams.filename)
    cd(saveDataDir_s)
    
    % Load some params struct. This may or may not include the necessary
    % fields based on how many times this script run.
    if ~isfield('DATAFLAGparams','matchObjPosDelta')
        cd(saveReconDir_s)
        if exist('DATAFLAGparams.matchObjPosDelta.mat','file')
            load('DATAFLAGparams.matchObjPosDelta.mat')
            DATAFLAGparams.matchObjPosDelta = matchObjPosDelta;
        end
    end
    if ~isfield('SMAPparams','filesize__I_smap_img_4d')
        cd(saveReconDir_s)
        if exist('SMAPparams.filesize__I_smap_img_4d.mat','file')
            load('SMAPparams.filesize__I_smap_img_4d.mat')
            SMAPparams.filesize__I_smap_img_4d = filesize__I_smap_img_4d;
        end
    end
    if ~isfield('SMAPparams','filesize__ref')
        cd(saveReconDir_s)
        if exist('SMAPparams.filesize__ref.mat','file')
            load('SMAPparams.filesize__ref.mat')
            SMAPparams.filesize__ref = filesize__ref;
        end
    end
    cd(saveDataDir_s)
            
    % Check if I_smap_img_4d is saved.
    % If I_smap_img_4d size is > 2GB, load this later.
    if exist('I_smap_img_4d.mat','file') || ...
            ~(SMAPparams.filesize__I_smap_img_4d > SMAPparams.filesize__ref)
        load  I_smap_img_4d
        load  I_mask_S_img_3d
    end
    cd(sharedDir_s)
    load  PSI
end

% Get parameters.
Ns = double(DWIparams.nSHOT);
R = double(DWIparams.SENSE_FACTOR);

% Match foldover direction.
if exist('I_smap_img_4d.mat','file')
    s_match_foldover_direction_IMG
end



%% Recon IMG
for dw_ori = 0:DWIparams.nDW_GRAD
    fprintf('\n\n')
    fprintf('RECON IMG DATA, ori[%d], RECONFLAGparams.psiMethod[%d]\n',dw_ori,RECONFLAGparams.psiMethod)
    
    % Load data.
    cd(saveDataDir_s)
    eval(sprintf('load  dw_data_ori%.2d__R%dS%d',dw_ori,R,Ns))
    
    % Get data size.
    % Take care of old format, K_alias_5d.
    if exist('K_alias_6d','var')
        if strcmpi(DWIparams.scan_mode,'3D')
            [ny,nx,nz,nchunk,nc,navg,ndw] = size(K_alias_6d);
            if DWIparams.nEC==2
                [nny,nnx,nnz,nnchunk,nnc,nnavg,nndw] = size(K_nav_6d);
                clear  K_nav_6d
            end
            % Report.
            fprintf('    K_alias_6d,        [ny,nx,nz,nchunk,nc,navg,ndw]          = [%d,%d,%d,%d,%d,%d,%d]\n', ...
                ny,nx,nz,nchunk,nc,navg,ndw)
        else
            [ny,nx,nz,nc,navg,ndw] = size(K_alias_6d);
            if DWIparams.nEC==2
                [nny,nnx,nnz,nnc,nnavg,nndw] = size(K_nav_6d);
                clear  K_nav_6d
            end
            nchunk = 1;
            nnchunk = 1;
            % Report.
            fprintf('    K_alias_6d,        [ny,nx,nz,nc,navg,ndw]          = [%d,%d,%d,%d,%d,%d]\n', ...
                ny,nx,nz,nc,navg,ndw)
        end
    else
        [ny,nx,nz,nc,navg] = size(K_alias_5d);
        [nny,nnx,nnz,nnc,nnavg] = size(K_nav_5d);
        ndw = 1;
        nndw = 1;
        K_alias_6d = reshape(K_alias_5d,[ny,nx,nz,nc,navg,ndw]);
        K_nav_6d = reshape(K_nav_5d,[nny,nnx,nnz,nnc,nnavg,nndw]);
        clear  K_alias_5d  K_nav_5d
        
        [ny,nx,nz,nc,navg,ndw] = size(K_alias_6d);
        if DWIparams.nEC==2
            [nny,nnx,nnz,nnc,nnavg,nndw] = size(K_nav_6d);
            clear  K_nav_6d
        end
        
        % Report.
        fprintf('    K_alias_6d,        [ny,nx,nz,nc,navg,ndw]          = [%d,%d,%d,%d,%d,%d]\n', ...
            ny,nx,nz,nc,navg,ndw)
    end
    pack
    
    
    %% Prepare IMG data
    if strcmpi(DWIparams.scan_mode,'3D')
        %I_alias_temp_7d = zeros(ny,nx,nz,nchunk,nc,navg,ndw,Ns,'single');
        H_alias_temp_7d = zeros(ny,nx,nz,nchunk,nc,navg,ndw,Ns,'single');
    else
        I_alias_temp_7d = zeros(ny,nx,nz,nc,navg,ndw,Ns,'single'); % each shot data
    end
    for ind_coil = 1:nc
        for ind_avg = 1:navg
            for ind_dw = 1:ndw
                for ind_chunk = 1:nchunk
                    
                    % For scan mode == 3D.
                    if strcmpi(DWIparams.scan_mode,'3D')
                        K_alias_temp_4d = zeros(ny,nx,nz,Ns,'single');
                    end
                    
                    for ind_sl = 1:nz
                        if strcmpi(DWIparams.scan_mode,'3D')
                            K_alias_m = K_alias_6d(:,:,ind_sl,ind_chunk,ind_coil,ind_avg,ind_dw);
                        else
                            K_alias_m = K_alias_6d(:,:,ind_sl,ind_coil,ind_avg,ind_dw);
                        end
                        
                        for ind_shot = 1:Ns
                            % Recon into original IMG space.
                            K_alias_shot_m = zeros(ny,nx);
                            
                            % Take each shot k-space data. From bottom to top of k-space image.
                            ky_shot_v = ny+1-ind_shot:-Ns:1;
                            K_alias_shot_m(ky_shot_v,:) = K_alias_m(ky_shot_v,:);
                            
                            %---------- START ----------
                            % Need this for 3T.
                            if GENparams.B0==30000
                                if ~strcmpi(GENparams.coilID,'SENSE-Breast-4')
                                    % [scan20100921_r2111__Smith_4117_3T] must multiply (-1).
                                    K_alias_shot_m(1:2:end,:) = K_alias_shot_m(1:2:end,:)*(-1); % original
                                    
                                    % [scan20101026_r2359__Smith_2402_3T] shouldn't multiply.
                                    % [scan20101207_r2359__Smith_4202_3T] should apply
                                end
                            elseif GENparams.B0==70000
                                if DWIparams.ZOOM==true
                                    K_alias_shot_m(1:2:end,:) = K_alias_shot_m(1:2:end,:)*(-1); % original
                                end
                                if strcmpi(GENparams.coilID,'RX-Intf-1_Quad-TR-1')
                                    K_alias_shot_m(1:2:end,:) = K_alias_shot_m(1:2:end,:)*(-1);
                                    % [scan20120420_r4140__Anderson_306871]
                                    % need *(-1)
                                end
                            end
                            %---------- END ----------
                            
                            % Just use original size data.
                            if strcmpi(DWIparams.scan_mode,'3D')
                                K_alias_temp_4d(:,:,ind_sl,ind_shot) = K_alias_shot_m;
                                
                                % No IFT: all dimensions in k-space. Later
                                % on this would be in hybrid-space as 2D
                                % in-plane image-space and k-space along
                                % kz.                                
                                %H_alias_temp_7d(:,:,ind_sl,ind_chunk,ind_coil,ind_avg,ind_dw,ind_shot) = ...
                                %    ift2(K_alias_shot_m);
                            else
                                I_alias_shot_m = ift2(K_alias_shot_m);
                                I_alias_temp_7d(:,:,ind_sl,ind_coil,ind_avg,ind_dw,ind_shot) = I_alias_shot_m;
                            end
                        end % ind_shot
                        
                    end % ind_sl
                    
                    % For 3D scan mode: Do 3D IFT.
                    if strcmpi(DWIparams.scan_mode,'3D')
                        for ind_shot = 1:Ns
                            
                            % 3D IFT: all dimensions in image-space.                            
                            %I_alias_temp_7d(:,:,:,ind_chunk,ind_coil,ind_avg,ind_dw,ind_shot) = ...
                            %    circshift(ift3(K_alias_temp_4d(:,:,:,ind_shot)),[0,0,-1]);
                            %I_alias_temp_7d(:,:,:,ind_chunk,ind_coil,ind_avg,ind_dw,ind_shot) = ...
                            %    ift3(K_alias_temp_4d(:,:,:,ind_shot));
                            
                            % kz shift of -1.
                            %H_alias_temp_7d(:,:,:,ind_chunk,ind_coil,ind_avg,ind_dw,ind_shot) = ...
                            %    circshift(H_alias_temp_7d(:,:,:,ind_chunk,ind_coil,ind_avg,ind_dw,ind_shot),[0,0,-1]);
                            
                            % image space shift then ft1 along z.
                            H_alias_temp_7d(:,:,:,ind_chunk,ind_coil,ind_avg,ind_dw,ind_shot) = ...
                                fftshift(ft1( circshift( ift3(K_alias_temp_4d(:,:,:,ind_shot)), [0 0 -1] ), 3 ), 3);
                            %*** use ift1(H_alias_temp_7d,3) for image
                            %*** space shifted data.
                            
                        end
                    end
                    
                end % ind_chunk
            end % ind_dw
        end % ind_avg
    end % ind_coil
    clear  K_alias_*  K_alias_shot_m  I_alias_shot_m  ky_shot_v
    
    
    % Rename.
    if strcmpi(DWIparams.scan_mode,'3D')
        H_alias_sense_7d = H_alias_temp_7d;
        clear  H_alias_temp_7d        
    else
        I_alias_sense_7d = I_alias_temp_7d;
        clear  I_alias_temp_7d
    end    
    
    % Get data size.
    if strcmpi(DWIparams.scan_mode,'3D')
        %[ny,nx,nz,nchunk,nc,navg,ndw,ns] = size(I_alias_sense_7d);
        [ny,nx,nz,nchunk,nc,navg,ndw,ns] = size(H_alias_sense_7d);
        % Report.
        fprintf('    H_alias_sense_7d,  [ny,nx,nz,nchunk,nc,navg,ndw,ns]       = [%d,%d,%d,%d,%d,%d,%d,%d]\n', ...
            ny,nx,nz,nchunk,nc,navg,ndw,ns)
    else
        [ny,nx,nz,nc,navg,ndw,ns] = size(I_alias_sense_7d);
        % Report.
        fprintf('    I_alias_sense_7d,  [ny,nx,nz,nc,navg,ndw,ns]       = [%d,%d,%d,%d,%d,%d,%d]\n', ...
            ny,nx,nz,nc,navg,ndw,ns)
    end
    
    
    % Pack.
    pack
    fprintf('    Memory packed\n')
    
    
    %% Save I_alias_sense_6d for all z
    fprintf('\n')
    cd(saveDataDir_s)
    
    % Slice by slice saving.
    % For 3D acquisition, add chunk.
    for ind_chunk = 1:nchunk
        for ind_slice = 1:nz
            if strcmpi(DWIparams.scan_mode,'3D')
%                 fname_s = sprintf('I_alias_sense_sl%.2d_chunk%.2d_ori%.2d__R%dS%d', ...
%                     ind_slice,ind_chunk,dw_ori,R,Ns);
%                 eval(sprintf('I_alias_sense_sl%.2d_chunk%.2d_7d = I_alias_sense_7d(:,:,ind_slice,ind_chunk,:,:,:,:);', ...
%                     ind_slice,ind_chunk))
%                 eval(sprintf('save  %s  I_alias_sense_sl%.2d_chunk%.2d_7d', ...
%                     fname_s,ind_slice,ind_chunk))

                fname_s = sprintf('H_alias_sense_sl%.2d_chunk%.2d_ori%.2d__R%dS%d', ...
                    ind_slice,ind_chunk,dw_ori,R,Ns);
                eval(sprintf('H_alias_sense_sl%.2d_chunk%.2d_7d = H_alias_sense_7d(:,:,ind_slice,ind_chunk,:,:,:,:);', ...
                    ind_slice,ind_chunk))
                eval(sprintf('save  %s  H_alias_sense_sl%.2d_chunk%.2d_7d', ...
                    fname_s,ind_slice,ind_chunk))
                
                fprintf('    [%s] is saved\n',fname_s)
                eval(sprintf('clear  H_alias_sense_sl%.2d*',ind_slice))
            else
                fname_s = sprintf('I_alias_sense_sl%.2d_ori%.2d__R%dS%d',ind_slice,dw_ori,R,Ns);
                eval(sprintf('I_alias_sense_sl%.2d_7d = I_alias_sense_7d(:,:,ind_slice,:,:,:,:);', ...
                    ind_slice))
                eval(sprintf('save  %s  I_alias_sense_sl%.2d_7d',fname_s,ind_slice))
                
                fprintf('    [%s] is saved\n',fname_s)
                eval(sprintf('clear  I_alias_sense_sl%.2d*',ind_slice))
            end            
        end
    end
    clear  I_alias_sense_7d  H_alias_sense_7d
    pack
    fprintf('    I(H)_alias_sense_7d is cleared and memory packed\n')
    
    
    %% Recon selectively for IMG
    if (dw_ori~=0) || (dw_ori~=1)
        fprintf('    dw_ori[%d], then pass\n',dw_ori)
        continue
    else
        fprintf('    dw_ori[%d], then recon\n',dw_ori)
    end
    
    
    %% Recon all slices IMG for all b and all z
    
    % IMG recon size.
    Ny = ny*R;
    Nx = nx;
    Nz = nz; % number of slices per chunk in 3D multichunk acq.
    Nc = nc;
    Nchunk = nchunk;
    
    % Construct noise correlation.
    %* This is done by the Cholesky factorization of the original
    %* noise correltation matrix to decorrelate coil sensitivity
    %* map and aliased image data.
    %* RECONFLAGparams.psiMethod of 1 and 2 generates the same
    %* results.
    
    %* Take care of Volume coil case (REFparams.filename is empty).
    if isempty(REFparams.filename)
        RECONFLAGparams.psiMethod = 0;
        PSI = 1;
    end
    
    if RECONFLAGparams.psiMethod==1
        %*** Method 1: E_prime(j) = PHI * E(j,:,s)
        L = chol(PSI,'lower');
        PP = kron(inv(L),eye(Ny/(R*Ns)));
        clear  L
    elseif RECONFLAGparams.psiMethod==2
        %*** Method 2: E_prime(j) = PHI * E(j)
        L = chol(PSI,'lower');
        PP = kron(inv(L),eye(Ny/R));
        clear  L
    elseif RECONFLAGparams.psiMethod==0
        %*** Method 0: Do nothing
        PP = kron(eye(size(PSI)),eye(Ny/(R*Ns)));
    else
        error('s_recon_IMG:main','Unknown RECONFLAGparams.psiMethod.')
    end
    
    
    % Reserve output.
    %* IMG is reconstructed for every 'avg', then keep recon data in 4D.
    if strcmpi(DWIparams.scan_mode,'3D')
        img_recon_5d = zeros(Ny,Nx,Nz,Nchunk,navg,ndw, 'single');
    else
        img_recon_5d = zeros(Ny,Nx,Nz,navg,ndw, 'single');
    end
    
    % Loop through slice.
    fprintf('\n')
    
    for ind_chunk = 1:nchunk 
        
        % $$$ Need to verify if SMAP slice position matches IMG kz position for each chunk $$$
        
        if DATAFLAGparams.matchObjPosSMAP_IMG==1
            slice_v = [floor(nz/2)+1,1:nz];
            slice_counter = -1;
        else
            slice_v = 1:nz;
            slice_counter = nan;
        end
        if ind_chunk==1
            DATAFLAGparams.matchObjPosDelta = []; % will be updated at below
        end
        
        for ind_slice = slice_v % slices/chunk for 3D case
            
            % Generate slice index for SMAP, since SMAP slices are at
            % 1:nkz*nchunk, while IMG data are at 1:nkz for each chunk.
            if strcmpi(DWIparams.scan_mode,'3D')
                %ind_slice1 = ind_slice + nz * (ind_chunk - 1);
                %ind_slice1 = (nz - ind_slice + 1) + nz * (ind_chunk - 1);
                
                % For [s_reslice_smap_v2.m] with unreversed SMAP slice
                % order prescribed in [s_get_SMAPparams.m] and
                % [f_coilsensemap.m]. 2012.07.05.
                %ind_slice1 = (nchunk-ind_chunk)*nz + ind_slice;
                
                % For [s_reslice_smap_v3.m].
                ind_slice1 = ind_slice;
            else
                ind_slice1 = ind_slice;
            end
            
            
            % Load I_smap_img_sl%.2d_4d for each slice.
            % If I_smap_img_4d size is > 2GB, load this.
            if (~exist('I_smap_img_4d.mat','file'))
                cd(saveDataDir_s)
                eval(sprintf('load  I_smap_img_chunk%.2d_4d',ind_chunk))
                eval(sprintf('load  I_mask_S_img_chunk%.2d_3d',ind_chunk))
                ind_slice2 = ind_slice1; % always 1
                
                s_match_foldover_direction_IMG
            else
                ind_slice2 = ind_slice1;
            end
            
                        
            % Control over ind_avg and ind_dw based on
            % DATAFLAGparams.matchObjPosSMAP_IMG for matching object position
            % between SMAP and IMG.
            if ~isnan(slice_counter)
                slice_counter = slice_counter+1;
                if slice_counter==0
                    avg_v = 1;
                    dw_v = 1;
                else
                    avg_v = 1:navg;
                    dw_v = 1:ndw;
                end
            else
                avg_v = 1:navg;
                dw_v = 1:ndw;
            end
            
            % Load I_alias_sense_sl%.2d_7d or
            % I_alias_sense_sl%.2d_chunk%.2d_7d.
            fprintf('\n')
            cd(saveDataDir_s)
            if strcmpi(DWIparams.scan_mode,'3D')
                
                % All dimensions are in image-space.  
%                 fname_s = sprintf('I_alias_sense_sl%.2d_chunk%.2d_ori%.2d__R%dS%d', ...
%                     ind_slice,ind_chunk,dw_ori,R,Ns);
%                 eval(sprintf('load  %s',fname_s))
%                 eval(sprintf('I_alias_sense_7d = I_alias_sense_sl%.2d_chunk%.2d_7d;', ...
%                     ind_slice,ind_chunk))
%                 eval(sprintf('clear  I_alias_sense_sl%.2d_chunk%.2d_7d',ind_slice,ind_chunk))

                % Rename Hybrid-space H_... data.
                fname_s = sprintf('H_alias_sense_sl%.2d_chunk%.2d_ori%.2d__R%dS%d', ...
                    ind_slice,ind_chunk,dw_ori,R,Ns);
                eval(sprintf('load  %s',fname_s))
                eval(sprintf('H_alias_sense_7d = H_alias_sense_sl%.2d_chunk%.2d_7d;', ...
                    ind_slice,ind_chunk))
                %I_alias_sense_7d = ift1(H_alias_sense_7d,3); % not right
                I_alias_sense_7d = H_alias_sense_7d;
                clear  H_alias_sense_7d
                eval(sprintf('clear  H_alias_sense_sl%.2d_chunk%.2d_7d',ind_slice,ind_chunk))
                
            else
                
                fname_s = sprintf('I_alias_sense_sl%.2d_ori%.2d__R%dS%d',ind_slice,dw_ori,R,Ns);
                eval(sprintf('load  %s',fname_s))
                eval(sprintf('I_alias_sense_7d = I_alias_sense_sl%.2d_7d;',ind_slice))
                eval(sprintf('clear  I_alias_sense_sl%.2d_7d',ind_slice))
                
            end
            fprintf('    [%s] is loaded and renamed for recon\n',fname_s)
            
            % Shift IMG data to match object position betweem SMAP and IMG.
            if ~isnan(slice_counter) && slice_counter~=0
                I_alias_sense_temp_7d = I_alias_sense_7d*0;
                if strcmpi(DWIparams.scan_mode,'3D')
                    [ny1,nx1,nz1,nchunk1,nc1,ndw1,nori1,nshot1] = size(I_alias_sense_7d);
                else
                    [ny1,nx1,nz1,nc1,ndw1,nori1,nshot1] = size(I_alias_sense_7d);
                    nchunk1 = 1; % unnecessary
                end
                
                % ===== slices and chunks are all 1 in I_alias_sense_7d =====
                
                for ind1=1:nc1
                    for ind2 = 1:ndw1
                        for ind3 = 1:nori1
                            for ind4 = 1:nshot1
                                if strcmpi(DWIparams.fat_shift_dir,'P')
                                    if strcmpi(DWIparams.scan_mode,'3D')
                                        I_alias_sense_temp_7d(:,:,:,:,ind1,ind2,ind3,ind4) = ...
                                            circshift(squeeze( ...
                                            I_alias_sense_7d(:,:,:,:,ind1,ind2,ind3,ind4)), ...
                                            [-round(delta),0,0,0]);
                                    else
                                        I_alias_sense_temp_7d(:,:,:,ind1,ind2,ind3,ind4) = ...
                                            circshift(squeeze( ...
                                            I_alias_sense_7d(:,:,:,ind1,ind2,ind3,ind4)), ...
                                            [-round(delta),0,0]);
                                    end
                                elseif strcmpi(DWIparams.fat_shift_dir,'L')
                                    if strcmpi(DWIparams.scan_mode,'3D')
                                        I_alias_sense_temp_7d(:,:,:,:,ind1,ind2,ind3,ind4) = ...
                                            circshift(squeeze( ...
                                            I_alias_sense_7d(:,:,:,:,ind1,ind2,ind3,ind4)), ...
                                            [0,-round(delta),0,0]);
                                    else
                                        I_alias_sense_temp_7d(:,:,:,ind1,ind2,ind3,ind4) = ...
                                            circshift(squeeze( ...
                                            I_alias_sense_7d(:,:,:,ind1,ind2,ind3,ind4)), ...
                                            [0,-round(delta),0]);
                                    end
                                else
                                    error('s_recon_IMG3:main','Unknown DWIparams.fat_shift_dir')
                                end
                            end
                        end
                    end
                end
                I_alias_sense_7d = I_alias_sense_temp_7d;
                clear  I_alias_sense_temp_7d
                fprintf('      I_alias_sense_7d is shifted by -delta=[%d]\n',-round(delta))
                pack
            end
            
            % Loop through avg, dw and col.
            for ind_avg = avg_v%1:navg
                for ind_dw = dw_v%1:ndw
                    
                    % Get mask.
                    if ~isempty(REFparams.filename)
                        I_mask_img_m = I_mask_S_img_3d(:,:,ind_slice1);
                    else
                        I_mask_img_m = ones(Ny,nx);
                    end
                    
                    
                    %-------------------- START --------------------
                    % Erode mask to eliminate very high intensity noisy voxel
                    % in SMAP.
                    
                    if ~isempty(REFparams.filename)
                        % This is volume coil case and erode must not be
                        % applied since the data is shifted (for Extremity-T/R
                        % coil) and exist at borders in the image.
                        
                        % This is default.
                        %erodeval = SMAPparams.vox_erode+3;
                        erodeval = SMAPparams.vox_erode;
                        I_mask_img_temp_m = I_mask_img_m;
                        %I_mask_img_temp_m([1:erodeval,Ny-erodeval+1:Ny],:) = 0;
                        I_mask_img_temp_m(:,[1:erodeval,Nx-erodeval+1:Nx]) = 0;
                        I_mask_img_temp_m = imerode(I_mask_img_temp_m,strel('disk',erodeval-1));
                        I_mask_img_m = I_mask_img_temp_m;
                        
                        if strcmpi(GENparams.coilID,'SENSE-Breast-4')
                            %erodeval = SMAPparams.vox_erode+3;
                            erodeval = round(SMAPparams.vox_dilate/1.5);
                            I_mask_img_temp_m = I_mask_img_m;
                            I_mask_img_temp_m(:,[1,end]) = 0;
                            I_mask_img_temp_m([1,end],:) = 0;
                            I_mask_img_temp_m = imerode(I_mask_img_temp_m,strel('disk',erodeval));
                            I_mask_img_m = I_mask_img_temp_m;
                        end
                        clear  I_mask_img_temp_m  erodeval
                    end
                    %-------------------- END --------------------
                    
                    
                    % Get image-space mixing matrix.
                    A_c = cell(1,Ns);
                    for ind_shot = (1:Ns)
                        A_m = [];
                        for ind_fold = 1:R*Ns
                            A_m = [A_m, eye(Ny/(R*Ns))*exp(1i*2*pi*(ind_fold-1)*(ind_shot-1)/Ns)];
                        end
                        A_c{1,ind_shot} = sparse(A_m);
                    end
                    
                    
                    % Generate aliased images for phase correction.
                    rstart_v = 1+Ny/(R*Ns)*(0:Nc*Ns-1);
                    g_m = zeros(Ny/(R*Ns)*Nc*Ns,Nx);
                    
                    
                    % Change aliased image order for coil and shots to incorporate
                    % noise correlation.
                    if RECONFLAGparams.psiMethod==0 || RECONFLAGparams.psiMethod==1
                        %*** Method 0 or 1.
                        for ind_shot = 1:Ns
                            for ind_coil = 1:Nc
                                ind = ind_coil+(ind_shot-1)*Nc;
                                r1 = rstart_v(ind);
                                r2 = r1+Ny/(R*Ns)-1;
                                
                                %g_m(r1:r2,:) = I_alias_sense_7d(1:Ny/(R*Ns), ...
                                %    :,ind_slice,ind_coil,ind_avg,ind_dw,ind_shot);
                                
                                % 2010.12.21.
                                % 2012.06.16: 3D
                                slice_idx = 1;
                                if strcmpi(DWIparams.scan_mode,'3D')
                                    chunk_idx = 1;
                                    g_m(r1:r2,:) = I_alias_sense_7d(1:Ny/(R*Ns), ...
                                        :,slice_idx,chunk_idx,ind_coil,ind_avg,ind_dw,ind_shot);
                                else
                                    g_m(r1:r2,:) = I_alias_sense_7d(1:Ny/(R*Ns), ...
                                        :,slice_idx,ind_coil,ind_avg,ind_dw,ind_shot);
                                end
                                
                                if ind_coil==1
                                    rr1 = r1;
                                end
                                if ind_coil==Nc
                                    rr2 = r2;
                                end
                            end
                            % Decorrelate g_m.
                            g_m(rr1:rr2,:) = PP * g_m(rr1:rr2,:);
                        end
                    elseif RECONFLAGparams.psiMethod==2
                        %*** Method 2.
                        for ind_coil = 1:Nc
                            for ind_shot = 1:Ns
                                ind = ind_shot+(ind_coil-1)*Ns;
                                r1 = rstart_v(ind);
                                r2 = r1+Ny/(R*Ns)-1;
                                
                                %g_m(r1:r2,:) = I_alias_sense_7d(1:Ny/(R*Ns), ...
                                %    :,ind_slice,ind_coil,ind_avg,ind_dw,ind_shot);
                                
                                % 2010.12.21.
                                slice_idx=1;
                                if strcmpi(DWIparams.scan_mode,'3D')
                                    chunk_idx = 1;
                                    g_m(r1:r2,:) = I_alias_sense_7d(1:Ny/(R*Ns), ...
                                        :,slice_idx,chunk_idx,ind_coil,ind_avg,ind_dw,ind_shot);
                                else
                                    g_m(r1:r2,:) = I_alias_sense_7d(1:Ny/(R*Ns), ...
                                        :,slice_idx,ind_coil,ind_avg,ind_dw,ind_shot);
                                end
                                
                            end
                        end
                        % Decorrelate g_m.
                        rr1 = rstart_v(1);
                        rr2 = r2;
                        g_m(rr1:rr2,:) = PP * g_m(rr1:rr2,:);
                    end
                    
                    
                    %-------------------- START --------------------
                    % SENSE recon without phase correction. Loop through columns.
                    
                    %aliasedMax = max(max(max(max(abs(I_alias_sense_7d(:,:,ind_slice,:,ind_avg,ind_dw,:))))));
                    
                    % 2010.12.21.
                    % 2012.06.16: 3D
                    slice_idx=1;
                    if strcmpi(DWIparams.scan_mode,'3D')
                        aliasedMax = max(max(max(max(abs(I_alias_sense_7d(:,:,slice_idx,:,:,ind_avg,ind_dw,:))))));
                    else
                        aliasedMax = max(max(max(max(abs(I_alias_sense_7d(:,:,slice_idx,:,ind_avg,ind_dw,:))))));
                    end
                    
                    % Reserve output.
                    img_recon_m = zeros(Ny,Nx);
                    
                    % Get EPI factor.
                    EPI_FACTOR = Ny/(R*Ns);
                    
                    fprintf('\n')
                    fprintf('    STARTING IMG recon DW_ORI[%d], AVG[%d], DW[%d], CHUNK[%d], SLICE[%d]\n', ...
                        dw_ori,ind_avg,ind_dw,ind_chunk,ind_slice)
                    
                    tic
                    for col = 1:nx
                        % Use g_m.
                        aliased_v = g_m(:,col);
                        
                        % Skip this column if there's little signal:
                        if (max(abs(aliased_v)) < 0.01*aliasedMax)
                            continue
                        end
                        if ~any(I_mask_img_m(:,col)==1)
                            continue
                        end
                        
                        % Get NAV phase in current column for each shot. Must be
                        % identity matrix here.
                        %for ind_shot = 1:ns
                        %    eval(sprintf('P%d_m = eye(Ny);',ind_shot,ind_shot))
                        %end
                        
                        % Calculate net mixing matrix. Method 1 and 2 work ok.
                        % Method 0 doesn' use PSI.
                        if RECONFLAGparams.psiMethod==0 || ...
                                RECONFLAGparams.psiMethod==1
                            
                            % Method 0 or 1.
                            E_m = zeros(Ny/(R*Ns)*Nc*Ns,Ny);
                            rstart_v = (1+Ny/(R*Ns)*(0:Nc*Ns-1));
                            for ind_shot = 1:Ns
                                P_v = ones(Ny,1);
                                A_m = A_c{ind_shot};
                                for ind_coil = 1:Nc
                                    ind = ind_coil+(ind_shot-1)*Nc;
                                    r1 = rstart_v(ind);
                                    r2 = r1+Ny/(R*Ns)-1;
                                    
                                    % Original
                                    %E_m(r1:r2,:) = A_m*diag(I_smap_img_4d(:,col,ind_slice,ind_coil))*P_m;
                                    
                                    % Modified
                                    if ~isempty(REFparams.filename)
                                        s_p_v = double(I_smap_img_4d(:,col,ind_slice2,ind_coil).*P_v);
                                    else
                                        s_p_v = double(ones(Ny,1).*P_v);
                                    end
                                    ASP_m = A_m;
                                    for indr = 1:EPI_FACTOR
                                        col_begin = indr;
                                        a_v = A_m(indr,:);
                                        ASP_m(indr,col_begin:EPI_FACTOR:end) = ...
                                            a_v(col_begin:EPI_FACTOR:end).*...
                                            s_p_v(col_begin:EPI_FACTOR:end).';
                                    end
                                    E_m(r1:r2,:) = ASP_m;
                                    
                                    % Reserve indices for PP.
                                    if ind_coil==1
                                        r00 = r1;
                                    end
                                    if ind_coil==Nc
                                        r11 = r2;
                                    end
                                end
                                
                                % Decorrelate E_m.
                                E_m(r00:r11,:) = sparse(PP)*sparse(E_m(r00:r11,:));
                            end
                            
                        elseif RECONFLAGparams.psiMethod==2
                            
                            % Method 2.
                            E_m = zeros(Ny/(R*Ns)*Nc*Ns,Ny);
                            rstart_v = (1+Ny/(R*Ns)*(0:Nc*Ns-1));
                            for ind_coil = 1:Nc
                                for ind_shot = 1:Ns
                                    P_v = ones(Ny,1);
                                    
                                    A_m = A_c{ind_shot};
                                    ind = ind_shot + (ind_coil-1)*Ns;
                                    r1 = rstart_v(ind);
                                    r2 = r1+Ny/(R*Ns)-1;
                                    
                                    % Original.
                                    %E_m(r1:r2,:) = A_m*diag(I_smap_img_4d(:,col,ind_slice,ind_coil))*P_m;
                                    
                                    % Modified.
                                    s_p_v = double(I_smap_img_4d(:,col,ind_slice2,ind_coil).*P_v);
                                    ASP_m = A_m;
                                    for indr = 1:EPI_FACTOR
                                        col_begin = indr;
                                        a_v = A_m(indr,:);
                                        ASP_m(indr,col_begin:EPI_FACTOR:end) = ...
                                            a_v(col_begin:EPI_FACTOR:end).*...
                                            s_p_v(col_begin:EPI_FACTOR:end).';
                                    end
                                    E_m(r1:r2,:) = ASP_m;
                                end
                            end
                            
                            % Decorrelate E_m.
                            %* PP is multiplied repeatedly over averages
                            E_m(rstart_v(1):r2,:) = sparse(PP)*sparse(E_m(rstart_v(1):r2,:));
                            
                        end % if RECONFLAGparams.psiMethod
                        
                        % Recon image.
                        %cImage_v = E_m \ aliased_v;
                        A = E_m'*E_m;
                        b = E_m'*aliased_v;
                        S = sparse(A);
                        cImage_v = S \ b;
                        cImage_v(isnan(cImage_v))=0;
                        
                        % Save recon image
                        img_recon_m(:, col) = cImage_v;
                        
                    end % for col
                    clear  cImage_v  aliased_v  E_m  P*_m  A  b  S  ASP_m  A_m
                    
                    % Mask recon image.
                    %m = img_recon_m.*imerode(I_mask_img_m,strel('disk',19));
                    m = img_recon_m;%.*imerode(I_mask_img_m,strel('disk',1));
                    %showimage(jet,imrotate(abs(m),90))
                    %showimage(jet,abs(m))
                    
                    % Match object position between SMAP and IMG data.
                    if slice_counter==0
                        
                        % Method 1 --------------------------------------
                        % Match centroid between SMAP and IMG (recon).
                        
                        % Calculate number of shifted voxels.
                        smap = squeeze(I_smap_img_4d(:,:,ind_slice2,:));
                        smap = sqrt(sum(abs(smap).^2,3));
                        mask_smap = abs(smap) > 0;
                        
                        bw = imerode(abs(m),strel('disk',4));
                        m1 = abs(m).*bw;
                        m1 = m1/max(m1(:));
                        th = graythresh(m1);
                        mask_img = imfill(m1 > th*0.1,'holes');
                        %mask_img = abs(m) > 0;
                        
                        stat_img = regionprops(mask_img,'Area','Centroid');
                        stat_smap = regionprops(mask_smap,'Area','Centroid');
                        
                        area_smap = [stat_smap(:).Area];
                        ind = find(area_smap==max(area_smap));
                        cent_smap = stat_smap(ind).Centroid; %[col,row]
                        
                        area_img = [stat_img(:).Area];
                        ind = find(area_img==max(area_img));
                        cent_img = stat_img(ind).Centroid;
                        
                        if strcmpi(DWIparams.fat_shift_dir,'P')
                            % Use centroid distance.
                            %delta = cent_img(2) - cent_smap(2); % (2) is along readout direction
                            
                            % Use slice OFFC.
                            delta = DWIparams.Offc.AP / RECONparams.ACQREC.ACQ_voxel_MPS(2);
                        elseif strcmpi(DWIparams.fat_shift_dir,'L')
                            % Use centroid distance.
                            %delta = cent_img(1) - cent_smap(1);
                            
                            % Use slice OFFC.
                            %delta = DWIparams.Offc.RL / RECONparams.ACQREC.ACQ_voxel_MPS(1);
                            delta = abs(DWIparams.Offc.AP) / RECONparams.ACQREC.ACQ_voxel_MPS(1);
                            %*** Not DWIPparams.Offc.AP ***
                        else
                            error('s_recon_IMG3:main','Unknown DWIparams.fat_shift_dir')
                        end
                        DATAFLAGparams.matchObjPosDelta(ind_chunk) = delta;
                        
                        fprintf('      shift distance calculated for chunk[%d] = %f\n',...
                            ind_chunk,delta)
                    end
                    
                    % Keep output.
                    % If slice_counter==1, this recon is for calculating
                    % centroid only, the don't keep it.
                    if slice_counter~=0
                        if strcmpi(DWIparams.scan_mode,'3D')
                            img_recon_5d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw) = m;
                        else
                            img_recon_5d(:,:,ind_slice,ind_avg,ind_dw) = m;
                        end
                    end
                    
                    %fprintf('\n')
                    t = toc;
                    fprintf('    Finished unfolding IMG in [%.3f]sec\n',t)
                    
                    % Show results.
%                     showimage(jet,abs(m),angle(m), ...
%                         sprintf('RECONFLAGparams.psiMethod[%d], DATAFLAGparams.epiCorrFitMethod[%d]', ...
%                         RECONFLAGparams.psiMethod,DATAFLAGparams.epiCorrFitMethod),' ')
                    
                end % for ind_dw
            end % for ind_avg
            
            % Clear used data.
            clear  I_alias_sense_7d
            
        end % for ind_slice
        
        % IFT along kz for 3D acquisition.
        if strcmpi(DWIparams.scan_mode,'3D')
            for ind1 = 1:ndw
                for ind2 = 1:navg
%                     img_recon_5d(:,:,:,ind_chunk,ind2,ind1) = ...
%                         fftshift(ift1(squeeze(img_recon_5d(:,:,:,ind_chunk,ind2,ind1)),3),3);
                    img_recon_5d(:,:,:,ind_chunk,ind2,ind1) = ...
                        ift1(squeeze(img_recon_5d(:,:,:,ind_chunk,ind2,ind1)),3);
                end
            end
        end
        fprintf('\n*** IFT along kz is performed for chunk [%d] ***\n',ind_chunk)
        
    end % for ind_chunk
    pack
    
    
    %% Save runtime parameters
    % This is for updated DATAFLAGparams.
    cd(saveReconDir_s)
    save('params','DATAFLAGparams','-append')
    fprintf('    DATAFLAGparams.matchObjPosDelta is updated\n')
    matchObjPosDelta = DATAFLAGparams.matchObjPosDelta;
    save DATAFLAGparams.matchObjPosDelta.mat  matchObjPosDelta
    
    
    %% Adjust reconstructed object position variable with foldover and fatshift.
    
    % Transverse:   A-L-P-R from North-East-South-West in image-plane.
    % Sagittal:     A-F-P-H from North-East-South-West in image-plane.
    % Coronal:      R-H-L-F from North-East-South-West in image-plane.
    
    % This is for no-phase correction data. So the object position here
    % actually doesn't matter for the final image reconstruction in
    % [s_recon_DWI.m] with phase correction. The adjustment made here is
    % for the comparison between with and without phase correction because
    % this data can also be resized in [s_resize_DWI.m] for comparison. Of
    % course the phase-corrected data goes with the same procedure in
    % [s_resize_DWI.m].
    
    %---------------------------------------------------------
    % Transverse.
    if strcmpi(DWIparams.slice_orientation,'transverse')
        if strcmpi(DWIparams.fold_over_dir,'AP')
            
            if strcmpi(DWIparams.fat_shift_dir,'P')
                % Do nothing. It is already in the right object position mentioned above
                
            elseif strcmpi(DWIparams.fat_shift_dir,'A')
                if strcmpi(DWIparams.scan_mode,'3D')
                    img_recon_temp_5d = zeros(Nx,Ny,Nz,Nchunk,navg,ndw, 'single');
                else
                    img_recon_temp_5d = zeros(Nx,Ny,Nz,navg,ndw, 'single');
                end
                for ind_avg = 1:navg
                    for ind_dw = 1:ndw
                        for ind_slice = 1:Nz
                            if strcmpi(DWIparams.scan_mode,'3D')
                                img_recon_temp_5d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw) = ...
                                    flipud(fliplr(img_recon_5d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw)));
                            else
                                img_recon_temp_5d(:,:,ind_slice,ind_avg,ind_dw) = ...
                                    flipud(fliplr(img_recon_5d(:,:,ind_slice,ind_avg,ind_dw)));
                            end
                        end
                    end
                end
                clear  img_recon_5d
                img_recon_5d = img_recon_temp_5d;
                clear  img_recon_temp_5d
            end
            
        elseif strcmpi(DWIparams.fold_over_dir,'RL')
            
            if strcmpi(DWIparams.fat_shift_dir,'L')
                if strcmpi(DWIparams.scan_mode,'3D')
                    img_recon_temp_5d = zeros(Nx,Ny,Nz,Nchunk,navg,ndw, 'single');
                else
                    img_recon_temp_5d = zeros(Nx,Ny,Nz,navg,ndw, 'single');
                end
                for ind_avg = 1:navg
                    for ind_dw = 1:ndw
                        for ind_slice = 1:Nz
                            if strcmpi(DWIparams.scan_mode,'3D')
                                img_recon_temp_5d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw) = ...
                                    flipud(permute(img_recon_5d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw),[2,1]));
                            else
                                img_recon_temp_5d(:,:,ind_slice,ind_avg,ind_dw) = ...
                                    flipud(permute(img_recon_5d(:,:,ind_slice,ind_avg,ind_dw),[2,1]));
                            end
                        end
                    end
                end
                clear  img_recon_5d
                img_recon_5d = img_recon_temp_5d;
                clear  img_recon_temp_5d
                
            elseif strcmpi(DWIparams.fat_shift_dir,'R')
                if strcmpi(DWIparams.scan_mode,'3D')
                    img_recon_temp_5d = zeros(Nx,Ny,Nz,Nchunk,navg,ndw, 'single');
                else
                    img_recon_temp_5d = zeros(Nx,Ny,Nz,navg,ndw, 'single');
                end
                for ind_avg = 1:navg
                    for ind_dw = 1:ndw
                        for ind_slice = 1:Nz
                            if strcmpi(DWIparams.scan_mode,'3D')
                                img_recon_temp_5d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw) = ...
                                    fliplr(permute(img_recon_5d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw),[2,1]));
                            else
                                img_recon_temp_5d(:,:,ind_slice,ind_avg,ind_dw) = ...
                                    fliplr(permute(img_recon_5d(:,:,ind_slice,ind_avg,ind_dw),[2,1]));
                            end
                        end
                    end
                end
                clear  img_recon_5d
                img_recon_5d = img_recon_temp_5d;
                clear  img_recon_temp_5d
                
            else
                error('s_recon_IMG:main','Unknown DWI fat_shift_dir')
            end
        else
            error('s_recon_IMG:main','Unknown DWI fold_over_dir')
        end
    end
    
    
    %---------------------------------------------------------
    % Sagittal.
    if strcmpi(DWIparams.slice_orientation,'sagittal')
        if strcmpi(DWIparams.fold_over_dir,'AP')
            
            if strcmpi(DWIparams.fat_shift_dir,'P')
                if strcmpi(DWIparams.scan_mode,'3D')
                    img_recon_temp_5d = zeros(Ny,Nx,Nz,Nchunk,navg,ndw, 'single');
                else
                    img_recon_temp_5d = zeros(Ny,Nx,Nz,navg,ndw, 'single');
                end
                for ind_avg = 1:navg
                    for ind_dw = 1:ndw
                        for ind_slice = 1:Nz
                            if strcmpi(DWIparams.scan_mode,'3D')
                                img_recon_temp_5d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw) = ...
                                    fliplr(img_recon_5d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw));
                            else
                                img_recon_temp_5d(:,:,ind_slice,ind_avg,ind_dw) = ...
                                    fliplr(img_recon_5d(:,:,ind_slice,ind_avg,ind_dw));
                            end
                        end
                    end
                end
                clear  img_recon_5d
                img_recon_5d = img_recon_temp_5d;
                clear  img_recon_temp_5d
                
            elseif strcmpi(DWIparams.fat_shift_dir,'A')
                img_recon_temp_5d = zeros(Ny,Nx,Nz,navg,ndw, 'single');
                for ind_avg = 1:navg
                    for ind_dw = 1:ndw
                        for ind_slice = 1:Nz
                            if strcmpi(DWIparams.scan_mode,'3D')
                                img_recon_temp_5d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw) = ...
                                    flipud(fliplr(img_recon_5d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw)));
                            else
                                img_recon_temp_5d(:,:,ind_slice,ind_avg,ind_dw) = ...
                                    flipud(fliplr(img_recon_5d(:,:,ind_slice,ind_avg,ind_dw)));
                            end
                        end
                    end
                end
                clear  img_recon_5d
                img_recon_5d = img_recon_temp_5d;
                clear  img_recon_temp_5d
                
            else
                error('s_recon_IMG:main','Unknown DWI fat_shift_dir')
            end
            
        elseif strcmpi(DWIparams.fold_over_dir,'FH') % this is not a good fold-over direction
            
            if strcmpi(DWIparams.fat_shift_dir,'F')
                warning('s_recon_IMG:main','SAGITTAL FH-F is not setup yet')
            elseif strcmpi(DWIparams.fat_shift_dir,'H')
                warning('s_recon_IMG:main','SAGITTAL FH-H is not setup yet')
            else
                error('s_recon_IMG:main','Unknown DWI fat_shift_dir')
            end
            
        else
            error('s_recon_IMG:main','Unknown DWI fold_over_dir')
        end
    end
    
    
    %---------------------------------------------------------
    % Coronal.
    if strcmpi(DWIparams.slice_orientation,'coronal')
        if strcmpi(DWIparams.fold_over_dir,'RL')
            
            if strcmpi(DWIparams.fat_shift_dir,'L')
                % Do nothing. It is already in the right object position mentioned above
                
            elseif strcmpi(DWIparams.fat_shift_dir,'R')
                warning('s_recon_IMG:main','CORONAL RL-R is not setup yet')
                
            else
                error('s_recon_IMG:main','Unknown DWI fat_shift_dir')
            end
            
        elseif strcmpi(DWIparams.fold_over_dir,'FH') % this is not a good fold-over direction
            
            if strcmpi(DWIparams.fat_shift_dir,'F')
                warning('s_recon_IMG:main','CORONAL FH-F is not setup yet')
                
            elseif strcmpi(DWIparams.fat_shift_dir,'H')
                warning('s_recon_IMG:main','CORONAL FH-H is not setup yet')
                
            else
                error('s_recon_IMG:main','Unknown DWI fat_shift_dir')
            end
            
        else
            error('s_recon_IMG:main','Unknown DWI fold_over_dir')
        end
    end
    
    
    %% Save IMG data.
    cd(saveReconDir_s)
    fname_s = sprintf('recon_img_ori%.2d__R%dS%d',dw_ori,R,Ns);
    eval(sprintf('save  %s  img_recon_5d',fname_s))
    fprintf('\n')
    fprintf('    [%s] is saved at \n',fname_s)
    fprintf('    [%s]\n',saveReconDir_s)
    
    % Clear data.
    clear  img_recon_m  g_m  aliasedMax  I_mask_img_m  I_alias_sense_7d
    clear  fname_s  m  img_recon_*  mask  PP  K_alias*  K_nav*
    pack
    
end % for dw_ori
clear  I_smap_img_4d  I_mask_S_img_3d

% Pack memory.
pack

fprintf('\n\n')



%% END





