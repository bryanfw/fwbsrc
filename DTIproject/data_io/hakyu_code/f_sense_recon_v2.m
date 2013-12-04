
function [recon_6d,gfactor_6d] = f_sense_recon_v2(sense_7d, smap_4d, psi_m, ...
    mask_3d, R, SMAPparams, GENparams, DATAFLAGparams, REFparams, DWIparams, dw_ori)
%[f_sense_recon_v2] performs conventional sense reconstruction.
%
% USAGE
%   [recon_6d,gfactor_6d] = f_sense_recon_v2( ...
%       sense_7d, smap_4d, psi_m, mask_3d, R, SMAPparams, GENparams, ...
%       DATAFLAGparams, REFparams, dw_ori)
%
% INPUT
%   sense_7d:   Aliased image data in [y,x,z,coil,avg,dw,shot]
%   smap_4d:    SMAP in [y,x,z,coil]
%   psi_m:      Coil noise correlation, [coil,coil]
%   mask_3d:    Mask in [y,x,z]
%   R:          SENSE factor
%   SMAPparams: SMAPparams
%   DATAFLAGparams: To match object position between SMAP and NAV
%   REFparams: To check if this is Volume Coil (REFparams.filename is
%              empty)
%   DWIparams: To check the fat_shift_dir.
%
% OUTPUT
%   recon_6d:   NAV recon image in [y,x,z,avg,dw,shot]
%   gfactor_6d: NAV g-factor image in [y,x,z,avg,dw,shot]
%
%
% Last modified
% 2009.11.24. Take input SMAPparams.
% 2010.08.10
%   This is generated from [sense_recon_v10.m].
% 2010.09.21.
%   Consider multiple non-zero b-values.
% 2010.09.30.
%   Erode mask only SMAPparams.vox_erode voxels.
% 2010.10.04.
%   Add '% Check if there is zero column in S_m.' to remove NaN rows in
%   recon image.
%   Erode mask by round(SMAPparams.vox_erode/1.5) voxels.
% 2010.10.06.
%   Take input of GENparams.
% 2010.12.07.
%   Adjust recon mask zeroing range.
% 2011.12.12.
%   Add dw_ori just to show info while reconstruction is going on.
% 2012.05.17.
%   Add DATAFLAGparams for matching object position betwen SMAP and NAV
%   (IMG) data.
% 2012.05.23.
%   Take care of Volume Coil case (REFparams.filename is empty).
% 2012.05.29.
%   DWIparams is added for input.
% 2012.06.16.
%   Take care of scan_mode = 3D case.
% 2012.07.05.
%   Keep Head-to-Foot slice order of SMAP by setting
%   SMAPparams.flip_slice_order == false for 3D case. Then regenerate
%   ind_slice1 accordingly.
% 2012.07.18.
%   gfactor_6d is not generated to reduce the memory and hdd space.
%
% Ha-Kyu



%% Check input
warning off

% TEST
% sense_7d = I_nav_sense_7d;
% smap_4d = I_smap_nav_4d;
% psi_m = PSI;
% mask_3d = I_mask_S_nav_3d;

% Check input
if length(size(sense_7d))<7
    error('Input [sense_7d] must be greater than 7-D.')
end
%if length(size(smap_4d))~=4
%    error('Input [smap_4d] must be 4-D.')
%end
%if length(size(mask_3d))~=3
%    error('Input [mask_3d] must be 3-D.')
%end

% Get size.
if strcmpi(DWIparams.scan_mode,'3D')
    [nny,nnx,nnz,nnchunk,nnc,nnavg,nndw,nns] = size(sense_7d);
else
    [nny,nnx,nnz,nnc,nnavg,nndw,nns] = size(sense_7d);
    nnchunk = 1;
end
if numel(smap_4d)==1 && all(smap_4d(:)==0) && ...
        (SMAPparams.filesize__I_smap_img_4d > SMAPparams.filesize__ref)
    Ny = nny * DWIparams.SENSE_FACTOR;
    Nx = nnx;
    %Nz = nnz*nnchunk;
    Nz = nnz; % [s_reslice_map_3d_v3.m]
    Nc = nnc;
else
    [Ny,Nx,Nz,Nc] = size(smap_4d);
end

% Calculate inverse of psi_m.
psi_inv_m = inv(psi_m);

% Shift NAV data to match object position between SMAP and NAV.
if DATAFLAGparams.matchObjPosSMAP_IMG==1
    %sense_temp_7d = sense_7d*0;
    delta_v = DATAFLAGparams.matchObjPosDelta;
    for ind0 = 1:nnchunk
        delta = delta_v(ind0);
    for ind1 = 1:nnc
        for ind2 = 1:nnavg
            for ind3 = 1:nndw
                for ind4 = 1:nns
                    if strcmpi(DWIparams.fat_shift_dir,'P')
                        if strcmpi(DWIparams.scan_mode,'3D')
                            sense_7d(:,:,:,ind0,ind1,ind2,ind3,ind4) = ...
                                circshift(squeeze(sense_7d(:,:,:,ind0,ind1,ind2,ind3,ind4)), ...
                                [-round(delta),0,0]);
                        else
                            sense_7d(:,:,:,ind1,ind2,ind3,ind4) = ...
                                circshift(squeeze(sense_7d(:,:,:,ind1,ind2,ind3,ind4)), ...
                                [-round(delta),0,0]);
                        end                        
                    elseif strcmpi(DWIparams.fat_shift_dir,'L')
                        if strcmpi(DWIparams.scan_mode,'3D')
                            sense_7d(:,:,:,ind0,ind1,ind2,ind3,ind4) = ...
                                circshift(squeeze(sense_7d(:,:,:,ind0,ind1,ind2,ind3,ind4)), ...
                                [0,-round(delta),0]);
                        else
                            sense_7d(:,:,:,ind1,ind2,ind3,ind4) = ...
                                circshift(squeeze(sense_7d(:,:,:,ind1,ind2,ind3,ind4)), ...
                                [0,-round(delta),0]);
                        end
                    else
                        error('f_sense_recon_v2:main','Unknown DWIparams.fat_shift_dir')
                    end
                end
            end
        end
    end
    end
    fprintf('        NAV data shifted\n')
    %sense_7d = sense_temp_7d;
    clear  sense_temp_7d
    pack
end



%% SENSE recon

% Reserve output.
recon_m = zeros(Ny,nnx);
%gmap_m = zeros(Ny,nnx);
if strcmpi(DWIparams.scan_mode,'3D')
    recon_6d = zeros(Ny,nnx,nnz,nnchunk,nnavg,nndw,nns);
    %gfactor_6d = zeros(Ny,nnx,nnz,nnchunk,nnavg,nndw,nns,'single');
    gfactor_6d = 0;
else
    recon_6d = zeros(Ny,nnx,nnz,nnavg,nndw,nns);
    %gfactor_6d = zeros(Ny,nnx,nnz,nnavg,nndw,nns,'single');
    gfactor_6d = 0;
end

fprintf('\n')
for ind_avg = 1:nnavg
    for ind_dw = 1:nndw
        for ind_chunk = 1:nnchunk
        for ind_slice = 1:nnz
            % Generate slice index for SMAP, since SMAP slices are at
            % 1:knz*nchunk, while NAV data are at 1:nz for each chunk.            
            if strcmpi(DWIparams.scan_mode,'3D')
                %ind_slice1 = ind_slice + nnz * (ind_chunk - 1);
                %ind_slice1 = (nnz - ind_slice + 1) + nnz * (ind_chunk - 1);
                
                % For [s_reslice_smap_v2.m] with unreversed SMAP slice
                % order prescribed in [s_get_SMAPparams.m] and
                % [f_coilsensemap.m]. 2012.07.05.
                %ind_slice1 = (nnchunk-ind_chunk)*nnz + ind_slice;
                
                % For [s_reslice_smap_v3.m].
                ind_slice1 = ind_slice;
            else
                ind_slice1 = ind_slice;
            end
            
            
            % Load I_smap_img_sl%.2d_4d for each slice.
            % If I_smap_img_4d size is > 2GB, load this.
            saveDataDir_s = DATAFLAGparams.saveDataDir_s;
            cd(saveDataDir_s)
            if (~exist('I_smap_img_4d.mat','file'))
                eval(sprintf('load  I_smap_img_chunk%.2d_4d',ind_chunk))
                eval(sprintf('load  I_mask_S_img_chunk%.2d_3d',ind_chunk))
                ind_slice2 = ind_slice1; % always 1
                
                s_match_foldover_direction
                
                smap_4d = I_smap_img_4d; clear I_smap_img_4d
                mask_3d = I_mask_S_img_3d; clear I_mask_S_img_3d
            else
                ind_slice2 = ind_slice1;
            end
            
            
            I_mask_nav_m = mask_3d(:,:,ind_slice1);
            
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
                I_mask_nav_temp_m = I_mask_nav_m;
                %I_mask_nav_temp_m([1:erodeval,Ny-erodeval+1:Ny],:) = 0;
                I_mask_nav_temp_m(:,[1:erodeval,Nx-erodeval+1:Nx]) = 0;
                I_mask_nav_temp_m = imerode(I_mask_nav_temp_m,strel('disk',erodeval-1));
                I_mask_nav_m = I_mask_nav_temp_m;
                
                if strcmpi(GENparams.coilID,'SENSE-Breast-4')
                    %erodeval = SMAPparams.vox_erode+3;
                    erodeval = round(SMAPparams.vox_dilate/1.5);
                    I_mask_nav_temp_m = I_mask_nav_m;
                    I_mask_nav_temp_m(:,[1,end]) = 0;
                    I_mask_nav_temp_m([1,end],:) = 0;
                    I_mask_nav_temp_m = imerode(I_mask_nav_temp_m,strel('disk',erodeval));
                    I_mask_nav_m = I_mask_nav_temp_m;
                end
                clear  I_mask_img_temp_m  erodeval
            end
            %-------------------- END --------------------
            
            for ind_shot = 1:nns
                fprintf('\n')
                fprintf('      STARTING NAV recon ORI[%d],AVG[%d],DW[%d],CHUNK[%d],SLICE[%d],SHOT[%d]\n', ...
                    dw_ori,ind_avg,ind_dw,ind_chunk,ind_slice,ind_shot)
                
                tic
                for indx = 1:nnx
                    for indy = 1:nny
                        rho_v = indy:nny:Ny;
                        
                        % Check coordinates.
                        if length(find(rho_v>0))~=R
                            fprintf('rho_v, '),disp(rho_v)
                            error('f_sense_recon:main','Check the number of rho_v.')
                        end
                        rho_v = rho_v(rho_v>0);
                        
                        % Calculate mask values.
                        if any(I_mask_nav_m(rho_v,indx) == 1)
                            
                            % Take rho_v which doesn't have zero values.
                            rho_v = rho_v(I_mask_nav_m(rho_v,indx) == 1);
                            
                            % Calculate unfolded voxel intensities.
                            if strcmpi(DWIparams.scan_mode,'3D')
                                m_v = squeeze(sense_7d(indy,indx,ind_slice,ind_chunk,:,ind_avg,ind_dw,ind_shot));
                            else
                                m_v = squeeze(sense_7d(indy,indx,ind_slice,:,ind_avg,ind_dw,ind_shot));
                            end
                            S_m = squeeze(smap_4d(rho_v,indx,ind_slice2,:));
                            if size(S_m,1)~=Nc
                                S_m = S_m.';
                            end
                            
                            % Check if there is zero column in S_m.
                            S_temp_m = [];
                            rho_temp_v = [];
                            for col_ind = 1:size(S_m,2)
                                s_v = S_m(:,col_ind);
                                if ~all(abs(s_v)==0)
                                    S_temp_m = [S_temp_m, s_v];
                                    rho_temp_v = [rho_temp_v,rho_v(col_ind)];
                                end
                            end
                            if isempty(S_temp_m) || isempty(rho_temp_v)
                                continue
                            end
                            S_m = S_temp_m;  clear  S_temp_m
                            rho_v = rho_temp_v;  clear  rho_temp_v
                            
                            % Calculate unfolded image values.
                            M = S_m'*psi_inv_m*S_m;
                            b_v = S_m'*psi_inv_m* m_v;
                            %M = S_m'*S_m;
                            %b_v = S_m'* m_v;
                            %M = S_m'*(psi_m\S_m);
                            %b_v = S_m'*(psi_m\m_v);
                            %a_v = M\b_v;
                            [q,r] = qr(M);
                            a_v = r \ (q'*b_v);
                            
                            % Output.
                            recon_m(rho_v,indx) = a_v;
                            %gmap_m(rho_v,indx) = diag( sqrt( inv(M) .* M ) );
                        end
                    end % for indy
                end % for indx
                %m = recon_m.*imerode(I_mask_nav_m,strel('disk',16));
                %showimage(jet,abs(m),angle(m))
                
                if strcmpi(DWIparams.scan_mode,'3D')
                    recon_6d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw,ind_shot) = recon_m.*I_mask_nav_m;
                    %gfactor_6d(:,:,ind_slice,ind_chunk,ind_avg,ind_dw,ind_shot) = gmap_m.*I_mask_nav_m;
                else
                    recon_6d(:,:,ind_slice,ind_avg,ind_dw,ind_shot) = recon_m.*I_mask_nav_m;
                    %gfactor_6d(:,:,ind_slice,ind_avg,ind_dw,ind_shot) = gmap_m.*I_mask_nav_m;
                end
                
                t = toc;
                fprintf('      Finished unfolding NAV in [%.3f]sec\n',t)
                %showimage(jet,abs(recon_m.*I_mask_nav_m))
                
            end % for ind_shot
        end % for ind_slice
        end % for ind_chunk
    end % for ind_dw
end % for ind_avg
clear  sense_7d  recon_m  gmap_m  M  b_v  a_v  m_v  S_m  q  r  I_mask_nav_m
pack

% imageviewxd(recon_6d)



%% END





