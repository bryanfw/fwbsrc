
%[s_check_smapSizeInfo] show the adjustment of image size of SMAP before
%reconstruction.
%
% USAGE:
%   s_check_smapSizeInfo
%
%
% Last modified
% 2010.08.06.
% 2010.09.01.
% 2010.09.02.
%   Due to non-integral sense factor, this generates error. Take care of
%   this as 'theoretical' and 'actual' and 'RECONparams' case. The
%   'theoretical' case use accurate sense factor which may not be integer.
%   The 'actual' case use integral sense factor as used for calculating
%   RECONparams.a_N_recon_m and RECONparams.a_N_recon_p.
% 2010.09.03.
%   Changed smap_matrix_v from double to uint16 class.
% 2010.09.16.
%   Adjust fov and acq_vox based on slice_orientation.
% 2012.05.22.
%   Take care of volume coil case.
%
% Ha-Kyu



%% Check SMAP and ALIASED matrix size

if ~isempty(REFparams.filename) % volume coil case
    
% REF.
if strcmpi(REFparams.slice_orientation,'transverse')
    if strcmpi(REFparams.fold_over_dir,'AP')
        ref_fov_m = REFparams.FOV.RL;
        ref_fov_p = REFparams.FOV.AP;
        ref_acq_vox_m = REFparams.acq_vox.RL;
        ref_acq_vox_p = REFparams.acq_vox.AP;
    elseif strcmpi(REFparams.fold_over_dir,'RL')
        ref_fov_m = REFparams.FOV.AP;
        ref_fov_p = REFparams.FOV.RL;
        ref_acq_vox_m = REFparams.acq_vox.AP;
        ref_acq_vox_p = REFparams.acq_vox.RL;
    else
        error('s_check_smapSizeInfo:main','Unknown REF fold-over direction')
    end
elseif strcmpi(REFparams.slice_orientation,'sagittal')
    if strcmpi(REFparams.fold_over_dir,'AP')
        ref_fov_m = REFparams.FOV.FH;
        ref_fov_p = REFparams.FOV.AP;
        ref_acq_vox_m = REFparams.acq_vox.FH;
        ref_acq_vox_p = REFparams.acq_vox.AP;
    elseif strcmpi(REFparams.fold_over_dir,'FH')
        ref_fov_m = REFparams.FOV.AP;
        ref_fov_p = REFparams.FOV.FH;
        ref_acq_vox_m = REFparams.acq_vox.AP;
        ref_acq_vox_p = REFparams.acq_vox.FH;
    else
        error('s_check_smapSizeInfo:main','Unknown REF fold-over direction')
    end
elseif strcmpi(REFparams.slice_orientation,'coronal')
    if strcmpi(REFparams.fold_over_dir,'RL')
        ref_fov_m = REFparams.FOV.FH;
        ref_fov_p = REFparams.FOV.RL;
        ref_acq_vox_m = REFparams.acq_vox.FH;
        ref_acq_vox_p = REFparams.acq_vox.RL;
    elseif strcmpi(REFparams.fold_over_dir,'FH')
        ref_fov_m = REFparams.FOV.RL;
        ref_fov_p = REFparams.FOV.FH;
        ref_acq_vox_m = REFparams.acq_vox.RL;
        ref_acq_vox_p = REFparams.acq_vox.FH;
    else
        error('s_check_smapSizeInfo:main','Unknown REF fold-over direction')
    end
else
    error('s_check_smapSizeInfo:main','Unknown REF slice orientation')
end
fprintf('SMAP    FOV                     : [M,P] = [%d,%d]\n', ref_fov_m,ref_fov_p)
fprintf('SMAP    acq_vox                 : [M,P] = [%.3f,%.3f]\n', ref_acq_vox_m,ref_acq_vox_p)
fprintf('SMAP    ovs_fac                 : [M,P] = [%.3f,%.3f]\n', ...
    REFparams.ovs_fac_m,REFparams.ovs_fac_p)
smap_matrix_v = uint16([ref_fov_m,ref_fov_p] ./ [ref_acq_vox_m,ref_acq_vox_p] .* ...
    [REFparams.ovs_fac_m,REFparams.ovs_fac_p]);
fprintf('SMAP    matrix                  : [M,P] = [%.3f,%.3f]\n', ...
    smap_matrix_v)
disp(' ')


% DWI.
if strcmpi(DWIparams.slice_orientation,'transverse')
    if strcmpi(DWIparams.fold_over_dir,'AP')
        dwi_fov_m = DWIparams.FOV.RL;
        dwi_fov_p = DWIparams.FOV.AP;
    elseif strcmpi(DWIparams.fold_over_dir,'RL')
        dwi_fov_m = DWIparams.FOV.AP;
        dwi_fov_p = DWIparams.FOV.RL;
    else
        error('s_check_smapSizeInfo:main','Unknown DWI fold-over direction')
    end
elseif strcmpi(DWIparams.slice_orientation,'sagittal')
    if strcmpi(DWIparams.fold_over_dir,'AP')
        dwi_fov_m = DWIparams.FOV.FH;
        dwi_fov_p = DWIparams.FOV.AP;
    elseif strcmpi(DWIparams.fold_over_dir,'FH')
        dwi_fov_m = DWIparams.FOV.AP;
        dwi_fov_p = DWIparams.FOV.FH;
    else
        error('s_check_smapSizeInfo:main','Unknown DWI fold-over direction')
    end
elseif strcmpi(DWIparams.slice_orientation,'coronal')
    if strcmpi(DWIparams.fold_over_dir,'FH')
        dwi_fov_m = DWIparams.FOV.RL;
        dwi_fov_p = DWIparams.FOV.FH;
    elseif strcmpi(DWIparams.fold_over_dir,'RL')
        dwi_fov_m = DWIparams.FOV.FH;
        dwi_fov_p = DWIparams.FOV.RL;
    else
        error('s_check_smapSizeInfo:main','Unknown DWI fold-over direction')
    end
else
    error('s_check_smapSizeInfo:main','Unknown DWI slice orientation')
end
fprintf('ALIASED FOV                     : [M,P] = [%d,%d]\n', dwi_fov_m,dwi_fov_p)
fprintf('ALIASED acq_vox                 : [M,P] = [%.3f,%.3f]\n', ...
    RECONparams.ACQREC.ACQ_voxel_MPS(1:2))
fprintf('ALIASED ovs_fac                 : [M,P] = [%.3f,%.3f]\n', ...
    RECONparams.ENCima.oversample_factors(1:2))


% Get theoretical and actual recon matrix for DWI.

% Here 'theoretical' means the sense_factor is non-integer. But actual
% initial recon matrix size is, E*Ns*R, here R must be integer. Then there
% is difference between 'theoretical' and 'actual' recon matrix for DWI.

alias_matrix_theoretical_v = [dwi_fov_m,dwi_fov_p] ./ ...
    RECONparams.ACQREC.ACQ_voxel_MPS(1:2) .* ...
    RECONparams.ENCima.oversample_factors(1:2);
% Above is the same as,
% (RECONparams.ENCima.epi_factors) * ...
%     (RECONparams.ACQREC.number_of_shots) * ...
%     (RECONparams.ENCima.sense_factors(2))
% Here, RECONparams.ENCima.sense_factor(2) can be non-integer.
        
alias_matrix_p_actual = (RECONparams.ENCima.epi_factors) * ...
    (RECONparams.ACQREC.number_of_shots) * ...
    (DWIparams.SENSE_FACTOR);
alias_matrix_actual_v = [alias_matrix_theoretical_v(1),alias_matrix_p_actual];
% Here, DWIparams.SENSE_FACTOR is integer as shown in ExamCard. This is the
% same as, [RECONparams.a_N_recon_m, RECONparams.a_N_recon_p].

alias_matrix_v = [RECONparams.a_N_recon_m, RECONparams.a_N_recon_p];

fprintf('ALIASED matrix (theoretical, after recon)    : [M,P] = [%.3f,%.3f]\n', ...
    alias_matrix_theoretical_v)
fprintf('ALIASED matrix (actual, after recon)         : [M,P] = [%.3f,%.3f]\n', ...
    alias_matrix_actual_v)
fprintf('ALIASED matrix (RECONparams, after recon)    : [M,P] = [%.3f,%.3f]\n', ...
    alias_matrix_v)
disp(' ')


% SMAP.
fprintf('SMAP    resize                  : [M,P] = [%.3f,%.3f]\n', ...
    RECONparams.s_N_res_m, RECONparams.s_N_res_p)
fprintf('SMAP    cut                     : [M,P] = [%.3f,%.3f]\n', ...
    RECONparams.s_cut_vox_m, RECONparams.s_cut_vox_p)
fprintf('SMAP    zeropad                 : [M,P] = [%.3f,%.3f]\n', ...
    RECONparams.s_zeropad_vox_m, RECONparams.s_zeropad_vox_p)
smap_final_matrix_v = [RECONparams.s_N_final_m, RECONparams.s_N_final_p];
fprintf('SMAP    after cut and zeropad   : [M,P] = [%.3f,%.3f]\n', ...
    smap_final_matrix_v)


% Check final matrix size.
if ~all(alias_matrix_v == smap_final_matrix_v)
    error('s_check_smapSizeInfo:main','ALIASED and SMAP matrix size is different for reconstruction')
end

end





