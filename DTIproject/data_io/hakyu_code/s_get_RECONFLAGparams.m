
%[s_get_RECONFLAGparams] get reconstruction related parameters,
%RECONFLAGparams. This parameter determines how the final IMG recon will be
%performed. For example,
%   phaseCorr: Do recon with NAV phase correction
%   navFmapCorr: Do fmap correction of NAV and use it
%   navDeform:  Do direction deformation of NAV and use it
%In all cases, windowed NAV recon image will be used.
%   
%
% USAGE:
%   s_get_RECONFLAGparams
%
% OUTPUT:
%   RECONFLAGparams: A struct with fields as,
%       phaseCorr, Phase corrected final IMG recon or not
%       navFmapCorr, Use fieldmap correction on nav recon image
%       navDeform, Use deformation of nav to img
%       $navFilt, Use low-pass filter on recon nav image
%       psiMethod, Method for incorporating noise correltation (PSI),
%                  [0,1,2], 1 and 2 generates same result
%       saveReconData, Save recon data
%       saveReconFig, Save recon figure
%
%
% Last modified
% 2010.08.09.
% 2010.08.24.
%   Clear some unused variables.
% 2010.09.16.
%   Clear unused parameters later.
% 2010.09.20.
%   Use psi_method and calc_diff for psiMethod and calcDiffProp fields.
% 2012.01.18.
%   Add 'sub_encoding_matrix_flag' and 'sub_encoding_matrix_siz_idx' and
%   'precalc_encoding_matrix_flag' for using sub-encoding matrix for faster
%   inversion of each image column.
% 2012.01.20.
%   Add 'run_pct' for running Parallel Computing Toolbox (PCT).
% 2012.02.13.
%   Add flag for resize to recon resolution or not.
%
% Ha-Kyu



%% Define RECONFLAGparams

% Define it.
RECONFLAGparams.phaseCorr = phase_corr;         % recon with phase correction
RECONFLAGparams.navFmapCorr = nav_fmap_corr;    % use fieldmap correction on nav recon image
RECONFLAGparams.navDeform = nav_deform;         % use deformation of nav to img
RECONFLAGparams.psiMethod = psi_method;         % method for incorporating noise correltation (PSI), [0,1,2]
RECONFLAGparams.calcDiffProp = calc_diff;       % calculate diffusion properties
RECONFLAGparams.coregDWI = [];                  % flag if DW images are co-registered to non-DW image, see [s_coreg_DWI.m]
RECONFLAGparams.precalc_encoding_matrix_flag = ...
    precalc_encoding_matrix_flag;               % precalculation of encoding matrix
RECONFLAGparams.sub_encoding_matrix_flag = ...
    sub_encoding_matrix_flag;                   % [true or false] for using sub-encoding matrix
RECONFLAGparams.sub_encoding_matrix_siz_idx = ...
    sub_encoding_matrix_siz_idx;                % [0,1,2,3]=[Ny,R,Ns,R*Ns] index of column size of sub-encoding matrix, Ny==R*Ns*EPI_factor
RECONFLAGparams.run_pct = run_pct;              % Run PCT
RECONFLAGparams.resize_DWI = flag_resize_DWI;   % resize to recon resolution or not

% Check error.
if RECONFLAGparams.navFmapCorr==1 && RECONFLAGparams.navDeform==1
    text1_s = 'RECONFLAGparams.navFmapCorr and RECONFLAGparams.navDeform are ';
    text2_s = 'mutually exclusive at current version EXCEPT for both zero.';
    error('s_get_RECONFLAGparams:main','%s%s',text1_s,text2_s)
end

if isempty(find(RECONFLAGparams.sub_encoding_matrix_siz_idx==[0,1,2,3]))
    text1_s = 'RECONFLAGparams.sub_encoding_matrix_siz_idx is ';
    text2_s = 'not one of [0,1,2,3].';
    error('s_get_RECONFLAGparams:main','%s%s',text1_s,text2_s)
end

% Uncomment this because sub-encoding matrix can be used with
% precalculation of encoding matrix.
% if RECONFLAGparams.precalc_encoding_matrix_flag==true && ...
%     RECONFLAGparams.sub_encoding_matrix_flag~=false
%     text1_s = 'RECONFLAGparams.sub_encoding_matrix_flag should be ';
%     text2_s = 'zero when RECONFLAGparams.precalc_encoding_matrix_flag is ';
%     text3_s = 'TRUE.';
%     error('s_get_RECONFLAGparams:main','%s%s%s',text1_s,text2_s,text3_s)
% end

% Clear.
% clear  phase_corr  nav_fmap_corr  nav_deform












