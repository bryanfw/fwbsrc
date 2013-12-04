
%[s_get_DATAFLAGparams] get data processing parameters, DATAFLAGparams.
%
% USAGE:
%   s_get_DATAFLAGparams
%
% OUTPUT:
%   DATAFLAGparams: A struct with fields as,
%       epiCorrFitMethod, Fitting method for EPI echo phase difference
%       readREF, Read REF data
%       readDWI, Read DWI data
%       genSMAP, Generate SMAP from coil reference data
%       genFMAP, Postprocess FMAP for fieldmap correction between NAV and IMG
%
%
% Last modified
% 2010.08.06.
% 2010.08.13.
%   Use values defined in [s_batch.m].
% 2010.08.18.
%   Removed fields in DATAFLAGparams,
%       saveREF, showREF, saveDWI, showDWI 
% 2010.09.16.
%   Clear unused parameters later.
% 2010.09.20.
%   Use epi_corr_method for epiCorrFitMethod field.
% 2012.05.17.
%   Use epi_corr_sel_total_var to select the best EPI ghost correction
%   results by using total variation in [f_epi_phase_corr_phc_v2.m].
%
%   Use match_obj_pos_smap_img to match object position between SMAP and
%   IMG data. This difference has been found and verified in
%   [s_batch__scan20120515_r4140__Jeong_999999.m] when 16 ch spine coil was
%   used at 7T.
%
%   Add another term at the end of e.g., raw_010_R2E9S8 such as
%   raw_010_R2E9S8_recon_again etc.
% 2012.07.12.
%   Add coreg_nav_to_img for use in [s_fmap_corr_NAV.m].
%
% Ha-Kyu



%% Define DATAFLAGparams

DATAFLAGparams.epiCorrFitMethod = epi_corr_fit_method; % default: 1, [1,2,4] can be used in [f_epi_phase_corr_phc.m]
DATAFLAGparams.readREF = read_ref;
DATAFLAGparams.readDWI = read_dwi;
DATAFLAGparams.genSMAP = gen_smap;   % generate coil sensitivity map
DATAFLAGparams.genFMAP = gen_fmap;
DATAFLAGparams.epiCorrSelTotalVar = epi_corr_sel_total_var;
DATAFLAGparams.matchObjPosSMAP_IMG = match_obj_pos_smap_img;
DATAFLAGparams.addTermAtDataDir = add_term_at_data_dir; % additional term at the end of dataDir_s
DATAFLAGparams.coreg_nav_to_img = coreg_nav_to_img; % coreg NAV to IMG in [s_fmap_corr_NAV.m]

% clear  read_ref  read_dwi  gen_smap  gen_fmap
















