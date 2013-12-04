
%[s_get_DWI_RECONparams] get general parameters, GENparams.
%
% USAGE:
%   s_get_DWI_RECONparams
%
% OUTPUT:
%   GENparams: A struct with fields as,
%       B0, in Gauss
%       gamma, in rad/(mT*s)
%
%
% Last modified
% 2010.08.06.
% 2010.08.11.
%   DWIparams.BW is directly acquired from ExamCard, removing,
%       DWIparams.AQbase_dur
%       DWIparams.irfe_sq_Tw_epi
%       DWIparams.AQME_dur
%       DWIparams.irfe_sq_epi_shift
% 2010.08.13.
%   Use [f_read_examcard_params.m].
% 2010.08.24.
%   Clear some unused parameters.
% 2010.09.01.
%   Add 32 channel coil (as coilType = 'RX-Intf-1').
%   Calculate `UGN1_ACQ_reduced_acq_time manually (for 3T) based on
%   ugeo0_aspect_ratio()-[mpugeo__g.c]. This is funtionalized as
%   [f_ugeo0_aspect_ratio.m].
% 2010.09.16.
%   '% FOV, Voxel size and Slice thickness.' and 'DWI parameters: Internal
%   parameters' are adjusted based on the slice_orientation.
% 2010.09.24.
%   Remove verbose input argument in [f_get_recon_params.m].
% 2010.10.22.
%   Add slice_orientation to input of [f_get_dw_orientation.m].
% 2010.11.15.
%   dwiparams.trb -> dwiparams.trb*100 for sec/cm^2. This change is made in
%   [f_read_examcard_params_7T.m] and [f_read_examcard_params_3T.m].
% 2010.12.09.
%   Add dwiparams.number_of_directions.
% 2011.01.19.
%   Get DWIparams.DW_GRAD in CRP coordinate system for <default display
%   orientation>.
% 2011.03.03.
%   Add Offc. parameters and ZOOM.
% 2011.03.31.
%   Add SENSE-Torso coil for 3T.
% 2011.06.20.
%   Use [f_get_dw_orientation2.m].
% 2011.09.13.
%   Modify for multiple SSH files.
% 2011.10.19.
%   Use isfield(dwiparams... statement for non-diffusion weighting case.
% 2012.02.03.
%   Modify for 3-D scan mode.
% 2012.03.22.
%   Consider Halfscan and halfscan factor.
% 2012.04.04.
%   Use [f_get_dw_orientation3.m] to consider different coilID and
%   corresponding diffusion gradient orientations.
% 2012.04.16.
%   Add RX-Intf-1_Quad-TR-1 for 7T 16 ch Spine coil.
% 2012.05.17.
%   Add diff_mode for diffusion mode.
% 2012.06.13.
%   Add 'Overcontiguous_slices' and 'Chunks'.
% 2012.06.20.
%   Take care for 3T ExamCard (some parameters do not exist).
%
% Ha-Kyu



%% Define DWIparams and RECONparams

% Clear.
clear  DWIparams



%----- DWI parameters: Raw filename -----
DWIparams.filename = DWIrawfilename;



%----- DWI parameters: Raw filename -----
DWIparams.ZOOM = zoom_s;



%----- DWI parameters: shot mode -----
DWIparams.shot_mode = dwiparams.shot_mode;



%----- DWI parameters: MOTION -----
DWIparams.nNSA = dwiparams.nsa;     % number of NSA



%----- DWI parameters: GEOMETRY -----
% [Coil selection], number of actual receiving coil channels.
if GENparams.B0==70000    
    if strcmpi(GENparams.coilID,'SENSE-Head-7T')
        DWIparams.nCOIL = 16;
    elseif strcmpi(GENparams.coilID,'RX-Intf-1')
        DWIparams.nCOIL = 32;
    elseif strcmpi(GENparams.coilID,'RX-Intf-1_Quad-TR-1')
        DWIparams.nCOIL = 16;
    elseif strcmpi(GENparams.coilID,'Extremity-T/R')
        DWIparams.nCOIL = 1;
    else
        error('s_get_DWI_RECONparams:main','Other coilType is not yet defined')
    end
elseif GENparams.B0==30000    
    if strcmpi(GENparams.coilID,'SENSE-Head-8')
        DWIparams.nCOIL = 8;
    elseif strcmpi(GENparams.coilID,'SENSE-NV-16')
        DWIparams.nCOIL = 16;
    elseif strcmpi(GENparams.coilID,'SENSE-Breast-4')
        DWIparams.nCOIL = 4;
    elseif strcmpi(GENparams.coilID,'SENSE-Torso')
        DWIparams.nCOIL = 6;
    elseif strcmpi(GENparams.coilID,'SENSE-Head-32')
        DWIparams.nCOIL = 32;
    else
        error('s_get_DWI_RECONparams:main','Other coilType is not yet defined')
    end
else
    error('s_get_DWI_RECONparams:main','Unknown GENparams.B0')
end


% FOV, Voxel size and Slice thickness.
DWIparams.FOV = dwiparams.fov; % dwiparams.fov.[AP,RL,FH]
if strcmpi(dwiparams.slice_orientation,'transverse')    
    if strcmpi(dwiparams.fold_over_dir,'AP')
        DWIparams.acq_vox_m = dwiparams.acq_vox.RL; % M_ORI = RL
    elseif strcmpi(dwiparams.fold_over_dir,'RL')
        DWIparams.acq_vox_m = dwiparams.acq_vox.AP; % M_ORI = AP
    else
        error('s_get_DWI_RECONparams:main','Unknown DWI fold-over direction')
    end
    DWIparams.slice_thickness = dwiparams.acq_vox.FH; % S_ORI = FH
elseif strcmpi(dwiparams.slice_orientation,'sagittal')
    if strcmpi(dwiparams.fold_over_dir,'AP')
        DWIparams.acq_vox_m = dwiparams.acq_vox.FH; % M_ORI = FH
    elseif strcmpi(dwiparams.fold_over_dir,'FH')
        DWIparams.acq_vox_m = dwiparams.acq_vox.AP; % M_ORI = AP
    else
        error('s_get_DWI_RECONparams:main','Unknown DWI fold-over direction')
    end
    DWIparams.slice_thickness = dwiparams.acq_vox.RL; % S_ORI = RL
elseif strcmpi(dwiparams.slice_orientation,'coronal')
    if strcmpi(dwiparams.fold_over_dir,'RL')
        DWIparams.acq_vox_m = dwiparams.acq_vox.FH; % M_ORI = FH
    elseif strcmpi(dwiparams.fold_over_dir,'FH')
        DWIparams.acq_vox_m = dwiparams.acq_vox.RL; % M_ORI = RL
    else
        error('s_get_DWI_RECONparams:main','Unknown DWI fold-over direction')
    end
    DWIparams.slice_thickness = dwiparams.acq_vox.AP; % S_ORI = AP
else
    error('s_get_DWI_RECONparams:main','Unknown DWI slice orientation')
end


% Others.
DWIparams.ACQ_RESOL = dwiparams.acq_matrix_mp(1);   % [Matrix scan], (M_ORI)
DWIparams.RECON_RESOL = dwiparams.recon_resol;  % [Recon - reconstruction]
if isfield(dwiparams,'sense_factor')
    DWIparams.SENSE_FACTOR = dwiparams.sense_factor;    % [SENSE - P reduction]
else
    DWIparams.SENSE_FACTOR = 1; % 1 or SENSE = 'no'
end
DWIparams.nSLICE = dwiparams.slices;    % [slices] -- must be updated for 3D acquisition
if isfield(dwiparams,'gap')
    DWIparams.slice_gap = dwiparams.gap;    % [slice gap]
else
    DWIparams.slice_gap = 0; % scan mode == 3D case
end
DWIparams.slice_orientation = dwiparams.slice_orientation;  % [slice orientation], (transverse, coronal, sagittal)
DWIparams.fold_over_dir = dwiparams.fold_over_dir;  % [fold-over direction]
DWIparams.fat_shift_dir = dwiparams.fat_shift_dir;  % [fat shift direction]
if isfield(dwiparams,'overcont_slices')
    DWIparams.overcont_slices = dwiparams.overcont_slices; % [Overcontiguous slices]
end
if isfield(dwiparams,'chunks') % there is no parameter from 3T examcard
    DWIparams.chunks = dwiparams.chunks; % [Chunks]
else
    DWIparams.chunks = 1; % make chunks one
end



%----- DWI parameters: OFFC/ANG -----
DWIparams.Offc.AP = dwiparams.offc(1);  % [Stack Offc.], slice offcenter AP
DWIparams.Offc.RL = dwiparams.offc(2);  % [Stack Offc.], slice offcenter RL
DWIparams.Offc.FH = dwiparams.offc(3);  % [Stack Offc.], slice offcenter FH
DWIparams.Ang.AP = dwiparams.ang(1);    % [Stack Ang.], slice angulation AP
DWIparams.Ang.RL = dwiparams.ang(2);    % [Stack Ang.], slice angulation RL
DWIparams.Ang.FH = dwiparams.ang(3);    % [Stack Ang.], slice angulation FH



%----- DWI parameters: CONTRAST -----
DWIparams.isEPI = true;     % [Fast Imaging mode]
DWIparams.scan_mode = dwiparams.scan_mode; %[Scan mode], (MS,M2D,2D,3D)
DWIparams.technique = dwiparams.technique; %[technique], (SE,...)
DWIparams.EPI_FACTOR = dwiparams.epi_factor;    % [EPI factor], (echo1, echo2)
DWIparams.nEC = dwiparams.echoes;   % [Echoes], number of echo
DWIparams.halfscan = dwiparams.halfscan; % [Halfscan] yes, no
if isfield(dwiparams,'halfscan_factor')
    DWIparams.halfscan_factor = dwiparams.halfscan_factor;
end
if isfield(dwiparams,'diff_mode')
    DWIparams.diff_mode = dwiparams.diff_mode;
end
if isfield(dwiparams,'diff_weightings')
    DWIparams.nROW = dwiparams.diff_weightings;     % [nr of b-factors] number of diffusion weighting, 0 and >0
else
    DWIparams.nROW = 1; % non-DW case
end
if isfield(dwiparams,'directional_resolution')
    DWIparams.directional_resolution = dwiparams.directional_resolution;    % [directional resolution]
end
if isfield(dwiparams,'number_of_directions')
    DWIparams.number_of_directions = dwiparams.number_of_directions;
end
if isfield(dwiparams,'gradient_overplus')
    DWIparams.gradient_overplus = dwiparams.gradient_overplus;
end

% Read diffusion encoding gradient orientations in 
if isfield(dwiparams,'directional_resolution')
    dw_grad = f_get_dw_orientation3(GENparams,DWIparams,user_defined_dw_ori);
end


%------------
% OLD METHOD.
% if isfield(dwiparams,'number_of_directions')
%     DWIparams.number_of_directions = dwiparams.number_of_directions;    
% %     dw_grad = f_get_dw_orientation(DWIparams.directional_resolution, ...
% %         DWIparams.Ang, DWIparams.fold_over_dir, DWIparams.fat_shift_dir, ...
% %         DWIparams.slice_orientation, DWIparams.number_of_directions);
%     dw_grad = f_get_dw_orientation(GENparams.B0,DWIparams.directional_resolution, ...
%         DWIparams.Ang, DWIparams.fold_over_dir, DWIparams.fat_shift_dir, ...
%         DWIparams.slice_orientation, GENparams.patient_orientation, ...
%         GENparams.patient_position, DWIparams.number_of_directions); % CRP coordinate system
% else
% %     dw_grad = f_get_dw_orientation(DWIparams.directional_resolution, ...
% %         DWIparams.Ang, DWIparams.fold_over_dir, DWIparams.fat_shift_dir, ...
% %         DWIparams.slice_orientation);
%     dw_grad = f_get_dw_orientation(GENparams.B0,DWIparams.directional_resolution, ...
%         DWIparams.Ang, DWIparams.fold_over_dir, DWIparams.fat_shift_dir, ...
%         DWIparams.slice_orientation, GENparams.patient_orientation, ...
%         GENparams.patient_position); % CRP coordinate system
% end
%------------

if exist('dw_grad','var')
    DWIparams.DW_GRAD = dw_grad;  % CRP coordinate in <default display orientation>
    clear  dw_grad;
    DWIparams.nDW_GRAD = size(DWIparams.DW_GRAD,1);
else
    if strcmpi(DWIparams.diff_mode,'DWI')
        DWIparams.nDW_GRAD = 3;
    end
end
if isfield(dwiparams,'trb')
    DWIparams.trb = dwiparams.trb*100;  % [max b-factor], sec/mm^2 -> sec/cm^2
end



%----- DWI parameters: Internal parameters -----
% Bandwidth
DWIparams.BW = dwiparams.bw; % [image-echo,navigator-echo]

% This is used for [f_get_ENCima.m].
% 7T examcard has `EX_ACQ_reduced_acq_time in GEOMETRY tab. But not in 3T.
% Then calculate it in the way it is calculated in
% [ugeo0_adjust_reduced_acq_time()-[mpugeo__g.c].
if GENparams.B0==70000
    DWIparams.UGN1_ACQ_reduced_acq_time = dwiparams.acq_reduced_acq_time;
elseif GENparams.B0==30000
    if strcmpi(dwiparams.slice_orientation,'transverse')
        if strcmpi(dwiparams.fold_over_dir,'AP')
            fov_m = dwiparams.fov.RL;
            fov_p = dwiparams.fov.AP;
            acq_vox_p = dwiparams.acq_vox.AP;
        elseif strcmpi(dwiparams.fold_over_dir,'RL')
            fov_m = dwiparams.fov.AP;
            fov_p = dwiparams.fov.RL;
            acq_vox_p = dwiparams.acq_vox.RL;
        else
            error('s_get_DWI_RECONparams:main','Unknown DWI fold-over direction')
        end
    elseif strcmpi(dwiparams.slice_orientation,'sagittal')
        if strcmpi(dwiparams.fold_over_dir,'AP')
            fov_m = dwiparams.fov.FH;
            fov_p = dwiparams.fov.AP;
            acq_vox_p = dwiparams.acq_vox.AP;
        elseif strcmpi(dwiparams.fold_over_dir,'FH')
            fov_m = dwiparams.fov.AP;
            fov_p = dwiparams.fov.FH;
            acq_vox_p = dwiparams.acq_vox.FH;
        else
            error('s_get_DWI_RECONparams:main','Unknown DWI fold-over direction')
        end
    elseif strcmpi(dwiparams.slice_orientation,'coronal')
        if strcmpi(dwiparams.fold_over_dir,'FH')
            fov_m = dwiparams.fov.RL;
            fov_p = dwiparams.fov.FH;
            acq_vox_p = dwiparams.acq_vox.FH;
        elseif strcmpi(dwiparams.fold_over_dir,'RL')
            fov_m = dwiparams.fov.FH;
            fov_p = dwiparams.fov.RL;
            acq_vox_p = dwiparams.acq_vox.RL;
        else
            error('s_get_DWI_RECONparams:main','Unknown DWI fold-over direction')
        end
    else
        error('s_get_DWI_RECONparams:main','Unknown slice orientation')
    end
    ovs_factor_p = 1.0; % this seems wrong, but it is 1.0 (maybe temporarily) in the source code
    EX_ACQ_scan_resol = dwiparams.acq_matrix_mp(1); % M_ORI is scan resolution
    EX_GEO_rect_fov_perc = min(fov_m,fov_p) / max(fov_m,fov_p) * 100;
    aspect_ratio = f_ugeo0_aspect_ratio(EX_ACQ_scan_resol,EX_GEO_rect_fov_perc);
    rounded_matrix_p = round(fov_p * ovs_factor_p / acq_vox_p);
    scan_percentage = 100*rounded_matrix_p / round(ovs_factor_p*EX_ACQ_scan_resol*aspect_ratio);
    DWIparams.UGN1_ACQ_reduced_acq_time = scan_percentage;
    clear  fov_p  fov_m  ovs_factor_p  acq_vox_p  EX_ACQ_scan_resol  EX_GEO_rect_fov_perc
    clear  rounded_matrix_p  scan_percentage
    
else
    error('s_get_DWI_RECONparams:main','Unknown GENparams.B0')    
end



%----- DWI parameters: Others 1 -----
% gradient starting polarity
grad_start_sign = f_get_grad_start_sign(fullfile(loadDir_s,DWIparams.filename));
DWIparams.STD_start_sign = grad_start_sign.STD_start_sign;   % gradient polarity of first k-space acquisition for STD data
if isfield(grad_start_sign,'PHC_start_sign')
    DWIparams.PHC_start_sign = grad_start_sign.PHC_start_sign;   % gradient polarity of first k-space acquisition for PHC data
end
if isfield(grad_start_sign,'PHX_start_sign')
    DWIparams.PHX_start_sign = grad_start_sign.PHX_start_sign;   % gradient polarity of first k-space acquisition for PHX data
end
clear  grad_start_sign


% Philips SSH DWI scan PAR filename.
if ~isempty(SSHparrecfilename)
    DWIparams.SSH.filename = SSHparrecfilename;
    for ind_ssh_data = 1:length(SSHparrecfilename)
        [data,info] = f_read_parrec(fullfile(loadParDir_s,SSHparrecfilename{ind_ssh_data}));
        %DWIparams.SSH.filename = SSHparrecfilename;
        DWIparams.SSH.trb{ind_ssh_data} = info.imgdef.diffusion_b_factor.uniq(2:end) * 100; % non-DW, in s/cm^2
        DWIparams.SSH.nSLICE(ind_ssh_data) = str2num(info.Max_number_of_slices);
        DWIparams.SSH.slice_thickness(ind_ssh_data) = info.imgdef.slice_thickness.uniq;
        DWIparams.SSH.nDW_GRAD(ind_ssh_data) = str2num(info.Max_number_of_gradient_orients)-1;
        DWIparams.SSH.REC_voxel_MPS{ind_ssh_data} = [info.imgdef.pixel_spacing.uniq,DWIparams.SSH.slice_thickness];
    end
else
    DWIparams.SSH.filename = {};
    DWIparams.SSH.trb = {}; % in s/cm^2
    DWIparams.SSH.nSLICE = [];
    DWIparams.SSH.slice_thickness = [];
    DWIparams.SSH.nDW_GRAD = [];
    DWIparams.SSH.REC_voxel_MPS = {};
end
clear  data  info



%----- STACKima -----
if isfield(dwiparams,'rfov')
    % Use rfov/100 as STACK`ima.aspect ratio.
    STACKima.aspect_ratio = double(uint8(dwiparams.rfov))/100;
else
    % Use aspect_ratio previously calculated as STACK`ima:aspect_ratio.
    STACKima.aspect_ratio = aspect_ratio;
end
clear aspect_ratio



%----- RECON parameters -----
DWIparams.loadDir_s = loadDir_s;
RECONparams = f_get_recon_params(REFparams,DWIparams,STACKima);

% 3-D scan mode related parameters
if strcmpi(DWIparams.scan_mode,'3D')
    scan_mode_3D_related_params = f_get_scan_mode_3D_params( ...
        fullfile(loadDir_s,DWIparams.filename));
    
    % Get scan mode = 3D related parameters in RECONparams above.
    max_min_enc_num = scan_mode_3D_related_params.kz_range;
    ovs_fac = scan_mode_3D_related_params.kz_oversample_factor;
    recon_res = scan_mode_3D_related_params.Z_resolution;
    sense_fac = scan_mode_3D_related_params.Z_direction_SENSE_factor;
    
    % Calculate kz directional parameters using the params got above.
    ovs_res = diff(max_min_enc_num)+1;
    scan_res = ovs_res/ovs_fac;
    profiles = ovs_res;
    interp_fac = recon_res/(sense_fac*ovs_res)*ovs_fac;
    ft_length = uint16(profiles*interp_fac);
    
    % Update RECONparams.ENCima.
    RECONparams.ENCima.oversample_resolutions(3) = ovs_res;
    RECONparams.ENCima.oversample_factors(3) = ovs_fac;
    RECONparams.ENCima.scan_resolutions(3) = scan_res;
    RECONparams.ENCima.profiles(3) = profiles;
    RECONparams.ENCima.interp_factors(3) = interp_fac;
    RECONparams.ENCima.ft_lengths(3) = ft_length;
    RECONparams.ENCima.recon_resolutions(3) = recon_res;
    RECONparams.ENCima.sense_factors(3) = sense_fac;
    RECONparams.ENCima.max_encoding_numbers(3) = max_min_enc_num(2);
    RECONparams.ENCima.min_encoding_numbers(3) = max_min_enc_num(1);
    clear  max_min_enc_num ovs_fac recon_res sense_fac ovs_res scan_res
    clear  profiles interp_fac ft_length
end

% Verify.
% ovs_res = RECONparams.ENCima.oversample_resolutions;
% ovs_fac = RECONparams.ENCima.oversample_factors;
% scan_res = RECONparams.ENCima.scan_resolutions;
% profiles = RECONparams.ENCima.profiles;
% interp_fac = RECONparams.ENCima.interp_factors;
% ft_length = RECONparams.ENCima.ft_lengths;
% recon_res = RECONparams.ENCima.recon_resolutions;
% sense_fac = RECONparams.ENCima.sense_factors;
% max_enc_num = RECONparams.ENCima.max_encoding_numbers;
% min_enc_num = RECONparams.ENCima.min_encoding_numbers;
% 
% disp([ovs_res./ovs_fac.*sense_fac.*interp_fac; recon_res])




%----- DWI parameters: Others 2 -----
DWIparams.acq_vox_p = RECONparams.ACQREC.ACQ_voxel_MPS(2);
DWIparams.nSHOT = RECONparams.ACQREC.number_of_shots;



%----- Take some parameters -----
Ns = double(DWIparams.nSHOT);
R = double(DWIparams.SENSE_FACTOR);


% Update actual number of slice for 3D acquisition.
if strcmpi(DWIparams.scan_mode,'3D')
    DWIparams.kz = RECONparams.ENCima.oversample_resolutions(3);
end


% TEMP ----- START
% This is temporary solution for BW(2) because in scan_mode=3D, no NAV BW
% is shown in INFO page.
DWIparams.BW(2) = DWIparams.BW(1) / DWIparams.nSHOT;
% TEMP ----- END


% Clear.
% clear  dwiparams  DWIrawfilename  SSHparrecfilename  data  info



