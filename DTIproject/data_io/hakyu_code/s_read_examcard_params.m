
%[s_read_examcard_params] reads ExamCard parameters.
%
%
% NOTE:
%   group index is as following,
%       1: INFO PAGE
%       2: GEOMETRY
%       3: CONTRAST
%       4: MOTION
%       5: DYN/ANG, angio
%       6: POST/PROC
%       7: OFFC/ANG, angulation
%
%
% Last modified
% 2010.08.13.



%% Simulate input
rawfilename_REF = 'raw_007'; % exported .LIST filename
examcardindex_REF = 7; % index of the file in the ExamCard
rawfilename_DWI = 'raw_008';
examcardindex_DWI = 8;
examcardindex_FMAP = 6; % FMAP is PAR/REC, then only examcardindex is required -> option


%% Read examcard
examcard = f_read_examcard(xmlfile);


%% Set group index - ExamCard tab
info_idx = 1;
geometry_idx = 2;
contrast_idx = 3;
motion_idx = 4;
dyn_ang_idx = 5;
post_proc_idx = 6;
offc_ang_idx = 7;


%% Read REF scan

% Set parameter index - ExamCard parameters
fov_ind = 9;
acq_vox_ind = 5;


% Scan name
scan_name = examcard.protocol(examcardindex_REF).ATTRIBUTE.id;

% FOV
fov = [0,0,0]; %[AP,RL,FH]
fov(2) = examcard.protocol(examcardindex_REF).group(geometry_idx). ...
    parameter(fov_ind).ATTRIBUTE.value; % RL
fov(1) = examcard.protocol(examcardindex_REF).group(geometry_idx). ...
    parameter(fov_ind).subparameter(1).ATTRIBUTE.value; % AP
fov(3) = examcard.protocol(examcardindex_REF).group(geometry_idx). ...
    parameter(fov_ind).subparameter(2).ATTRIBUTE.value; % FH

% acq_vox_(m,p,s)
acq_vox_mps = [0,0,0];
acq_vox_s = examcard.protocol(examcardindex_REF).group(info_idx). ...
    parameter(acq_vox_ind).ATTRIBUTE.value; % [m,p,s]
ind = regexp(acq_vox_s,'/');
acq_vox_mps(1) = str2num(acq_vox_s(1:ind(1)-1));
acq_vox_mps(2) = str2num(acq_vox_s(ind(1)+1:ind(2)-1));
acq_vox_mps(3) = str2num(acq_vox_s(ind(2)+1:end));

% Output for REFparams.
clear  refparams
refparams.scan_name = scan_name;
refparams.fov = fov;
refparams.acq_vox_mps = acq_vox_mps;

% Report.
fprintf('\n')
fprintf('  REF params\n');
disp(refparams)

% Clear.
clear  fov  acq_vox_mps  acq_vox_s  scan_name


%% Read DWI params

% Set parameter index - ExamCard parameters
fov_ind = 8;
acq_vox_ind = 6;
acq_matrix_ind = 5;
recon_resol_ind = 13;
sense_ind = 6;
stacks_ind = 15;
epi_factor_ind = 6;
echoes_ind = 7;
diffusion_mode_ind = 20;
stack_offc_ind = 2;
wfs_bw_echo1_ind = 12;
wfs_bw_echo2_ind = 13;
scan_percentage_ind = 14;
rfov_ind = 11;

% Scan name
scan_name = examcard.protocol(examcardindex_DWI).ATTRIBUTE.id;

% FOV
fov = [0,0,0]; %[AP,RL,FH]
fov(1) = examcard.protocol(examcardindex_DWI).group(geometry_idx). ...
    parameter(fov_ind).ATTRIBUTE.value; % AP
fov(2) = examcard.protocol(examcardindex_DWI).group(geometry_idx). ...
    parameter(fov_ind).subparameter(1).ATTRIBUTE.value; % RL
fov(3) = examcard.protocol(examcardindex_DWI).group(geometry_idx). ...
    parameter(fov_ind).subparameter(2).ATTRIBUTE.value; % FH

% acq matrix (m,p)
acq_matrix_mp = [0,0];
acq_matrix_s = examcard.protocol(examcardindex_DWI).group(info_idx). ...
    parameter(acq_matrix_ind).ATTRIBUTE.value; % [M,P]
ind = regexp(acq_matrix_s,'x');
acq_matrix_mp(1) = str2num(acq_matrix_s(1:ind-1));
acq_matrix_mp(2) = str2num(acq_matrix_s(ind+1:end));
clear  acq_matrix_s  ind

% acq_vox_(m,p,s)
acq_vox_mps = [0,0,0];
acq_vox_s = examcard.protocol(examcardindex_DWI).group(info_idx). ...
    parameter(acq_vox_ind).ATTRIBUTE.value; % [m,p,s]
ind = regexp(acq_vox_s,'/');
acq_vox_mps(1) = str2num(acq_vox_s(1:ind(1)-1));
acq_vox_mps(2) = str2num(acq_vox_s(ind(1)+1:ind(2)-1));
acq_vox_mps(3) = str2num(acq_vox_s(ind(2)+1:end));
clear  acq_vox_s  ind

% reconstruction resolution
recon_resol = examcard.protocol(examcardindex_DWI).group(geometry_idx). ...
    parameter(recon_resol_ind).subparameter.ATTRIBUTE.value; % [m,p,s]

% sense factor
sense_factor = examcard.protocol(examcardindex_DWI).group(geometry_idx). ...
    parameter(sense_ind).subparameter(1).ATTRIBUTE.value; % P reduction

% slices
slices = examcard.protocol(examcardindex_DWI).group(geometry_idx). ...
    parameter(stacks_ind).subparameter(2).ATTRIBUTE.value; % number of slices

% slice gap
gap = examcard.protocol(examcardindex_DWI).group(geometry_idx). ...
    parameter(stacks_ind).subparameter(5).ATTRIBUTE.value; % slice gap

% slice orientation
slice_orientation = examcard.protocol(examcardindex_DWI).group(geometry_idx). ...
    parameter(stacks_ind).subparameter(6).ATTRIBUTE.value; % slice orientation

% fold over and fat shift
fold_over = examcard.protocol(examcardindex_DWI).group(geometry_idx). ...
    parameter(stacks_ind).subparameter(7).ATTRIBUTE.value; % fold-over direction
fat_shift = examcard.protocol(examcardindex_DWI).group(geometry_idx). ...
    parameter(stacks_ind).subparameter(8).ATTRIBUTE.value; % fat shift direction

% EPI factors
epi_factor_v = [0,0];
epi_factor_v(1) = examcard.protocol(examcardindex_DWI).group(contrast_idx). ...
    parameter(epi_factor_ind).ATTRIBUTE.value; % image-echo only
epi_factor_v(2) = 17; % navigator-echo epi factor is fixed to 17

% echoes
echoes = examcard.protocol(examcardindex_DWI).group(contrast_idx). ...
    parameter(echoes_ind).ATTRIBUTE.value; % number of echoes

% nr of b factors
diff_weightings = examcard.protocol(examcardindex_DWI).group(contrast_idx). ...
    parameter(diffusion_mode_ind).subparameter(4).ATTRIBUTE.value; % number of diffusion weightings

% directional resolution
directional_resolution = examcard.protocol(examcardindex_DWI).group(contrast_idx). ...
    parameter(diffusion_mode_ind).subparameter(3).ATTRIBUTE.value; % number of diffusion gradient orientations

% max b-factor
trb = examcard.protocol(examcardindex_DWI).group(contrast_idx). ...
    parameter(diffusion_mode_ind).subparameter(6).ATTRIBUTE.value; % diffusion b-value, s/mm^2
trb = trb*100; % s/cm^2

% angulation
ang_v = [0,0,0]; % [AP,RL,FH]
ap = examcard.protocol(examcardindex_DWI).group(offc_ang_idx). ...
    parameter(stack_offc_ind).subparameter(3).ATTRIBUTE.value; % AP
rl = examcard.protocol(examcardindex_DWI).group(offc_ang_idx). ...
    parameter(stack_offc_ind).subparameter(4).ATTRIBUTE.value; % RL
fh = examcard.protocol(examcardindex_DWI).group(offc_ang_idx). ...
    parameter(stack_offc_ind).subparameter(5).ATTRIBUTE.value; % FH
ang_v = [ap,rl,fh];
clear  ap  rl  fh

% bandwidth
bw_v = [0,0]; % image-echo and navigator-echo
bw1_s = examcard.protocol(examcardindex_DWI).group(info_idx). ...
    parameter(wfs_bw_echo1_ind).ATTRIBUTE.value;
ind = regexp(bw1_s,'/');
bw1 = str2num(bw1_s(ind+1:end));
bw2_s = examcard.protocol(examcardindex_DWI).group(info_idx). ...
    parameter(wfs_bw_echo2_ind).ATTRIBUTE.value;
ind = regexp(bw2_s,'/');
bw2 = str2num(bw2_s(ind+1:end));
bw_v = [bw1,bw2];
clear  bw1_s  bw2_s  ind

% UGN1_ACQ_reduced_acq_time
%* There are 'IF_scan_percentage' in INFO and 'EX_ACQ_reduced_acq_time' in
%* GEOMETRY both termed as 'Scan percentage(%s)'. Here,
%* EX_ACQ_reduced_acq_time will be used for the UGN1 substitute. But it must
%* be verified in 3T.
acq_reduced_acq_time = examcard.protocol(examcardindex_DWI).group(geometry_idx). ...
    parameter(scan_percentage_ind).ATTRIBUTE.value;

% RFOV
rfov = examcard.protocol(examcardindex_DWI).group(geometry_idx). ...
    parameter(rfov_ind).ATTRIBUTE.value;

% Output for DWIparams.
clear  dwiparams
dwiparams.scan_name = scan_name;
dwiparams.fov = fov;
dwiparams.acq_matrix_mp = acq_matrix_mp;
dwiparams.acq_vox_mps = acq_vox_mps;
dwiparams.recon_resol = recon_resol;
dwiparams.sense_factor = sense_factor;
dwiparams.slice_orientation = slice_orientation;
dwiparams.fold_over = fold_over;
dwiparams.fat_shift = fat_shift;
dwiparams.epi_factor = epi_factor_v;
dwiparams.echoes = echoes;
dwiparams.diff_weightings = diff_weightings;
dwiparams.directional_resolution = directional_resolution;
dwiparams.trb = trb;
dwiparams.ang = ang_v;
dwiparams.bw = bw_v;
dwiparams.acq_reduced_acq_time = acq_reduced_acq_time;
dwiparams.rfov = rfov;

% Report.
fprintf('\n')
fprintf('  DWI params\n');
disp(dwiparams)

% Clear.
clear  scan_name  fov  acq_matrix_mp  acq-vox_mps  recon_resol  sense_factor
clear  slice_orientation  fold_over  fat_shift  epi_factor_v  echoes  
clear  diff_weightings  directional_resolution  trb  ang_v  bw_v
clear  acq_reduced_acq_time  rfov










