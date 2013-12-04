
function [refparams,dwiparams] = f_read_examcard_params_7T(...
    xmlfile, REFrawfilename, DWIrawfilename, ...
    REFexamcardindex, DWIexamcardindex, coilType)
%[f_read_examcard_params] reads ExamCard parameters for 7T
%
% USAGE:
%   [refparams,dwiparams] = f_read_examcard_params_7T(...
%     xmlfile, REFrawfilename, DWIrawfilename, ...
%     REFexamcardindex, DWIexamcardindex, coilType)
%
% INPUT:
%   xmlfile:            XML format ExamCard
%   REFrawfilename:     Reference raw data filename (list/data data)
%   DWIrawfilename:     DWI raw data filename (list/data data)
%   REFexamcardindex:   Index of the REF scan in the ExamCard, e.g., if the
%                       REF scan was the 6th scan, then it is 6
%   REFexamcardindex:   Index of the DWI scan in the ExamCard
%   coilType:           'Coil selection' parameter in GEOMETRY tab in ExamCard
%
% OUTPUT:
%   refparams:  REF scan parameters acquired from ExamCard
%   dwiparams:  DWI scan parameters acquired from ExamCard
%
% NOTE:
%   dwiparams is only for custom-made msh-EPI SENSE scan, not for any ssh
%   or msh scans using Philips protocol which can be exported in PAR/REC
%   format.
%
%   REF scan must be exported as raw .LIST format to process SMAP.
%
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
% 2010.09.01.
%   Use GEOMETRY tab for reading acq_vox_mps. When INFO PAGE tab is used,
%   the acquisition voxel size is rounded.
% 2010.09.02.
%   This is saved for 7T scan.
% 2010.09.03.
%   Modify for 32 channel Rx coil.
% 2010.09.16.
%   Modify trb.
% 2010.09.18.
%   Add slice_orientation.
% 2010.10.08.
%   Add Fold-over suppression.
% 2010.10.21.
%   Check CONTRAST - Fast Imaging mode - shot mode which is [single-shot, multishot]
%   then read single-shot DWI data using examcard.
% 2010.11.11.
%   Take care coronal slice orientation.
% 2010.11.15.
%   trb from s/cm^2 to s/mm^2.
% 2010.11.17.
%   Major modifications are made for reading ExamCard.
% 2010.12.09.
%   Add diffusion 'nr of directions'.
% 2011.03.02.
%   Add 'direction' in 'Diffusion mode' under Contrast tab for DWI case.
% 2011.03.03.
%   Add Offc. parameters.
% 2011.06.16.
%   Add 'gradient overplus'.
% 2011.07.28.
%   Add 'Matrix scan' for 'reconstruction' not under 'Recon voxel size
%   (mm)'. This is for 7T R3.2.1.0.
% 2011.10.31.
%   Add 'if ~isempty(REFrawfilename)' and 'if ~isempty(DWIrawfilename)' for
%   reading raw data.
%   Add REFrawfilename and DWIrawfilename as input.
% 2012.02.01.
%   Modified for 3-D scan mode.
% 2012.03.22.
%   Consider Halfscan and halfscan factor.
% 2012.05.17.
%   Add diff_mode for diffusion_mode.
% 2012.06.13.
%   Add 'Overcontiguous_slices' and 'Chunks'. 
%
% Ha-Kyu



%% Simulate input
examcardindex_REF = REFexamcardindex; % index of the file in the ExamCard
examcardindex_DWI = DWIexamcardindex;



%% Read examcard
examcard = f_read_examcard(xmlfile);



%% Set group index - ExamCard tab
info_idx = 1;
info_name_s = 'InfoPage';

geometry_idx = 2;
geometry_name_s = 'geometry';

contrast_idx = 3;
contrast_name_s = 'contrast';

motion_idx = 4;
motion_name_s = 'motion';

dyn_ang_idx = 5;
dyn_ang_name_s = 'dyn/ang';

post_proc_idx = 6;
post_proc_name_s = 'postproc';

offc_ang_idx = 7;
offc_ang_name_s = 'offc/ang';





%% Read REF scan ExamCard
if ~isempty(REFrawfilename)
    
    scan_name = examcard.protocol(examcardindex_REF).ATTRIBUTE.id;
    ngroup = length(examcard.protocol(examcardindex_REF).group);
    for ind_group = 1:ngroup
        tab_name_s = examcard.protocol(examcardindex_REF).group(ind_group).ATTRIBUTE.name;
        
        
        
        % Process each of the ExamCard tabs.
        %-------------------- InfoPage -----------------------------------------------------------
        if strcmpi(tab_name_s,info_name_s)
            % do nothing
            
            
            
            %-------------------- GEOMETRY -----------------------------------------------------------
        elseif strcmpi(tab_name_s,geometry_name_s)
            nparam = length(examcard.protocol(examcardindex_REF).group(ind_group).parameter);
            
            for ind_param = 1:nparam
                param_name_s = examcard.protocol(examcardindex_REF). ...
                    group(ind_group).parameter(ind_param).ATTRIBUTE.name;
                
                switch param_name_s
                    case 'Stacks'
                        struct_temp = examcard.protocol(examcardindex_REF). ...
                            group(ind_group).parameter(ind_param);
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'slices'
                                        slices = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'slice thickness (mm)'
                                        slice_thickness = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'slice orientation'
                                        slice_orientation = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'fold-over direction'
                                        fold_over_dir = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'fat shift direction'
                                        fat_shift_dir = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''Stacks'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                        
                    case 'Fold-over suppression'
                        struct_temp = examcard.protocol(examcardindex_REF). ...
                            group(ind_group).parameter(ind_param);
                        fold_over_suppression = struct_temp.ATTRIBUTE.value;
                        clear  struct_temp
                        
                    case 'FOV          FH (mm)' % coronal RL-F or RL-H | sagital AP-F or AP-H
                        struct_temp = examcard.protocol(examcardindex_REF). ...
                            group(ind_group).parameter(ind_param);
                        fov.FH = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'RL (mm)' % phase
                                        fov.RL = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'stack    AP (mm)' % slice
                                        fov.AP = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'AP (mm)' % phase
                                        fov.AP = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'stack    RL (mm)' % slice
                                        fov.RL = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''FOV'' parameter must have ''subparameter''')
                        end
                        
                    case 'FOV          RL (mm)' % coronal FH-R or FH-L | transverse AP-R or AP-L
                        struct_temp = examcard.protocol(examcardindex_REF). ...
                            group(ind_group).parameter(ind_param);
                        fov.RL = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'FH (mm)'
                                        fov.FH = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'stack    AP (mm)'
                                        fov.AP = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'AP (mm)'
                                        fov.AP = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'stack    FH (mm)'
                                        fov.FH = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''FOV'' parameter must have ''subparameter''')
                        end
                        
                    case 'FOV          AP (mm)' % transverse RL-A or RL-P | sagittal FH-A or FH-P
                        struct_temp = examcard.protocol(examcardindex_REF). ...
                            group(ind_group).parameter(ind_param);
                        fov.AP = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'RL (mm)'
                                        fov.RL = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'stack    FH (mm)'
                                        fov.FH = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'FH (mm)'
                                        fov.FH = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'stack    RL (mm)'
                                        fov.RL = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''FOV'' parameter must have ''subparameter''')
                        end
                        
                    case 'Voxel size   FH (mm)' % coronal RL-F or RL-H | sagittal AP-F or AP-H
                        struct_temp = examcard.protocol(examcardindex_REF). ...
                            group(ind_group).parameter(ind_param);
                        acq_vox.FH = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'RL (mm)'
                                        acq_vox.RL = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'AP (mm)'
                                        acq_vox.AP = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''Voxel size'' parameter must have ''subparameter''')
                        end
                        
                    case 'Voxel size   RL (mm)' % coronal FH-R or FH-L | transverse AP-R or AP-L
                        struct_temp = examcard.protocol(examcardindex_REF). ...
                            group(ind_group).parameter(ind_param);
                        acq_vox.RL = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'FH (mm)'
                                        acq_vox.FH = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'AP (mm)'
                                        acq_vox.AP = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''Voxel size'' parameter must have ''subparameter''')
                        end
                        
                    case 'Voxel size   AP (mm)' % transverse RL-A or RL-P | sagittal FH-A or FA-P
                        struct_temp = examcard.protocol(examcardindex_REF). ...
                            group(ind_group).parameter(ind_param);
                        acq_vox.AP = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'RL (mm)'
                                        acq_vox.RL = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'FH (mm)'
                                        acq_vox.FH = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''Voxel size'' parameter must have ''subparameter''')
                        end
                        
                end % switch param_name_s
            end % for ind_param
            
        elseif strcmpi(tab_name_s,contrast_name_s)
            % do nothing
            
        elseif strcmpi(tab_name_s,motion_name_s)
            % do nothing
            
        elseif strcmpi(tab_name_s,dyn_ang_name_s)
            % do nothing
            
        elseif strcmpi(tab_name_s,post_proc_name_s)
            % do nothing
            
            
            
            %-------------------- OFFC/ANG -----------------------------------------------------------
        elseif strcmpi(tab_name_s,offc_ang_name_s)
            nparam = length(examcard.protocol(examcardindex_REF).group(ind_group).parameter);
            
            for ind_param = 1:nparam
                param_name_s = examcard.protocol(examcardindex_REF). ...
                    group(ind_group).parameter(ind_param).ATTRIBUTE.name;
                
                switch param_name_s
                    case 'Stack Offc. AP (P=+mm)'
                        struct_temp = examcard.protocol(examcardindex_REF). ...
                            group(ind_group).parameter(ind_param);
                        ap_offc = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'RL (L=+mm)'
                                        rl_offc = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'FH (H=+mm)'
                                        fh_offc = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'Ang.  AP (deg)'
                                        ap_ang = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'RL (deg)'
                                        rl_ang = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'FH (deg)'
                                        fh_ang = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                            offc_v = [ap_offc,rl_offc,fh_offc];
                            ang_v = [ap_ang,rl_ang,fh_ang];
                            clear  ap_offc  rl_offc  fh_offc  ap_ang  rl_ang  fh_ang
                        end
                        clear  struct_temp
                end % switch param_name_s
            end % for ind_param
            
        end % if strcmpi(tab_name_s,
    end % for ind_group
    
    
    % Output for REFparams.
    clear  refparams
    refparams.scan_name = scan_name;
    refparams.fov = fov;
    refparams.acq_vox = acq_vox;
    refparams.slices = slices;
    refparams.slice_thickness = slice_thickness;
    refparams.slice_orientation = slice_orientation;
    refparams.fold_over_dir = fold_over_dir;
    refparams.fat_shift_dir = fat_shift_dir;
    refparams.fold_over_suppression = fold_over_suppression;
    refparams.offc = offc_v;
    refparams.ang = ang_v;
    
    
    % Report.
    fprintf('\n')
    fprintf('  refparams\n');
    disp(refparams)
    
    
    % Clear.
    clear  fov  acq_vox  acq_vox_s  scan_name  fold_over_dir  slices  slice_thickness
    clear  slice_orientation  fold_over_suppression
        
end





%% Read DWI scan ExamCard
if ~isempty(DWIrawfilename)
    
    scan_name = examcard.protocol(examcardindex_DWI).ATTRIBUTE.id;
    % shot_mode = examcard.protocol(examcardindex_DWI).group(contrast_idx). ...
    %     parameter(fast_imaging_mode_ind). ...
    %     subparameter(fast_imaging_mode_shot_mode_sub_ind).ATTRIBUTE.value; % [single-shot, multishot]
    
    ngroup = length(examcard.protocol(examcardindex_DWI).group);
    for ind_group = 1:ngroup
        tab_name_s = examcard.protocol(examcardindex_DWI).group(ind_group).ATTRIBUTE.name;
        
        
        
        % Process each of the ExamCard tabs.
        %-------------------- CONTRAST -----------------------------------------------------------
        if strcmpi(tab_name_s,contrast_name_s)
            nparam = length(examcard.protocol(examcardindex_DWI).group(ind_group).parameter);
            
            for ind_param = 1:nparam
                param_name_s = examcard.protocol(examcardindex_DWI). ...
                    group(ind_group).parameter(ind_param).ATTRIBUTE.name;
                
                switch param_name_s
                    case 'Scan mode'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        
                        scan_mode = struct_temp.ATTRIBUTE.value; % scan mode: MS,2D,3D,M2D
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'technique'
                                        technique = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''Scan mode'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp                        
                        
                    case 'Fast Imaging mode'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'shot mode'
                                        shot_mode = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''Fast Imaging mode'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                        
                    case 'EPI factor' % multishot case
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        
                        epi_factor_v = [0,0];
                        epi_factor_v(1) = struct_temp.ATTRIBUTE.value; % image-echo
                        epi_factor_v(2) = 17; % navigator-echo
                        clear  struct_temp
                        
                    case 'Echoes'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        echoes = struct_temp.ATTRIBUTE.value;
                        clear  struct_temp
                        
                    case 'Halfscan'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        halfscan = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'factor'
                                        halfscan_factor = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        end                                        
                        clear  struct_temp
                        
                    case 'Diffusion mode'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        diff_mode = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'direction' % DWI case
                                        val = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                        a = regexp(val,'[MPS]');
                                        directional_resolution = val(a);
                                        clear  val  a
                                    case 'directional resolution' % DTI case
                                        directional_resolution = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'nr of b-factors'
                                        diff_weightings = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'max b-factor'
                                        maxtrb = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value; % s/mm^2
                                        trb = linspace(0,maxtrb,diff_weightings); % s/mm^2
                                    case 'nr of directions'
                                        number_of_directions = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value; % number
                                    case 'gradient overplus'
                                        gradient_overplus = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value; % yes, no
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''Diffusion mode'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                end % switch param_name_s
            end % for ind_param
            
            
            
            %-------------------- GEOMETRY -----------------------------------------------------------
        elseif strcmpi(tab_name_s,geometry_name_s)
            nparam = length(examcard.protocol(examcardindex_DWI).group(ind_group).parameter);
            
            for ind_param = 1:nparam
                param_name_s = examcard.protocol(examcardindex_DWI). ...
                    group(ind_group).parameter(ind_param).ATTRIBUTE.name;
                
                switch param_name_s
                    case 'SENSE'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'P reduction (RL)' % transverse or coronal
                                        sense_factor = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'P reduction (AP)' % transverse or sagittal
                                        sense_factor = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'P reduction (FH)' % coronal or sagittal
                                        sense_factor = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''SENSE'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                        
                    case 'RFOV (%)'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        rfov = struct_temp.ATTRIBUTE.value;
                        clear  struct_temp
                        
                    case 'Recon voxel size (mm)'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'reconstruction' % 7T R2.5.3.4
                                        recon_resol = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''SENSE'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                        
                    case 'Matrix scan'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'reconstruction' % 7T R3.2.1.0
                                        recon_resol = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_3T:main','''Matrix scan'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                        
                    case 'Scan percentage (%)'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        
                        % UGN1_ACQ_reduced_acq_time
                        %* There are 'IF_scan_percentage' in INFO and
                        %* 'EX_ACQ_reduced_acq_time' in * GEOMETRY both termed as
                        %* 'Scan percentage(%s)'. Here, * EX_ACQ_reduced_acq_time
                        %* will be used for the UGN1 substitute. But it must * be
                        %* verified in 3T.
                        acq_reduced_acq_time = struct_temp.ATTRIBUTE.value;
                        
                    case 'Overcontiguous slices'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        overcont_slices = struct_temp.ATTRIBUTE.value;
                        clear  struct_temp
                        
                    case 'Stacks'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'slices'
                                        slices = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'gap (mm)'
                                        slice_gap = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'slice orientation'
                                        slice_orientation = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'fold-over direction'
                                        fold_over_dir = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'fat shift direction'
                                        fat_shift_dir = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''Stacks'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                        
                    case 'Chunks'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        chunks = struct_temp.ATTRIBUTE.value;
                        clear  struct_temp
                        
                    case 'FOV          FH (mm)' % coronal RL-F or RL-H | sagittal AP-F or AP-H
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        fov.FH = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'RL (mm)'
                                        fov.RL = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'AP (mm)'
                                        fov.AP = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''FOV'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                        
                    case 'FOV          RL (mm)' % coronal FH-R or FH-L | transverse AP-R or AP-L
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        fov.RL = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'FH (mm)'
                                        fov.FH = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'AP (mm)'
                                        fov.AP = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''FOV'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                        
                    case 'FOV          AP (mm)' % transverse RL-A or RL-P | sagittal FH-A or FH-P
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        fov.AP = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'RL (mm)'
                                        fov.RL = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'FH (mm)'
                                        fov.FH = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''FOV'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                        
                    case 'Voxel size   FH (mm)' % coronal RL-F or RL-H | sagittal AP-F or AP-H
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        acq_vox.FH = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'RL (mm)'
                                        acq_vox.RL = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'AP (mm)'
                                        acq_vox.AP = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''Voxel size'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                        
                    case 'Voxel size   RL (mm)' % coronal FH-R or FH-L | transverse AP-R or AP-L
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        acq_vox.RL = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'FH (mm)'
                                        acq_vox.FH = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'AP (mm)'
                                        acq_vox.AP = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''Voxel size'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                        
                    case 'Voxel size   AP (mm)' % transverse RL-A or RL-P | sagittal FH-A or FH-P
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        acq_vox.AP = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'RL (mm)'
                                        acq_vox.RL = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'FH (mm)'
                                        acq_vox.FH = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                        else
                            error('f_read_examcard_params_7T:main','''Voxel size'' parameter must have ''subparameter''')
                        end
                        clear  struct_temp
                        
                    case 'Slice thickness (mm)'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        slice_thickness = struct_temp.ATTRIBUTE.value;
                        if isfield(acq_vox,'FH') && isfield(acq_vox,'RL') % coronal RL-F or RL-H
                            acq_vox.AP = slice_thickness;
                        elseif isfield(acq_vox,'FH') && isfield(acq_vox,'AP') % sagittal AP-F or AP-H
                            acq_vox.RL = slice_thickness;
                        elseif isfield(acq_vox,'RL') && isfield(acq_vox,'FH') % coronal FH-R or FH-L
                            acq_vox.AP = slice_thickness;
                        elseif isfield(acq_vox,'RL') && isfield(acq_vox,'AP') % transverse AP-R or AP-L
                            acq_vox.FH = slice_thickness;
                        elseif isfield(acq_vox,'AP') && isfield(acq_vox,'RL') % transverse RL-A or RL-R
                            acq_vox.FH = slice_thickness;
                        elseif isfield(acq_vox,'AP') && isfield(acq_vox,'FH') % sagittal FH-A or FH-P
                            acq_vox.RL = slice_thickness;
                        else
                            error('f_read_examcard_params_7T:main','''Slice thickness (mm)'' parameter has unknown direction of acq_vox')
                        end
                        clear  slice_thickness  struct_temp
                        
                end % switch param_name_s
            end % for ind_param
            
            
            
            %-------------------- InfoPage -----------------------------------------------------------
        elseif strcmpi(tab_name_s,info_name_s)
            nparam = length(examcard.protocol(examcardindex_DWI).group(ind_group).parameter);
            
            for ind_param = 1:nparam
                param_name_s = examcard.protocol(examcardindex_DWI). ...
                    group(ind_group).parameter(ind_param).ATTRIBUTE.name;
                
                switch param_name_s
                    case 'ACQ matrix M x P'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        acq_matrix_mp = [0,0];
                        acq_matrix_s = struct_temp.ATTRIBUTE.value; % [M,P]
                        ind = regexp(acq_matrix_s,'x');
                        acq_matrix_mp(1) = str2num(acq_matrix_s(1:ind-1));
                        acq_matrix_mp(2) = str2num(acq_matrix_s(ind+1:end));
                        clear  acq_matrix_s  ind  struct_temp
                        
                    case 'WFS (pix) / BW (Hz) [echo1]'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        
                        bw1_s = struct_temp.ATTRIBUTE.value;
                        ind = regexp(bw1_s,'/');
                        bw1 = str2num(bw1_s(ind+1:end));
                        clear  bw1_s  ind  struct_temp
                        
                    case 'WFS (pix) / BW (Hz) [echo2]'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        
                        bw2_s = struct_temp.ATTRIBUTE.value;
                        ind = regexp(bw2_s,'/');
                        bw2 = str2num(bw2_s(ind+1:end));
                        
                        % Combine BW for echo1 and echo2.
                        bw_v = [bw1,bw2]; % image-echo and navigator-echo
                        clear  bw2_s  bw1  bw2  ind  struct_temp
                        
                    case 'WFS (pix) / BW (Hz)' % single-shot case
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        
                        bw1_s = struct_temp.ATTRIBUTE.value;
                        ind = regexp(bw1_s,'/');
                        bw1 = str2num(bw1_s(ind+1:end));
                        bw_v = [bw1];
                        clear  bw1_s  bw1  ind  struct_temp
                        
                    case 'EPI factor' % single-shot case
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        
                        epi_factor_v = struct_temp.ATTRIBUTE.value; % image-echo
                        clear  struct_temp
                end % switch param_name_s
            end % for ind_param
            
            
            
            %-------------------- MOTION -----------------------------------------------------------
        elseif strcmpi(tab_name_s,motion_name_s)
            nparam = length(examcard.protocol(examcardindex_DWI).group(ind_group).parameter);
            
            for ind_param = 1:nparam
                param_name_s = examcard.protocol(examcardindex_DWI). ...
                    group(ind_group).parameter(ind_param).ATTRIBUTE.name;
                
                switch param_name_s
                    case 'NSA'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        nsa = struct_temp.ATTRIBUTE.value;
                end
            end
            
            
            
            %-------------------- OFFC/ANG -----------------------------------------------------------
        elseif strcmpi(tab_name_s,offc_ang_name_s)
            nparam = length(examcard.protocol(examcardindex_DWI).group(ind_group).parameter);
            
            for ind_param = 1:nparam
                param_name_s = examcard.protocol(examcardindex_DWI). ...
                    group(ind_group).parameter(ind_param).ATTRIBUTE.name;
                
                switch param_name_s
                    case 'Stack Offc. AP (P=+mm)'
                        struct_temp = examcard.protocol(examcardindex_DWI). ...
                            group(ind_group).parameter(ind_param);
                        ap_offc = struct_temp.ATTRIBUTE.value;
                        
                        if isfield(struct_temp,'subparameter')
                            nsubparam = length(struct_temp.subparameter);
                            for ind_subparam = 1:nsubparam
                                subparam_name_s = struct_temp.subparameter(ind_subparam).ATTRIBUTE.name;
                                switch subparam_name_s
                                    case 'RL (L=+mm)'
                                        rl_offc = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'FH (H=+mm)'
                                        fh_offc = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'Ang.  AP (deg)'
                                        ap_ang = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'RL (deg)'
                                        rl_ang = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                    case 'FH (deg)'
                                        fh_ang = struct_temp.subparameter(ind_subparam).ATTRIBUTE.value;
                                end
                            end
                            offc_v = [ap_offc,rl_offc,fh_offc];
                            ang_v = [ap_ang,rl_ang,fh_ang];
                            clear  ap_offc  rl_offc  fh_offc  ap_ang  rl_ang  fh_ang
                        end
                        clear  struct_temp
                end % switch param_name_s
            end % for ind_param
            
        elseif strcmpi(tab_name_s,dyn_ang_name_s)
            % do nothing
            
        elseif strcmpi(tab_name_s,post_proc_name_s)
            % do nothing
            
        end % if strcmpi(tab_name_s, ...
    end % for ind_group
    
    
    % Output for DWIparams.
    clear  dwiparams
    dwiparams.scan_name = scan_name;
    dwiparams.scan_mode = scan_mode;
    dwiparams.technique = technique;
    dwiparams.shot_mode = shot_mode;
    if exist('diff_mode','var')
        dwiparams.diff_mode = diff_mode;
    end
    dwiparams.acq_matrix_mp = acq_matrix_mp;
    dwiparams.acq_vox = acq_vox;
    dwiparams.bw = bw_v;
    if exist('sense_factor','var')
        dwiparams.sense_factor = sense_factor; % SENSE can be 'no'
    end
    dwiparams.fov = fov;
    dwiparams.rfov = rfov;
    dwiparams.recon_resol = recon_resol;
    dwiparams.acq_reduced_acq_time = acq_reduced_acq_time;
    dwiparams.slices = slices;
    if exist('slice_gap','var')
        dwiparams.gap = slice_gap;
    end
    dwiparams.slice_orientation = slice_orientation;
    dwiparams.fold_over_dir = fold_over_dir;
    dwiparams.fat_shift_dir = fat_shift_dir;
    if exist('overcont_slices','var')
        dwiparams.overcont_slices = overcont_slices;
    end
    if exist('chunks','var')
        dwiparams.chunks = chunks;
    end
    dwiparams.epi_factor = epi_factor_v;
    dwiparams.echoes = echoes;
    if exist('directional_resolution','var')
        dwiparams.directional_resolution = directional_resolution;
    end
    dwiparams.halfscan = halfscan;
    if exist('halfscan_factor','var')
        dwiparams.halfscan_factor = halfscan_factor;
    end
    dwiparams.diff_weightings = diff_weightings;
    dwiparams.trb = trb;
    if exist('number_of_directions','var')
        dwiparams.number_of_directions = number_of_directions;
    end
    if exist('gradient_overplus','var')
        dwiparams.gradient_overplus = gradient_overplus;
    end
    dwiparams.nsa = nsa;
    dwiparams.offc = offc_v;
    dwiparams.ang = ang_v;
    
    % Report.
    fprintf('\n')
    fprintf('  dwiparams\n');
    disp(dwiparams)
    
    % Clear.
    clear  scan_name  fov  acq_matrix_mp  acq_vox  recon_resol  sense_factor
    clear  slice_orientation  fold_over_dir  fat_shift_dir  epi_factor_v  echoes
    clear  diff_weightings  directional_resolution  trb  offc_v  ang_v  bw_v
    clear  acq_reduced_acq_time  rfov  nsa  shot_mode  number_of_directions
    clear  gradient_overplus
    
end





%% Output
if isempty(REFrawfilename)
    refparams = [];
end
if isempty(DWIrawfilename)
    dwiparams = [];
end





function examcard = f_read_examcard(xmlfile)
%[f_read_examcard] reads examcard and returns examcard as a Matlab struct.
%
% USAGE:
%   examcard = f_read_examcard(xmlfile)
%
% INPUT:
%   xmlfile:    ExamCard in XML format
%
% OUTPUT:
%   examcard:   A struct array including all ExamCard parameters for all
%               scans in the ExamCard
%
% NOTE:
%   This function uses [xml_read.m] by Tuszynski.
%
%
% Last modified
% 2010.08.12.



%% Read ExamCard in xml format
fprintf('Reading ExamCard\n')
tic
tree = xml_read(xmlfile);
t = toc;
fprintf('Examcard is read in %.3f sec\n',t)


%% Take teh whole ExamCard

% Declare examcard as a struct array.
clear examcard
examcard = tree.folder.examcard;









