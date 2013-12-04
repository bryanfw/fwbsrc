
%[s_read_DWI6] read, save and show DWI data in LIST/DATA format.
%
%
% Last modified
% 2011.04.09.
%   This script is generated and modified for using [f_read_listdata2.m].
% 2011.04.17.
%   This script is generated and modified for using [f_read_listdata3.m].
% 2011.04.18.
%   This script is generated and modified for using [f_read_listdata4.m].
% 2011.06.08.
%   Use 'if filesize_GB < filesize_ref_GB' for cutting zero k-space regions
%   when [f_read_listdata.m] is used.
% 2011.10.28.
%   filesize_ref_GB was set to 0.5GB.
% 2011.12.09.
%   [s_read_DWI5.m] is generated from [s_read_DWI4.m]. This script uses
%   [f_read_listdata5.m].
% 2012.02.01.
%   Use [f_read_listdata5.m] for filesize_GB < filesize_ref_GB also.
%   Take care of scan mode = 3D, then use kz than loc.
% 2012.02.03.
%   [f_read_listdata5.m] is modified for saving *_attr_m.mat files.
% 2012.03.22.
%   Single-shot with NSA=1 for ind_avg==2 for non-DW data is taken care.
% 2012.05.16.
%   Add total variation measure to determine which method work best for EPI
%   ghost correction using [f_epi_phase_corr_phc_v2.m].
% 2012.06.11.
%   Take care of 'scan_mode' = '3D'. Use [f_epi_phase_corr_phc_v3.m] for
%   this.
%
% Ha-Kyu



%% Read DWI data

fprintf('Read DWI data using [%s.m].\n',mfilename)

if DATAFLAGparams.readDWI==1
    %% Read data
    
    % Use flag_cell__read_data for using or skipping 'Read data' cell.
    tic
    flag_cell__read_data = 1;
    
    if flag_cell__read_data == 1
        
        % Get DWI raw data filesize and choose {selective reading, whole
        % reading} based on the filesize.
        cd(loadDir_s)
        dirinfo = dir(loadDir_s);
        index = strcmpi({dirinfo.name},[DWIparams.filename '.data']);
        filesize = dirinfo(index).bytes;
        filesize_GB = filesize/(10^9); % in gigabyte
        filesize_ref_GB = 0.5; % reference file size in gigabyte        
        
        
        %if ~strcmpi(GENparams.coilID,'RX-Intf-1') % not a 32 channel coil
        if filesize_GB < filesize_ref_GB % read whole raw data at once            
            
            % Read DATA data using [f_read_listdata5.m].
            cd(loadDir_s)
            [noi_data,phc_data,std_data,nav_data,frc_data] = ...
                f_read_listdata5(DWIparams.filename,saveDataDir_s);
            
            % PHC data.
            [nx,ny,nz,nsl,ndyn,ncard,nec,nmix,ncoil,navg,nsign,ngrad] = size(phc_data);
            
            % STD data.
            [kx,ky,kz,loc,dyn,card,echo,mix,aver,chan,extr1,extr2] = size(std_data);
            
            % Extract data.
            % Don't permute to save memory.
            data_nd = squeeze( std_data );
            clear  std_data
            
            % Reshape data to keep singleton dimension.
            if strcmpi(DWIparams.scan_mode,'3D')
                data_nd = reshape(data_nd,[kx,ky,kz,echo,aver,chan,extr1,extr2]);
                [Nx,Ny,Nz,Nec,Navg,Nc,Ndw,Nori] = size(data_nd);
                fprintf('    Size of data_nd for [nKx,nKy,nKz,nECHO,nAVG,nCOIL,nROW,nDW_GRAD]\n')
            else
                data_nd = reshape(data_nd,[kx,ky,loc,echo,aver,chan,extr1,extr2]);
                [Nx,Ny,Nsl,Nec,Navg,Nc,Ndw,Nori] = size(data_nd);
                fprintf('    Size of data_nd for [nKx,nKy,nSLICE,nECHO,nAVG,nCOIL,nROW,nDW_GRAD]\n')
            end
            fprintf('       ')
            disp(size(data_nd))
            
            
            % Read and save data_nd for each echo, dw and ori.
            fprintf('\n')
            cd(saveDataDir_s)
            for ind_ec = 1:Nec
                for ind_dw = 1:Ndw % For ind_dw=2, all nonzero b-value data are saved
                    if ind_dw==1
                        ori_v = 1;
                        dw_idx = 1;
                    else
                        ori_v = 1:Nori;
                        dw_idx = 2:Ndw;
                    end
                    for ind_ori = ori_v
                        if ind_dw==1
                            ori_idx = 0;
                        else
                            ori_idx = ind_ori;
                        end
                        
                        % Read each echo and dw orientation data.
                        % Keep this as 8-D.
                        eval(sprintf('data_nd__ec%.2d_ori%.2d = data_nd(:,:,:,ind_ec,:,:,dw_idx,ind_ori);', ...
                            ind_ec,ori_idx));
                        
                        % Save each echo and slice data.
                        eval(sprintf('save  data_nd__ec%.2d_ori%.2d  data_nd__ec%.2d_ori%.2d', ...
                            ind_ec,ori_idx,ind_ec,ori_idx))
                        
                        % Clear the data.
                        eval(sprintf('clear  data_nd__ec%.2d_ori%.2d',ind_ec,ori_idx))
                        
                        % Report.
                        fprintf('    data_nd__ec%.2d_ori%.2d is saved.\n', ind_ec,ori_idx)
                    end
                end
            end
            save  phc_data  phc_data
            clear  phc_data
            fprintf('    phc_data is saved and cleared.\n')
            clear  data_nd
            pack
            fprintf('    data_nd is cleared and memory packed.\n')
            
            
        else % GENparams.coilID = 'RX-Intf-1' case --------------------------------------
            
            %**********
            % Use [f_read_listdata2.m] here.
            % This function reads the DATA data as a whole.
            %
            % Use [f_read_listdata3.m] here.
            % This function reads the DATA data as a part.
            %
            % Use [f_read_listdata4.m] here.
            % This function reads the DATA data as a part.
            %
            % Use [f_read_listdata5.m] here.
            % This function takes care of reading NAV data that is
            % interleaved in .LIST file. This function is the same as
            % [f_read_listdata4_v2.m].
            %**********
            
            % Change directory to where all the raw data and (TYP)_attr_m files are saved.
            cd(saveDataDir_s) % for saving (TYP_attr_m files separate from other raw data
            
            % LIST/DATA filename.
            listdataname_s = fullfile(loadDir_s,DWIparams.filename);
            
            % Read raw data.
            clear  loadopts
            for ind_ec = 1:DWIparams.nEC
                loadopts.echo = ind_ec-1; % index of loadopts starts from zero
                dw_v = 0:1; % diffusion weighting
                
                for ind_dw = dw_v % ind_dw=1 can include multiple non-zero b-values
                    if ind_dw==0
                        ori_v = 0;
                    else
                        ori_v = 1:DWIparams.nDW_GRAD;
                    end
                    for ind_ori = ori_v
                        if ind_dw==0
                            loadopts.extr1 = 0; % starts from zero
                            loadopts.extr2 = 0;
                        else
                            loadopts.extr1 = 1:DWIparams.nROW-1;
                            loadopts.extr2 = ind_ori-1; % 0 ~ DWIparams.nDW_GRAD-1
                        end
                        
                        t0=clock;
                        fprintf('-----------------------------------------\n')
                        fprintf('    Reading ECHO[%d], DW[%d], ORI[%d]\n', ...
                            loadopts.echo,loadopts.extr1,loadopts.extr2)
                        disp('      loadopts is')
                        disp(loadopts)
                        
                        cd(saveDataDir_s)
                        
                        % Read whole PHC data.
                        if (ind_ec==1) && (ind_ori==0)
                            selected_TYP = {'phc'};
                            [phc_data] = f_read_listdata5(listdataname_s,saveDataDir_s,selected_TYP);
                            clear  noi_data
                            
                            % Change order of 'sign' in PHC data.
                            % If [f_read_listdata5.m] is used, sign(-1) = sing index 1,
                            % sign(+1) = sign index 2 in the 11th dimension. With this
                            % change of order, remained code may not be modified.
                            phc_data_tmp = phc_data*0;
                            phc_data_tmp(:,:,:,:,:,:,:,:,:,:,1,:) = phc_data(:,:,:,:,:,:,:,:,:,:,2,:);
                            phc_data_tmp(:,:,:,:,:,:,:,:,:,:,2,:) = phc_data(:,:,:,:,:,:,:,:,:,:,1,:);
                            phc_data = phc_data_tmp;
                            clear  phc_data_tmp
                            fprintf('    Order of ''sign'' attr. in phc_data is changed\n')
                            
                            % Save PHC data.
                            [nx,ny,nz,nsl,ndyn,ncard,nec,nmix,ncoil,navg,nsign,ngrad] = size(phc_data);
                            cd(saveDataDir_s)
                            save  phc_data  phc_data
                            clear  phc_data
                            fprintf('    phc_data is saved and cleared.\n')
                        end
                        
                        % Read selected STD data.
                        selected_TYP = {'std'};
                        [sense_data] = f_read_listdata5(listdataname_s,saveDataDir_s,selected_TYP,loadopts);
                        
                        t1 = etime(clock,t0);
                        fprintf('    DW[%d], ORI[%d] is read in %.3f sec\n', ...
                            ind_dw,ind_ori,t1)
                        
                        % STD data.
                        [kx,ky,kz,loc,dyn,card,echo,mix,aver,chan,extr1,extr2] = size(sense_data);
                        
                        % Extract data.
                        % Don't permute to save memory.
                        data_nd = squeeze( sense_data );
                        clear  sense_data
                        
                        % Reshape data to keep singleton dimension. 
                        % 8-D data (MS) or 9-D data (3D).
                        if strcmpi(DWIparams.scan_mode,'MS')
                            data_nd = reshape(data_nd,[kx,ky,loc,echo,aver,chan,extr1,extr2]);
                            [Nx,Ny,Nsl,Nec,Navg,Nc,Ndw,Nori] = size(data_nd);
                            fprintf('    Size of data_nd for [nKx,nKy,nSLICE,nECHO,nAVG,nCOIL,nROW,nDW_GRAD]\n')
                            fprintf('       ')
                        elseif strcmpi(DWIparams.scan_mode,'3D')
                            data_nd = reshape(data_nd,[kx,ky,kz,loc,echo,aver,chan,extr1,extr2]);
                            [Nx,Ny,Nz,Nloc,Nec,Navg,Nc,Ndw,Nori] = size(data_nd);
                            fprintf('    Size of data_nd for [nKx,nKy,nKz,loc,nECHO,nAVG,nCOIL,nROW,nDW_GRAD]\n')
                            fprintf('       ')
                        else
                            error('s_read_DWI6:main','Unknown DWIparams.scan_mode')
                        end                        
                        disp(size(data_nd))
                        
                        % Keep each echo and dw orientation data.
                        eval(sprintf('data_nd__ec%.2d_ori%.2d = data_nd;', ...
                            ind_ec,ind_ori));
                        
                        % Save each echo and slice data.
                        cd(saveDataDir_s)
                        eval(sprintf('save  data_nd__ec%.2d_ori%.2d  data_nd__ec%.2d_ori%.2d', ...
                            ind_ec,ind_ori,ind_ec,ind_ori))
                        fprintf('    data_nd__ec%.2d_ori%.2d is saved.\n',ind_ec,ind_ori)
                        
                        % Clear and pack the data.
                        eval(sprintf('clear  data_nd__ec%.2d_ori%.2d  data_nd',ind_ec,ind_ori))
                        pack
                        fprintf('    data_nd* is cleared and memory packed.\n')
                        fprintf('\n')
                        
                    end % for ind_ori
                end % for ind_dw
                
                % Clear loadopts.
                clear  loadopts
                
            end % for ind_ec
        end % if filesize_GB < filesize_ref_GB 
        
        
        % Cut off zero data space
        % This is due to the difference of EPI_FACTOR for IMG and NAV
        % acquisitions.
        Ns = double(DWIparams.nSHOT);
        ny = DWIparams.EPI_FACTOR(1)*Ns;      % echo 1 ky size
        if strcmpi(DWIparams.shot_mode,'multishot') && Nec==2
            nny = DWIparams.EPI_FACTOR(2)*Ns;     % echo 2 ky size
        else
            nny = ny;
        end
        
        % Below is not necessary when [f_read_listdata5.m] is used.
        if filesize_GB < filesize_ref_GB
            ycutoff = abs(ny-nny)/2;
            if ny > nny
                ec_cutoff = 2; % echo to be cut off
                if ny~=Ny
                    error('s_read_DWI:main','Check EPI_FACTOR and Ns for data size.')
                end
            else
                ec_cutoff = 1;
                if nny~=Ny
                    error('s_read_DWI:main','Check EPI_FACTOR and Ns for data size.')
                end
            end
        end
        
    end % if flag_cell__read_data
    
    
    
    %% Safety zone
    % This is to further process data reading and correction steps
    % when there is an error and this script is stopped after
    % the time consuming raw data reading. What we need after reading the
    % raw data is,
    % 1. size of data_nd
    % 2. phc_data
    % When this cell is to be used, above 'Read data' cell must be commented.
    
    flag_cell__safety_zone = 1;
    
    if flag_cell__safety_zone == 1
        cd(saveDataDir_s)
        load  phc_data
        load  data_nd__ec01_ori00
        
        % Get size of data_nd.
        
        % For single-shot NSA==1 case, Navg measured from ori00 would be
        % 2, but Navg for DW data would be 1. This need to be considered
        % when EPI phase correction is applied below.
        
        if strcmpi(DWIparams.scan_mode,'3D')
            [Nx,Ny,Nz,Nloc,Nec,Navg,Nc,Ndw,Nori] = size(data_nd__ec01_ori00);
            Nsl = Nz;
            Nchunk = Nloc;
        else
            [Nx,Ny,Nsl,Nec,Navg,Nc,Ndw,Nori] = size(data_nd__ec01_ori00);
            Nchunk = 1;
        end
        Nec = DWIparams.nEC;
        Ndw = DWIparams.nROW;
        if ~isfield(DWIparams,'nDW_GRAD')
            DWIparams.nDW_GRAD = extr2;
        end
        Nori = DWIparams.nDW_GRAD; % no DW directions in 3-D DWI case
        
        % Get cutoff and ny, nny.
        Ns = double(DWIparams.nSHOT);
        ny = DWIparams.EPI_FACTOR(1)*Ns;      % echo 1 ky size
        %nny = DWIparams.EPI_FACTOR(2)*Ns;     % echo 2 ky size
        if strcmpi(DWIparams.shot_mode,'multishot') && Nec==2
            nny = DWIparams.EPI_FACTOR(2)*Ns;     % echo 2 ky size
        else
            nny = ny;
        end
        
        % Below is not necessary when [f_read_listdata5.m] is used.
        if filesize_GB < filesize_ref_GB
            ycutoff = abs(ny-nny)/2;
            if ny > nny
                ec_cutoff = 2; % echo to be cut off
                if ny~=Ny
                    error('s_read_DWI:main','Check EPI_FACTOR and Ns for data size.')
                end
            else
                ec_cutoff = 1;
                if nny~=Ny
                    error('s_read_DWI:main','Check EPI_FACTOR and Ns for data size.')
                end
            end
        end        
        clear  data_nd__ec01_ori00
        
    end % if flag_cell__safety_zone
    
    
    
    %% Read and correct k-space data
    
    % Load PHC data.
    cd(saveDataDir_s)
    load  phc_data
    
    % Correct raw data.
    for ind_dw = 0:1
        if ind_dw==0
            ori_v = 1;
        else
            ori_v = 1:Nori;
        end
        
        % Consider different Navg between ind_dw of 0 and 1 when NSA==1.
        if strcmpi(DWIparams.shot_mode,'single-shot') && (DWIparams.nNSA==1)
            if ind_dw==0;
                Navg_v = 1:Navg; % non-DW case, Navg=2 by default
            else
                Navg_v = 1; % DW case, Navg==NSA
            end
        else
            Navg_v = 1:Navg; % multi-shot, single-shot with NSA>1 cases
        end
        
        for ind_ori = ori_v
            % Get orientation index starting from 0 for b=0.
            if ind_dw==0
                ori_idx = 0;
                nDW = 1;
            else
                ori_idx = ind_ori;
                nDW = DWIparams.nROW-1;
            end
            
            for ind_ec = 1:Nec
                
                % Read raw data for each echo and ori: Method 2.
                cd(saveDataDir_s)
                eval(sprintf('load  data_nd__ec%.2d_ori%.2d',ind_ec,ori_idx))
                
                fprintf('--------------------------------------------------------\n')
                fprintf('    [data_nd__ec%.2d_ori%.2d] is read for correction\n', ...
                    ind_ec,ori_idx)
                eval(sprintf('data_nd = data_nd__ec%.2d_ori%.2d;',ind_ec,ori_idx))
                eval(sprintf('clear  data_nd__ec%.2d_ori%.2d',ind_ec,ori_idx))
                
                
                % Generate k_img_corr_ori##_6d and k_nav_corr_ori##_6d.
                % Generate only once for each ind_ec and for all ind_sl,
                % ind_coil, ind_avg.
                if strcmpi(DWIparams.scan_mode,'3D')
                    if ind_ec==1
                        k_data_corr_name_s = sprintf('k_img_corr_ori%.2d_6d',ori_idx);
                        eval(sprintf('%s = zeros(ny,Nx,Nsl,Nchunk,Nc,length(Navg_v),nDW,''single'');', ...
                            k_data_corr_name_s))
                    elseif ind_ec==2
                        k_data_corr_name_s = sprintf('k_nav_corr_ori%.2d_6d',ori_idx);
                        eval(sprintf('%s = zeros(nny,Nx,Nsl,Nchunk,Nc,length(Navg_v),nDW,''single'');', ...
                            k_data_corr_name_s))
                    else
                        error('s_read_DWI:main','Unknown ind_ec,ind_sl,ind_coil,ind_avg.')
                    end
                else                    
                    if ind_ec==1
                        k_data_corr_name_s = sprintf('k_img_corr_ori%.2d_6d',ori_idx);
                        eval(sprintf('%s = zeros(ny,Nx,Nsl,Nc,length(Navg_v),nDW,''single'');', ...
                            k_data_corr_name_s))
                    elseif ind_ec==2
                        k_data_corr_name_s = sprintf('k_nav_corr_ori%.2d_6d',ori_idx);
                        eval(sprintf('%s = zeros(nny,Nx,Nsl,Nc,length(Navg_v),nDW,''single'');', ...
                            k_data_corr_name_s))
                    else
                        error('s_read_DWI:main','Unknown ind_ec,ind_sl,ind_coil,ind_avg.')
                    end
                end
                
                % Keep PHC and STD start sign.
                PHC_start_sign_keep = DWIparams.PHC_start_sign;
                STD_start_sign_keep = DWIparams.STD_start_sign;
                
                for ind_mult_dw = 1:nDW % multiple DW
                    for ind_sl = 1:Nsl % slice or kz
                        for ind_coil = 1:Nc                            
                            for ind_avg = Navg_v
                                for ind_chunk = 1:Nchunk % chunk for 3D multi-chunk acq
                                
                                    % Adjust PHC and STD starting gradient polarity.
                                    %* This has to be done only for SSH case.
%                                     if DWIparams.nSHOT==1
%                                         if ind_avg==1
%                                             DWIparams.PHC_start_sign = DWIparams.PHC_start_sign;
%                                             DWIparams.STD_start_sign = DWIparams.STD_start_sign;
%                                         else
%                                             DWIparams.PHC_start_sign = -DWIparams.PHC_start_sign;
%                                             DWIparams.STD_start_sign = -DWIparams.PHC_start_sign;
%                                         end
%                                     end
                                    
                                    % Read k_m from the raw data.
                                    % data dimension = [kx,ky,sl,ec,avg,coil,dw,ori]
                                    idx_ec = 1;
                                    %idx_dw = 1;
                                    idx_dw = ind_mult_dw;
                                    idx_ori = 1;
                                    
                                    if strcmpi(DWIparams.scan_mode,'3D')
                                        k_m = squeeze(data_nd(:,:,ind_sl,ind_chunk,idx_ec, ...
                                            ind_avg,ind_coil,idx_dw,idx_ori));
                                        k_m = permute(k_m,[2,1]); % [x,y]->[y,x]
                                    else
                                        k_m = squeeze(data_nd(:,:,ind_sl,idx_ec, ...
                                            ind_avg,ind_coil,idx_dw,idx_ori));
                                        k_m = permute(k_m,[2,1]); % [x,y]->[y,x]
                                    end
                                    
                                    % Cutoff unused k-space.
                                    % Below is not necessary when [f_read_listdata5.m] is used.
                                    if filesize_GB < filesize_ref_GB
                                        if ec_cutoff==1 && ind_ec==1
                                            k_m = k_m(ycutoff+1:ycutoff+ny,:);
                                        elseif ec_cutoff==2 && ind_ec==2
                                            k_m = k_m(ycutoff+1:ycutoff+nny,:);
                                        end
                                    end
                                    
                                    % Get each shot ky index.
                                    ky_shot_odd_ind_m = [];
                                    ky_shot_even_ind_m = [];
                                    for ind_shot = 1:Ns
                                        if ind_ec==1
                                            eval(sprintf('ky_shot%d_v = ind_shot:Ns:ny;',ind_shot))
                                        elseif ind_ec==2
                                            eval(sprintf('ky_shot%d_v = ind_shot:Ns:nny;',ind_shot))
                                        end
                                        
                                        eval(sprintf('ky_shot_odd_ind_m = [ky_shot_odd_ind_m; ky_shot%d_v(1:2:end)];', ...
                                            ind_shot))  % odd ky lines in each shot
                                        eval(sprintf('ky_shot_even_ind_m = [ky_shot_even_ind_m; ky_shot%d_v(2:2:end)];', ...
                                            ind_shot))  % even ky lines in each shot
                                    end
                                    
                                    % Reserve corrected k-space data.
                                    k_corr_m = k_m*0;
                                    
                                    
                                    %-------------------------------------------------------------
                                    DATAFLAGparams.epiCorrMethod = 2; % 2010.08.05.
                                    
                                    if DATAFLAGparams.epiCorrMethod==1 % kxshift
                                        
                                        % Correct k-space data.
                                        for ind_shot = 1:Ns
                                            % Set kxshift.
                                            kxshift = DATAFLAGparams.kxshift;
                                            %kxshift = 8.7;
                                            
                                            % Correction based on starting gradient polarity.
                                            if DWIparams.PHC_start_sign == DWIparams.STD_start_sign
                                                
                                                % Generate phase gradient
                                                phase_grad_m = repmat( ...
                                                    exp(1i*2*pi*(-Nx/2:(Nx/2)-1)/Nx*kxshift),...
                                                    length(ky_shot_even_ind_m(ind_shot,:)),1);
                                                
                                                % Keep ODD (+ve) ky lines. 'keeping ky lines'
                                                k_corr_m(ky_shot_odd_ind_m(ind_shot,:),:) = ...
                                                    k_m(ky_shot_odd_ind_m(ind_shot,:),:);
                                                
                                                % Correct EVEN (-ve) ky lines.
                                                h_m_temp = ift1( ...
                                                    k_m(ky_shot_even_ind_m(ind_shot,:),:) * ...
                                                    exp(1i*pi),2) .* conj(phase_grad_m);
                                                
                                                % Take corrected ky lines.
                                                k_corr_m(ky_shot_even_ind_m(ind_shot,:),:) = ...
                                                    ft1(h_m_temp,2);
                                            else
                                                
                                                % Generate phase gradient
                                                phase_grad_m = repmat( ...
                                                    exp(1i*2*pi*(-Nx/2:(Nx/2)-1)/Nx*kxshift),...
                                                    length(ky_shot_odd_ind_m(ind_shot,:)),1);
                                                
                                                % Keep EVEN (+ve) ky lines. 'keeping ky lines'
                                                k_corr_m(ky_shot_even_ind_m(ind_shot,:),:) = ...
                                                    k_m(ky_shot_even_ind_m(ind_shot,:),:);
                                                
                                                % Correct ODD (-ve) ky lines.
                                                h_m_temp = ift1( ...
                                                    k_m(ky_shot_odd_ind_m(ind_shot,:),:) * ...
                                                    exp(1i*pi),2) .* conj(phase_grad_m);
                                                
                                                % Take corrected ky lines.
                                                k_corr_m(ky_shot_odd_ind_m(ind_shot,:),:) = ...
                                                    ft1(h_m_temp,2);
                                            end
                                        end % for ind_shot
                                        
                                        %-------------------------------------------------------------
                                    elseif DATAFLAGparams.epiCorrMethod==2 % phc
                                        
                                        % Set null kxshift.
                                        kxshift = [];
                                        
                                        % Take corrected ky lines.
                                        %k_corr_m = f_epi_phase_corr_phc(k_m,phc_data, ...
                                        %    DATAFLAGparams,DWIparams,ind_sl,ind_ec,ind_coil);
                                        
                                        % phc_data: {'kx','ky','kz','loca' or chunk,
                                        % 'dyn','card','echo','mix','chan','aver','sign','grad'};
                                        if strcmpi(DWIparams.scan_mode,'3D')
                                            %k_corr_m = f_epi_phase_corr_phc_v2(k_m,phc_data, ...
                                            %    DATAFLAGparams,DWIparams,1,ind_ec,ind_coil,ind_avg,ind_dw);
                                            k_corr_m = f_epi_phase_corr_phc_v3(k_m,phc_data, ...
                                                DATAFLAGparams,DWIparams,1,ind_chunk,ind_ec, ...
                                                ind_coil,ind_avg,ind_dw);
                                        else
                                            %k_corr_m = f_epi_phase_corr_phc(k_m,phc_data, ...
                                            %    DATAFLAGparams,DWIparams,ind_sl,ind_ec,ind_coil,ind_avg,ind_dw);
                                            
                                            % Use total variation.
                                            % Then
                                            % DATAFLAGparams.epiCorrFitMethod
                                            % doesn't mean anything in this
                                            % version.
                                            %k_corr_m = f_epi_phase_corr_phc_v2(k_m,phc_data, ...
                                            %    DATAFLAGparams,DWIparams,ind_sl,ind_ec,ind_coil,ind_avg,ind_dw);
                                            k_corr_m = f_epi_phase_corr_phc_v3(k_m,phc_data, ...
                                                DATAFLAGparams,DWIparams,ind_sl,ind_chunk,ind_ec, ...
                                                ind_coil,ind_avg,ind_dw);
                                            %showimage(jet,abs(ift2(k_m)),abs(ift2(k_corr_m)),'tot var')
                                        end
                                        
                                        % TEMP--START
%                                         for ind=[1,4]
%                                             DATAFLAGparams.epiCorrFitMethod = ind;
%                                             %k_corr_m = f_epi_phase_corr_phc_v2(k_m,phc_data, ...
%                                             %    DATAFLAGparams,DWIparams,ind_sl, ...
%                                             %    ind_ec,ind_coil,ind_avg,ind_dw);
%                                             k_corr_m = f_epi_phase_corr_phc_v3(k_m,phc_data, ...
%                                                 DATAFLAGparams,DWIparams,ind_sl, ...
%                                                 ind_chunk,ind_ec,ind_coil,ind_avg,ind_dw);
%                                             showimage(jet,abs(k_m),abs(ift2(k_m)),...
%                                                 abs(k_corr_m),abs(ift2(k_corr_m)), ...
%                                                 'k','',sprintf('k corr %d,coil %d',ind,ind_coil))
%                                             eval(sprintf('i_corr%d_m = abs(ift2(k_corr_m));',ind))
%                                         end
%                                         input('Press enter to close figures...')
%                                         close all
                                        % TEMP--END
                                    else
                                        error('s_read_DWI:main','Unknown DATAFLAGparams.epiCorrMethod')
                                    end                                    
                                    
                                    
                                    % TEST show folded image.
%                                     if ind_ec==1
%                                         showimage(jet,abs(k_m),abs(ift2(k_m)),abs(k_corr_m),abs(ift2(k_corr_m)), ...
%                                             sprintf('EC[%d]',ind_ec),'Original',sprintf('EC[%d]',ind_ec), ...
%                                             sprintf('epiCorrMethod[%d],kxshift[%.3f]', ...
%                                             DATAFLAGparams.epiCorrFitMethod,kxshift))
%                                     else
%                                         idx_shot = 1;
%                                         idx_shot_ky = idx_shot:DWIparams.nSHOT:DWIparams.EPI_FACTOR(2)*DWIparams.nSHOT;
%                                         ny_res = DWIparams.EPI_FACTOR(2)*DWIparams.nSHOT;
%                                         
%                                         showimage(jet,imresize(abs(ift2(k_m(idx_shot_ky,:))),[ny_res,size(k_m,2)]), ...
%                                             imresize(abs(ift2(k_corr_m(idx_shot_ky,:))),[ny_res,size(k_corr_m,2)]), ...
%                                             sprintf('EC[%d],original',ind_ec), ...
%                                             sprintf('EC[%d],epiCorrMethod[%d],kxshift[%.3f]', ...
%                                             ind_ec,DATAFLAGparams.epiCorrMethod,kxshift))
%                                     end
                                    
                                    
                                    % Output corrected k-space data.
%                                     if strcmpi(DWIparams.fold_over_dir,'AP')
%                                         if strcmpi(DWIparams.fat_shift_dir,'A')
%                                             %k_ec_corr_7d(:,:,ind_sl,ind_coil,ind_avg,ind_dw,ind_ori) = ...
%                                             %            flipud(fliplr(k_corr_m));
%                                             eval(sprintf('%s(:,:,ind_sl,ind_coil,ind_avg,ind_mult_dw) = flipud(fliplr(k_corr_m));', ...
%                                                 k_data_corr_name_s))
%                                         else
%                                             %k_ec_corr_7d(:,:,ind_sl,ind_coil,ind_avg,ind_dw,ind_ori) = ...
%                                             %    k_corr_m;
%                                             eval(sprintf('%s(:,:,ind_sl,ind_coil,ind_avg,ind_mult_dw) = k_corr_m;', ...
%                                                 k_data_corr_name_s))
%                                         end
%                                     elseif strcmpi(DWIparams.fold_over_dir,'RL')
%                                         if strcmpi(DWIparams.fat_shift_dir,'L')
%                                             % 2010.07.09: Don't forget to permute(...,[2,1]) and flipud(...)
%                                             % after recon.
%                                             %k_ec_corr_7d(:,:,ind_sl,ind_coil,ind_avg,ind_dw,ind_ori) = k_corr_m;
%                                             eval(sprintf('%s(:,:,ind_sl,ind_coil,ind_avg,ind_mult_dw) = k_corr_m;', ...
%                                                 k_data_corr_name_s))
%                                         else
%                                             % 2010.07.09: Check if this is correct. And don't forget to
%                                             % permute(...,[2,1]) and fliplr(...) after recon.
%                                             %k_ec_corr_7d(:,:,ind_sl,ind_coil,ind_avg,ind_dw,ind_ori) = k_corr_m;
%                                             eval(sprintf('%s(:,:,ind_sl,ind_coil,ind_avg,ind_mult_dw) = k_corr_m;', ...
%                                                 k_data_corr_name_s))
%                                         end
%                                     end
                                    
                                    if strcmpi(DWIparams.scan_mode,'3D')
                                        eval(sprintf('%s(:,:,ind_sl,ind_chunk,ind_coil,ind_avg,ind_mult_dw) = k_corr_m;', ...
                                            k_data_corr_name_s))
                                    else
                                        eval(sprintf('%s(:,:,ind_sl,ind_coil,ind_avg,ind_mult_dw) = k_corr_m;', ...
                                            k_data_corr_name_s))
                                    end
                                    clear  k_m  k_corr_m
                                
                                end % for ind_chunk                                
                            end % for ind_avg
                            %input(' ')
                            %close all
                        end % for ind_coil                        
                    end % for ind_sl (ind_kz)
                end % for ind_mult_dw
                
                % Clear data_nd.
                clear  data_nd
                
                % Delete data_nd__echo##_ori##.
                filename_s = sprintf('data_nd__ec%.2d_ori%.2d.mat',ind_ec,ori_idx);
                eval(sprintf('delete  %s',filename_s))
                fprintf('      [%s] is deleted\n',filename_s)
                
                % Report read raw data name.
                fprintf('\n')
                fprintf('      SENSE data [%s] is read for echo[%d].\n', ...
                    k_data_corr_name_s, ind_ec)
                
                % Save DWI data.
                fprintf('      Save [%s] data for echo[%d].\n', ...
                    k_data_corr_name_s,ind_ec)
                
                % Rename data.
                if ind_ec==1
                    eval(sprintf('K_alias_6d = %s;', k_data_corr_name_s))
                elseif ind_ec==2
                    eval(sprintf('K_nav_6d = %s;', k_data_corr_name_s))
                else
                    error('s_read_DWI:main','Unknown ind_ec')
                end
                
                % Clear data.
                eval(sprintf('clear  %s',k_data_corr_name_s))
                
                % Save each orientation data.
                cd(saveDataDir_s)
                if DWIparams.nEC==1
                    eval(sprintf('save  dw_data_ori%.2d__R%dS%d  K_alias_6d', ...
                        ori_idx,DWIparams.SENSE_FACTOR,DWIparams.nSHOT));
                    if strcmpi(DWIparams.scan_mode,'3D')
                        [ny,nx,nsl,nchunk,ncoil,navg,ndw] = size(K_alias_6d);
                    else
                        [ny,nx,nsl,ncoil,navg,ndw] = size(K_alias_6d);
                    end
                    clear  K_alias_6d
                    
                    fprintf('      dw_data_ori%.2d__R%dS%d is saved with [ny,nx,nsl,ncoil,navg,ndw]=[%d,%d,%d,%d,%d,%d]\n', ...
                        ori_idx,DWIparams.SENSE_FACTOR,DWIparams.nSHOT, ...
                        ny,nx,nsl,ncoil,navg,ndw)
                elseif DWIparams.nEC==2
                    if ind_ec==1
                        eval(sprintf('save  dw_data_ori%.2d__R%dS%d  K_alias_6d', ...
                            ori_idx,DWIparams.SENSE_FACTOR,DWIparams.nSHOT));
                        if strcmpi(DWIparams.scan_mode,'3D')
                            [ny_al,nx_al,nsl_al,nchunk_al,ncoil_al,navg_al,ndw_al] = size(K_alias_6d); % keep this over ind_ec loop                            
                        else
                            [ny_al,nx_al,nsl_al,ncoil_al,navg_al,ndw_al] = size(K_alias_6d);
                        end
                        clear  K_alias_6d
                        fprintf('      K_alias_6d is saved\n')
                    elseif ind_ec==2
                        eval(sprintf('save  dw_data_ori%.2d__R%dS%d  K_nav_6d  -append', ...
                            ori_idx,DWIparams.SENSE_FACTOR,DWIparams.nSHOT));
                        if strcmpi(DWIparams.scan_mode,'3D')
                            [ny_na,nx_na,nsl_na,nchunk_na,ncoil_na,navg_na,ndw_na] = size(K_nav_6d);
                        else
                            [ny_na,nx_na,nsl_na,ncoil_na,navg_na,ndw_na] = size(K_nav_6d);
                        end
                        clear  K_nav_6d
                        fprintf('      K_nav_6d is saved\n')
                        
                        % Report after saving K_alias_6d and K_nav_6d.
                        fprintf('\n')
                        fprintf('    dw_data_ori%.2d__R%dS%d is saved with\n', ...
                            ori_idx,DWIparams.SENSE_FACTOR,DWIparams.nSHOT)
                        if strcmpi(DWIparams.scan_mode,'3D')
                            fprintf('      [ny,nx,nsl,nchunk,ncoil,navg,ndw]: K_alias_6d[%d,%d,%d,%d,%d,%d,%d], K_nav_6d[%d,%d,%d,%d,%d,%d,%d]\n', ...
                                ny_al,nx_al,nsl_al,nchunk_al,ncoil_al,navg_al,ndw_al, ...
                                ny_na,nx_na,nsl_na,nchunk_na,ncoil_na,navg_na,ndw_na)
                        else
                            fprintf('      [ny,nx,nsl,ncoil,navg,ndw]: K_alias_6d[%d,%d,%d,%d,%d,%d], K_nav_6d[%d,%d,%d,%d,%d,%d]\n', ...
                                ny_al,nx_al,nsl_al,ncoil_al,navg_al,ndw_al, ...
                                ny_na,nx_na,nsl_na,ncoil_na,navg_na,ndw_na)
                        end
                    else
                        error('s_read_DWI:main','Unknown ind_ec')
                    end
                else
                    error('s_read_DWI:main','Unknown DWIparams.nEC')
                end % if DWIparams.nEC
                fprintf('\n')
                
            end % for ind_ec
        end % for ind_ori
    end % for ind_dw
    clear  k_ec*  phx_data  data_nd  data_epiref_nd  info  h*  ky_shot*  kp0*  kp1* p0* p1*
    clear  k_m  k_corr_m  k_m_odd*  k_m_even*  k_corr_m_odd*  k_corr_m_even*  phc_data
    
    % Report.
    tt = toc;
    fprintf('    SENSE raw data is read and extracted in [%f]min\n',tt/60)
    
    
end % if DATAFLAGparams.readDWI==1

% Pack memory.
pack

fprintf('\n\n')



%% Save DWI data



%% Show DWI data



%% END




