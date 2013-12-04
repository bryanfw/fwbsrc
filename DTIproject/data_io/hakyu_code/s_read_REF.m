
%[s_read_REF] read, save and show REF data in LIST/DATA format.
%
%
% Last modified
% 2010.07.26. 
%   Change psi to PSI to remove a conflict with Matlab function 'psi()'.
% 2010.08.09.
%   This is generated from [REF_readsaveshow.m]. This only contains the
%   last comment from the file.
% 2010.08.18.
%   Remove cells for 'Save ...' and 'Show ...' part.
%   Pack memory.
% 2010.09.06.
%   Define bodycoil position in 32 channel coil at 7T.
% 2010.09.09.
%   Consider actual body-coil position, 26th(0~) at 7T and 27th(0~) at 3T.
%   But [test_f_read_list_data.m] considers this (as of 2010.09.09) to reduce
%   the size of data. Then keep using nCHAN for reading body-coil data.
% 2010.09.16.
%   Consider 'if nCHAN < body_coil_position' to correctly read bodycoil.
% 2010.09.17.
%   Modify the method for reading bodycoil. Divided into field strength.
% 2010.09.18.
%   Use REFparams.nCOIL instead of DWIparams.nCOIL.
% 2010.09.20.
%   Use [test_f_read_listdata.m] saved as [f_read_listdata.m].
% 2010.09.22.
%   Modify reading data_body_3d.
% 2010.12.21.
%   Generate info_list out of [f_read_listdata.m] for fast reading as done
%   in [s_read_DWI.m].
% 2011.08.29.
%   Adjust body_coil_position for 3T NV16 coil.
% 2011.11.15.
%   Save noi_data. It can be used for SNR calculation.
%
% Ha-Kyu



%% Read REFERENCE and BODYCOIL data

fprintf('Read REFERENCE and BODYCOIL data per each NSA.\n')

% Read data.
if DATAFLAGparams.readREF==1
    
    %----------------------------------------
    % Generate info struct for .LIST file.
    
    % Get .LIST file, listtext.
    name = REFparams.filename;
    toks = regexp(name,'^(.*?)(\.list|\.data)?$','tokens');
    prefix = toks{1}{1};
    listname = sprintf('%s.list',prefix);
    
    % Read list file, listtext.
    cd(loadDir_s)
    fid = fopen(listname,'r');
    if fid~=-1,
        listtext = fread(fid,inf,'uint8=>char')';
        fclose(fid);
    else
        error( sprintf('cannot open %s for reading', listname) );
    end
    clear  fid
    fprintf('    LIST file [%s] is read\n',listname)
    
    % Set attributes.
    attr = {'mix','dyn','card','echo','loca','chan','extr1','extr2','ky', ...
        'kz','n.a.','aver','sign','rf','grad','enc','rtop','rr','size','offset'};
    
    % Clean attr names (will be used as variablenames and fieldnames)
    for k=1:length(attr),
        attr{k} = f_cleanFieldname( attr{k} );
    end
    
    % Set pattern for tags.
    pattern = '(?<typ>\w+)';
    for k=1:length(attr),
        pattern = sprintf('%s\\s+(?<%s>-?\\d+)',pattern,attr{k});
    end
    
    %----- Check point #2 ----- start
    tic
    info_list = regexp(listtext, pattern, 'names');
    t = toc;
    fprintf('    CP#2: %.3f sec\n',t)
    %----- Check point #2 ----- end
    
    %----------------------------------------
    
    
    % Read data.
    cd(loadDir_s)
    [ref_info, ref_data,phx_data,frx_data,noi_data] = f_read_listdata(REFparams.filename,info_list);
    clear  phx_data  frx_data  ref_info  info_list
    [nKx,nKy,nKz,nLOC,nDYN,nCARD,nECHO,nMIX,nAVG,nCHAN] = size(ref_data);
    
    
    % Calculate noise correlation matrix.
    m = squeeze(noi_data(:,1,1:REFparams.nCOIL));    %[kx,loc,chan]->[kx,chan]
    cd(sharedDir_s)
    save  noi_data  noi_data
    clear  noi_data
    PSI = zeros(REFparams.nCOIL,REFparams.nCOIL);
    for ind1 = 1:REFparams.nCOIL
        for ind2= 1:REFparams.nCOIL
            PSI(ind1,ind2) = sum(m(:,ind1).*conj(m(:,ind2))) / ...
                sqrt(sum(m(:,ind1).*conj(m(:,ind1)))*sum(m(:,ind2).*conj(m(:,ind2))));
        end
    end
    clear  m
    
    
    % Average over NSA and extract valid coil data.
    data_5d = squeeze( mean(ref_data,9) );
    clear  ref_data
    fprintf('    ref_data averaged and extracted for valid coil.\n')
    
    
    % Divide into coil reference and bodycoil images. Bodycoil data is in
    % 26th (7T) and 27th (3T) coil position starting from 0. 
    data_ref_4d = squeeze(data_5d(:,:,:,1,1:REFparams.nCOIL));
    if GENparams.B0==70000
        body_coil_position = REFparams.bodycoil_idx+1;
        if nCHAN < body_coil_position
            body_coil_position = nCHAN;
            data_body_3d = squeeze(data_5d(:,:,:,2,body_coil_position));
        else            
            data_body_3d = squeeze(data_5d(:,:,:,2,body_coil_position));
        end
        clear data_5d
    elseif GENparams.B0==30000
        coil_idx_offset = 16; % coil idx starts from 16 when nCOIL is less than 16
        if REFparams.nCOIL < coil_idx_offset
            offset = coil_idx_offset;
            coil_idx_max = offset+REFparams.nCOIL-1;
        else
            offset = 0;
            coil_idx_max = offset+REFparams.nCOIL-1;
        end
            
        if coil_idx_max > REFparams.bodycoil_idx
            body_coil_position = REFparams.bodycoil_idx-offset+1;
            data_body_3d = squeeze(data_5d(:,:,:,2,body_coil_position));
        else
            sl_temp = floor(size(data_5d,3)/2);
            body_coil_position_temp = zeros(1,size(data_5d,5));
            for ind_temp = 1:size(data_5d,5)
                m_temp = squeeze(data_5d(:,:,sl_temp,2,ind_temp));
                if isempty(find(m_temp(:)))==1
                    body_coil_position_temp(ind_temp) = 0;
                else
                    body_coil_position_temp(ind_temp) = 1;
                end
            end
            if length(find(body_coil_position_temp))>1
                error('body_coil_position is > 1')
            end
            body_coil_position = find(body_coil_position_temp);                
            %data_body_3d = squeeze(data_5d(:,:,:,2,nCHAN));
            data_body_3d = squeeze(data_5d(:,:,:,2,body_coil_position));            
        end
        clear  data_5d  body_coil_position_temp
    else
        error('s_read_REF:main','Unknown GENparams.B0')
    end
    fprintf('      body coil position is %d\n',body_coil_position)
    
    
    % Revert RF phase cycling.
    for ind_kz = 1:size(data_ref_4d,3)
        for ind_c = 1:size(data_ref_4d,4)
            data_ref_4d(:,1:2:end,ind_kz,ind_c) = data_ref_4d(:,1:2:end,ind_kz,ind_c)*exp(1i*pi);
            if ind_c==1     % run once for bodycoil image
                data_body_3d(:,1:2:end,ind_kz) = data_body_3d(:,1:2:end,ind_kz)*exp(1i*pi);
            end
        end
    end    
    fprintf('    RF phase cycling is reverted.\n')
    
    
    % Output data.
    ref_raw_extr = permute(data_ref_4d,[2,1,3,4]);
    clear  data_ref_4d
    fprintf('    ref_raw_extr is generated.\n')
    body_raw_extr = permute(data_body_3d,[2,1,3]);
    clear  data_body_3d
    fprintf('    body_raw_extr is generated.\n')
    
    % Report.
    fprintf('\n')
    fprintf('    ref_raw_extr,\n')
    fprintf('        [nKy,nKx,nKz,nCOIL] = [%d,%d,%d,%d]\n',size(ref_raw_extr))
    fprintf('    body_raw_extr,\n')
    fprintf('        [nKy,nKx,nKz] = [%d,%d,%d]\n',size(body_raw_extr))
    fprintf('    PSI,\n')
    fprintf('        [nCOIL,nCOIL] = [%d,%d]\n',size(PSI))
    
    
    % Save reference and bodycoil data.
    K_ref = ref_raw_extr;   % [nKy,nKx,nKz,nCOIL]
    K_body = body_raw_extr(:,:,:,1);  % only the first coil has bodycoil image
    clear  ref_raw_extr  body_raw_extr
    
    % Save reference and bodycoil data.
    cd(sharedDir_s)
    save    K_ref   K_ref
    save    K_body  K_body
    save    PSI     PSI
    
    % Report.
    fprintf('    K_ref, K_body are saved.\n')
    
    % Clear data.
    clear  K_ref  K_body  PSI
    
end % if DATAFLAGparams.readREF

% Pack memory.
pack

fprintf('\n\n')



%% Save REFERENCE and BODYCOIL data



%% Show REFERENCE and BODYCOIL data



%% END




