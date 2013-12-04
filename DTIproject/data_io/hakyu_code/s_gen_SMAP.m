
%[s_gen_SMAP] generates coil sensitivity map.
%
%
% Last modified
% 2010.04.07. Save individual slice of SMAP (3D data) so that this can be 
%             run parallel.
% 2010.08.09
%   This is generated from [SMAP.m].
% 2010.08.18.
%   Comment saving I_smap_3d__kz%.2d
%   Pack memory.
% 2010.12.09.
%   Add 'REFparams.filename' for saving 'I_smap'. This is to distinguish
%   the REF files for different slice orientations in 'shared' directory.
% 2012.05.14.
%   Add [s_gen_bodycoil_from_ref.m] to generate bodycoil equivalent using
%   reference coil images.
% 
% Ha-Kyu



%% Generate coil sensitivity map.

fprintf('Generate coil sensitivity map.\n')

% Generate and save or load coil sensitivity map.
if DATAFLAGparams.genSMAP==1
    t0 = clock;
    
    % Load REF data.
    cd(sharedDir_s)
    load  K_ref
    
    % Select bodycoil equivalent.
    if strcmpi(SMAPparams.gen_bodycoil,'body')
        % do nothing
    elseif strcmpi(SMAPparams.gen_bodycoil,'ref')
        % Run the script.
        s_gen_bodycoil_from_ref
        
        % Save the bodycoil equivalent.
        cd(sharedDir_s)
        save  K_body_from_ref  K_body
    else
        error('s_gen_SMAP:main','Unknown SMAPparams.gen_bodycoil')
    end
    
    % Load BODY data.
    cd(sharedDir_s)
    if strcmpi(SMAPparams.gen_bodycoil,'body')
        load  K_body
        fprintf('  Load K_body\n')
    elseif strcmpi(SMAPparams.gen_bodycoil,'ref')
        load  K_body_from_ref
        fprintf('  Load K_body_from_ref\n')
    else
        error('s_gen_SMAP:main','Unknown SMAPparams.gen_bodycoil')
    end
    
    % Generate and save coil sensitivity map.
    K_ref_4d = permute(K_ref,[2,1,3,4]);   % [Ky,Kx,Kz,COIL]->[Kx,Ky,Kz,COIL]
    K_body_3d = permute(K_body,[2,1,3]);   % [Ky,Kx,Kz]->[Kx,Ky,Kz]
    [nKx,nKy,nKz,nCOIL] = size(K_ref_4d);
    
    % output matrix size.
    I_smap_4d = zeros(nKy,nKx,nKz,nCOIL,'single');
    
    % generate coil sensitivity map.
    kz_v = 1:nKz;
    fprintf('    Processing SMAP for kz_REF=[%d ~ %d] of [%d].\n', ...
        kz_v(1),kz_v(end),nKz)
    for ind_kz = kz_v   % kz_REF_v comes from initial cell
        %fprintf('    Coil sensitivity map for Kz=[%d] being processed\n',ind_kz)
        SMAPparams.kz_REF = ind_kz;
        [I_smap_3d,mask_m] = f_coilsensemap(K_ref_4d,K_body_3d,SMAPparams);
        I_smap_3d = permute(I_smap_3d,[2,1,3]);
        mask_m = permute(mask_m,[2,1]);
        I_smap_4d(:,:,ind_kz,:) = I_smap_3d;
        
        % Save individual slice of SMAP.
        % This is for running separate cpu.
        %cd(sharedDir_s)
        %eval(sprintf('save  I_smap_3d__kz%.2d  I_smap_3d',ind_kz))
    end
    fprintf('\n\n')
    
    % Save smap.
    I_smap.smap = I_smap_4d;
    
    cd(sharedDir_s)
    %save  I_smap  I_smap
    eval(sprintf('save  I_smap__%s      I_smap',REFparams.filename))
    
    % Report.    
    fprintf('    Coil sensitivity map is generated and saved in [%f]sec.\n',etime(clock,t0))
end
clear  I_smap*  K_ref  K_body  mask_m  kz_v

% Pack memory.
pack

fprintf('\n\n')



%% END





