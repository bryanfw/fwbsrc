
%[s_gen_bodycoil_from_ref] generates bodycoil image and k-space data from
%reference coil images. The final bodycoil image and its k-space data is
%already removed for object phase.
%
% This is useful when no bodycoil images are acquired.
%
%
% Last modified
% 2011.06.06.
% 2012.05.14.
%   Use ift3() and ft3().
%
% Ha-Kyu



%% Generate bodycoil image data from reference image data

% Generate image-space reference data.
I_ref = zeros(size(K_ref));
for ind_coil = 1:size(K_ref,4)
    %I_ref(:,:,:,ind_coil) = ifftshift(ifftshift(ifftn( ...
    %    fftshift(fftshift(K_ref(:,:,:,ind_coil),2),1) ),2),1);
    I_ref(:,:,:,ind_coil) = ift3(K_ref(:,:,:,ind_coil));
end
% I_ref = circshift(I_ref,[0,0,-1,0]); % this will be performed in [f_coilsensemap.m]
%[Ny,Nx,Nz,Nc]=size(I_ref);

% Generate image-space bodycoil data from reference data.
I_body = sqrt(sum((I_ref.*conj(I_ref)),4));
%K_body = fftshift(fftshift(fftn( ifftshift(ifftshift(I_body,2),1) ),2),1);
%K_body = ft3(I_body);
%I_body1 = ifftshift(ifftshift(ifftn( fftshift(fftshift(K_body,2),1) ),2),1); % same as I_body



%% Calculate object phase -- this comes from [smap_HKJ.m]

[Ny,Nx,Nz,Nc] = size(I_ref);

row = floor(Ny/2)+1;
col = floor(Nx/2)+1;
% I_body1 = I_body*0;
s = zeros(Ny,Nx);
for ind_sl = 1:Nz
    I_body_m = I_body(:,:,ind_sl);
    for ind_coil = 1:Nc
        p = squeeze(I_ref(:,:,ind_sl,ind_coil));
        a = angle(p);
        
        psi_c = a(row,col);
        s = s + abs(p).*p.*exp(-1i*psi_c);
    end
    psi = angle(s);
    I_body_m = I_body_m.*exp(1i*psi);
    I_body(:,:,ind_sl) = I_body_m; % object phase removed bodycoil image data
end
%K_body = fftshift(fftshift(fftn( ifftshift(ifftshift(I_body,2),1) ),2),1); % object phase removed k-space bodycoil data
K_body = ft3(I_body);


%% NOTE
% Save K_body as k-space bodycoil data. This must be put in [s_gen_SMAP.m] before,
%
% ----
% Generate and save coil sensitivity map.
%     K_ref_4d = permute(K_ref,[2,1,3,4]);   % [Ky,Kx,Kz,COIL]->[Kx,Ky,Kz,COIL]
%     K_body_3d = permute(K_body,[2,1,3]);   % [Ky,Kx,Kz]->[Kx,Ky,Kz]
% ----
%
% above procedure. Then proceed as they are for the remainders.















