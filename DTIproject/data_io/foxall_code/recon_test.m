%
% Basic IEPI Reconstruction Script
%
%

%---------------------------------------------------------------
% Read the data description in the .LIST file
%---------------------------------------------------------------
filename              = 'raw_201';
plotname              = strrep(filename,'_','-');
listfile              = [filename,'.list'];
datafile              = [filename,'.data'];
[data_des,error_flag] = read_list(listfile,'SSEPI');


%---------------------------------------------------------------
% Read the raw data from the .DATA file
%---------------------------------------------------------------
data_raw              = read_data(data_des,datafile);


%-------------------------------------------------------------
% Data Description
%
% Find the size of the data description and set up the 
% data description keys
%-------------------------------------------------------------
ddi_sizes = size(data_des.data_attributes);
ddi_nkeys = ddi_sizes(2);
Natt      = ddi_sizes(1);

switch (ddi_nkeys)
     case 19
          ddi_19keys
     case 20
          ddi_20keys
     otherwise
          error('Check the number of data attributes in the .LIST file');
end


%-------------------------------------------------------------
% Data Sizes
%
% The STD data lines will each have Nx data points in them
% There will be Ns*number_of_gradient_echoes or Ny lines per
% kspace. Where Ns is the number of shots taken to fill kspace.
% 
%-------------------------------------------------------------
Nx   = round(data_des.kx_range(2) - data_des.kx_range(1) + 1);
Nx   = Nx*data_des.kx_oversample_factor;
Npts = Nx/2;

Ny   = round(data_des.ky_range(2) - data_des.ky_range(1) + 1);
Ny   = Ny*data_des.ky_oversample_factor;
Ns   = round(Ny/data_des.number_of_gradient_echoes);

Nstd = data_raw.std_count;
Ndyn = data_des.number_of_dynamic_scans;



%-------------------------------------------------------------
% Data Sort
%
% Sort the STD data lines by dynamic, echo bin, and points
%-------------------------------------------------------------
data_sort = zeros(Ndyn,Ny,Npts);
data_sign = zeros(Ndyn,Ny);
std_I     = 0;

for I = 1:Natt
    if (char(data_des.data_type(I)) == 'STD')
       std_I                 = std_I + 1;   
       ky                    = data_des.data_attributes(I,ddi.ky);
       echo                  = round(ky - data_des.ky_range(1) + 1);
       dyn                   = data_des.data_attributes(I,ddi.dyn) + 1;
       data_sort(dyn,echo,:) = complex_unshuffle(squeeze(data_raw.std(std_I,:))); 
       data_sign(dyn,echo)   = data_des.data_attributes(I,ddi.sign);
    end
end





%--------------------------------------------------------------------------------------
% INPUT Management
%--------------------------------------------------------------------------------------
data_iepi = squeeze(data_sort(1,:,:));
echo_sign = squeeze(data_sign(1,:));



%--------------------------------------------------------------------------------------
% Create blank kspace for each epi shot and sort the input data into it.
% Create blank pspace for the positive echoes and sort them into it.
% Create blank nspace for the negative echoes and sort them into it.
%
% sum of kspace(1,:,:) + kspace(2,:,:) ... should give back the input data
% sum of pspace(1,:,:) + pspace(2,:,:) ... should give all the postive echoes
% sum of nspace(1,:,:) + nspace(2,:,:) ... should give all the negative echoes
% sum of pspace sum and the nspace sum should give back the input data. 
%
%--------------------------------------------------------------------------------------

for S = 1:Ns
    kspace(S,:,:) = zeros(size(data_iepi));
    pspace(S,:,:) = zeros(size(data_iepi));
    nspace(S,:,:) = zeros(size(data_iepi));
end


for S = 1:Ns
for Y = 1:Ny
    shot_index = 1 + rem(Y-1,Ns);

    if (shot_index == S)
        kspace(S,Y,:) = squeeze(data_iepi(Y,:));

        if (echo_sign(Y) > 0)
           pspace(S,Y,:) = squeeze(data_iepi(Y,:));
        else
           nspace(S,Y,:) = squeeze(data_iepi(Y,:));
        end
    end

end
end


%----------------------------------------------------------------------------------------
% ZERO FILL, FOURIER TRANSFORM & DECIMATE OVERSAMPLED REGION
%----------------------------------------------------------------------------------------
x_matrix  =  data_des.X_resolution;
x_osample =  data_des.kx_oversample_factor;
x_ZF      =  round((x_matrix*x_osample - Npts)/2);
x_OS      =  round(((x_matrix*x_osample) - x_matrix)/2);
x_beg     =  x_ZF;
x_end     =  x_ZF + Npts - 1;
x_OS_beg  =  x_OS;
x_OS_end  =  x_OS + x_matrix - 1;


y_matrix  =  data_des.Y_resolution;
y_osample =  data_des.ky_oversample_factor;
y_ZF      =  round((y_matrix*y_osample - Ny)/2);
y_beg     =  y_ZF;
y_end     =  y_ZF + Ny - 1;



for S=1:Ns
    PK_space                          =  zeros(y_matrix*y_osample,x_matrix*x_osample);
    PK_space(y_beg:y_end,x_beg:x_end) =  squeeze(pspace(S,:,:));
    NK_space                          =  zeros(y_matrix*y_osample,x_matrix*x_osample);
    NK_space(y_beg:y_end,x_beg:x_end) =  squeeze(nspace(S,:,:));  
    P_image                           =  fft2(fftshift(PK_space));
    N_image                           =  fft2(fftshift(NK_space));
    pimage(S,:,:)                     =  P_image(:,x_OS_beg:x_OS_end);
    nimage(S,:,:)                     =  N_image(:,x_OS_beg:x_OS_end);
end


%----------------------------------------------------------------------------------------
% LINEAR PHASE CORRECTION ESTIMATES 
%
% Made using the Ahn & Cho method applied to the whole matrix in the read direction to 
% compensate for the echo shift between positive and negative echoes.
% 
% Linear phase corrected copies of the sub images are made.
%----------------------------------------------------------------------------------------

for S = 1:Ns

    pim          = squeeze(pimage(S,:,:));
    nim          = squeeze(nimage(S,:,:));
    auto_cor_pim = pim.*conj(scrollY(pim,1));
    auto_cor_nim = nim.*conj(scrollY(nim,1));

    for Y = 1:y_matrix
        fpc_pim(Y) = 2*x_matrix*angle(sum(auto_cor_pim(Y,:)));
        fpc_nim(Y) = 2*x_matrix*angle(sum(auto_cor_nim(Y,:)));
    end

    mask_pim     = sum(abs(transpose(pim)));
    mask_nim     = sum(abs(transpose(nim)));
    Fp(S)        = sum(fpc_pim.*mask_pim)/sum(mask_pim);
    Fn(S)        = sum(fpc_nim.*mask_nim)/sum(mask_nim);

    xval         = 1:x_matrix;
    pim_pc       = (xval - x_matrix/2 - 0.5)*Fp(S)/(2*x_matrix);
    nim_pc       = (xval - x_matrix/2 - 0.5)*Fn(S)/(2*x_matrix);
    pexp         = exp(i*pim_pc);
    nexp         = exp(i*nim_pc);
    
    for Y = 1:y_matrix
        fpc_pimg(Y,:) = pim(Y,:).*pexp;
        fpc_nimg(Y,:) = nim(Y,:).*nexp;
    end

    fpc_pimage(S,:,:) = fpc_pimg;
    fpc_nimage(S,:,:) = fpc_nimg;

end

%----------------------------------------------------------------------------------------------
% RECONSTRUCTION FILTERS
%----------------------------------------------------------------------------------------------
pn_filter                           = zeros(y_matrix,x_matrix);
c_filter                            = zeros(y_matrix,x_matrix);
filter_1                            = round(y_matrix/8.0);
filter_2                            = round(y_matrix/4.0);

f_beg_1                             = filter_1;
f_end_1                             = filter_1 + filter_2 - 1;
f_beg_2                             = filter_1 + 2*filter_2;
f_end_2                             = filter_1 + 3*filter_2 - 1;
pn_filter(f_beg_1:f_end_1,:)        = ones(filter_2,x_matrix);
pn_filter(f_beg_2:f_end_2,:)        = ones(filter_2,x_matrix);

c_beg_1                             = 1;
c_end_1                             = filter_1;
c_beg_2                             = y_matrix - filter_1 + 1;
c_end_2                             = y_matrix;
c_filter(c_beg_1:c_end_1,:)         = ones(filter_1,x_matrix);
c_filter(c_beg_2:c_end_2,:)         = ones(filter_1,x_matrix);



%----------------------------------------------------------------------------------------------
% ZERO PHASE CORRECTION ESTIMATES
%
% Made by combining the postive and negative echo images for each shot at different zero 
% order phase corrections. The figure of merit: Q total sum of the absval of the combined
% data will be minimum at the desired phase correction. The minimum is found by brute force.
%----------------------------------------------------------------------------------------------
Nz        = 180;
Zp        = 2*pi*((1:Nz) - 1)/(Nz - 1);


for S = 1:Ns

    pim    = squeeze(fpc_pimage(S,:,:));
    nim    = squeeze(fpc_nimage(S,:,:));

for Z = 1:Nz
    zpc    = exp(i*Zp(Z))*ones(size(pim));
    zpc_im = abs(pim + nim.*zpc);
    Q(S,Z) = sum_Q(zpc_im.*pn_filter);
end
    

    [min_val,min_index] = min(Q(S,:));
    Z_min(S)            = Zp(min_index);
    zpc                 = exp(i*Z_min(S))*ones(size(pim));
    zpc_image(S,:,:)    = pim + nim.*zpc;

end





%-------------------------------------------------------------------------------------------------
% SHOT COMBINATION CORRECTION
%-------------------------------------------------------------------------------------------------
cpc_image = squeeze(zpc_image(1,:,:));
C_min     = zeros(1,Ns);
Qc        = zeros(Ns,Nz);
Zc        = 2*pi*((1:Nz) - 1)/(Nz - 1);

for S = 2:Ns
    sim = squeeze(zpc_image(S,:,:));
   
for Z = 1:Nz
    cpc     = exp(i*Zc(Z))*ones(size(sim));
    cpc_im  = abs(cpc_image + sim.*cpc);
    Qc(S,Z) = sum_Q(cpc_im.*c_filter);
end
    
    [min_val,min_index] = min(squeeze(Qc(S,:)));
    C_min(S)            = Zc(min_index);
    cpc                 = exp(i*C_min(S))*ones(size(sim));
    cpc_image           = cpc_image + sim.*cpc;
end



