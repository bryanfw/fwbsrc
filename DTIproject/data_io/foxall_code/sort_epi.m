%
% SORT_EPI - Calculate some basic properties of the input data 
%            set from the data description. Create a matrix of
%            all the STD data lines for onward reconstruction
%

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

