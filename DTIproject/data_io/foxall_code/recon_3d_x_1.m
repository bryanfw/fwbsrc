%
% RECON 3D X 1
%
% Reconstruct a 3D raw data set held in raw data format files from 
% a Philips NV scanner. 
%
% FORMAT
%
% raw_image = recon_3d_x_1(list_filename,data_filename,conversion_spec)
%
% INPUT
%
% list_filename        (string)    Pathname of the .list file that holds the
%                                  data format description.
% data_filename        (string)    Pathname of the .data file that holds the
%                                  binary data from the scanner.
% conversion_spec      (string)    The M-FILE holding linekeys and string
%                                  conversion specifications used to interpet
%                                  the .list file by read_list.m. 
%
% OUTPUT
%
% raw_image            (float[])   A three dimensional array holding the full
%                                  3D complex image formed from the binary data. 
%                                  Oversampled regions are not stripped. 
%
% M_FILE DEPENDENCY
%
% read_list.m           Used to read the data description from the  .list file
% read_data.m           Used to read the raw data from the .data file
% complex_unshuffle.m   Used to unshuffle the raw_data lines and apply the complex
%                       phase factor needed to center the data along the read out
%                       (X) direction.
%
% $Id: recon_3d_x_1.m,v 1.2 2005/04/21 15:00:08 dfoxall Exp $
%

function raw_image = recon_3d_x_1(list_filename,data_filename,conversion_spec)

%--------------------------------------------------------------------
% CONSTANTS
%
% ddi                 (structure)    Data description indices
%
%--------------------------------------------------------------------
ddi.mix    =   1;
ddi.dyn    =   2;
ddi.card   =   3;
ddi.echo   =   4;
ddi.loca   =   5;
ddi.chan   =   6;
ddi.extr1  =   7;
ddi.extr2  =   8;
ddi.ky     =   9;
ddi.kz     =   10;
ddi.aver   =   11;
ddi.sign   =   12;
ddi.rf     =   13;
ddi.grad   =   14;
ddi.enc    =   15;
ddi.rtop   =   16;
ddi.rr     =   17;
ddi.size   =   18;
ddi.offset =   19;   

FLOAT32_BYTES = 4;


%--------------------------------------------------------------------
% READ THE DATA DESCRIPTION PARAMETERS
% READ RAW DATA
%--------------------------------------------------------------------
[data_des,error_flag] = read_list(list_filename,conversion_spec);
raw_data              = read_data(data_des,data_filename);


%--------------------------------------------------------------------
% CREATE ONE 3D K SPACE
%--------------------------------------------------------------------
N_d    = raw_data.std_count;
N_x    = round( data_des.X_resolution * data_des.kx_oversample_factor );
N_y    = round( data_des.Y_resolution * data_des.ky_oversample_factor );
N_z    = round( data_des.Z_resolution * data_des.kz_oversample_factor );
kspace = zeros(N_x,N_y,N_z);


for I = 1:N_d
    ky                = raw_data.std_att(I,ddi.ky);
    kz                = raw_data.std_att(I,ddi.kz);
    y_i               = ky - data_des.ky_range(1) + 1;
    z_i               = kz - data_des.kz_range(1) + 1;
    kspace(:,y_i,z_i) = complex_unshuffle(raw_data.std(I,:));
end


%--------------------------------------------------------------------
% 3D FFT
%--------------------------------------------------------------------
raw_image = fftn(fftshift(kspace));


%--------------
% END
%--------------

