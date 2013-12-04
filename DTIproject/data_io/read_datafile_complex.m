function [complex_img] = read_datafile_complex(file_name, lstparms, DVattribs)
% Read '.data' file for complex data, modified from Amol's codes for raw data
% Inputs: file_name: cpx_XXX.data
%         lstparms, DVattribs: outputs from .list file
% Outputs: complex_img:  a 6D array with complex number for 2D scan,
%                        (Nx, Ny, Nphase, Nlocation, Nchannel, Ndyn)
% Example:
% [complex_img] = read_datafile_complex('cpx_040(1028).data', lstparms, DVattribs);
%==========================================================================
% Author        Hui Wang, Ph.D student
%               University of Louisville, KY
%               12/20/2010
% Updated: 6/20/2011
% Updated: 7/11/2011  added dynamic scan dimension in complex_img
% Updated: 8/15/2011  read by image, instead of by lines, make the program
%                     faster.
%==========================================================================
fprintf('Reading .DATA file ...');

%--------------------------------------------------------------------
% CALCULATE DATA FILE SIZE FROM THE DATA DESCRIPTION
% READ THE .data FILE (IEEE FLOATS)
%--------------------------------------------------------------------
FLOAT32_BYTES = 4;
% index column numbers of data attributes
% these index numbers are used in lstparms.attrib array
col=0;
for i=1:size(DVattribs{1},1)
    temp = char(DVattribs{1}(i));
    if (temp(1)~='#')
        temp(~isstrprop(temp, 'alphanum'))='_'; % all non-alphanum becomes '_'
        DVidx.(temp) = col;col=col+1;
    end
end
% (The last vector offset + last vector size)/FLOAT32_BYTES, hui
total_calc_floats = (lstparms.attrib(size(lstparms.attrib,1),DVidx.offset)+...
    lstparms.attrib(size(lstparms.attrib,1),DVidx.size))/FLOAT32_BYTES;
% Open and read binary complex data
fid = fopen(file_name,'r');
[raw_data,raw_data_size] = fread(fid,inf,'float32');
fclose(fid);
% Validate file size
if (raw_data_size ~= total_calc_floats)
    error([ 'File size: expected ',num2str(total_calc_floats), ...
        ', reading ', num2str(raw_data_size)]);
end
%--------------------------------------------------------------------------
% SORT DATA N_y*LINES BY TYPE INTO OUTPUT STRUCTURE
% This will accelarate the sorting process dramatically
% sort_data is a 2*1 structure, where sort_data.std_data is the data vector
% which corresponds to one image, std_cnt is the total number of images,
% which is equal to N_chan*N_card*N_cyn*..., hui
%--------------------------------------------------------------------------
N_x = lstparms.X_resolution;
N_y = lstparms.Y_resolution;
N_z = lstparms.Z_resolution;
N_card = lstparms.number_of_cardiac_phases;
N_loca = lstparms.number_of_locations;
N_dyn = lstparms.number_of_dynamic_scans;

std_cnt = 0;
% Assumes from the first STD vector, all are STD vectors
% makes the sort_data struct in which the raw data is better organized
beg = min(strmatch('STD', lstparms.lable)); % find 1st STD label
for i = beg :N_y :size(lstparms.attrib,1) 
    copy_beg = 1+lstparms.attrib(i,DVidx.offset)/FLOAT32_BYTES; % get starting byte from correct column (index given by DVidx.offset)
    copy_end = (lstparms.attrib(i+N_y-1,DVidx.size)+...
        lstparms.attrib(i+N_y-1,DVidx.offset))/FLOAT32_BYTES;
    std_cnt = std_cnt + 1;
    sort_data.std_data(std_cnt,:) = raw_data(copy_beg:copy_end);
    sort_data.std_attr(std_cnt,:) = lstparms.attrib(i,:);
end

% Initialize complex image matrix
if (isempty(lstparms.number_of_coil_channels))
    N_chan = 1;   % for Q-body, the there is no channel # in .list file
else
    N_chan = lstparms.number_of_coil_channels;
end

% For the case like two channels with number 16, 17 
% (FWB - moved above initialization )
[unique_chan, m_chan, n_chan] = unique(sort_data.std_attr(:,DVidx.chan));

if (lstparms.number_of_encoding_dimensions==2)  % 2D imaging
    complex_img = zeros(N_x, N_y, N_card, N_loca, N_chan, N_dyn);
else % 3D imaging
    complex_img = zeros(N_x, N_y, N_z, N_card, N_loca, N_chan, N_dyn);
end

% FIXME: In the case of SENSE ref scan, the body coil is another "location" wiht
% a unique coil number so end up with  matrix that is filled with lots of zeros
% for loc1 coil 9 and  location 2, coils 1:8 (for 8chan head coil)

%--------------------------------------------------------------------------
% Arrange the vector into image matrix
%--------------------------------------------------------------------------
for I = 1:std_cnt
    cardiacphase = sort_data.std_attr(I,DVidx.card);
    loca = sort_data.std_attr(I,DVidx.loca)+1;
    dyn_i = sort_data.std_attr(I,DVidx.dyn)+1;
    z_loc = sort_data.std_attr(I,DVidx.z)+1; % (FWB - this line was missing)
    complex_img(:,:,z_loc,cardiacphase+1,loca,n_chan(I),dyn_i) = reshape(complex_unshuffle(sort_data.std_data(I,:)), N_x, N_y);
%     figure, imshow(abs(complex_img(:,:,cardiacphase+1,loca,n_chan(I),dyn_i)),[]);
end

fprintf(' complete.\n');

end

%%
function complex_data = complex_unshuffle(real_data)
%--------------------------------------------------------------------------
% The odd number of real_data is the real component, while the even number
% of real_data is the imaginary component.
% Hui Wang, 12/20/2010
%--------------------------------------------------------------------------
real_size = max(size(real_data));
if ( rem(real_size,2) ~= 0 )
    error(['complex unshuffle: odd input data size: ', num2str(real_size)]);
else
    half_size = real_size / 2;
end
odd = 1 + 2*((1:half_size) - 1);
even = odd + 1;
complex_data = real_data(odd) + 1i*real_data(even);

%to shift it one half period on Fourier transformation.
%x = 1 + (1:half_size);
%complex_data = (real_data(odd) + 1i*real_data(even)).*exp(1i*pi*x);
end


