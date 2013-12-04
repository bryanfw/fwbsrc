function [listfile,error_flag] = ReadLIST_withPhaseCorrection(lstfile,conversion_spec)
% Read the data description parameters in an exported .list file obtained from a Philips scanner.
%
% Inputs:
%           If no inputs, via UI, select name of .LIST file
%           assumes that the corresponding .DATA file is in the same folder
%           With inputs:
%                       filename = filename (with path) of .LIST file
%
% Output:
%
% Steps:
% 1: Read LIST/DATA data in a matrix
%==========================================================================
% Author        Amol Pednekar
%               MR Clinical Science
%               Philips medical systems
%               12/29/2009
%==========================================================================
if nargin < 1
    [fname,pname] = uigetfile('*.LIST','Select *.LIST file');
    lstfile=[pname fname];
end

[lstparms, DVattribs] = read_listfile_either_fwb(lstfile);
[raw_data] = read_datafile([lstfile(1:(end-5)) '.data'], lstparms, DVattribs);
raw_data.std_count=find(cell2mat(strfind(cellstr(lstparms.lable),'STD')), 1, 'last' );


function [sort_data] = read_datafile(file_name, lstparms, DVattribs)
%--------------------------------------------------------------------
% CALCULATE DATA FILE SIZE FROM THE DATA DESCRIPTION
% READ THE .data FILE (IEEE FLOATS)
%--------------------------------------------------------------------
FLOAT32_BYTES = 4;
% index column numbers of data attributes
col=0;
for i=1:size(DVattribs{1},1)
    temp = char(DVattribs{1}(i));
    if (temp(1)~='#')
        temp(~isstrprop(temp, 'alphanum'))='_';
        DVidx.(temp) = col;col=col+1;
    end
end
total_calc_floats = (lstparms.attrib(size(lstparms.attrib,1),DVidx.offset)+...
    lstparms.attrib(size(lstparms.attrib,1),DVidx.size))/FLOAT32_BYTES;
% Open binary raw data
fid = fopen(file_name,'r');
[raw_data,raw_data_size] = fread(fid,inf,'float32');
fclose(fid);
% Validate file size
if (raw_data_size ~= total_calc_floats)
    error([ 'File size: expected ',num2str(total_calc_floats), ...
        ', reading ', num2str(raw_data_size)]);
end
%--------------------------------------------------------------------
% SORT DATA LINES BY TYPE INTO OUTPUT STRUCTURE
%--------------------------------------------------------------------
disp ('Sort');
std_cnt = 0;
copy_beg = 1; copy_end = 0; copy_off = 0;
for i = 1:size(lstparms.attrib,1)
    copy_beg = copy_beg + copy_off;
    copy_off = lstparms.attrib(i,DVidx.size)/FLOAT32_BYTES;
    copy_end = (copy_beg - 1) + copy_off;
    copy_typ = char(DVidx.typ);
    % sort by label, channel
    switch (char(lstparms.lable(i,:)))
        case 'STD'
%             disp(lstparms.attrib(i,DVidx.offset))
            std_cnt = std_cnt + 1;
            sort_data.std_data(std_cnt,:) = raw_data(copy_beg:copy_end);
            sort_data.std_attr(std_cnt,:) = lstparms.attrib(i,:);
        case 'REJ'
            
        case 'PHX'
            
        case 'FRX'
            
        case 'NAV'
            
        case 'NOI'
            
        otherwise
            error('.LIST file has invalid lable ',char(lstparms.lable(i,:)));
    end
end
% The Fourier transformation lengths for STD, NAV and REJ data of a SENSE
% scan are (by definition) equal to
% reconstruction matrix length * k-space oversample sample factor / SENSE factor
if (isempty(lstparms.X_direction_SENSE_factor))
    lstparms.X_direction_SENSE_factor = 1;
end
if (isempty(lstparms.Y_direction_SENSE_factor))
    lstparms.Y_direction_SENSE_factor = 1;
end
N_x = round( lstparms.X_resolution * lstparms.kx_oversample_factor / lstparms.X_direction_SENSE_factor );
N_y = round( lstparms.Y_resolution * lstparms.ky_oversample_factor / lstparms.Y_direction_SENSE_factor );
N_z = round( lstparms.Z_resolution * lstparms.kz_oversample_factor / lstparms.Z_direction_SENSE_factor );
if (lstparms.number_of_encoding_dimensions==2)
    kspace = zeros(max(sort_data.std_attr(:,DVidx.chan))+1,N_x,N_y);  % FIXME: bug here will be too big if channel numbers dont start at 0
else
    kspace = zeros(max(sort_data.std_attr(:,DVidx.chan))+1,N_x,N_y,N_z);
end
% Acquired kspace
K_x = lstparms.kx_range(2)-lstparms.kx_range(1)+1; KxOffset = ceil((N_x-K_x)/2);
K_y = lstparms.ky_range(2)-lstparms.ky_range(1)+1; KyOffset = ceil((N_y-K_y)/2);
if (~isempty(lstparms.kz_range))
    K_z = lstparms.kz_range(2)-lstparms.kz_range(1)+1; KzOffset = ceil((N_z-K_z)/2);
end

for I = 1:std_cnt
    ky = sort_data.std_attr(I,DVidx.ky);
    y_i = (ky - lstparms.ky_range(1) + 1)+ KyOffset;
    kspace(sort_data.std_attr(I,DVidx.chan)+1,KxOffset+1:end-KxOffset,y_i) = complex_unshuffle(sort_data.std_data(I,:));
    %     kspace(sort_data.std_attr(I,DVidx.chan),:,y_i) = sort_data.std_data(I,:);
end
% Recon only 1 channel
disp_chnum = round((min(sort_data.std_attr(:,DVidx.chan))+max(sort_data.std_attr(:,DVidx.chan)))/2);
k_dispCH = squeeze(kspace(disp_chnum+1,:,:));
% % Direct 2d/3D recon
raw_image = fftn(fftshift(k_dispCH),[N_x,N_y]);
figure, subplot(1,2,1),imshow(abs(raw_image),[]);
subplot(1,2,2),imshow(abs(log(k_dispCH)),[]);
% FFT in readout direction
fftx=fft(k_dispCH,N_x,1);
% 1D linear phase correction
phi_0 = lstparms.zero_order_phase_error_X;
phi_1 = lstparms.first_order_phase_error_X;
mult = [1:N_x];
phi = ((mult-N_x)*phi_1)+phi_0;
for iy = 1:N_y
    fftpx(:,iy) = fftx(:,iy).*exp(1i*phi)';
end
% FFT in phase direction
fftxy=fft(fftx,N_y,2);
fftpxy=fft(fftpx,N_y,2);
figure, 
subplot(2,2,1), imshow(abs(fftxy)',[]);
subplot(2,2,2), imshow(abs(fftpxy)',[]);
subplot(2,2,3), imshow(abs(fftpxy)'-abs(fftxy)',[]);
% subplot(2,3,1),imshow(abs((k_dispCH))',[]);
% subplot(2,3,2),imshow(abs(fftx)',[]);
% subplot(2,3,3),imshow(abs(fftxy)',[]);
% subplot(2,3,4),imshow(abs(fftpx)',[]);
% subplot(2,3,5),imshow(abs(fftpxy)',[]);
% subplot(2,3,6),imshow((abs(fftpxy)'-abs(fftxy)'),[]);
% subplot(2,3,5),imshow(abs(raw_image)'-abs(fftxy)',[]);

function complex_data = complex_unshuffle(real_data)
%--------------------------------------------------------------------------
% The odd number of real_data is the real component, while the even number
% of real_data is the imaginary component.
% 
%--------------------------------------------------------------------------
real_size = max(size(real_data));
if ( rem(real_size,2) ~= 0 )
    error(['complex unshuffle: odd input data size: ', num2str(real_size)]);
else
    half_size = real_size / 2;
end
x = 1 + (1:half_size);
odd = 1 + 2*((1:half_size) - 1);
even = odd + 1;
%to shift it one half period on Fourier transformation.
complex_data = (real_data(odd) + 1i*real_data(even)).*exp(1i*pi*x);


