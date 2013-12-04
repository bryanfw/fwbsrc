%-----------------------------------------------------------------------------
% TEST SCRIPT - For navigator echoes
%-----------------------------------------------------------------------------


%-----------------------------------------------------------------------------
% Load Dynamic series
%-----------------------------------------------------------------------------
load data_raw_200

%-----------------------------------------------------------------------------
% Sort Series In Image Blocks
%-----------------------------------------------------------------------------
sort_epi;
sort_nav;

%----------------------------------------------------------------------------------------------------
% FFT   -  Perform a 1D FFT on all the navigator echoes to create profiles
%----------------------------------------------------------------------------------------------------

for D = 1:Ndyn
for S = 1:Ns
    nav_prof           = fft(fftshift(squeeze(nav_data(D,S,:))));
    nav_profile(D,S,:) = nav_prof;
    nav_A(D,S)         = abs(nav_profile(D,S,round(nav_len/2)));
    nav_P(D,S)         = angle(nav_profile(D,S,round(nav_len/2)));
end
end



figure

subplot(2,1,1)
plot(1:Ndyn,nav_A,'.')
xlabel('Dynamics')
ylabel(['Signal Amplitudes - ',int2str(Ns),' shots'])
title([plotname,' - Navigator Signal Amplitudes'])
grid


subplot(2,1,2);
plot(1:Ndyn,nav_P,'.')

xlabel('Dynamics')
ylabel(['Signal Phase (rads) - ',int2str(Ns),' shots'])
title([plotname,' - Navigator Signal Phases'])
grid

