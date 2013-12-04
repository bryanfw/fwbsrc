% 
% RECON_NAV - Reconstruct the raw navigator in data_nav() compute the amplitude and phase differences
%             between navigators for each dynamic
%


%----------------------------------------------------------------------------------------------------
% FFT   -  Perform a 1D FFT on all the navigator echoes to create profiles
%
%          Compute the navigator amplitude as the amplitude of the central ordinate of the profile
%          Compute the navigator phase as the phase of the central ordinate of the profile
%----------------------------------------------------------------------------------------------------

for D = 1:Ndyn
for S = 1:Ns
    nav_prof           = fft(fftshift(squeeze(nav_data(D,S,:))));
    nav_profile(D,S,:) = nav_prof;
    nav_A(D,S)         = abs(nav_profile(D,S,round(nav_len/2)));
    nav_P(D,S)         = angle(nav_profile(D,S,round(nav_len/2)));
end
end



