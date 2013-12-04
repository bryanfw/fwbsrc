 % MR Motion correction based on k-space data
% Based on Automatic Compensation of Motion Artifacts in MRI
% By: Kevin Wilson and Tuhin Sinha
% Date: June 12, 2007
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

im = load('Brain4.dat');
fft_im = fft2(im);
shifted_fft_im = fftshift(fft_im);
kStart = 160;
kStop = 256;
shift = exp(-i*2*pi/256 .* [0:256-1] .* -0.022 * 256);
shift = repmat(shift,[kStop-kStart+1 1]);
shifted_fft_im(kStart:kStop,:) = shifted_fft_im(kStart:kStop,:).*shift;

corr_kSpace = translationMotionCorrection(shifted_fft_im,4,16);