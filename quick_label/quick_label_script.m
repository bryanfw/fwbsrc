% find(obs.havetruth==3,1)
ind = 1; %


% 
% file = '~/Code/fwbsrc/quick_label/sample_data/p10s2_rawimg.nii.gz';% home
file = obs.intensity{ind};

% out = quick_label(file);

%% TESTING - coomment out later
% get some testing data in

vol = load_nii_gz(file);
invol = vol.img; % cast to double

% saturate btwen .5 and 99.5 percentile
invol = double(invol);
lowval = prctile(invol(:), .5);
highval = prctile(invol(:), 99.5);
invol(invol<lowval) = lowval;
invol(invol>highval) = highval;
% scale to one;
invol = (invol-lowval)/(highval-lowval);
clear lowval highval

% flip flop dimensions until happy
invol = flipdim(permute(invol,[2 1 3]),1); 

% pxdim stuff
pixdim = vol.hdr.dime.pixdim(2:4);
pixdim = pixdim([2 1 3]);
axial_aspect = pixdim.^-1;

% other stuff 
I = imread('eight.tif');

%% real work
% out = quick_label(file);

