% find(obs.havetruth==3,1)
ind = 8; % have 3 truths

file = obs.intensity{8};

out = quick_label(file);

% home
file = '~/Code/fwbsrc/quick_label/sample_data/p10s2_rawimg.nii.gz';