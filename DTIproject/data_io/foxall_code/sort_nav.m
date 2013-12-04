%
% SORT_NAV
%
% Extract and sort navigator NAV vectors for use in reconstruction
%

nav_size  = data_raw.nav_size;
nav_len   = round(nav_size/2);
nav_data  = zeros(Ndyn,Ns,nav_len);


for D = 1:Ndyn
for S = 1:Ns
    nav_I           = (Ndyn - 1)*Ns + S;
    profile         = squeeze(data_raw.nav(nav_I,:));
    nav_data(D,S,:) = complex_unshuffle(profile);
end
end