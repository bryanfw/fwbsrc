function t2map = T2mapping(data, TE, idx);
%T2MAPPING Performs T2 mapping of multi-echo MR data.
%   T2MAPPING calculates apparent T2 values for MR data using a linear
%   fit of multi-echo T2-weighted images.
%
%   'data' is the multi-echo MR data organized as an M-by-N matrix, where
%   M corresponds to the number of echos and N corresponds to the number of
%   pixels ones wishes to map.
%
%   'TE' is a 1-by-M vector containing the echo times for each row in 'data'.
%
%   'idx' is a 1-by-P (1<P<N) vector listing the column indices of 'data'
%   that should be fit.  All other columns are ignored.
%
%   Example
%   -------
%	data1 = load_MR('/some/sems_data1');
%       data2 = load_MR('/some/sems_data2');
%       TE = [12, 24]*1e-3; %milliseconds
%       data = [data1(:); data2(:)];
%       idx = find(data1(:)>0.005);
%	T2map = T2Mapping(data, TE, idx);
%
% $Id: T2mapping.m 16 2006-08-18 16:20:03Z vuiis-svn $
% $Log:$

num_TE = size(data,1);
data_size = size(data,2);

t2map = zeros(1,data_size);
XX=waitbar(0.01, 'Curve Fitting Progress . . .');
for ii = 1:length(idx),
    A = [TE', ones(num_TE,1)];
    b = log(data(:,idx(ii)));
    x = linsolve(A,b);
    t2map(idx(ii)) = -1./x(1);
    waitbar(ii/length(idx));
end
close(XX);
msgbox('I am done now');
