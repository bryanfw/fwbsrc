function [t1map, residual] = T1FlipMapping(data, fliplist, idx);
%T1FLIPMAPPING Performs T1 mapping of multi-flip MR data
%   T1FLIPMAPPING calculates apparent T1 values for MR data using a non-linear
%   fit of T1-weighted images collected with different flip angles.
%
%   'data' is the MR data organized as an M-by-N matrix, where M
%   corresponds to the number of flip angle images, and N
%   corresponds to the number of pixels ones wishes to map.
%
%   'fliplist' is a 1-by-M vector containing the flip angles (in
%   radians) for each row in 'data'.
%
%   'idx' is a 1-by-P (1<P<N) vector listing the column indices of 'data'
%   that should be fit.  All other columns are ignored.
%
%   Notes
%   -----
%       Due to some constraints in MATLAB you will need to specify
%       a global variable 'TR' initialized to (you guessed it) TR.
%       
%       We can fix this, but I'm too lazy to... :P
%
%   Example
%   -------
%       global TR;
%       TR = 200e-3; %milliseconds
%	data1 = load_MR('/some/gems_data1');
%       data2 = load_MR('/some/gems_data2);
%       data3 = load_MR('/some/gems_data3');
%       data4 = load_MR('/some/gems_data4');
%       data5 = load_MR('/some/gems_data5');
%       data6 = load_MR('/some/gems_data6');
%       fliplist = [15, 30, 45, 60, 75, 90]*pi/180;
%       data = [data1(:);data2(:);data3(:);data4(:);data5(:);data6(:)];
%       idx = find(data1(:)>0.015);
%	T1map = T1FlipMapping(data, fliplist, idx);
%
% $Id: T1FlipMapping.m 6 2006-06-28 19:42:29Z vuiis-svn $

warning off;
num_angles = size(data,1);
data_size = size(data,2);

t1map = zeros(1,data_size);
residual = zeros(1,data_size);

mean_decay = mean(data(:,idx),2);
stupid_guess = [0.1,0.1,0.1];
[init_cond, R]=nlinfit(fliplist, mean_decay, @t1flip, stupid_guess);

XX=waitbar(0.01, 'Curve Fitting Progress . . .');
for ii = 1:length(idx),
    decay=data(:,idx(ii));
    [beta, R]=nlinfit(fliplist, decay, @t1flip, init_cond);
    t1map(idx(ii))=beta(2);
    residual(idx(ii))=sqrt(sum(R.*R));
    waitbar(ii/length(idx));
end
close(XX)
warning on;



function Y = t1flip(beta, X)
global TR
if beta(2) <= 0,
    Y = beta(1)*sin(X)+beta(3);
else,
    Y = (beta(1)*(1 - exp(-TR/beta(2)))*sin(X)./(1-cos(X)*exp(-TR/beta(2))))+beta(3);
end
Y=Y';
