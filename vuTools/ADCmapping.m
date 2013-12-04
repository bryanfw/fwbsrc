function [adcmap] = ADCmapping(data, bvals, idx);
%ADCMAPPING Performs ADC mapping of MR data
%   ADCMAPPING calculates apparent ADC values for MR data using a linear
%   fit of T2-weighted images collected with different diffusion gradients.
%
%   'data' is the MR data organized as an M-by-N matrix, where M
%   corresponds to the number of different diffusion gradient
%   images, and N corresponds to the number of pixels one wishes to map.
%
%   'bvals' is a 1-by-M vector containing the b-values each row in 'data'.
%
%   'idx' is a 1-by-P (1<P<N) vector listing the column indices of 'data'
%   that should be fit.  All other columns are ignored.
%
%   Example
%   -------
%	data1 = load_MR('/some/semsdw_data1');
%       data2 = load_MR('/some/semsdw_data2');
%       data3 = load_MR('/some/semsdw_data3');
%       data4 = load_MR('/some/semsdw_data4');
%       bvals = [0, 100, 400, 800];
%       data = [data1(:); data2(:);data3(:);data4(:)];
%       idx = find(data1(:)>0.01);
%	ADCmap = ADCMapping(data, bvals, idx);
%
% $Id: ADCmapping.m 13 2006-07-17 23:14:31Z vuiis-svn $
num_bvals = size(data,1);
data_size = size(data,2);

adcmap = zeros(1,data_size);

XX=waitbar(0.01, 'Curve Fitting Progress . . .');
for ii = 1:length(idx),
    A = [bvals', ones(num_bvals,1)];
    b = log(data(:,idx(ii)));
    x = linsolve(A,b);
    adcmap(idx(ii)) = -x(1);
    waitbar(ii/length(idx));
end
close(XX);
