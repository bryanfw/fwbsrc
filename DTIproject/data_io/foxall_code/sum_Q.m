%
% SUM_Q - Return the value for the sum of the magnitude of the input data matrix 
%         as the quality factor for the positive to negative echo correction.

function Q = sum_Q(data)

Q = sum(sum(abs(data)));