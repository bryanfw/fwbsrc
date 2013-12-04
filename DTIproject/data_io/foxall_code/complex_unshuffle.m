%
% COMPLEX UNSHUFFLE
%
% Convert an array of real numbers stored as (real,imag) pairs
% into an array of MATLAB complex numbers. Apply a exp(i*pi*x)
% factor to the data to shift it one half period on Fourier 
% transformation.
%
% $Id: complex_unshuffle.m,v 1.2 2005/04/21 14:19:29 dfoxall Exp $
%

function complex_data = complex_unshuffle(real_data);

real_size = max(size(real_data));

if ( rem(real_size,2) ~= 0 )
   error(['complex unshuffle: odd input data size: ', num2str(real_size)]);
else
   half_size = real_size / 2;
end

x            = 1 + (1:half_size);
odd          = 1 + 2*((1:half_size) - 1);
even         = odd + 1;
complex_data = (real_data(odd) + i*real_data(even)).*exp(i*pi*x);

%--------------------
% END
%--------------------

