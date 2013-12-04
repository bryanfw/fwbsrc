function Q = ifft_Q(data)

projection = sum(abs(transpose(data)));
ft_proj    = fftshift(fft(fftshift(projection)));
Q          = -1.0*sum(abs(ft_proj));