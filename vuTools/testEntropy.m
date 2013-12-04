im = vuOpenImage('Brain4.dat');
clear e
kSpace = fftshift(fft2(im));
kStart = 127;
kStop = size(kSpace,2);
count = 1;
for t = -0.1:0.01:0.1
    e(count) = EntropyMetric(t,kSpace,kStart,kStop,1,size(kSpace,2),2);
    count = count +1;
end
figure,plot(e,'.')

%%
clear e
clear x
im = load('Brain4.dat');
fft_im = fft2(im);
shifted_fft_im = fftshift(fft_im);
kStart = 100;
shift = exp(-i*2*pi/256 .* [0:256-1] .* 0.01 * 256);
shift = repmat(shift,[256-kStart+1 1]);
shifted_fft_im(kStart:256,:) = shifted_fft_im(kStart:256,:).*shift;
count = 1;
for t = -0.1:0.01:0.1
    e(count) = EntropyMetric(t,shifted_fft_im,kStart,256,1,256,1);
    x(count) = t;
    count = count +1;
end
figure,plot(x,e,'.')
%%
clear e
clear x
im = load('Brain1.dat');
fft_im = fft2(im);
countx = 1;
for s = -0.2:0.01:0.2
    shifted_fft_im = fftshift(fft_im);
    kStart = 64;
    shift = exp(-i*2*pi/256 .* [0:256-1] .* s * 256);
    shift = repmat(shift,[32 1]);
    shifted_fft_im(kStart:kStart+31,:) = shifted_fft_im(kStart:kStart+31,:).*shift;
    county = 1;
    for t = -0.2:0.01:0.2
        e(county,countx) = EntropyMetric(t,shifted_fft_im,kStart,kStart+31,1,256,2);
        x(county,countx) = t;
        county = county +1;
    end
    countx = countx +1;
end
%figure,plot(x,e,'.')