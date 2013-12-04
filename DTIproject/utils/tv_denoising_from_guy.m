addpath toolbox_image/
addpath toolbox_image/toolbox/


% Removing noises from the SENSE-produced images
disp(' smoothing the phase maps ... ');
PDmap1r = real(PDmap1);
PDmap1i = imag(PDmap1);
PDmap2r = real(PDmap2);
PDmap2i = imag(PDmap2);
PDmap3r = real(PDmap3);
PDmap3i = imag(PDmap3);
PDmap4r = real(PDmap4);
PDmap4i = imag(PDmap4);
mr1 = max(PDmap1r(:));
if option4 == 2
    hanning20 = hanning(floor(xsize/16))*hanning(floor(xsize/16))';
    for cnt = 1:number_of_images
        tmp = conv2(PDmap1r(:,:,cnt),hanning20,'same');
        PDmap1r(:,:,cnt)=tmp;
        tmp = conv2(PDmap1i(:,:,cnt),hanning20,'same');
        PDmap1i(:,:,cnt)=tmp;
        tmp = conv2(PDmap2r(:,:,cnt),hanning20,'same');
        PDmap2r(:,:,cnt)=tmp;
        tmp = conv2(PDmap2i(:,:,cnt),hanning20,'same');
        PDmap2i(:,:,cnt)=tmp;
        tmp = conv2(PDmap3r(:,:,cnt),hanning20,'same');
        PDmap3r(:,:,cnt)=tmp;
        tmp = conv2(PDmap3i(:,:,cnt),hanning20,'same');
        PDmap3i(:,:,cnt)=tmp;
        tmp = conv2(PDmap4r(:,:,cnt),hanning20,'same');
        PDmap4r(:,:,cnt)=tmp;
        tmp = conv2(PDmap4i(:,:,cnt),hanning20,'same');
        PDmap4i(:,:,cnt)=tmp;
    end
end
if option4 == 1
    options.lambda = 0.3*mr1;
    options.niter_inner = 5;
    options.niter =  50;
    for cnt = 1:number_of_images
        tmp = perform_tv_denoising(PDmap1r(:,:,cnt),options);
        PDmap1r(:,:,cnt)=tmp;
        tmp = perform_tv_denoising(PDmap1i(:,:,cnt),options);
        PDmap1i(:,:,cnt)=tmp;
        tmp = perform_tv_denoising(PDmap2r(:,:,cnt),options);
        PDmap2r(:,:,cnt)=tmp;
        tmp = perform_tv_denoising(PDmap2i(:,:,cnt),options);
        PDmap2i(:,:,cnt)=tmp;
        tmp = perform_tv_denoising(PDmap3r(:,:,cnt),options);
        PDmap3r(:,:,cnt)=tmp;
        tmp = perform_tv_denoising(PDmap3i(:,:,cnt),options);
        PDmap3i(:,:,cnt)=tmp;
        tmp = perform_tv_denoising(PDmap4r(:,:,cnt),options);
        PDmap4r(:,:,cnt)=tmp;
        tmp = perform_tv_denoising(PDmap4i(:,:,cnt),options);
        PDmap4i(:,:,cnt)=tmp;
    end
end
PDmap1sm = PDmap1r + 1i*PDmap1i;
PDmap2sm = PDmap2r + 1i*PDmap2i;
PDmap3sm = PDmap3r + 1i*PDmap3i;
PDmap4sm = PDmap4r + 1i*PDmap4i;
phasemap1_s = PDmap1sm ./ abs(PDmap1sm);
phasemap2_s = PDmap2sm ./ abs(PDmap2sm);
phasemap3_s = PDmap3sm ./ abs(PDmap3sm);
phasemap4_s = PDmap4sm ./ abs(PDmap4sm);
