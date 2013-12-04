function mi = mi(im1,im2,nBins)
    joint_pdf = zeros(nBins,nBins);
    for i = 1:length(im1(:))
        joint_pdf(im1(i)+1,im2(i)+1) = joint_pdf(im1(i)+1,im2(i)+1)+1;
    end
    joint_pdf = joint_pdf./length(im1(:));
    
    ln_joint_pdf = joint_pdf;
    ln_joint_pdf(ln_joint_pdf==0) = 1;
    ln_joint_pdf = log(ln_joint_pdf);
    
    fixed_pdf = sum(joint_pdf,2);
    ln_fixed_pdf = fixed_pdf;
    ln_fixed_pdf(ln_fixed_pdf==0) = 1;
    ln_fixed_pdf = log(ln_fixed_pdf);
    
    moving_pdf = sum(joint_pdf,1);
    ln_moving_pdf = moving_pdf;
    ln_moving_pdf(ln_moving_pdf==0) = 1;
    ln_moving_pdf = log(ln_moving_pdf);
    
    fixed_entropy = sum(fixed_pdf.*ln_fixed_pdf);
    moving_entropy = sum(moving_pdf.*ln_moving_pdf);
    joint_entropy = sum(joint_pdf(:).*ln_joint_pdf(:));
    
    mi = -(fixed_entropy+moving_entropy) ./joint_entropy;