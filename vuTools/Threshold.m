orgIm = im.Data(:,:,64);
for i = 1:256
    for j = 1:256
        if orgIm(i,j,1) > 135
            orgIm(i,j,1) = 256;
        else
            orgIm(i,j,1) = 0;
        end
    end
end
corrIm = corr_im(:,:,1);
for i = 1:256
    for j = 1:256
        if corrIm(i,j,1) > 135
            corrIm(i,j,1) = 256;
        else
            corrIm(i,j,1) = 0;
        end
    end
end

subIm = im.Data(:,:,1)-corr_im(:,:,1);
