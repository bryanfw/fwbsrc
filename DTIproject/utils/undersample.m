function imgu = undersample(img,Rx,Ry)

imgu=ft2(img);

maskx = zeros(size(imgu));
masky = zeros(size(imgu));
maskx(:,1:Rx:end) = 1; 
masky(1:Ry:end,:) = 1;

imgu = imgu .* maskx .* masky;

imgu=ift2(imgu);
