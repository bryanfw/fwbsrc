%
% RECON MOVIE - Make and display movies of dynamic image series
%
%
%

GRAY_VALUES = 256;
GRAY_RANGE  = 96;
GRAY_CENTER = 128;

x_matrix  =  data_des.X_resolution;
y_matrix  =  data_des.Y_resolution;

sum_image = zeros(y_matrix,x_matrix);

for D = 1:Ndyn
    sum_image = sum_image + squeeze(abs(pic(D,:,:)));
end

mean_image = sum_image/Ndyn;

pic_data_scale = GRAY_RANGE/(max(max(mean_image)));

for D = 1:Ndyn
    ds_pic = GRAY_CENTER + pic_data_scale*squeeze(abs(pic(D,:,:)));
    image(ds_pic);
    axis('image')
    title(['magnitude image: ',int2str(D)]);
    colormap(gray(GRAY_VALUES));
    movie_image(D) = getframe(gca); 
end

movie2avi(movie_image,'data_recon_200_mag','FPS',5);