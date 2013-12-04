%
% PARREC MOVIE - Make and display movies of a dynamic image series with window and centering to make
%                EPI ghosting visible with the image. 
%
% Preconditions - The times series has been read into MATLAB with the read_in_data.m script
%
% Output        - A .AVI file is created from stored screen shots with the movie2avi command.
%


GRAY_VALUES = 256;
GRAY_RANGE  = 96;
GRAY_CENTER = 128;

output_filename = 'matlab_movie';

Ndyn     = parfile.max_dynamics;
x_matrix = parfile.scan_tags(10);
y_matrix = parfile.scan_tags(11);

sum_image = zeros(y_matrix,x_matrix);

for D = 1:Ndyn
    sum_image = sum_image + transpose(squeeze(abs(rs_data(D,:,:))));
end

mean_image    = sum_image/Ndyn;
rs_data_scale = GRAY_RANGE/(max(max(mean_image)));


for D = 1:Ndyn
    ds_pic = GRAY_CENTER + rs_data_scale*transpose(squeeze(abs(rs_data(D,:,:))));
    image(ds_pic);
    axis('image')
    title(['magnitude image: ',int2str(D)]);
    colormap(gray(GRAY_VALUES));
    movie_image(D) = getframe(gca); 
end

movie2avi(movie_image,output_filename,'FPS',5);