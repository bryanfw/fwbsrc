%-----------------------------------------------------------------------------
% TEST SCRIPT
%-----------------------------------------------------------------------------


%-----------------------------------------------------------------------------
% Load Dynamic series
%-----------------------------------------------------------------------------
load data_raw_200

%-----------------------------------------------------------------------------
% Sort Series In Image Blocks
%-----------------------------------------------------------------------------
sort_epi;
sort_nav;


%-----------------------------------------------------------------------------
% Reconstruct Each Image In The Time Series
%-----------------------------------------------------------------------------
figure


for D = 1:Ndyn

       data_iepi     = squeeze(data_sort(D,:,:));
       echo_sign     = squeeze(data_sign(D,:));
       recon_results = recon_2_shot(data_iepi,echo_sign,data_des);
       Fp(D,:)       = recon_results.Fp;
       Fn(D,:)       = recon_results.Fn;
       Z(D,:)        = recon_results.Z;
       C(D,:)        = recon_results.C;
       pic(D,:,:)    = recon_results.combined_image;

       subplot(3,1,1)
       imagesc(abs(recon_results.combined_image));
       colormap(gray(256));
       axis('image');
       ylabel('Encode');
       xlabel('Read');
       title(['Dynamic Image = ',int2str(D)]);

       subplot(3,1,2)
       plot(recon_results.Zp,recon_results.Qp,'.');
       xlabel('radians');
       title(['Z = ',num2str(recon_results.Z),'    Fp = ',num2str(recon_results.Fp),'    Fn = ',num2str(recon_results.Fn)]);
       grid;

       subplot(3,1,3)
       plot(recon_results.Zc,recon_results.Qc,'.');
       xlabel('radians');
       title(['C = ',num2str(recon_results.C)]);
       grid;

       pause(2);

end

%-----------------------------------------------------------------------------
% Reconstruct Each Navigator In The Time Series
%-----------------------------------------------------------------------------
recon_nav;


