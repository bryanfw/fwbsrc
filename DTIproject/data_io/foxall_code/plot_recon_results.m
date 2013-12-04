%
% PLOT_RECON_RESULTS - Plot results from a time series reconstruction
%
%

figure
subplot(4,2,1)
plot(1:Ndyn,Fp(:,1),'r.');
grid
xlabel('Dynamic')
ylabel('Phase Correction (rads)')
title('First Order Phase Correction - Positive Echoes - Shot 1')

subplot(4,2,2)
plot(1:Ndyn,Fp(:,2),'b.')
grid
xlabel('Dynamic')
ylabel('Phase Correction (rads)')
title('First Order Phase Correction - Positive Echoes - Shot 2')

subplot(4,2,3)
plot(1:Ndyn,Fn(:,1),'r.')
grid
xlabel('Dynamic')
ylabel('Phase Correction (rads)')
title('First Order Phase Correction - Negative Echoes - Shot 1')

subplot(4,2,4)
plot(1:Ndyn,Fn(:,2),'b.')
grid
xlabel('Dynamic')
ylabel('Phase Correction (rads)')
title('First Order Phase Correction - Negative Echoes - Shot 2')


subplot(4,2,5)
plot(1:Ndyn,Z(:,1),'r.')
grid
xlabel('Dynamic');
ylabel('Phase Correction (rads)');
title('Zero Order Phase Correction - Shot 1');


subplot(4,2,6)
plot(1:Ndyn,Z(:,2),'b.')
grid
xlabel('Dynamic');
ylabel('Phase Correction (rads)');
title('Zero Order Phase Correction - Shot 2');


subplot(4,2,7)
plot(1:Ndyn,C(:,2),'g.')
grid
xlabel('Dynamic');
ylabel('Phase Correction (rads)');
title('Combination Zero Order Phase Correction - Shots 1 and 2');

subplot(4,2,8)
plot(1:Ndyn,nav_P(:,1),'r.',1:Ndyn,nav_P(:,2),'b.')
grid
xlabel('Dynamic');
ylabel('Navigator Phase (rads)');
title('Navigator Phase - Shots 1 and 2');
