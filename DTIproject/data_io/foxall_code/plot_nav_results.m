figure

subplot(2,1,1)
plot(1:Ndyn,nav_A,'.')
xlabel('Dynamics')
ylabel(['Signal Amplitudes - ',int2str(Ns),' shots'])
title([plotname,' - Navigator Signal Amplitudes'])
grid


subplot(2,1,2);
plot(1:Ndyn,nav_P,'.')

xlabel('Dynamics')
ylabel(['Signal Phase (rads) - ',int2str(Ns),' shots'])
title([plotname,' - Navigator Signal Phases'])
grid

