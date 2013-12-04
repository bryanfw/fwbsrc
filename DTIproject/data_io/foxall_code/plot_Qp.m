%
% PLOT_QP  - Plot the Q parameter curves for the positive/negative echo correction from IEPI recon test 
%


figure

for S = 1:Ns
    subplot(Ns,1,S);
    plot(Zp,Q(S,:),'r.')
    xlabel('radians')
    title(['Qpn Plot - Shot ',int2str(S),':   Z min = ',num2str(Z_min(S)),'   Fp = ',num2str(Fp(S)),'   Fn = ',num2str(Fn(S))]);
    grid
end