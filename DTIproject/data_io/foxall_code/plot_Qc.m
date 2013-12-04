%
% Plot Qc  - Plot the Qc parameter for the shot combination phase correction
%


figure

for S=1:Ns

    subplot(Ns,1,S);
    plot(Zc,Qc(S,:),'b.')
    xlabel('radians')
    title([ 'Qc Plot - Shot',int2str(S),'     C min = ',num2str(C_min(S))   ]);
    grid

end