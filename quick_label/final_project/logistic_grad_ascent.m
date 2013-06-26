function B = logistic_grad_ascent(X, Y)

close all;

fore = find(Y ==1);
back = find(Y ==0);
figure(4); hold on; 
plot(X(back,2), X(back,3),'ro');
plot(X(fore,2), X(fore,3),'bx');
hold off
xl = xlim;
yl= ylim;

[n,m] = size(X); 

B = [-8;.15;-.5]; %initial estimate
eta = [.1; .01; .01]; % gradient change parameter
maxiters = 1000; 

oldexpy = -ones(size(Y));
iter = 0;
for i = 1:maxiters
    
    iter = iter + 1; 
    
    expy = 1 ./ (1 + exp(-X*B)); % expected value of y = sigmoid funct/ logistic funct

    figure(5);imagesc(expy); colorbar; caxis([0 1]);
    

   
    grad = X' * (Y - expy) / n; %grad(1) =0; 
    
    B = B + eta .* grad;    
    y = -xl * B(2)/B(3) - B(1)/B(3);
    likelihood= sum(Y .* log (expy) + (1-Y) .* log(1-(expy-.000000000001)));
    
    figure(4); hold on; plot(xl,y,'k');axis([xl yl]); hold off;
    
    fprintf('%3d: [', iter);
    fprintf(' %g', B(1:end)); 
    fprintf(' ] loglikelihood = %4.3f\n',  likelihood);
    
    if sum(abs(expy - oldexpy)) < 10e-10; 
        fprintf('Converged\n')
        return; 
    end
%     fprintf('\nBhat using "hill climber"   = '); fprintf('%4.2f\t',B); fprintf('B2/B3 = %4.2f, B1/B3 = %4.2f', B(2)/B(3),B(1)/B(3));
    
    oldexpy = expy; 
end

 
 fprintf('\nBhat using "hill climber"   = '); fprintf('%4.2f\t',B); fprintf('B2/B3 = %4.2f, B1/B3 = %4.2f', B(2)/B(3),B(1)/B(3));
 
  keyboard;
 
% expy = 1 ./ (1 + exp(-X*B)
%{    
Bhat2_1 = logistic(X,Y,[],[],[]); 
fprintf('\nBhat using "logistic" = '); fprintf('%4.2f\t',Bhat2_1); 
fprintf('B2/B3 = %4.2f, B1/B3 = %4.2f', Bhat2_1(2)/Bhat2_1(3),Bhat2_1(1)/Bhat2_1(3));
%}

B1 = -8.16; % close to optimum (-8.1643)
% B2 = -2:.02:2; % includes optimum (.1633)
B2 = 0.1633;
B3 = -2:.02:2; % includes optimum (-.5582)

likelihoodmat = zeros(length(B2), length(B3));

for i = 1:length(B2);
    for j = 1:length(B3);
        B = [B1;B2(i);B3(j)]; 
    
        expy = 1 ./ (1 + exp(-X*B));
        likelihoodmat(i,j) = sum(Y .* log (expy) + (1-Y) .* log(1-(expy-.000000000001)));
    end
    fprintf('.');
end
fprintf('\n');

surf(B2, B3, likelihoodmat);

