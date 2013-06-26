function loglike = loglikelihood(B, Y, X, robust, g00, g11)

% gamma 0,0: P(obs labe =0 | true label = 0) aka specificity (about .7)
g01 = 1-g00;  

% g11 = 1;      % P(obs label =1 | true_label =1) aka sensitivity (about 1)
g10 = 1-g11; 


Yhat = 1 ./ (1 + exp(- X* B)); 

if ~robust
    loglike = -1 * sum( Y.*log(Yhat) + (1-Y).*log(1-Yhat)); 
elseif robust
    loglike = - 1 * sum(Y.*log(g11*Yhat+g01*(1-Yhat)) + (1-Y).*log(g00*(1-Yhat)+g10*Yhat));
else
    error('Need to enter robust parameter as 0/1');
end
