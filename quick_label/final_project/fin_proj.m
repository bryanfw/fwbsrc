% simulation example - to get the robustness thing to work

% 2 Distributions: 2-dimensional multi-variate guassian 
mu1 = [10 16 ];
var1 = [100 150 ];
corr1 = rand(2); corr1 = triu(corr1) + triu(corr1)'; corr1(logical(diag(ones(2,1)))) = ones(2,1); % a valid random corr matrix
sig1 = sqrt(var1)'* sqrt(var1) .* corr1;

mu2 = [-10 2];
var2 = [60 59];
corr2 = rand(2); corr2 = triu(corr2) + triu(corr2)'; corr2(logical(diag(ones(2,1)))) = ones(2,1); % a valid random corr matrix
sig2 = sqrt(var2)'* sqrt(var2) .* corr2;

% make some samples using Multi Variate Normal RaNDom numbers function
T = mvnrnd(mu1,sig1,1000); 
F = mvnrnd(mu2,sig2,1000);
% figure; hold on; 
% plot(T(:,1), T(:,2),'x');
% plot(F(:,1), F(:,2),'ro');
% axis([-10 120 -4 10]);
% hold off

% Good data! ------ Y1
True_Y = [ones(1,1000), zeros(1,1000)]';
Y1 = True_Y;
X = horzcat( ones(2000,1),[T;F]); % column of ones, stuck to T stacked on top of F

% Bad data! ------- Y2
percent_T = sum(Y1 == 1)/length(Y1);
percent_F = 1-percent_T;
false_pos = 0.250;
Y2 = [ones(1, floor(length(Y1)*(percent_T+false_pos*percent_T))) , zeros(1, length(Y1)-floor(length(Y1)*(percent_T+false_pos*percent_T)))]';

%using funct found online
Bhat1_1 = logistic(X,Y1,[],[],[]); % make sure this includes column of 1's
Bhat2_1 = logistic2(X,Y2,[],[],[]);

%using glmfit - cant feed this thing the whole X or it flips out (no 1's column)
Bhat1_2 = glmfit(X(:,2:3),Y1,'binomial','link','logit');  
Bhat2_2 = glmfit(X(:,2:3),Y2,'binomial','link','logit'); 

%using fminsearch
[Bhat1_3, ML1_3] = fminsearch(@loglikelihood, [0 0 0]', optimset('TolFun', 0.0001), Y1, X, 0,[],[]);
[Bhat2_3, ML2_3] = fminsearch(@loglikelihood, [0 0 0]', optimset('TolFun', 0.0001), Y2, X, 0, [],[]); %#ok<NASGU>

% ROBUST METHOD
fprintf('\nROBUST ITERATIVE METHOD'); 
g00_init = .7; g00 = g00_init; old_g00= 0;
g11_init = .99; g11 = g11_init;
max_its = 25; iter = 0;

[Bhat3_3, ML2_3] = fminsearch(@loglikelihood, [0 0 0]', optimset('TolFun', 0.0001), Y2, X, 1, g00_init, g11_init); % once and then iterate
while (abs(old_g00 - g00) > .001 && iter < max_its)  
    old_g00 = g00;
    iter = iter+1;
    fprintf('\n %u. g00 = %.3f \t g11 = %.3f ', iter, g00, g11);
    Yhat = 1./(1+exp(-X*Bhat3_3));
    g00 = sum(Yhat < .5 & Y2==0)/sum(Yhat<.5); g01 = 1-g00;
    g11 = sum(Yhat >= .5 & Y2==1)/sum(Yhat>.5); g10 = 1-g11; 
    [Bhat3_3, ML2_3] = fminsearch(@loglikelihood, [0 0 0]', optimset('TolFun', 0.0001), Y2, X, 1, g00, g11); % once and then iterate
end
iter = iter+1;
fprintf('\n %u. g00 = %.3f \t g11 = %.3f ', iter, g00, g11);
fprintf('\n');



% print results
fprintf('\nUSING GOOD DATA');
fprintf('\nBhat using "logistic"     = '); fprintf('%4.2f\t',Bhat1_1); fprintf('B2/B3 = %4.2f, B1/B3 = %4.2f', Bhat1_1(2)/Bhat1_1(3),Bhat1_1(1)/Bhat1_1(3));
fprintf('\nBhat using "glmfit"       = '); fprintf('%4.2f\t',Bhat1_2); fprintf('B2/B3 = %4.2f, B1/B3 = %4.2f', Bhat1_2(2)/Bhat1_2(3),Bhat1_2(1)/Bhat1_2(3));
fprintf('\nBhat using "fminsearch"   = '); fprintf('%4.2f\t',Bhat1_3); fprintf('B2/B3 = %4.2f, B1/B3 = %4.2f', Bhat1_3(2)/Bhat1_3(3),Bhat1_3(1)/Bhat1_3(3));
fprintf('\n');
% 
fprintf('\nUSING DATA WITH HIGH FALSE POS RATE');
fprintf('\nBhat using "logistic"     = '); fprintf('%4.2f\t',Bhat2_1); fprintf('B2/B3 = %4.2f, B1/B3 = %4.2f', Bhat2_1(2)/Bhat2_1(3),Bhat2_1(1)/Bhat2_1(3));
fprintf('\nBhat using "glmfit"       = '); fprintf('%4.2f\t',Bhat2_2); fprintf('B2/B3 = %4.2f, B1/B3 = %4.2f', Bhat2_2(2)/Bhat2_2(3),Bhat2_2(1)/Bhat2_2(3));
fprintf('\nBhat using "fminsearch"   = '); fprintf('%4.2f\t',Bhat2_3); fprintf('B2/B3 = %4.2f, B1/B3 = %4.2f', Bhat2_3(2)/Bhat2_3(3),Bhat2_3(1)/Bhat2_3(3));
fprintf('\n');
% 
fprintf('\nUSING DATA WITH HIGH FALSE POS RATE WITH ROBUSTNESS ACTIVATED');
fprintf('\nBhat using "fminsearch"   = '); fprintf('%4.2f\t',Bhat3_3); fprintf('B2/B3 = %4.2f, B1/B3 = %4.2f', Bhat3_3(2)/Bhat3_3(3),Bhat3_3(1)/Bhat3_3(3));
fprintf('\n');

% plot results
figure; 
hold on; 
fore = find(Y2 ==1);
back = find(Y2 ==0);
plot(X(back,2), X(back,3),'ro');
plot(X(fore,2), X(fore,3),'bx');
xl = xlim;
yl= ylim;
y1_1 = -xl * Bhat1_1(2)/Bhat1_1(3) - Bhat1_1(1)/Bhat1_1(3); % good data - logistic
y1_2 = -xl * Bhat1_2(2)/Bhat1_2(3) - Bhat1_2(1)/Bhat1_2(3); % good data - glmfit
y1_3 = -xl * Bhat1_3(2)/Bhat1_3(3) - Bhat1_3(1)/Bhat1_3(3); % good data - fminsearch
%
y2_1 = -xl * Bhat2_1(2)/Bhat2_1(3) - Bhat2_1(1)/Bhat2_1(3); % bad data - logistic
y2_2 = -xl * Bhat2_2(2)/Bhat2_2(3) - Bhat2_2(1)/Bhat2_2(3); % bad data - glmfit
y2_3 = -xl * Bhat2_3(2)/Bhat2_3(3) - Bhat2_3(1)/Bhat2_3(3); % bad data - fminsearch
%
y3_3 = -xl * Bhat3_3(2)/Bhat3_3(3) - Bhat3_3(1)/Bhat3_3(3); % bad data - fminsearch
plot(xl,y1_2,'b', xl,y2_2,'b','LineWidth',6);
plot(xl,y1_3,'c', xl,y2_3,'c','LineWidth',4);
plot(xl,y1_1,'r', xl,y2_1,'r','LineWidth',2);
plot(xl,y3_3,'k--', 'LineWidth',4);
axis([xl yl]);
hold off

Yhat(:,1) = 1 ./ (1 + exp(-X * Bhat1_3 ));  
Yhat(:,2) = 1 ./ (1 + exp(-X * Bhat2_3 ));  
Yhat(:,3) = 1 ./ (1 + exp(-X * Bhat3_3 )); 
spec = sum(Yhat < .5 & repmat(Y1,1,3)==0,1)./sum(Yhat<.5,1);
sens = sum(Yhat >= .5 & repmat(Y1,1,3)==1,1)./sum(Yhat>.5,1);
fprintf('\nSPECIFICITY AND SENSITIVITY FOR THE THREE METHODS\n')
fprintf(' %.3f ',spec); fprintf('\n');
fprintf(' %.3f ',sens); fprintf('\n');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END EASY EXAMPLE BEGIN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

monte_carlo_its = 10;
spec_mat = zeros( monte_carlo_its,5);
sens_mat = zeros( monte_carlo_its,5);



% % Image size is 256x256
% 
% % truth is a circle (maybe try a ring later)
var = 50;
u = 0;
x = -127:1:128;
y = 1/sqrt(var * 2*pi) * exp(-(x-u).^2/(2*var));
img = y'*y;
thresh = img(75,75);
ring_thresh = img(90,100);
circ_img = img>thresh;
ring_img = img>thresh & img<ring_thresh;
figure; imagesc(circ_img); axis equal
figure; imagesc(ring_img); axis equal

truth = ring_img(:); % vectorized

for i=1:monte_carlo_its
    
% create fake distributions - RGB - in the 0/255 scale
Background_Sigma = [50 0 50; ...
                    0 200 0; ...
                    50 0 100];  
Foreground_Sigma = [150 0 50; ...
                    0 150 0;...
                    50 0 80];
Background_Mu = [150; 130; 250];
Foreground_Mu = [125; 150; 200];

Background_data = min(255, max(0, mvnrnd(Background_Mu,Background_Sigma,sum(truth==0)))) / 255; % data is ensured to be 0-255
Foreground_data = min(255, max(0, mvnrnd(Foreground_Mu,Foreground_Sigma,sum(truth==1)))) / 255;


X = zeros(length(truth), 4);
X(:,1) = 1;
X(truth==0,2:4) = Background_data;
X(truth==1,2:4) = Foreground_data; 

% VISUALIZE
% color_img = reshape(X(:,2:4), [size(circ_img) 3]);
% figure; image(color_img); axis equal;
% figure; imshow(rgb2gray(color_img)); axis equal;

% GOOD DATA - Y1
Y1 = truth;

% DATA WITH HIGH LEVEL OF FALSE POSITIVES = POORLY LABELED DATA
% original guassian was thresholded at the (75,75) point
% call the label boundary thresholded at 10-50 pixels to the left of that point; 
% labeled = img > img(75,75-floor(10+50*rand)); FOR CIRCLE IMAGE
labeled = (img > img(75,75-floor(10+45*rand)) ) & img < img(90,100+floor(8+45*rand)); % FOR RING IMAGE
% imagesc(ring_img + labeled); axis equal; pause(2); 
Y2 = labeled(:);

%using glmfit - cant feed this thing the whole X or it flips out (no 1's column)
Bhat1_2 = glmfit(X(:,2:4),Y1,'binomial','link','logit');  Yhat1_2 = 1 ./(1+exp(-X*Bhat1_2));
Bhat2_2 = glmfit(X(:,2:4),Y2,'binomial','link','logit');  Yhat2_2 = 1 ./(1+exp(-X*Bhat2_2));

%using fminsearch
[Bhat1_3, ~] = fminsearch(@loglikelihood, [0 0 0 0]', optimset('TolFun', 0.0001), Y1, X, 0,[],[]);  Yhat1_3 = 1 ./(1+exp(-X*Bhat1_3));
[Bhat2_3, ~] = fminsearch(@loglikelihood, [0 0 0 0]', optimset('TolFun', 0.0001), Y2, X, 0, [],[]); Yhat2_3 = 1 ./(1+exp(-X*Bhat2_3));

% ROBUST METHOD
fprintf('\nROBUST ITERATIVE METHOD'); 
g00_init = .8; g00 = g00_init; old_g00= 0;  % specificity
g11_init = .8; g11 = g11_init;              % sensitivity  (approx 1)
max_its = 25; iter = 0;

[Bhat3_3, ~] = fminsearch(@loglikelihood, [0 0 0 0]', optimset('TolFun', 0.0001), Y2, X, 1, g00_init, g11_init); % once and then iterate
while (abs(old_g00 - g00) > .0005 && iter < max_its)  
    old_g00 = g00;
    iter = iter+1;
    fprintf('\n %u. g00 = %.3f \t g11 = %.3f ', iter, g00, g11);
    Yhat = 1./(1+exp(-X*Bhat3_3));
    g00 = sum(Yhat < .5 & Y2==0)/sum(Yhat<.5); g01 = 1-g00;
    g11 = sum(Yhat >= .5 & Y2==1)/sum(Yhat>.5); g10 = 1-g11; 
    [Bhat3_3, ~] = fminsearch(@loglikelihood, Bhat3_3, optimset('TolFun', 0.0001), Y2, X, 1, g00, g11); % once and then iterate
end
Yhat3_3 = 1 ./(1+exp(-X*Bhat3_3));
fprintf('\n');

% % print results
% fprintf('\nUSING GOOD DATA');
% fprintf('\nBhat using "glmfit"       = '); fprintf('%4.2f   ',Bhat1_2); fprintf('\tB2/B3 = %4.2f, B1/B3 = %4.2f', Bhat1_2(2)/Bhat1_2(3),Bhat1_2(1)/Bhat1_2(3));
% fprintf('\nBhat using "fminsearch"   = '); fprintf('%4.2f   ',Bhat1_3); fprintf('\tB2/B3 = %4.2f, B1/B3 = %4.2f', Bhat1_3(2)/Bhat1_3(3),Bhat1_3(1)/Bhat1_3(3));
% fprintf('\n');
% % 
% fprintf('\nUSING DATA WITH HIGH FALSE POS RATE');
% fprintf('\nBhat using "glmfit"       = '); fprintf('%4.2f   ',Bhat2_2); fprintf('\tB2/B3 = %4.2f, B1/B3 = %4.2f', Bhat2_2(2)/Bhat2_2(3),Bhat2_2(1)/Bhat2_2(3));
% fprintf('\nBhat using "fminsearch"   = '); fprintf('%4.2f   ',Bhat2_3); fprintf('\tB2/B3 = %4.2f, B1/B3 = %4.2f', Bhat2_3(2)/Bhat2_3(3),Bhat2_3(1)/Bhat2_3(3));
% fprintf('\n');
% % 
% fprintf('\nUSING DATA WITH HIGH FALSE POS RATE WITH ROBUSTNESS ACTIVATED');
% fprintf('\nBhat using "fminsearch"   = '); fprintf('%4.2f   ',Bhat3_3); fprintf('\tB2/B3 = %4.2f, B1/B3 = %4.2f', Bhat3_3(2)/Bhat3_3(3),Bhat3_3(1)/Bhat3_3(3));
% fprintf('\n');


Yhat_all = [Yhat1_2, Yhat1_3, Yhat2_2, Yhat2_3, Yhat3_3];
spec = sum(Yhat_all <  .5 & repmat(Y1,1,size(Yhat_all,2))==0,1)./sum(Y1==0);
sens = sum(Yhat_all >= .5 & repmat(Y1,1,size(Yhat_all,2))==1,1)./sum(Y1==1);

fprintf('\nSPECIFICITY AND SENSITIVITY FOR THE THREE METHODS\n')
fprintf(' %.3f ',spec); fprintf('\n');
fprintf(' %.3f ',sens); fprintf('\n');

spec_mat(i,:) = spec;
sens_mat(i,:) = sens;

figure; 
subplot(2,2,1); imagesc(ring_img + labeled); axis equal; colormap(gray); title('manually segmented area')
subplot(2,2,2); imagesc(reshape(Yhat1_3,size(circ_img))); axis equal; colormap(gray); title('LR trained on GOOD data')
subplot(2,2,3); imagesc(reshape(Yhat2_3,size(circ_img))); axis equal; colormap(gray); title('LR trained on BAD data')
subplot(2,2,4); imagesc(reshape(Yhat3_3,size(circ_img))); axis equal; colormap(gray); title('rLR trainted on BAD data')
pause(2);

end


big_mat = [spec_mat, sens_mat];
figure; boxplot(big_mat)




figure; boxplot(spec_mat(:,[1 3 5]));title('Specificity')
figure; boxplot(sens_mat(:,[1 3 5]));title('Sensitivity')













