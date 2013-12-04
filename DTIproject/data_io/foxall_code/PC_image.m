%
% PC_IMAGE - Create Phase Corrected EPI Image From Positive Echo and Negative Echo Sub matrices
%

function data = PC_image(Z,F,pimage,nimage,x_matrix,y_matrix)

x          = 1:x_matrix;
pc         = Z + (F/(2*x_matrix))*(x - x_matrix/2 -0.5);
cexp       = exp(i*pc);

for Y = 1:y_matrix
    pic(Y,:) = pimage(Y,:) + nimage(Y,:).*cexp;
end

data = pic;