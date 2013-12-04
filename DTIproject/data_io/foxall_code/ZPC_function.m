%
% FPC_FUNCTION - Return the minimization Q value For EPI reconstruction
%

function Q = FPC_function(ZPC,FPC,rfilter,pimage,nimage,x_matrix,y_matrix)

pic = PC_image(ZPC,FPC,pimage,nimage,x_matrix,y_matrix);
Q   = sum_Q(pic.*rfilter);