function SMOO_IM = vuDWIAnisotropicDiffusion(X, varargin)
% vuDWIAnisotropicDiffusion(...) performs anisotropic diffusion on DWI image X,
% using methods presented by Xu, Ding, et. al.
%
%   SYNTAX:
%       SMOO_IM = vuDWIAnisotropicDiffusion(X)
%       SMOO_IM = vuDWIAnisotropicDiffusion(X, options)       
%
%   OPTIONS:
%       sigma     : standard deviation of noise in the image
%       multiTimeStep : break up the time step into several iterations resulting
%                       is a slightly better smoothing, at a computational cost.
%
%   OPTIONS & DEFAULTS:
%       sigma = 0.05;
%       multiTimeStep = false;
%
%   OUTPUTS:
%       SMOO_IM is the smoothed filtered image
%
%   EXAMPLE:
%       im = vuOpenImage('DWI_Image');
%       smooth_im = vuDWIAnisotropicDiffusion(im);
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

if nargin < 1
  error('MATLAB:vuDWIAnisotropicDiffusion:NotEnoughInputs', 'Not enough input arguments.');
end

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
    disp('X is a valid MetaImage struct.')
  else
      error('MATLAB:vuDWIAnisotropicDiffusion:InvalidStruct', 'The input image structure is not valid.');
  end
  isStruct = true;
else
  X = vuGenerateMetaImage(single(X),ones(1,3),zeros(1,3));
  disp('X is a matrix.')
  isStruct = false;
end

if (length(X.Dims) ~= 4)
  error('MATLAB:vuDWIAnisotropicDiffusion:UnknownDims', 'vuDWIAnisotropicDiffusion can only handle images of 4 dimensions.');
end

% Get and parse options
p = inputParser;
p.addParamValue('sigma',0.05,@(x) isa(x,'double'));
p.addParamValue('multiTimeStep',false,@(x) isa(x,'double')||isa(x,'logical'));
p.FunctionName='vuDWIAnisotropicDiffusion';
p.parse(varargin{:});


% other variable dependent
rho = p.Results.sigma*2;
if(p.Results.multiTimeStep)
    iterations = ceil(p.Results.sigma/0.05);
    remainder = mod(p.Results.sigma,0.05)*100;
    timeStep(1:iterations-1)=(3/44)*5;
    timeStep(iterations) = (3/44)*remainder;
else
    iterations = 1;
    timeStep = (3/44)*p.Results.sigma*100;
end
% Left in other code just in case!
type = 'implicit_ADI';
isFastEig = true;

if ~(strcmp(type,'explicit')||strcmp(type,'implicit_multisplitting')||strcmp(type,'implicit_AOS')||strcmp(type,'implicit_ADI'))
    error('MATLAB:vuDWIAnisotropicDiffusion:TypeMismatch','The defined numerical scheme is not reconized');
end

% Call Smoothing Function
if (isStruct)
    SMOO_IM.Data = X.Data;
    % Free some memory
    X.Data = [];
    for i = 1:iterations
        SMOO_IM.Data = aniso4D_smoothing_oneiteration(SMOO_IM.Data,p.Results.sigma,rho,timeStep(i),X.Spc,type,isFastEig);
    end
    SMOO_IM.Dims = X.Dims;
    SMOO_IM.Spc = X.Spc;
    SMOO_IM.Origin = X.Origin;
    % Copy Orientation if it exists
    if (isfield(X,'Orientation'))
        REG_IM.Orientation = X.Orientation;
    end
    % Copy parmeters if they exist
    if (isfield(X,'Parms'))
        REG_IM.Parms = X.Parms;
    end
else
    SMOO_IM = X.Data;
    % Free some memory
    X.Data = [];
    for i = 1:iterations
        SMOO_IM = aniso4D_smoothing_oneiteration(SMOO_IM,p.Results.sigma,rho,timeStep(i),X.Spc,type,isFastEig);
    end
end

return;

function img_U_s = aniso4D_smoothing_oneiteration(img_U,sigma,rho,delta_t,sr,type,isfasteig)
% DWI Anisotropic Diffusion function
% Original Code by: Qing Xu, Zhaohua Ding
% Adapted for vuTools May 16, Kevin Wilson
%
% Copyright (c) 2007 - Vanderbilt University Institute of Imaging Science

[img_height, img_width, img_slice, img_type] = size(img_U);
% Smooth image iteratively
Gx = zeros(img_height, img_width, img_slice);
Gy = zeros(img_height, img_width, img_slice);
Gz = zeros(img_height, img_width, img_slice);
J_p = zeros(img_height, img_width, img_slice, 3, 3);
NeVal = zeros(3, 1);
G_sigma = GaussianKernel([1 1 1],[sigma sigma sigma],sr);
G_rho   = GaussianKernel([1 1 1],[rho rho rho],sr);

% Convolve G_sigma with raw image
for t=1:img_type
    U_sigma(:, :, :) = convn(img_U(:,:,:,t), G_sigma, 'same');

    % Get gradient info
    [Gx(:, :, :), Gy(:, :, :), Gz(:, :, :)] = gradient(U_sigma(:,:,:), sr(1), sr(2), sr(3));
    % Construct structure tensor
    J_p(:, :, :, 1, 1) = J_p(:, :, :, 1, 1) + Gx(:, :, :).^2;
    J_p(:, :, :, 1, 2) = J_p(:, :, :, 1, 2) + Gx(:, :, :).*Gy(:, :, :);
    J_p(:, :, :, 1, 3) = J_p(:, :, :, 1, 3) + Gx(:, :, :).*Gz(:, :, :);
    J_p(:, :, :, 2, 2) = J_p(:, :, :, 2, 2) + Gy(:, :, :).^2;
    J_p(:, :, :, 2, 3) = J_p(:, :, :, 2, 3) + Gy(:, :, :).*Gz(:, :, :);
    J_p(:, :, :, 3, 3) = J_p(:, :, :, 3, 3) + Gz(:, :, :).^2;
end
J_p(:,:,:,2,1) = J_p(:,:,:,1,2);
J_p(:,:,:,3,1) = J_p(:,:,:,1,3);
J_p(:,:,:,3,2) = J_p(:,:,:,2,3);

clear Gx Gy Gz U_sigma

% Smooth structure tensor with G_rho
for m=1:3
    for n=1:3
        J_p(:, :, :, m, n) =  convn(J_p(:, :, :, m, n), G_rho, 'same');
    end
end

D = zeros(img_height, img_width, img_slice, 3, 3);
[x, y, z] = meshgrid((1:img_width)*sr(1), (1:img_height)*sr(2), (1:img_slice)*sr(3));

if isfasteig 
[eVec,eVal,mask] = eigss_fast(J_p);
idx = find(mask==0);
for m = idx'
    [iy,ix,iz] = ind2sub(size(mask),m);
    [eVec(iy,ix,iz,:,:),eValTmp] = eigs(squeeze(J_p(iy, ix, iz, 1:3, 1:3)), 3);
    eVal(iy,ix,iz,:) = diag(eValTmp);
end


mask1 =  eVal(:,:,:,1)<0.1^30;
mask2 = (eVal(:,:,:,2)<0.1^30)&(~mask1);
mask3 = (eVal(:,:,:,3)<0.1^30)&(~mask2)&(~mask1);
mask4 = (~mask1)&(~mask2)&(~mask3);
idx1 = find(mask1==1);
idx2 = find(mask2==1);
idx3 = find(mask3==1);
idx4 = find(mask4==1);
eVal = reshape(eVal,[img_height*img_width*img_slice 3]);
eVal(idx1,:) = 1;
eVal(idx2,1) = 0;
eVal(idx2,2:3) = 1.5;
eVal(idx3,1:2) = 0;
eVal(idx3,3) = 3;
eVal(idx4,:) = 1./eVal(idx4,:);
eVal(idx4,:) = 3*eVal(idx4,:)./repmat(sum(eVal(idx4,:),2),[1 3]);
eVal = reshape(eVal,[img_height img_width img_slice 3]);

eVal = reshape(eVal,[img_height img_width img_slice 1 3]);
eVal = repmat(eVal,[1 1 1 3 1]);
S = eVec.*eVal;
D(:,:,:,1,1) = squeeze(sum(S(:,:,:,1,:).*eVec(:,:,:,1,:),5));
D(:,:,:,1,2) = squeeze(sum(S(:,:,:,1,:).*eVec(:,:,:,2,:),5));
D(:,:,:,1,3) = squeeze(sum(S(:,:,:,1,:).*eVec(:,:,:,3,:),5));
D(:,:,:,2,2) = squeeze(sum(S(:,:,:,2,:).*eVec(:,:,:,2,:),5));
D(:,:,:,2,3) = squeeze(sum(S(:,:,:,2,:).*eVec(:,:,:,3,:),5));
D(:,:,:,3,3) = squeeze(sum(S(:,:,:,3,:).*eVec(:,:,:,3,:),5));
D(:,:,:,2,1) = D(:,:,:,1,2);
D(:,:,:,3,1) = D(:,:,:,1,3);
D(:,:,:,3,2) = D(:,:,:,2,3);
else
% construct new structure tensors used in the smoothing
for i=1:img_height
    for j=1:img_width
        for k=1:img_slice
            no = 3;
            
            [eVec, eVal] = eigs(squeeze(J_p(i, j, k, 1:3, 1:3)), 3);
            eVal = diag(eVal);
            if  eVal(1) == 0
                NeVal(1:3) = [no/3 no/3 no/3];
            elseif  eVal(2) == 0
                NeVal(1:3) = [0 no/2 no/2];
            elseif  eVal(3) == 0
                NeVal(1:3) = [0 0 no];
            else
                NeVal(1:3) = 1./eVal;
                NeVal = no*NeVal/sum(NeVal);
            end
            D(i, j, k, 1:3, 1:3) = eVec*[NeVal(1) 0 0; 0 NeVal(2) 0; 0 0 NeVal(3)]*eVec';
        end
    end
end
    
end

%explicit scheme
switch type 
    case 'explicit'
        for t=1:img_type
            U_sigma(:, :, :) = convn(img_U(:,:,:,t), G_sigma, 'same');
            [Gx(:, :, :), Gy(:, :, :), Gz(:, :, :)] = gradient(U_sigma(:,:,:), sr(1), sr(2), sr(3));
            img_U(:,:,:,t) = img_U(:,:,:,t) + ...
                delta_t*divergence(x,y,z, ...
                squeeze(D(:,:,:,1,1)).*Gx(:,:,:) + squeeze(D(:,:,:,1,2)).*Gy(:,:,:) + squeeze(D(:,:,:,1,3)).*Gz(:,:,:),...
                squeeze(D(:,:,:,2,1)).*Gx(:,:,:) + squeeze(D(:,:,:,2,2)).*Gy(:,:,:) + squeeze(D(:,:,:,2,3)).*Gz(:,:,:),...
                squeeze(D(:,:,:,3,1)).*Gx(:,:,:) + squeeze(D(:,:,:,3,2)).*Gy(:,:,:) + squeeze(D(:,:,:,3,3)).*Gz(:,:,:));
        end
    
    case 'implicit_multisplitting'
        img_U_tmp = img_U;
        for t=1:img_type
            U_sigma(:, :, :) = convn(img_U(:,:,:,t), G_sigma, 'same');
            [Gx(:, :, :), Gy(:, :, :), Gz(:, :, :)] = gradient(U_sigma(:,:,:), sr(1), sr(2), sr(3));
            %explicit
            img_U(:,:,:,t) = img_U(:,:,:,t) + ...
                delta_t*divergence(x,y,z, ...
                squeeze(D(:,:,:,1,2)).*Gy(:,:,:) + squeeze(D(:,:,:,1,3)).*Gz(:,:,:),...
                squeeze(D(:,:,:,2,1)).*Gx(:,:,:)                                      + squeeze(D(:,:,:,2,3)).*Gz(:,:,:),...
                squeeze(D(:,:,:,3,1)).*Gx(:,:,:) + squeeze(D(:,:,:,3,2)).*Gy(:,:,:)                                    );
            %implicit AOS
            img_U(:,:,:,t) = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,1,1)),delta_t,1,2);
            img_U(:,:,:,t) = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,2,2)),delta_t,1,1);
            img_U(:,:,:,t) = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,3,3)),delta_t,1,3);
        end

    case 'implicit_AOS'
        img_U_tmp = img_U;
        for t=1:img_type
            U_sigma(:, :, :) = convn(img_U(:,:,:,t), G_sigma, 'same');
            [Gx(:, :, :), Gy(:, :, :), Gz(:, :, :)] = gradient(U_sigma(:,:,:), sr(1), sr(2), sr(3));
            %explicit
            img_U(:,:,:,t) = img_U(:,:,:,t) + ...
                delta_t*divergence(x,y,z, ...
                squeeze(D(:,:,:,1,2)).*Gy(:,:,:) + squeeze(D(:,:,:,1,3)).*Gz(:,:,:),...
                squeeze(D(:,:,:,2,1)).*Gx(:,:,:)                                      + squeeze(D(:,:,:,2,3)).*Gz(:,:,:),...
                squeeze(D(:,:,:,3,1)).*Gx(:,:,:) + squeeze(D(:,:,:,3,2)).*Gy(:,:,:)                                    );
            %implicit AOS
            img_Uz = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,3,3)),3*delta_t,1,3);
            img_Uy = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,2,2)),3*delta_t,1,1);
            img_Ux = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,1,1)),3*delta_t,1,2);
            img_U(:,:,:,t) = (img_Ux + img_Uy + img_Uz)/3;
        end
    
    case 'implicit_ADI'
        for t=1:img_type 
            img_U(:,:,:,t) = ParobolicEqnSolver3dOnce_ADI(squeeze(img_U(:,:,:,t)),D,delta_t,1,1,1);
        end
end
img_U_s = img_U;

return;

function H = GaussianKernel(sz,sigma,sr)
% This function generates a Gaussian kernel
% Input: 
%       sz : [nx, ny, nz] the size of the kernel is 2nx + 1; 2ny + 1...
%       sigma : the sigma for each dimenstion, note that it could be
%       vector
% Output:
%       H : the kernel

if ~(length(sz) == length(sigma))
    error('the size of sigma and sz does not match');
end
H = ones(2*sz+1);
for i = 1 : length(sz)
    x = (-sz(i):sz(i))*sr(i);
    kernel_1d = (1/(2*pi*sigma(i)^2)^(1/2))*exp(-x.^2/(2*sigma(i)^2));
    kernel_1d_col = kernel_1d';
    sz_col = (2*sz+1);
    sz_col = circshift(sz_col,[1 i-1]);
    kernel_3d_outorder = repmat(kernel_1d_col,[1 sz_col(2:end)]);
    kernel_3d_inorder = shiftdim(kernel_3d_outorder,length(sz)+1-i);
    H = H.*kernel_3d_inorder;
end
H = H/sum(H(:));

function U_next = SolveTriangleTomasBatch(U,G,dt,h,dim)
% Do the triangular batching job
sz = size(U);
U = shiftdim(U,dim-1);
G = shiftdim(G,dim-1);
sz = size(U);
U = reshape(U,[sz(1) sz(2)*sz(3)]);
G = reshape(G,[sz(1) sz(2)*sz(3)]);
U = U';
G = G';
for i=1:sz(2)*sz(3)
     U(i,:) = SolveTriangleTomas(squeeze(U(i,:)),squeeze(G(i,:)),dt,h);
end
U = U';
G = G';
U = reshape(U,[sz(1) sz(2) sz(3)]);
G = reshape(G,[sz(1) sz(2) sz(3)]);
U = shiftdim(U,4-dim);
G = shiftdim(G,4-dim);
U_next = U;
return

function u_next = SolveTriangleTomas(u,g,dt,h)
% Solve (I-dt*A(g))*u_next = u,where u is just a vector with boundary
% condition on two ends

sz = size(u);
d = u;
g = g*dt;
%%%%%%%%%%% Solve the linear system using Tomas %%%%%%%%%%%%%%%%%%%%%%
alpha2 = ((g(1)+g(2))/(2*h^2))/(1+(g(1)+g(2))/(2*h^2));
alphan = ((g(end)+g(end-1))/(2*h^2))/(1+(g(end)+g(end-1))/(2*h^2));
beta2 = d(1)/(1+(g(1)+g(2))/(2*h^2));
betan = d(end)/(1+(g(end)+g(end-1))/(2*h^2));

n = length(d);
% Solving the p and q
p = alpha2;
q = beta2;
for i=2:n-1
    a = (g(i-1)+g(i))/(2*h^2);
    c = (g(i)+g(i+1))/(2*h^2);
    b = -(a+c);
    a = -a;
    c = -c;
    b = 1-b;
    p_next = -c/(a*p(i-1)+b);
    q_next = (d(i)-a*q(i-1))/(a*p(i-1)+b);
    p = [p p_next];
    q = [q q_next];
end

% Solving the w
w = (alphan*q(end)+betan)/(1-alphan*p(end));

for i=n-1:-1:1
    w_prev = w(1)*p(i)+q(i);
    w = [w_prev w];
end

u_next = w;

return;

function I_next = ParobolicEqnSolver3dOnce_ADI(I,D,dt,dx,dy,dz)
% Solve 2d parabolic PDE It = D(:,:,:,1,1)Ixx + D(:,:,:,2,2)Iyy +
% D(:,:,:,3,3)Izz + 2*D(:,:,:,1,2)Ixy + 2*D(:,:,:,1,3)Ixz +
% 2*D(:,:,:,2,3)Iyz
%
% Note that reflecting boundary condition is used here, so there is no b.c.
% input
%
% Input: I--3d data on the current time step
%       D--Diffusion coefficients, which is a matrix 
%       dt,dx,dy,dz--the time and space step
%
% Output: I_next-- the 3d data on next time step
sz =     size(I);
[X, Y, Z] = meshgrid((1:sz(2))*dy, (1:sz(1))*dx, (1:sz(3))*dz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equation:(I-dt/2A11)(I-dt/2A22)(I-dt/2A33)Uk+1=(I+dt/2A11)(I+dt/2A22)(I+dt/2A33)Uk + explicit part 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Gx,Gy,Gz] = gradient(I);
R1 = dt*divergence(X,Y,Z, ...
            squeeze(D(:,:,:,1,2)).*Gy + squeeze(D(:,:,:,1,3)).*Gz,...
            squeeze(D(:,:,:,2,1)).*Gx + squeeze(D(:,:,:,2,3)).*Gz,...
            squeeze(D(:,:,:,3,1)).*Gx + squeeze(D(:,:,:,3,2)).*Gy);
R2 = I;

[junk,junk,Gz] = gradient(R2);
[junk,junk,Gzz] = gradient(squeeze(D(:,:,:,3,3)).*Gz);
R2 = R2 + dt/2*Gzz;

[junk,Gy,junk] = gradient(R2);
[junk,Gyy,junk] = gradient(squeeze(D(:,:,:,2,2)).*Gy);
R2 = R2 + dt/2*Gyy;

[Gx,junk,junk] = gradient(R2);
[Gxx,junk,junk] = gradient(squeeze(D(:,:,:,1,1)).*Gx);
R2 = R2 + dt/2*Gxx;

R = R1 + R2;

I_next = SolveTriangleTomasBatch(R,squeeze(D(:,:,:,1,1)),dt/2,1,2);
I_next = SolveTriangleTomasBatch(I_next,squeeze(D(:,:,:,2,2)),dt/2,1,1);
I_next = SolveTriangleTomasBatch(I_next,squeeze(D(:,:,:,3,3)),dt/2,1,3);

R3 = I_next-I;
[Gx,Gy,Gz] = gradient(R3);
R3 = .5*dt*divergence(X,Y,Z, ...
            squeeze(D(:,:,:,1,2)).*Gy + squeeze(D(:,:,:,1,3)).*Gz,...
            squeeze(D(:,:,:,2,1)).*Gx + squeeze(D(:,:,:,2,3)).*Gz,...
            squeeze(D(:,:,:,3,1)).*Gx + squeeze(D(:,:,:,3,2)).*Gy);
R = R1 + R2 + R3;
I_next = SolveTriangleTomasBatch(R,squeeze(D(:,:,:,1,1)),dt/2,1,2);
I_next = SolveTriangleTomasBatch(I_next,squeeze(D(:,:,:,2,2)),dt/2,1,1);
I_next = SolveTriangleTomasBatch(I_next,squeeze(D(:,:,:,3,3)),dt/2,1,3);

return;

function [evector,evalue,mask] = eigss_fast(D)
% Optimization based on Matlab analytical eigenvalue&eigenvector solution "Analytical Computation of the Eigenvalues and Eigenvectors in DT-MRI"
%
% Input: D 3Dx3x3 matrix, if you wanna input a 2D image,please input it as a
% 3D array.
%
% Output: evector: the column vectors are eigenvectors
%        evalue:  a row vector, eigenvalue
sz = size(D);
sx = sz(length(sz)-1);
sy = sz(length(sz));
fakezero = 0.1^30;

if (sx==3)
    if (sy==3)
    % I did not test symmetr here,because I do not wanna introduce any
    % additoanl computation.

        % Determination of eigenvalue
        I1 = D(:,:,:,1,1)+D(:,:,:,2,2)+D(:,:,:,3,3);
        I2 = (D(:,:,:,1,1).*D(:,:,:,2,2)+D(:,:,:,1,1).*D(:,:,:,3,3)+D(:,:,:,2,2).*D(:,:,:,3,3))-(D(:,:,:,1,2).^2+D(:,:,:,1,3).^2+D(:,:,:,2,3).^2);
        I3 = D(:,:,:,1,1).*D(:,:,:,2,2).*D(:,:,:,3,3)+D(:,:,:,1,2).*D(:,:,:,2,3).*D(:,:,:,3,1)+D(:,:,:,1,3).*D(:,:,:,2,1).*D(:,:,:,3,2)...
            -D(:,:,:,1,3).*D(:,:,:,2,2).*D(:,:,:,3,1)-D(:,:,:,1,2).*D(:,:,:,2,1).*D(:,:,:,3,3)-D(:,:,:,1,1).*D(:,:,:,2,3).*D(:,:,:,3,2);

        v = (I1/3).^2-I2/3+fakezero;
        s = (I1/3).^3-I1.*I2/6+I3/2;
        sita = acos(s./v.*sqrt(1./v))/3;
        evalue = zeros([sz(1:3) 3]);
        evalue(:,:,:,1) = I1/3+2*sqrt(v).*cos(sita);
        evalue(:,:,:,2) = I1/3-2*sqrt(v).*cos(pi/3*ones(sz(1:3))+sita);
        evalue(:,:,:,3) = I1-evalue(:,:,:,1)-evalue(:,:,:,2);
        % Build the mask for degenerate and non-pos case
        mask = ones(sz(1:3)); 
        mask = mask.*(I3>fakezero).*(squeeze(D(:,:,:,1,1))>=fakezero).*(squeeze(D(:,:,:,2,2))>=fakezero).*(squeeze(D(:,:,:,3,3))>=fakezero)...
            .*((squeeze(D(:,:,:,1,1)).*squeeze(D(:,:,:,2,2))-squeeze(D(:,:,:,1,2)).^2)>=fakezero)...
            .*((squeeze(D(:,:,:,1,1)).*squeeze(D(:,:,:,3,3))-squeeze(D(:,:,:,1,3)).^2)>=fakezero)...
            .*((squeeze(D(:,:,:,2,2)).*squeeze(D(:,:,:,3,3))-squeeze(D(:,:,:,2,3)).^2)>=fakezero)...
            .*~((abs(squeeze(D(:,:,:,1,3)))<0.1^19).*(abs(squeeze(D(:,:,:,2,3)))<0.1^19))...
            .*~((abs(squeeze(D(:,:,:,1,2)))<0.1^19).*(abs(squeeze(D(:,:,:,1,3)))<0.1^19))...
            .*~((abs(squeeze(D(:,:,:,1,2)))<0.1^19).*(abs(squeeze(D(:,:,:,2,3)))<0.1^19))...
            .*(abs(s./v.*sqrt(1./v))<=1).*(v>0);

        % Is it probabatily a degenrate case
        degentype = abs(evalue(:,:,:,1)-evalue(:,:,:,2))<0.1^5;
        % Determine the eigenvector
        A = repmat(D(:,:,:,1,1),[1 1 1 2]) - evalue(:,:,:,1:2);
        B = repmat(D(:,:,:,2,2),[1 1 1 2]) - evalue(:,:,:,1:2);
        C = repmat(D(:,:,:,3,3),[1 1 1 2]) - evalue(:,:,:,1:2);
        evector = zeros(sz);
        evector(:,:,:,1,1:2) =((repmat(D(:,:,:,1,2).*D(:,:,:,2,3),[1 1 1 2])-B.*repmat(D(:,:,:,1,3),[1 1 1 2]))...
                .*(repmat(D(:,:,:,1,3).*D(:,:,:,2,3),[1 1 1 2])-C.*repmat(D(:,:,:,1,2),[1 1 1 2])));

        evector(:,:,:,2,1:2) =((repmat(D(:,:,:,1,3).*D(:,:,:,2,3),[1 1 1 2])-C.*repmat(D(:,:,:,1,2),[1 1 1 2]))...
                .*(repmat(D(:,:,:,1,3).*D(:,:,:,1,2),[1 1 1 2])-A.*repmat(D(:,:,:,2,3),[1 1 1 2])));

        evector(:,:,:,3,1:2) =((repmat(D(:,:,:,1,2).*D(:,:,:,2,3),[1 1 1 2])-B.*repmat(D(:,:,:,1,3),[1 1 1 2]))...
                .*(repmat(D(:,:,:,1,3).*D(:,:,:,1,2),[1 1 1 2])-A.*repmat(D(:,:,:,2,3),[1 1 1 2])));
        % Norm of the first 2 eigenvectors
        evector(:,:,:,:,1) = evector(:,:,:,:,1)./(repmat(sqrt(sum(evector(:,:,:,:,1).^2,4)),[1 1 1 3])+fakezero);
        evector(:,:,:,:,2) = evector(:,:,:,:,2)./(repmat(sqrt(sum(evector(:,:,:,:,2).^2,4)),[1 1 1 3])+fakezero);
        % Use the tensor product of the first 2 eigenvectors to generate the
        % 3rd eigenvector
        evector(:,:,:,1,3) = evector(:,:,:,2,1).*evector(:,:,:,3,2)-evector(:,:,:,2,2).*evector(:,:,:,3,1); 
        evector(:,:,:,2,3) = evector(:,:,:,3,1).*evector(:,:,:,1,2)-evector(:,:,:,3,2).*evector(:,:,:,1,1);
        evector(:,:,:,3,3) = evector(:,:,:,1,1).*evector(:,:,:,2,2)-evector(:,:,:,1,2).*evector(:,:,:,2,1);


        A = repmat(D(:,:,:,1,1),[1 1 1 2]) - evalue(:,:,:,[1 3]);
        B = repmat(D(:,:,:,2,2),[1 1 1 2]) - evalue(:,:,:,[1 3]);
        C = repmat(D(:,:,:,3,3),[1 1 1 2]) - evalue(:,:,:,[1 3]);

        evector1 = zeros(sz);
        evector1(:,:,:,1,[1 3]) =((repmat(D(:,:,:,1,2).*D(:,:,:,2,3),[1 1 1 2])-B.*repmat(D(:,:,:,1,3),[1 1 1 2]))...
                .*(repmat(D(:,:,:,1,3).*D(:,:,:,2,3),[1 1 1 2])-C.*repmat(D(:,:,:,1,2),[1 1 1 2])));    
        evector1(:,:,:,2,[1 3]) =((repmat(D(:,:,:,1,3).*D(:,:,:,2,3),[1 1 1 2])-C.*repmat(D(:,:,:,1,2),[1 1 1 2]))...
                .*(repmat(D(:,:,:,1,3).*D(:,:,:,1,2),[1 1 1 2])-A.*repmat(D(:,:,:,2,3),[1 1 1 2])));    
        evector1(:,:,:,3,[1 3]) =((repmat(D(:,:,:,1,2).*D(:,:,:,2,3),[1 1 1 2])-B.*repmat(D(:,:,:,1,3),[1 1 1 2]))...
                .*(repmat(D(:,:,:,1,3).*D(:,:,:,1,2),[1 1 1 2])-A.*repmat(D(:,:,:,2,3),[1 1 1 2])));

        evector1(:,:,:,:,1) = evector1(:,:,:,:,1)./(repmat(sqrt(sum(evector1(:,:,:,:,1).^2,4)),[1 1 1 3])+fakezero);
        evector1(:,:,:,:,3) = evector1(:,:,:,:,3)./(repmat(sqrt(sum(evector1(:,:,:,:,3).^2,4)),[1 1 1 3])+fakezero);


        evector1(:,:,:,1,2) = evector1(:,:,:,2,1).*evector1(:,:,:,3,3)-evector1(:,:,:,2,3).*evector1(:,:,:,3,1); 
        evector1(:,:,:,2,2) = evector1(:,:,:,3,1).*evector1(:,:,:,1,3)-evector1(:,:,:,3,3).*evector1(:,:,:,1,1);
        evector1(:,:,:,3,2) = evector1(:,:,:,1,1).*evector1(:,:,:,2,3)-evector1(:,:,:,1,3).*evector1(:,:,:,2,1);

        evector = evector.*repmat(~degentype,[1 1 1 3 3])+evector1.*repmat(degentype,[1 1 1 3 3]);

        return;

    end
end

return;