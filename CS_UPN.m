%% This program is used for image denoising
% using compressed sensing (CS) method

%-------------Description------------------

% solve the following convex optimization problem

% min f(x)= min 0.5*||A*x-b||^2+alpha*||x||_TV

% A =system or projection matrix
% b= original undersampled projection data
% x= reconstructed image
% ||=l2-norm, ||_TV= total-variation (TV) norm
% alpha = control parameter 

% The minimization problem is solved using unknown parameter Nesterov (UPN)

% Author: Sovan Mukherjee, April 2015

% Ref: 
% Chen et al.,"Prior image constrained compressed sensing (PICCS):
% A method to accurately reconstruct dynamic CT images from highly
% undersampled projection data sets," Phys. Med. Biol.,35(2),
% 660-661(2008).

% Jensen et al., "Implementation of an optimal first-order method for 
%strongly convex total variation regularization," BIT Numer. Methods, 51,
% 1–28 (2011).

% Niu et al.,"Accelarated barrier optimization compressed sensing (ABOCS)
% for CT reconstruction with improved convergence," Phys. Med. Biol.,59,
% 1801-1814 (2014).


clear all;
close all;
clc;

%% Shepp-Logan Phantom
                                                        
Im= phantom(512);                                                          % generate Shepp-Logan phantom 

figure; imshow(Im); title('shepp-logan image');                            % show the original phantom

%get the projection for initial image

angle_proj= (1:4:360);                                                     % angle of projection

Proj=myRadon(Im,angle_proj);                                               % get the projection

figure;

imagesc(Proj);                                                             % show the projection for initial image

colormap('gray'); title('projection data');

% back-projection

x0= myFilteredBackprojectionSpatialDomain...
    (Proj,angle_proj);                                                     % back projection (FBP)

% x0= iradon...
%     (Proj(:,1:length(angle_proj)),angle_proj);                           % back projection

figure; imshow(mat2gray(x0));                                                        % show the projection image

n_d= size(x0,1);                                                           % number of detectors

x0=x0(:);                                                                  % resize x0

b=Proj(:);                                                                 % resize Projections                                                           % vectored projection data

%% Calculate the gradient operator

N=size(Im,1); 

[D1,D2]= grad_operator_new(N);

D= [D1;D2];

%% -- Get System or Projection Matrix -----------------------------------

A=get_system_matrix(n_d,angle_proj);

%% Unknown Parameter Nesterov's (UPN) algorithm

alpha=1;                                                                   % Regularization parameter

% alpha=0.1;                                                               % Regularization parameter

nu=0.01;

% UPN starts with backtracking line search (BT)

y = x0;                                                                    % initialized image

mu= 20;                                                                    % strong convexity parameter

L0= 1e3;                                                                   % Lipschitz constant

pl=1.3;

Nmax= 100;                                                                 % Maximum iterations (can be changed)

x1=x0;

mu0=mu;

theta0=sqrt(mu0/L0);

tic;

% error= zeros(1,Nmax);

for k= 1:Nmax
    
    fprintf('Iteration=%d\t',k)
    
    [x, L0]=BT(A,b,y,D1,D2,L0,alpha,pl,N);
    
    % calculate the objective functions

    fy = obj_func(D1,D2,y,A,b,alpha);                                      % f(y)
    fx1 = obj_func(D1,D2,x1,A,b,alpha);                                    % f(x1)
    
    % calculate the gradient of the objective function
    
    grad_fy= Grad_fy(A,b,y,alpha,D1,D2,N);
        
    mu0=min(mu0,(fx1-fy-grad_fy'*(x1-y))/(0.5*norm(x1-y)^2));              % update mu
    
    
    thetaNew= 0.5*((mu0/L0)-theta0^2+sqrt(((mu0/L0)-theta0^2)^2 ...        % update theta
   +4*theta0^2));
    
    beta= theta0*(1-theta0)/(theta0^2+thetaNew);                           % update beta
    
    if beta>1
        
        break;
        
    end
    
    beta(isnan(beta))=0;
    
    theta0=thetaNew;                                                       % update theta
    
    y= x+beta*(x-x1);                                                      % update y
    
    x1=x;
    
    grad_fx= Grad_fy(A,b,x,alpha,D1,D2,N);
    
    G_nu_x= nu*(x-max(0,x-(1/nu)*grad_fx));
    
    error= norm(G_nu_x)/(N*N);
    
    Error(:,k)=error;
    
    fprintf('Error=%d\t',error);
    
    fprintf('Beta=%d\n',beta);
    
    % stopping criterion
    
    if norm(G_nu_x)/(N*N)<5e-6
        
        break;
        
    end
    
end

%% UPN ends here

%% show the reconstructed image

x= reshape(x,N,N);
x=imrotate(x,270);
% 
figure;imshow(mat2gray(x)); title('Reconstructed image');
    
toc;


