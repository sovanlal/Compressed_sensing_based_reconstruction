% backtracking line search to find Lipschitz constant L
% author: Sovan Mukherjee, April 2015

function [x1, L0]= BT(A,b,y,D1,D2,L0,alpha,pl,N)

% calculate the gradient of the objective function

grad_fy= Grad_fy(A,b,y,alpha,D1,D2,N);

x1= max(0,y-(1/L0)*grad_fy);

% calculate the objective functions

fy = obj_func(D1,D2,y,A,b,alpha);                                          % f(y)
                                       
fx1 = obj_func(D1,D2,x1,A,b,alpha);                                        % f(x1)

% backtracking line search (BT)

while fx1>fy+grad_fy'*(x1-y)+0.5*L0*norm(x1-y)^2
    
    L0=L0*pl;
    
    x1= max(0,y-(1/L0)*grad_fy);
    
    fx1 = obj_func(D1,D2,x1,A,b,alpha); 
end