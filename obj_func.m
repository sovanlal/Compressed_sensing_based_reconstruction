% this function calculates the objective function

% f(x)= 0.5* ||A*x-b||^2+alpha*||x||_TV

% Author: Sovan Mukherjee, Aproil 2015

function fx  = obj_func(D1,D2,x,A,b,alpha)

D1_x= D1*x;

D2_x= D2*x;

TV_norm_x= sum(sqrt(abs(D1_x).^2+abs(D2_x).^2));

fx= 0.5*norm(A*x-b).^2+alpha*TV_norm_x;