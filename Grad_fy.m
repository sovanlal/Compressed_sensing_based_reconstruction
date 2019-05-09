% this function calculates the gradient of the objective function

% f(x)= 0.5* ||A*x-b||^2+alpha*||x||_TV

% Author: Sovan Mukherjee, April 2015

function grad_fx  = Grad_fy(A,b,x,alpha,D1,D2,N)

D=[D1;D2];

D1_x= D1*x;

D2_x= D2*x;

TV= [D1_x;D2_x];

dTV=alpha*D'*TV;

dData= A'*(A*x-b);

grad_fx= dData+dTV;