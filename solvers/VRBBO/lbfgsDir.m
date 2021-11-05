%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% lbfgsDir.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,par,info]=lbfgsDir(fun,point,step,par,tune,info);
%
% computes the approximate inverse Hessian multiplied by the gradient
%
% inputs: fun,point,step,par,tune,info
%
% outputs: point,par,info

% Ref.
% Nocedal, J., Updating Quasi-Newton Matrices with Limited Storage,
% MATHEMATICS OF COMPUTATION. 35(151), 773--782 (1980).

function step = lbfgsDir(point,step)

S     = point.S; 
Y     = point.Y;
g     = point.grad;
hdiag = point.hdiag;

[n,k] = size(S);
for i = 1:k, ro(i,1) = 1/(Y(:,i)'*S(:,i)); end

q = zeros(n,k+1);
r = zeros(n,1);
alpha =zeros(k,1);
beta =zeros(k,1);

q(:,k+1) = g;

for i = k:-1:1
    alpha(i) = ro(i)*S(:,i)'*q(:,i+1);
    q(:,i)   = q(:,i+1)-alpha(i)*Y(:,i);
end

r = hdiag*q(:,1);

for i = 1:k
    beta(i) = ro(i)*Y(:,i)'*r;
    r       = r + S(:,i)*(alpha(i)-beta(i));
end
step.p = -r;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%