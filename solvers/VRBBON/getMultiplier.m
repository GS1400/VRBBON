%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% getMultiplier.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [step]=getMultiplier(point,step)
% get multiplier to estimate the gradient vector and Hessian matrix
%
% for details of input and output structures see VRBBON.m
%


function [step]=getMultiplier(point,step)
X   = point.XX; 
xm  = point.xxm;
ind = point.ind; 
K   = point.K;
A   = [X(ind,:) - ones(K,1)*xm 0.5*(X(ind,:) - ones(K,1)*xm).^2];
n=size(X,2);
for i=1:n-1
  B = (X(ind,i)-xm(i))*ones(1,n-i);
  A = [A B.*(X(ind,i+1:n)-ones(K,1)*xm(i+1:n))];
end
A = A./(step.sc*ones(1,size(A,2)));

warning off
step.y = A\step.b;
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%