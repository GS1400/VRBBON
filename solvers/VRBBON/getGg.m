
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% getGg.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step]=getGg(point,step)
% estimate the gradient vector and Hessian matrix
%
% for details of input and output structures see VRBBON.m
%

function [point,step]=getGg(point,step)
n = point.msim;
point.g = step.y(1:n);
point.randg = all(step.y(n+1:2*n))==0;
if ~point.randg
    
     step.vec     = step.y;
    [step]        = adjustVec(step);
     step.y       = step.vec;
     
     for i=1:n, point.G(i,i) = step.y(n+i); end

     l = 2*n+1;
     for i = 1:n-1
       for j=i+1:n
         point.G(i,j) = step.y(l);
         point.G(j,i) = step.y(l);
         l = l + 1;
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%