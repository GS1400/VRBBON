
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updateSY.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point] = updateSY(point)
% update both S and Y
%
% for details of input and output structures see VSBBO.m
%
function [point] = updateSY(point,tune)

s0 = point.xm - point.xm0; 
y0 = point.grad-point.grad0;
point.hdiag = s0'*y0/(y0'*y0);
if s0'*y0>sqrt(eps)*max(eps,norm(s0)*norm(y0))
    
   if point.mq<tune.mqmax, point.mq=point.mq+1;
   else, point.mq=1; 
   end;
   
   point.S(:,point.mq)=s0;
   point.Y(:,point.mq)=y0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
