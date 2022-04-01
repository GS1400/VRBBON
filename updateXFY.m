
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updateXFY.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point] = updateXFY(point,step,tune)
% update both X, Y and F
%
% for details of input and output structures see VRBBON.m
%
function [point] = updateXFY(point,step,tune)

if all(point.X(:,1)==0)
   point.X(:,1)=[]; point.F(1)=[];  point.Y=[]; point.m= point.m-1;
   point.b=point.b-1;
end

if ~all(point.xm==0)

    if point.m>=tune.mmax
        [~,iw] = max(point.F); % find worst point
    else
        point.m=point.m+1;
        iw = point.m;
    end
  

    point.X(:,iw) = point.xm; point.F(iw) = point.fm;
    point.Y(iw)   = step.alphaE;
    point.b = iw;
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
