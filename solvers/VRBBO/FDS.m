
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,info] = FDS(fun,point,step,par,tune,info)
% repeatedly impove function value by Delta
%
% for details of input and output structures see VRBBO.m
%
function [point,info] = FDS(fun,point,step,par,tune,info)

ngood = 0;
sub = (tune.alg >=3); % random subspace is used

QN = (tune.alg==2 | tune.alg==5); % if QN is true, LBFGS is used
while 1
       if QN
          point.grad0 = point.grad; point.xm0 = point.xm;
       end
       if info.prt>=0,
          disp('=============')
          disp('start MLS')
          disp('=============')
       end
       [point,par,step,info] = MLS(fun,point,step,par,tune,info);
       if info.prt>=0,
          disp('=============')
          disp('end MLS')
          disp('=============')
       end
       % update X, F, and so on
       if sub, [point] = updateXF(point,tune); end
   
       if info.prt>=0, nfFDS=info.nf, end;
       if ~ par.good, break; end
       ngood      = ngood+1;
       
       if info.done, break; end;
       % update S, Y, and so on
       if QN, [point] = updateSY(point,tune); end
end
info.ngood = ngood;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
