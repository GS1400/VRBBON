
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MLS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,par,info]=MLS(fun,point,step,par,tune,info);
%
% try to significantly improve function value by gammamin*Delta
%
% inputs: fun,point,step,par,tune

% outputs: point,good,lam,info

function [point,par,step,info] = MLS(fun,point,step,par,tune,info)

if tune.alg >=3
  point.xm=point.X(:,point.b); point.fm = point.F(point.b);
end
point.xr = point.xm; point.finit=point.fm; 
info.stateTest = ''; step.r = 0;  step.q = zeros(point.n,1); par.good=0;
par.t=0; prt=(info.prt>=1);

QN = (tune.alg==2 | tune.alg==5); % the gradient is estimated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while 1
    
    par.t=par.t+1;   
   
    % identify dir
    [par,step,info]  = identifyDir(point,step,par,tune,info);
    
    % compute direction
    [step,par] = direction(point,step,par,tune);
    
    
    % perform extrapolation step
    [point,step,par,info]=extrapolatStep(fun,point,step,par,tune,info);
    
   
    if info.done, break;   end
    
    if info.nE == 0 
        % estimate the tth component of the gradient
        if par.dir==1 & QN
           point.grad(par.t,:) = (point.fr-point.finit)/(par.A(par.t));
        end

         if par.dir==1 
            point.xr(par.t) = point.xm(par.t);
         end
         % reduce the step size
         par.A(par.t) = max(tune.alphamin,par.A(par.t)/tune.gammaE); 
    else
        % estimate the tth component of the gradient
        if par.dir==1 & QN
           point.grad(par.t,:) = (point.fm-point.finit)/(par.A(par.t));
        end  
    end
  
    if par.t==par.T, break; end
    
end  % of while T>0


if prt
    info.stateTest = info.stateTest(1:end-1);
    info.stateTest=[info.stateTest,']'];
    stateTest = info.stateTest
end

info.laststate=info.stateTest;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
