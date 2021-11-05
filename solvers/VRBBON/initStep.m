%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initStep.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [step]=initStep(step,par,tune)
% generate initial step size
%
% for details of input and output structures see VSBBON.m
%
function [step]=initStep(step,par,tune)

if par.itc==1
   step.delta=tune.deltamax; step.alphaE = step.delta;
else

    ok = (step.alow~=0 & ~isinf(step.aupp));

    if ok
        step.alp1=step.alow; step.alp3=step.aupp;
        [step]=bisection(step);
        step.alphaE  = max(step.delta,step.alp);
    end
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%