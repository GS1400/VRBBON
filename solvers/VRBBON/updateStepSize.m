%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% updateStepSize.m %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [step]=updateStepSize(step,par,tune)
% update the step size
%
% for details of input and output structures see VRBBON.m
%

function [step]=updateStepSize(step,par,tune)

if par.totalngood == 0
   step.delta = step.delta/tune.Q;
else
    
    ok = (step.alow~=0 & ~isinf(step.aupp));

    if ok
        step.alp1=step.alow; step.alp3=step.aupp;
        [step]=bisection(step);
        step.delta  = max(step.delta,step.alp);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      