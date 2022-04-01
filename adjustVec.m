%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adjustVec.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [step]=adjustVec(step,tune)
% remove NaN or inf from a vector
%
% for details of input and output structures see VRBBON.m
%

function [step]=adjustVec(step,tune)

for i=1:length(step.vec)
    Inan  = isnan(step.vec);
    if any(Inan), step.vec(Inan) = tune.gammav;end
    Iinf= (step.vec==inf);
    if any(Iinf), step.vec(Iinf) = tune.gammav;end
    Iinf= (step.vec==-inf);
    if any(Iinf), step.vec(Iinf) = -tune.gammav;end  
end  

Iz = (step.vec==0);
if any(Iz)
   maxsc=max(step.vec);
   step.vec(Iz)=maxsc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%