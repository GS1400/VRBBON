
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updatecumu.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [step] = updatecum(point,step,tune)
% update r and q for generating cumulative step
%
% for details of input and output structures see VSBBO.m
%
function [step] = updatecum(point,step,tune)

fields = fieldnames(point);
for i=1:length(fields),
    field=fields{i};
    eval([field,'=point.',field,';']);
end;

a=tune.a; r=step.r; q=step.q; p=step.p;

% state depends on whether extrapolation or opposote direction
if fr<fm, d = 4*fm-3*fr-fl; else, d = fl-fr; end
h = fr+fl-2*fm; 



% finding alpha depends on whether h is positive or not
if h<=0, 
   if d>=0, alpha = a; else, alpha = -a; end
else
   if d>=0, alpha = min(a,d/(2*h)); else, alpha = max(-a,d/(2*h)); end 
end   


% update r and q
step.q = q +alpha*p; step.r = r + 0.5*alpha*(d-alpha*h);
step.d=d;step.h=h;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
