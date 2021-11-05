%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% initInterval.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step]=initInterval(point,step)
% Get initail interval for step size by using the information obtained
% by extrapolation stage and initail tuning inteval 
%
% for details of input and output structures see VRBBON.m
%

function [point,step]=initInterval(point,step)

ind1=find(point.dlist<0); 

if ~isempty(ind1), step.alow=min(step.alow,max(step.alist(ind1))); end

ind2=find(point.dlist>=0 & step.alist>step.alow);

if ~isempty(ind2), step.aupp=max(step.aupp,min(step.alist(ind2))); end;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
