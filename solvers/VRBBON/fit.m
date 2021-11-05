
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step]=fit(point,step,tune)
% construct a reduced quadratic model
%
% for details of input and output structures see VRBBON.m
%

function [point,step]=fit(point,step,tune)

[point]      = adjustX(point,tune); % adjust X
[point,step] = scaleModel(point,step,tune); % get sc


step.vec     = step.sc;
[step]       = adjustVec(step,tune);
step.sc      = step.vec;

step.vec     = (point.F(point.ind)-point.fm)'./step.sc;
[step]       = adjustVec(step,tune);
step.b       = step.vec;

[step]       = getMultiplier(point,step); % get y

step.vec     = step.y;
[step]       = adjustVec(step,tune);
step.y       = step.vec;

[point,step] = getGg(point,step);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%