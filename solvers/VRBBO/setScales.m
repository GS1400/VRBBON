
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setScales.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step,par,info]=setScales(fun,point,step,par,tune,info)
% estimate Deltamax, lam and s
%
% for details of input and output structures see VSBBO.m
%
function [point,step,par,info] = setScales(fun,point,step,par,tune,info)
 
% get and regularize first function value
info.nf=1;
point.x = max(point.low,min(point.upp,point.x));
fnew = fun(point.x);
if isnan(fnew), fnew=inf; end;
point.f=max(-1e50,min(fnew,1e50));

if ~isfield(info,'finit'), info.finit=point.f; end

% get initial values for subspace 

if tune.alg>=3
    point.F = point.f; point.b = 1; point.m=1; point.X = point.x;
end

point.xm = point.x; point.fm = point.f; 

% get initial values for T
switch tune.alg
    case 2 % random subspace directions are not used
      par.T =tune.C+tune.R+2; 
    case 3 %  random, random subspace, and cumulative are used
      par.T = tune.R+tune.S+1; 
    case 4 % LBFGS is not used
      par.T = tune.C+tune.R+tune.S+1; 
    case 5 % all directions are used
       par.T = tune.C+tune.R+tune.S+2;  
    otherwise  % basic version; only random directions are used    
      par.T = tune.R;
end

par.lam = 1;


step.deltamin = tune.deltamin;
step.deltamax = tune.deltamax;

step.delta    = tune.deltamin;  
step.Delta    = step.Deltamax;
   
info.nf    = 1;
par.A      = ones(par.T,1);
step.Delta = tune.Deltamax;	
point.fext = point.fm;

step.deltamin = tune.deltamin;
step.deltamax = tune.deltamax;
step.s        = ones(point.n,1);

sub = (tune.alg >=3); % random subspace is used

QN = (tune.alg==2 | tune.alg==5); % if QN is true, LBFGS is used

if QN
    point.S = []; % m previous parameter values
    point.Y = []; % m previous gradient values
    point.ms =0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             get X in order to estimate Deltamax                % 
%             MLS doesn't use subspace direction                 %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tt=1:tune.T0
    point.tt=tt;
    QN1 = (QN & tt>=2);
    if QN1, 
       point.grad0 = point.grad; point.xm0 = point.xm;
    end
    if info.prt>=0,
      disp('=============')
      disp('start MLS')
      disp('=============')
    end
   [point,par,step,info] = MLS(fun,point,step,par,tune,info);
   if info.prt>=0,
     disp('==============')
     disp('end MLS')
     disp('==============')
   end
   
   if sub, [point] = updateXF(point,tune); end
   
   if info.done, break; end;
   
   
   if info.prt>=0, nfSS=info.nf, end;
   
   if QN1, [point] = updateSY(point,tune); end
   
end

if ~sub, estim =  0;
    
    % (~info.done & sub & point.ms>=3)
else, estim =  (sub & point.ms>=3); 
end

if estim % estimate Deltamax, s, and lam
    
    II  = setdiff(1:point.ms,point.b);
    XX  = point.X(:,II);

    for i=1:point.ms-1
        dX(:,i)  = XX(:,i)-point.X(:,point.b);
    end

    step.s = max(dX')'; step.s(step.s==0) = 1;

    dF = median(abs(point.F-point.F(point.b))); 

    if dF==0
        step.Deltamax = tune.gammamax;
        step.Delta    = step.Deltamax;
        par.lam       = (tune.gammalambda)/sqrt(point.n);
    else
        step.Deltamax = tune.gammamax*min(dF,1);
        step.Delta    = step.Deltamax;
        par.lam       = (tune.gammalambda)*dF/sqrt(point.n);
    end
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
