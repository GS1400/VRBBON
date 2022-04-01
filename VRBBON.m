
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VRBBON.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,info] = VRBBON(fun,x,st,tune);
%
% solve unconstrained noisy black box optimization problem 
%    min f(x) 
%  
% fun      % function handle for f(.)
% x        % starting point (must be specified)
% st       % structure with stop and print criteria
%          % (indefinite run if no stopping criterion is given)
%  .secmax       %   stop if sec>=secmax (default: inf)
%  .nfmax        %   stop if nf>=nfmax   (default: inf)
%  .qf           %   stop if qf<=accf    (default: 1e-4)
%  .fbest        %   optimal function value (default: 0)
%  .prt          %   printlevel (default: -1)
%                %   -2: nothing, -1: litte, >=0: more and more
% tune     % optional structure containing tuning parameters
%          %   for details see below
%
% x        % best point found 
% f        % function value at best point found 
% info     % structure with performance info
%          %   for details see below
% 
function [x,f,info] = VRBBON(fun,x,st,tune)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info.error='';
% check function handle
if isempty(fun)
    
  message = 'VRBBON needs the function handle fun to be defined';
  info.error= message; 
  return
elseif ~isa(fun,'function_handle')
  message = 'fun should be a function handle';
  info.error= message;
  return
end

% initialize structure with information about points and function values
%
% point  % structure with information about points and function values
%  .x        %   point with best function value found
%  .f        %   best function value at x
%  .X        %   a list of the m best point so far
%  .F        %   the function values of X
%  .xr       %   newest point 
%  .fr       %   function value at xr
%  .xm       %   point at xr-p
%  .fm       %   function value at xm 
%  .xinit    %   initial point
%  .finit    %   initial function value
%  .b        %   index of best point
%  .m        %   number of best point kept
%
% starting point
if isempty(x)
  message = 'starting point must be defined';
  info.error= message;
  return      
elseif ~isa(x,'numeric')
  message = 'x should be a numeric vector'; 
  info.error= message;
  return       
end
% dimension
point.n = length(x);
% artificial lower and upper bounds to help prevent overflow
point.low=x-1e20;
point.upp=x+1e20;
point.x=x;
% other point info
point.X     = point.x; point.F  = NaN;   point.xinit = point.x; 
point.finit = NaN;     point.fr = NaN;   point.fm    = NaN;
point.b     = NaN;     point.m  = 0; 
point.xr    = NaN+ones(point.n,1); 
point.xm    = NaN+ones(point.n,1); 



% initialize structure containing all tuning parameters 
%
% tune % structure containing all tuning parameters 
%      % all parameters have a default that can be overwritten 
%      % by specifying it as an input
%  .mmax      % maximum number of best points kept
%  .T0        % maximal number of calls to RLS in DS
%  .R         % maximal number of random directions in RLS
%  .E         % maximal number of extrapolations in RLS
%  .deltamin  % minimum threshold for stopping VRBBON
%  .deltamax  % initial step size
%  .gamma     % factor for adjusting forcing function
%  .gammaE    % factor for extrapolation test
%  .gammaX    % factor for adjusting X
%  .gammaa    % tiny factor for finding a robustified interval
%  .gammar    % factor for adjusting the random direction
%  .gammad1   % factors for updating the trust region radius
%  .gammad2   % factors for updating the trust region radius
%  .gammav    % factor for adjusting sc and y
%  .alpmin    % minimum step size in RLS
%  .Q         % factor for adjusting Delta
%  .dmin      % minimum value for the trust region radius
%  .dmax      % maximum value for the trust region radius
%  .alow      % lower bound for initial interval
%  .aupp      % upper bound for initial interval  
%  .gammap    % factor for adjusting the trust region direction
%  .minqeps   % minimum threshild for stopping minq8
%  .minqmax   % maxiumu number of iterations in minq8
%  .model     % linear or quadratic model
%  .mbis      % parameter for bisection
%  .comBound  % complexity bound verified?


if ~exist('tune'), tune=[]; end;

% maximum number of best points kept
if ~isfield(tune,'mmax'), 
    tune.mmax = min(0.5*point.n*(point.n+3),230);
end;

% maximum number of calls to robustified line search
% using random direction in DS
if ~isfield(tune,'T0'), tune.T0 = point.n; end;


% complexity bound verified?
if ~isfield(tune,'comBound'), tune.comBound = 2; end;

switch  tune.comBound
    case 2
        if ~isfield(tune,'R'), tune.R = point.n; end;
    case 1
        if ~isfield(tune,'R'), tune.R = max(1,fix(point.n/2)); end;
        if ~isfield(tune,'C'), tune.C = max(1,fix(point.n/2)); end;
    case 0
        if ~isfield(tune,'C'), tune.C = point.n; end;
    otherwise 
        error('comBound should be one of 0,1,2')
end

if isfield(tune,'C')
  tune.C = min(tune.C,point.n);
end
if isfield(tune,'R')
    tune.R = min(tune.R,point.n);
end



% maximal number of random trials
if ~isfield(tune,'E'), tune.E= Inf; end;

% minimum threshold for stopping VRBBON
if ~isfield(tune,'deltamin'), tune.deltamin = 0; end;

% initial step size 
if ~isfield(tune,'deltamax'), tune.deltamax = 1; end;

% factor for adjusting forcing function
if ~isfield(tune,'gamma'), tune.gamma = 1e-6; end;

% factor for extrapolation test
if ~isfield(tune,'gammaE'), tune.gammaE = 3; end;

% factor for adjusting X
if ~isfield(tune,'gammaX'), tune.gammaX = 100; end;

% tiny factor for finding a robustified interval
if ~isfield(tune,'gammaa'), tune.gammaa = 1e-5; end;

% factor for adjusting the random direction
if ~isfield(tune,'gammar'), tune.gammar = 1e-30; end;

% factors for updating the trust region radius
if ~isfield(tune,'gammad1'), tune.gammad1 = 0.1; end;
if ~isfield(tune,'gammad2'), tune.gammad2 = 2; end;

% factor for adjusting sc and y
if ~isfield(tune,'gammav'), tune.gammav = 100; end;

% minimum step size in RLS
if ~isfield(tune,'alpmin'), tune.alpmin = 1e-3*rand; end;

% factor for adjusting delta
if ~isfield(tune,'Q'), tune.Q = 1.5; end;

% minimum value for the trust region radius
if ~isfield(tune,'dmin'), tune.dmin = 1e-4; end;

% maximum value for the trust region radius
if ~isfield(tune,'dmax'), tune.dmax = 1e+3; end;

% lower bound for initial interval
if ~isfield(tune,'alow'), tune.alow = 1e-2; end;

% upper bound for initial interval
if ~isfield(tune,'aupp'), tune.aupp = 0.99; end;

% factor for adjusting the trust region direction
if ~isfield(tune,'gammap'), tune.gammap = 0.25; end;

% minimum threshild for stopping minq8
if ~isfield(tune,'minqeps'), tune.minqeps = 1e-8; end;

% maxiumu number of iterations in minq8
if ~isfield(tune,'minqmax'), tune.minqmax = 10000; end;

% linear or quadratic model
if ~isfield(tune,'model'), tune.model = 1; end;

% parameter for perturbed random direction
if ~isfield(tune,'gammak'), tune.gammak = 0.85; end;



% initialize structure for step management
%

% step   % structure for step management
%  .p        % random search direction
%  .deltamin % 
%  .deltamax % 
%  .delta    % steps size
%


% others
step.delta = NaN; 
step.p     = NaN+ones(point.n,1); 
step.alphaE = NaN;
step.alow   = 0;
step.aupp   = Inf;


% par  % structure containing parameters modified during the search
%  .T        % maximal number of random directions
%  .df       % the difference of best point and the new point
%  .ngood    % number of extrapolation stages
%  .dir      % select direction type
%  .state    % state of cumulative step
%  .desent   % is there the initial interval or not ?
%
par.state = NaN;  par.dir = NaN; par.T = NaN; 
par.df= 1; par.T = NaN; par.ngood =0;
par.state=-1; par.desent=1;

% info  % performance information for VRBBON
%  .prt          % printlevel 
%                %   -1: nothing, 0: litte, >=1: more and more
%  .secmax       % stop if sec>=secmax 
%  .nfmax        % stop if nf>=nfmax 
%  .finit        % the initial f
%  .fbest        % the optimal f
%  .acc          % minimum thresold for qf 
%  .qf           % stop if qf <= acc
%  .initTime     % inital cputime
%  .done         % done with the search?
%
if ~exist('st'), st=[]; end;
% print level
if isfield(st,'prt'), info.prt = st.prt; 
else, info.prt = -2; 
end;
% stopping criteria
if isfield(st,'secmax'), info.secmax=st.secmax;
else, info.secmax=inf;
end;
if isfield(st,'nfmax'), info.nfmax=st.nfmax;
else, info.nfmax=500*point.n;
end;
if isfield(st,'fbest'), info.fbest=st.fbest;
else, info.fbest=0;
end;

if isfield(st,'finit'), info.finit=st.finit;
else, info.finit=fun(x);
end;

if isfield(st,'accf'), info.accf=st.accf;
else, info.accf=1e-4;
end;



if isfield(st,'reallife'), info.reallife=st.reallife;
else, info.reallife=0;
end;


info.initTime=cputime; info.nquadmodel=0; info.nlinearmodel=0;
info.model = 0;

if info.reallife, info.arrayqf =1; info.arraynf = 0; end
     


if info.prt>=0,
   disp(' ')
   disp(' ')
   disp('======================================================')
   disp('======================================================')
   disp('start VRBBON')
   disp('======================================================')
   disp('======================================================')
end




for round=1
    
  % one round only, needed to be ableto jump to the end
  % when a stopping test is satisfied

   % get and regularize first function value
    info.nf=1;
    point.x = max(point.low,min(point.upp,point.x));
    fnew = fun(point.x);
    if isnan(fnew), fnew=inf; end;
    point.f=max(-1e50,min(fnew,1e50));
    
    
    % initial function value for stopping test
    info.finit = point.f; 

    % get initial values for subspace 
    point.F = point.f; point.xm = point.x; point.fm = point.f; 
    point.b = 1; point.m =1; point.Y(1) = 1; par.itc = 0;
    
    step.aupp = tune.aupp; step.alow = tune.alow;
    
    %%%%%%%%%%%%%
    % main loop %
    %%%%%%%%%%%%%
    if info.prt>=0
      disp(' ')
      disp(' ') 
      disp('start of main loop:')
      disp(' ')
      disp(' ')
    end
   
    
    while 1
        
        % update the information
        point.xm=point.X(:,point.b); point.fm = point.F(point.b);
        point.xr = point.xm; point.finit=point.fm; par.stateTest = '';
        par.ngoodi=par.ngood; par.ngood = 0; par.ngoodold=0;
        point.xinit=point.xm; 
        
        par.itc = par.itc + 1;
        
        if info.prt>=0
           disp(['**************  DS at iter = ', ...
               num2str(par.itc),'  **************'])
        end
        
       
        
        % try to get a good decrease on the function value
        [point,par,step,info]=DS(fun,point,step,par,tune,info);
        
        
        if info.done, break; end  
        if step.delta <= tune.deltamin, break; end
        % expand or reduce the step size depending on
        % whether there exits any improvement on f or not
        
        if info.prt>=0
           disp(['delta = ',num2str(step.delta),';'])
        end
        
        [step]=updateStepSize(step,par,tune);
        
      
    end
    
    if info.done, break; end  
    
    
    if info.prt>=0
       disp(' ')
       disp('======================================================')
       disp('======================================================')
       disp('end  of main loop')
       disp('======================================================')
       disp('======================================================')
       disp(' ')
    end
   
end


% assign results

x = point.X(:,point.b);
f = point.F(point.b);


% status detection 

if info.qf<=info.accf
    info.status = 'accuracy reached';
elseif info.nf>=info.nfmax
     info.status = 'nfmax reached';
elseif info.sec>=info.secmax
    info.status = 'secmax reached';
elseif step.delta<=tune.deltamin
    info.status = 'deltamin reached';
else
    info.status = 'unknown';
end

if info.prt>=0,
  disp('==============================================================')
  disp('end VRBBON')
  disp('==============================================================')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

