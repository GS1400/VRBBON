
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VSBBO.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,info] = VRBBO(fun,x,st,tune);
%
% solve global unconstrained black box optimization problem 
%    min f(x) 
%  
% fun      % function handle for f(.)
% x        % starting point (must be specified)
% st       % structure with stop and print criteria
%          % (indefinite run if no stopping criterion is given)
%  .secmax       %   stop if sec>=secmax (default: inf)
%  .nfmax        %   stop if nf>=nfmax   (default: inf)
%  .ftarget      %   stop if f<=ftarget  (default: -inf)
%  .prt          %   printlevel (default: -1)
%                %   -1: nothing, 0: litte, >=1: more and more
% tune     % optional structure containing tuning parameters
%          %   for details see below
%
% x        % best point found 
% f        % function value at best point found 
% info     % structure with performance info
%          %   for details see below
% 
function [x,f,info] = VRBBO(fun,x,st,tune);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check function handle
if isempty(fun)
  message = 'VRBBO needs the function handle fun to be defined';
  disp(message)  
  return
elseif ~isa(fun,'function_handle')
  message = 'fun should be a function handle';
  disp(message)
  return
end


% initialize structure with information about points and function values
%
% point  % structure with information about points and function values
%  .x        %   point with best function value found
%  .f        %   best function value at x
%  .X        %   a list of the m best point so far
%  .F        %   the function values of X
%  .xl       %   point at xr-2p
%  .fl       %   function value at xr-2p
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
  disp(message) 
  return      
elseif ~isa(x,'numeric')
  message = 'x should be a numeric vector'; 
  disp(message ) 
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
point.fl    = NaN;     point.b  = NaN;   point.ms     = 0; 
point.xr    = NaN+ones(point.n,1); point.xm = NaN+ones(point.n,1); 
point.mq    = 0;


% initialize structure containing all tuning parameters 
%
% tune % structure containing all tuning parameters 
%      % all parameters have a default that can be overwritten 
%      % by specifying it as an input
%  .mmax        % maximum number of best points kept
%  .T0          % maximal number of MLS in setScale
%  .C           % maximal number of coordinate direction
%  .S           % maximal number of random subspace directions
%  .R           % maximal number of random directions in MLS
%  .E           % maximal number of extrapolations in MLS
%  .secSub      % scale subspace direction ?
%  .secCum      % scale cumulative direction ?
%  .cum         % select cumulative step type
%  .a           % upper bound to estimate cumulative step-size
%  .Deltamin    % minimium threshold for good impovement
%  .Deltamax    % initial threshold for good impovement
%  .deltamin    % minimum norm of trial step
%  .deltamax    % minimum norm of trial step
%  .gammadelta  % factor for finding delta
%  .gammamax    % factor for adjusting Deltamax
%  .gammaE      % factor for extrapolation test
%  .gammalambda % factor for finding initial lam
%  .Q           % factor for adjusting Delta
%  .alg         % identify the type of algorithm
%  .alphamin    % minimum threshold for extrapolation step sizes 

%
if ~exist('tune'), tune=[]; end;
% maximum number of best points kept
if ~isfield(tune,'msmax'), tune.msmax = min(point.n+1,5); end;
% maximum number of memory for L-BFGS
if ~isfield(tune,'mqmax'), tune.mqmax = min(point.n,5); end;

% maximal number of multi-line searches in setScales
if ~isfield(tune,'T0'), tune.T0 = 2*tune.msmax; end;


% identify the type of algorithm
if ~isfield(tune,'alg'), tune.alg = 5; end;

switch  tune.alg
    
    case 0 % basic version, R<n
        % number of random direction
        if ~isfield(tune,'R'),tune.R = fix(point.n/2)+1; end;
    case  1 % basic version R>=n
        % number of random direction 
        if ~isfield(tune,'R'),tune.R =point.n; end;
        
        
    case 2 % no random subspace direction
           
         % maximal number of coordinate trials
        if ~isfield(tune,'C'), 
           tune.C =point.n;
        end;
        
        % number of random direction  
        if ~isfield(tune,'R'), 
          tune.R = min(fix(point.n/10)+1,20); 
        end; 
                
    case 3 % no coordinate direction
        
        % maximal number of subspace directions in multi-line search
        if ~isfield(tune,'S'), tune.S=fix(point.n/5); end; 
        % number of random direction
        if ~isfield(tune,'R'),tune.R = fix(point.n/2)+1; end;
        
        
    otherwise % alg=4 (no quasi Newton direction)
              % alg=5 (all directions)
        % maximal number of coordinate trials
        if ~isfield(tune,'C'), 
           tune.C =point.n;
        end;
        % maximal number of subspace directions in multi-line search
         if ~isfield(tune,'S'), tune.S = min(fix(point.n/10)+1,5); end; 
        % number of random direction
        if ~isfield(tune,'R'), 
          tune.R = min(fix(point.n/10)+1,20); 
        end;
end

% maximal number of random trials
if ~isfield(tune,'E'), tune.E= 10; end;
% scale for subspace direction
if ~isfield(tune,'scSub'), tune.scSub = 0; end;
% scale for cumulative direction
if ~isfield(tune,'scCum'), tune.scCum = 0; end;
% scale for heuristic direction
if ~isfield(tune,'scCor'), tune.scCor = 0; end;
% select cumulative step type
if ~isfield(tune,'cum'), tune.cum = 1; end;
% upper bound for estimating cumulative step-size
if tune.cum==2
   if ~isfield(tune,'a'), tune.a = 3; end;
end
% minimium threshold for good impovement
if ~isfield(tune,'Deltamin'), tune.Deltamin = 0;end;
% initial threshold for good impovement
if ~isfield(tune,'Deltamax'), tune.Deltamax = 1e-3; end; 
% factor for finding delta
if ~isfield(tune,'gammadelta'), tune.gammadelta = 1e6; end;
% factor for adjusting Deltamax
if ~isfield(tune,'gammamax'), tune.gammamax = 1e-3; end;
% factor for extrapolation test
if ~isfield(tune,'gammaE'), tune.gammaE = 4; end;
% factor for finding initial lam
if ~isfield(tune,'gammalambda'), tune.gammalambda = 1e-6; end;
% factor for angle condition
if ~isfield(tune,'gammaw'), tune.gammaw = eps; end;
% factor for adjusting Delta
if ~isfield(tune,'Q'), tune.Q = 2; end;
% factor for angle condition
if ~isfield(tune,'gammaangle'), tune.gammaangle = 1e-20; end;
% minimum norm of trial step
if ~isfield(tune,'deltamin'), tune.deltamin = 1e-4*sqrt(point.n); end;
% maximum norm of trial step
if ~isfield(tune,'deltamax'), tune.deltamax = 0.1*sqrt(point.n); end;
% minimium threshold for extrapolation step sizes
if ~isfield(tune,'alphamin'), tune.alphamin = 1e-50; end;
% parameter for extrapolation test
if ~isfield(tune,'gammamin'), tune.gammamin = 1e-6; end;

if isfield(st,'reallife'), info.reallife=st.reallife;
else, info.reallife=0;
end;

% initialize structure for step management
%

% step   % structure for step management
%  .s        % scaling vector
%  .p        % random search direction
%  .dp       % scaled length of p
%  .deltamin % minimum norm of trial steps
%  .deltamax % maximum norm of tial steps
%  .delta    % norm of trial steps 
%  .Delta    % threshold for good impovement
%  .Deltamin % minimium threshold for good impovement
%  .Deltamax % initial threshold for good impovement
%  .q        % cumulative step
%  .r        % cumulative gain
%
% initial gain required for a good step
step.Deltamin = tune.Deltamin;
% maximal gain required for a good step
step.Deltamax = tune.Deltamax;
% others
onesn = ones(point.n,1); 
step.r = NaN; step.dp = NaN; step.Delta = NaN;
step.delta = NaN; step.s = onesn; 
step.q= NaN+onesn; step.p = NaN+onesn; 


% par  % structure containing parameters modified during the search
%  .T        % maximal number of random directions
%  .lam      % lower bound on the Lipschitz constant
%  .good     % indicator for improvement
%  .ss       % are we in setScale?
%  .dir      % select direction type
%  .state    % state of cumulative step
%  .hs       % parameter for heuristic direction
%  .nmax     % parameter for heuristic direction
%
par.lam = NaN;   par.state = NaN;  par.dir = NaN; par.T = NaN; 
par.ss  = NaN;   par.good  = NaN; par.df= NaN; 


% info  % performance information for VRBBO
%  .prt          % printlevel 
%                %   -1: nothing, 0: litte, >=1: more and more
%  .secmax       % stop if sec>=secmax 
%  .nfmax        % stop if nf>=nfmax 
%  .ftarget      % stop if f<=ftarget 
%  .initTime     % inital cputime
%  .done         % done with the search?
%  .nf           % number of function evaluations used
%  .nstep        % structure counting steps of a particular kind
%    .nCoor      % number of coordinate direction   
%    .nRand      % number of scaled Random direction   
%    .nExt       % number of extrapolations
%    .nOrt       % number of opposite direction
%    .nCum       % number of cumulative direction
%    .nSub       % number of random subspace direction
%    .nQN        % number of quasi Newton direction
%  .laststate    % state of last iteration
%  .ngood        % number of times that good holds in FDS
%  .nDelta       % number of times that Delta reduces in VRBBO
%
% print level

if ~exist('st'), st=[]; end;


if isfield(st,'prt'), info.prt = st.prt; 
else,  info.prt = -1; 
end
% stopping criteria
if isfield(st,'secmax'), info.secmax=st.secmax;
else, info.secmax=200*point.n;
end
if isfield(st,'nfmax'), info.nfmax=st.nfmax;
else, info.nfmax=200*point.n;
end
if isfield(st,'fbest'), info.fbest=st.fbest;
else, info.fbest=-inf; % 0.01*fun(x);
end
if isfield(st,'accf'), info.accf=st.accf;
else, info.accf=-inf;
end

if isfield(st,'finit'), info.finit=st.finit;
else, info.finit=fun(x);
end;



info.initTime=cputime;
% other
info.nstep.nCoo  = 0;  % number of coordinate direction
info.nstep.nNew  = 0;  % number of quasi Newton direction
info.nstep.nSub  = 0;  % number of subspace direction
info.nstep.nRan  = 0;  % number of scaled Random direction
info.nstep.nExt  = 0;  % number of extrapolate
info.nstep.nOpp  = 0;  % number of opposite direction
info.nstep.nCum  = 0;  % number of cumulative step
info.nDelta      = 0;  % number of times that Delta reduces in VRBBO

if info.reallife, info.arrayqf = 1; info.arraynf = 0; end


if info.prt>=0,
   disp('======================')
   disp('start VRBBO')
   disp('======================')
end


for round=1, 
  % one round only, needed to be ableto jump to the end
  % when a stopping test is satisfied

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                      get scales
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if info.prt>=0,
    disp('======================')
    disp('start setScale')
    disp('======================')
  end
  [point,step,par,info] = setScales(fun,point,step,par,tune,info);
  if info.done, break; end;
  if info.prt>=0,
  disp('======================')
  disp('end setScale')
  disp('======================')
  end
  
  if info.done, break; end;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                      fixed decrease search
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if info.prt>=0,
    disp('======================')
    disp('start FDS')
    disp('======================')
  end
  while 1
    [point,info] = FDS(fun,point,step,par,tune,info);
     % check stopping test
    ok = (info.done | step.Delta<=tune.Deltamin);
    if ok, break; end;
    % gain of Delta unlikely with step of lenghth delta
    step.Delta = step.Delta/tune.Q; % reduce Delta
    info.nDelta = info.nDelta+1;
  end
  if ok, break; end;
end

if info.prt>=0,
    disp('======================')
    disp('end FDS')
    disp('======================')
end

% assign results
if tune.alg>=3
    x = point.X(:,point.b);
    f = point.F(point.b);
else
    x=point.xm; f=point.fm;  
end


if info.prt>=0,
  disp('======================')
  disp('end VRBBO')
  disp('======================')
end

if info.qf<=info.accf,  
    info.status='accuarcy reached';
elseif info.nf>=info.nfmax, 
    info.status='nfmax reached';
elseif info.sec>=info.secmax,    
    info.status='secmax reached';
else
    info.status='unknown';
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

