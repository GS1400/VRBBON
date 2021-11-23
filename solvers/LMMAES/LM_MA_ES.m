%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Limited Memory LM-MA-ES with gamma z-paths (stored in p_matrix), 
%% cumulative step size adaptation (CSA), and weighted recombination
%% 
%% Note, this implementation is not optimized for vector operations in Matlab.
%% For large search space dimensionalities n > 500, it is also recommended to 
%% code this in C, C++, or Fortran
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I. Loshchilov, T. Glasmachers, and H.-G. Beyer
%% "Large Scale Black-Box Optimization by Limited-Memory Matrix Adaptation" 
%% IEEE Transactions on Evolutionary Computation, vol. 23, no. 2, pp. 353-358
%% DOI: 10.1109/TEVC.2018.2855049 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% n - search space dimensionality of y-vector in f(y) -> opt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% returns:
%% y_opt - approximation of the optimizer
%% f_dyn - the dynamics of the parental f-values
%% sigma_dyn - dynamics of the mutation strength sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% n -             search space dimensionality of y-vector in f(y) -> opt.
%% lambda -        offspring population size: lambda = 4 + floor(3*log(n))
%% mu -            parental population size: default size: lambda/2
%% gamma -         number of path vectors used: default size gamma = lambda
%% goal_f_name -   name of objective function f to be optimized
%% y_init -        initial y-vector (start vector)
%% sigma_init -    initial mutation stregth sigma
%% stepsize_stop - if parental change in search space is smaller then stop
%% f_stop -        stop if f gets smaller (minimize) or larger (maximization)
%% g_stop -        maximum number of generations until termination
%% opt -           "maximization" or "minimization"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y_opt, f_dyn, sigma_dyn,info] = ...
  LM_MA_ES( fun,y_init,st,tune)

 n=length(y_init);
 
 if ~exist('tune'), tune=[]; end;

if ~isfield(tune,'lambda')
   tune.lambda = 4 + floor(3*log(n));
end

lambda = tune.lambda;

if ~isfield(tune,'mu')
   tune.mu = floor(tune.lambda/2);
end

mu=tune.mu;

if ~isfield(tune,'gamma')
   tune.gamma = tune.lambda;
end

gamma=tune.gamma;

if ~isfield(tune,'sigma_init')
   tune.sigma_init = 1;
end

sigma_init=tune.sigma_init;

if ~isfield(tune,'stepsize_stop')
   tune.stepsize_stop = 1e-10;
end
stepsize_stop=tune.stepsize_stop;
if ~isfield(tune,'opt')
   tune.opt = 'minimization';
end

opt=tune.opt;
  

       
if ~exist('st'), st=[]; end;


if isfield(st,'prt'), info.prt = st.prt; 
else,  info.prt = -1; 
end;
% stopping criteria
if isfield(st,'secmax'), info.secmax=st.secmax;
else, info.secmax=200*n;
end;
if isfield(st,'nfmax'), info.nfmax=st.nfmax;
else, info.nfmax=200*n;
end;
if isfield(st,'fbest'), info.fbest=st.fbest;
else, info.fbest=-inf; 
end;

if isfield(st,'accf'), info.accf=st.accf;
else, info.accf=-inf;
end;
if isfield(st,'finit'), info.finit=st.finit;
else, info.finit=fun(y_init);
end;



if isfield(st,'reallife'), info.reallife=st.reallife;
else, info.reallife=0;
end;




info.initTime=cputime;  info.nf=0;


if info.reallife, info.arrayqf = 1; info.arraynf = 0; end
     


  dim = length(y_init);
  wi_raw = log(lambda/2 + 0.5) - log((1:mu));
  wi = wi_raw/sum(wi_raw);
  mu_eff = 1/sum(wi .^2);
  c_s = 2*lambda/dim;
  sqrt_s = sqrt(c_s*(2-c_s)*mu_eff);
  cp = lambda/dim*4.^-(0:gamma-1);
  cd = 1.5.^-(0:gamma-1) / dim;
  sqrt_cp = sqrt(cp .* (2-cp) * mu_eff);
  Parent.y = y_init;
  sigma = sigma_init;
  s = ones(dim, 1);
  if strcmp(opt, 'maximization')
    ordering = 'descend';
    maxi = 1;
    f_bsf = -1e300;
  else
    ordering = 'ascend';
    maxi = 0;
    f_bsf = 1e300;
  end
  f_dyn = [];
  sigma_dyn = [];
  g = 0;
  p_matrix = zeros(dim, gamma);
  while(1)
    clear OffspringPop;
    for l=1:lambda
      Offspring.z = randn(dim, 1);
      d = Offspring.z;
      for j = 1:min(g, gamma)
        d = (1-cd(j))*d + (cd(j)*(p_matrix(:, j)'*d)) * p_matrix(:, j);
      end
      Offspring.d = d;
      [Offspring.f,ftrue] = fun( Parent.y + sigma*Offspring.d);
      OffspringPop(l) = Offspring;
      info.nf=info.nf+1;
      % check stopping test
      sec       = (cputime-info.initTime);
      info.done = (sec>info.secmax)|(info.nf>=info.nfmax);
      info.qf   = abs((ftrue-info.fbest)/(info.finit-info.fbest));
      info.done = (info.done|info.qf<=info.accf);
      info.sec  = sec;
      
      if info.reallife
           if  info.qf>0
              info.arrayqf = [info.arrayqf   min(1,info.qf)];
              info.arraynf = [info.arraynf   info.nf];
         else
              info.arrayqf = [info.arrayqf   1];
              info.arraynf = [info.arraynf   info.nf];
         end
      end
      
      if info.done, break; end
      
    end
    if info.done, break; end
    ranks = RankPop(OffspringPop, ordering);
    sum_z = zeros(dim, 1);
    sum_d = zeros(dim, 1);
    for m = 1:mu
      sum_z = sum_z + wi(m)*OffspringPop(ranks(m)).z;
      sum_d = sum_d + wi(m)*OffspringPop(ranks(m)).d;
    end
    Parent.z = sum_z;
    Parent.d = sum_d;
    Parent.y = Parent.y + sigma*Parent.d;
    s = (1-c_s)*s + sqrt_s*Parent.z;
    for i=1:gamma
      p_matrix(:, i) = (1-cp(i))*p_matrix(:, i) + sqrt_cp(i)*Parent.z;
    end
    sigma = sigma*exp(0.5*c_s*(norm(s)^2/dim - 1));
    g = g+1;
% Statistics:
    [Parent.f, ftrue]= fun(Parent.y);
    f_dyn(g) = Parent.f; 
    sigma_dyn(g) = sigma;
    info.nf=info.nf+1;
       % check stopping test
      sec       = (cputime-info.initTime);
      info.done = (sec>info.secmax)|(info.nf>=info.nfmax);
      info.qf   = abs((ftrue-info.fbest)/(info.finit-info.fbest));
      info.done = (info.done|info.qf<=info.accf);
      info.sec  = sec;
      
      if info.reallife
         if  info.qf>0
              info.arrayqf = [info.arrayqf   min(1,info.qf)];
              info.arraynf = [info.arraynf   info.nf];
         else
              info.arrayqf = [info.arrayqf   1];
              info.arraynf = [info.arraynf   info.nf];
         end
      end
      
      if info.done, break; end
% Termination conditions:
    if ( maxi )
      if (Parent.f > f_bsf)
        f_bsf = Parent.f;
        y_opt = Parent.y;
      end
      if ( Parent.f > f_stop ...
        || g > info.nfmax ...
        || norm(sigma*Parent.d) < stepsize_stop ...
        || ( (g > 10) && abs(f_dyn(g)-f_dyn(g-10)) <= 100*eps(abs(f_dyn(g))) ))
       break;
      end
    else
      if (Parent.f < f_bsf)
        f_bsf = Parent.f;
        y_opt = Parent.y;
      end
      if (Parent.f < info.fbest ...
        || g > info.nfmax ...
        || norm(sigma*Parent.d) < stepsize_stop ...
        || ( (g > 10) && abs(f_dyn(g)-f_dyn(g-10)) <= 100*eps(abs(f_dyn(g))) )) 
        break;
      end
    end
  end
  
end
