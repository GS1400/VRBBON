
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% direction.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step,par] = direction(point,step,par,tune,info)
% generate random direction, random subspace direction, 
% trust-region or random scaled gradient direction
%
% for details of input and output structures see VRBBON.m
%
function [point,step,par] = direction(point,step,par,tune,info)
switch par.dir
    
    
    case 0 % scaled random direction
        
        step.p = rand(point.n,1)-0.5;
        step.p = step.p/norm(step.p);
    
    case 1 % random approximate coordinate directions
        
        step.p = tune.gammar*(rand(point.n,1)-0.5);
        if par.t>point.n
           t  = mod( par.t,point.n);
        else
            t = par.t;
        end
        step.p(par.Iperm(t)) = 1 ;  
        step.p = step.p/norm(step.p);
        
   case 2 % random subspace direction
       
        alpha = rand(point.m,1)-0.5;
        alpha=alpha/norm(alpha);
        for i=1:point.m
            dX(:,i) = alpha(i)*(point.X(:,i)-point.xm);
        end
        step.p = sum(dX')';  
     
   case 3 % direction obtained by solving the trust-region subproblem

           [A,D]    = eig(point.G);
           diagD     = diag(D);
           ind = diagD==0;
           if ~isempty(ind)
               [Dmax]=max(diagD);
               diagD(ind)=Dmax;
           end

           data.gam = point.fm;
           data.c   = point.g;
           data.D   = diagD'; % positive D
           data.b   = zeros(point.msim,1);
           data.A   = A;
           xl       = point.xxm'-step.d;
           xu       = point.xxm'+step.d;
           
           
           % call minq8
           warning off
           [zeta,~,~,~,~] = minq8 ...
                (data,xl,xu,point.xxm',tune.minqmax,tune.minqeps,0);
           warning on
           
           p              = tune.gammap*(zeta-point.xxm');
           step.p         = zeros(point.n,1); 
           I              = randperm(point.n,point.msim);
           step.p(I)      = p;
           xc             = mean(point.X')';
           step.p         = xc-point.xm+step.p;
           
    case 4 % random scaled gradient direction in subspace
        

           beta      = 1/(1+info.nf)^tune.gammak ;
           p         = rand(point.msim,1)-0.5;       
           alp       = (1+beta*point.g'*p)/(point.g'*point.g);
           p         = beta*p-alp*point.g;
           step.p    = zeros(point.n,1); 
           I         = randperm(point.n,point.msim);
           step.p(I) = p;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
