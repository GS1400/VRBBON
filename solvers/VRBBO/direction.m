
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% direction.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [step,par] = direction(point,step,par,tune,info)
% generate coordinate, quasi-Newton, subspace, random and
% cumulative directions
%
% for details of input and output structures see VRBBO.m
%
function [step,par] = direction(point,step,par,tune)
switch par.dir
            
   case 1 % coordinate direction
        step.p  = 1;  
        scale   = 0;  
   case 2 % quasi Newton direction
        if point.mq==0
           step.p = -point.grad;
           scale   = 1;
        else
           step     = lbfgsDir(point,step);
           sigma    = point.grad'*step.p;
          if sigma>0
             % move away from maximizer or saddle point
            step.p=-step.p;sigma=-sigma;
          end;
          sigma1  =  point.grad'* point.grad;
          sigma2  = step.p'*step.p;
          sigma12 = sigma1*sigma2;
          c       = sigma/sqrt(sigma12);
          if c<-tune.gammaangle &~isnan(c) & isfinite(c) & isreal(c)
              scale   = 0;
          else
            w=(sigma12*max(tune.gammaw,1-c^2))/(1-tune.gammaangle^2);
            t=(sigma+tune.gammaangle*sqrt(w))/sigma1;
            ok = (w>0 & isfinite(t)& t>0 & ~isnan(w)&~isnan(t));
            if ok
              step.p  =  step.p-t* point.grad;
              scale   =  0;
            else 
              step.p  = - point.grad/norm(point.grad);
              scale   = 0;
            end;
         end; 

        end
   case 3 % random subspace direction
        alpha = rand(point.ms,1)-0.5;
        alpha=alpha/norm(alpha);
        for i=1:point.ms
            dX(:,i) = alpha(i)*(point.X(:,i)-point.xm);
        end
        step.p=sum(dX')'; 
        scale = tune.scSub; 
  
   case 4 % random direction
        step.p = rand(point.n,1)-0.5;
        step.p(step.p==0) = 0.5; 
        scale = 1; 
   case 5 % cumulative direction
        step.p=step.q; scale = tune.scCum;
end
% if step is null or contaminated by nan or inf,
% it is replaced by a random direction
ok = (norm(step.p)==0|isnan(norm(step.p))|isinf(norm(step.p)));
if ok, step.p=rand(point.n,1)-0.5; end
if scale % scale to ||p||=delta
   
     step.delta = max(step.deltamin,...
     min(sqrt(par.A(par.t)*tune.gammadelta*step.Delta/par.lam),...
     step.deltamax))/par.A(par.t);
    
     step.p  = step.s.*step.p; 
     step.p  = step.p*(step.delta/norm(step.p));  
     step.dp = step.delta;  
else
    step.dp = norm(step.p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
