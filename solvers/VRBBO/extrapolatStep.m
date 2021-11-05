%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% extrapolatStep %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform an extrapolation step to hopefully get a good decrease
% in the function value
%
% for details of input and output structures see VRBBO.m
%

function [point,step,par,info]=extrapolatStep(fun,point,step,par,...
                                   tune,info)

info.nE=0; par.alphaE = par.A(par.t); opp=0; prt=(info.prt>=1); 

    
while 1
     
       
     if prt, info.stateTest=[info.stateTest,info.stepType]; end;
     if par.dir==1
        point.xr(par.t) = point.xm(par.t)+par.alphaE*step.p;
     else
        point.xr = point.xm+par.alphaE*step.p;
     end
     
      [point.fr,info.ftrue] = fun(point.xr);
      
      if isnan(point.fr), point.fr=inf; end;
      info.nf   = info.nf+1; 
      point.fr  = max(-1e50,min(point.fr,1e50));
      
      
      
      % check stopping test
      sec       = (cputime-info.initTime);
      info.done = (sec>info.secmax)|(info.nf>=info.nfmax);
      info.qf   = (info.ftrue-info.fbest)/(info.finit-info.fbest);
      info.done = info.done|(info.qf<=info.accf);
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
      
     if info.done, point.fm = point.fr; point.xm=point.xr;break; end

      df = point.fm-point.fr;
     
     if (opp==1)
        % estimate the Lipschitz constant
      
        par.lam = max(par.lam,abs(point.fl+point.fr-2*point.fm)...
                  /(step.dp)^2);
              
       % update r and q
       if tune.cum ==2, [step] = updatecum(point,step,tune); end
       
     end
     
     if info.nE==1
         fm = point.fm;
         point.fm=point.fext;
         par.lam = max(par.lam,abs(point.fl+point.fr-2*point.fm)...
                  /(step.dp)^2);
         % update r and q
         if tune.cum ==2, [step] = updatecum(point,step,tune); end
        point.fm=fm;
            
     end
     
    ext = (df>tune.gammamin*step.Delta & info.nE<tune.E );
   
    
    if ext % extrapolation is tried
        info.nE=info.nE+1;
        point.fext = point.fr; 
        if info.nE==1
           point.fl = point.fm;
        else
             % expand step size
             par.alphaE = tune.gammaE*par.alphaE;  
        end
                         
        if prt
          info.stepType  = 'Ext|'; % extrapolate
          info.nstep.nExt = info.nstep.nExt+1;
        end 
    elseif (opp==0&info.nE==0) % opposite direction is tried
     point.fl   = point.fr;
     step.p     = -step.p;
     opp = 1;
     if prt
        info.stepType  = 'Opp|'; % opposite direction
        info.nstep.nOpp = info.nstep.nOpp+1;
     end
    elseif info.nE>=1
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % the end of extrapolation %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           par.alphaE    = par.alphaE/tune.gammaE;
           par.A(par.t) = par.alphaE;
           
           if par.dir==1
              point.xm(par.t) = point.xm(par.t) + par.alphaE*step.p;
              point.xr(par.t) = point.xm(par.t);
           else
              point.xm = point.xm+ par.alphaE*step.p;
           end
           
           point.fm = point.fext; par.good=1;
           
            if prt
               disp(['function value improved at nf=',...
               num2str(info.nf),' to f=',num2str(point.fm)]) 
            end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           break; % extrapolation step ends up
    else % no decrease in f; both p and -p were tried
        break; % extrapolation step ends up
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%