
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MLS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,par,step,info] = MLS(fun,point,step,par,tune,info)
%
% try to significantly improve function value
%
% for details of input and output structures see VRBBON.m

function [point,par,step,info] = MLS(fun,point,step,par,tune,info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step.alist=[]; point.dlist=[]; 


if isfield(tune,'C')
  par.Iperm = randperm(point.n);
end


par.ngood=0;

% find a good initial step size
[step]=initStep(step,par,tune);

while 1
  switch par.state
  case -1 % new direction
    par.t=par.t+1; nE=0;  
        
    [point,step,par] = direction(point,step,par,tune,info);
    
  case 0  % extrapolate stage
      
    nE=nE+1;
    
    point.fext = point.fr; 
    
    
    extF(nE)   = point.fr; extalp(nE) =  step.alphaE;
      
    step.alphaE = step.alphaE*tune.gammaE;
    
   
  case 1 % opposite direction

     step.p = -step.p; nE = 0;
  end
      point.xr = point.xm+step.alphaE*step.p;
      [point.fr,info.ftrue] = fun(point.xr);
      if isnan(point.fr), point.fr=inf; end;
      info.nf   = info.nf+1; 
      point.fr  = max(-1e50,min(point.fr,1e50));
     
      % check stopping test
      sec       = (cputime-info.initTime);
      info.done = (sec>info.secmax)|(info.nf>=info.nfmax);
      info.qf   = abs((info.ftrue-info.fbest)/(info.finit-info.fbest));
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
     
      
      if info.done
          point.fm = point.fr; point.xm=point.xr;
          par.state = -1;
          break; 
      end

    par.df = point.fm-point.fr;

    ext = (par.df>=tune.gamma*step.alphaE^2 & ...
           nE <= tune.E); 
       
    if ext, par.state=0; continue;
    else
        opp = (par.state<0);
        if opp, par.state=1; 
            if info.prt>=3, disp('opposite stage used'); end  

            continue;
        end
        if par.state==0 
            
           if info.prt>=3
                disp(['extrapolation stage performed at nf=',...
                      num2str(info.nf)])
           end       
            
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % the end of extrapolation %
          %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
          
          [fb,ind] = min(extF); ind =ind(1);
                    
           step.alphaE  = extalp(ind); 
           point.xm     =  point.xm+step.alphaE*step.p;        
           point.fm     =  fb;   
                      
           % update X F Y
           [point] = updateXFY(point,step,tune); 
           
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if par.desent==1
              [point,step]=initInterval(point,step); 
              if info.prt>=3, disp('initial interval found'); end  
              par.desent =0;
           else
              [step]=updateInterval(step,par);
              if info.prt>=3, disp('interval updated'); end  
           end
           
           par.ngood = par.ngood+1;
            
        else
            
           if info.prt>=3
               disp('there was no decease in f; reduce step size')
           end  
           
           if par.df>0 
              point.xm  = point.xr;  
              point.fm  = point.fr;
              [point]   = updateXFY(point,step,tune);
              par.ngood = par.ngood+1;
           end
           
          [step]=reducedStepSize(step,tune); % reduce the step size
          
          [step]=updateInterval(step,par); % update interval
            
        end
                
        par.state = -1;
        
        if par.t==par.T, break; end
    
    end
    
   
   
end  % of while T>0



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
