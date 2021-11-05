%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% identifyDir.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identify kind of direction
%
% for details of input and output structures see VRBBO.m
%
function [par,step,info]  = identifyDir(point,step,par,tune,info)
switch    tune.alg
             
    case 2  % no subspace
       
        if par.t<=tune.C
          par.dir       = 1;   % pick coordinate direction
          if info.prt
            if par.t==1
               info.stepType      = '[Coo|'; % coordinate direction
            else
               info.stepType      = 'Coo|'; % coordinate direction 
            end
            info.nstep.nCoo = info.nstep.nCoo+1; 
          end
        elseif par.t== tune.C+1 % pick quasi Newton direction
              par.dir = 2;   % pick quasi Newton direction
              if info.prt
                 info.stepType      = 'New|'; % quasi Newton direction
                 info.nstep.nNew = info.nstep.nNew+1; 
              end
        elseif par.t<=tune.C +tune.R+1 
            
              par.dir   = 4;   % pick random direction
              if info.prt
              info.stepType  = 'Ran|'; % Random direction
              info.nstep.nRan = info.nstep.nRan+1;
              end
        else
          if tune.cum==2 & step.r>=step.Delta % second cum step
          elseif tune.cum==1, step.q=point.xm-point.xinit; % first cum step
          end
          if tune.cum>0 & all(step.q==0)
              par.dir=4; % pick random step
             if info.prt 
               info.stepType  = 'Sub|'; % subspace direction
               info.nstep.nSub= info.nstep.nSub+1; 
             end
          else
            par.dir = 5; % pick cumulative step

                           % use initial value
             if  info.prt             
                info.stepType  = 'Cum|'; % cumulstive direction
                info.nstep.nCum = info.nstep.nCum+1; 
             end
          end
        end
       
    case 3 % none of quasi Newton and  coordinate diections
       
        if par.t<=tune.R
         par.dir   = 4;   % pick random direction
          if info.prt
          if par.t==1    
              info.stepType  = '[Ran|'; % Random direction
           else
              info.stepType  = 'Ran|'; % Random direction
           end
          info.nstep.nRan = info.nstep.nRan+1;
          end
        
        elseif par.t<=tune.R+tune.S & point.ms>=3
           par.dir       = 3;   % pick random subspace direction
          if info.prt
            info.stepType      = 'Sub|'; % random subspace direction
            info.nstep.nSub = info.nstep.nSub+1; 
          end  
    
        else
              if tune.cum==2 & step.r>=step.Delta % second cum step
              elseif tune.cum==1, 
                  step.q=point.xm-point.xinit; % first cum step
              end
              if tune.cum>0 & all(step.q==0)
                  par.dir=4; % pick random  step
                 if info.prt 
                   info.stepType  = 'Sub|'; % subspace direction
                   info.nstep.nSub= info.nstep.nSub+1; 
                 end
              else
                par.dir = 5; % pick cumulative step

                               % use initial value
                 if  info.prt             
                    info.stepType  = 'Cum|'; % cumulstive direction
                    info.nstep.nCum = info.nstep.nCum+1; 
                 end
              end
        end
    case 4 % no quasi Newton
    
        if par.t<=tune.C
          par.dir       = 1;   % pick coordinate direction
          if info.prt
            if par.t==1
               info.stepType      = '[Coo|'; % coordinate direction
            else
               info.stepType      = 'Coo|'; % coordinate direction 
            end
            info.nstep.nCoo = info.nstep.nCoo+1; 
          end
       
        elseif par.t<=tune.C +tune.S & point.ms>=3
           par.dir       = 3;   % pick random subspace direction
          if info.prt
            info.stepType      = 'Sub|'; % random subspace direction
            info.nstep.nSub = info.nstep.nSub+1; 
          end  
        elseif par.t<par.T 
              par.dir   = 4;   % pick random direction
              if info.prt
              info.stepType  = 'Ran|'; % Random direction
              info.nstep.nRan = info.nstep.nRan+1;
              end
        else
          if tune.cum==2 & step.r>=step.Delta % second cum step
          elseif tune.cum==1, step.q=point.xm-point.xinit; % first cum step
          end
          if tune.cum>0 & all(step.q==0)
              par.dir=4; % pick random step
             if info.prt 
               info.stepType  = 'Sub|'; % subspace direction
               info.nstep.nSub= info.nstep.nSub+1; 
             end
          else
            par.dir = 5; % pick cumulative step

                           % use initial value
             if  info.prt             
                info.stepType  = 'Cum|'; % cumulstive direction
                info.nstep.nCum = info.nstep.nCum+1; 
             end
          end
        end
       
    
    case 5 % all directions
      
   if par.t<=tune.C
      par.dir       = 1;   % pick coordinate direction
      if info.prt
        if par.t==1
           info.stepType      = '[Coo|'; % coordinate direction
        else
           info.stepType      = 'Coo|'; % coordinate direction 
        end
        info.nstep.nCoo = info.nstep.nCoo+1; 
      end
    elseif par.t== tune.C+1 % pick quasi Newton direction
          par.dir = 2;   % pick quasi Newton direction
          if info.prt
             info.stepType      = 'New|'; % quasi Newton direction
             info.nstep.nNew = info.nstep.nNew+1; 
          end
    elseif par.t<=tune.C +tune.S+1 & point.ms>=3
       par.dir       = 3;   % pick random subspace direction
      if info.prt
        info.stepType      = 'Sub|'; % random subspace direction
        info.nstep.nSub = info.nstep.nSub+1; 
      end  
    elseif par.t<par.T 
          par.dir   = 4;   % pick random direction
          if info.prt
          info.stepType  = 'Ran|'; % Random direction
          info.nstep.nRan = info.nstep.nRan+1;
          end
    else
      if tune.cum==2 & step.r>=step.Delta % second cum step
      elseif tune.cum==1, step.q=point.xm-point.xinit; % first cum step
      end
      if tune.cum>0 & all(step.q==0)
          par.dir=4; % pick random step
         if info.prt 
           info.stepType  = 'Sub|'; % subspace direction
           info.nstep.nSub= info.nstep.nSub+1; 
         end
      else
        par.dir = 5; % pick cumulative step

                       % use initial value
         if  info.prt             
            info.stepType  = 'Cum|'; % cumulstive direction
            info.nstep.nCum = info.nstep.nCum+1; 
         end
      end
   end
    otherwise
          par.dir   = 4;   % pick random direction
          if info.prt
           if par.t==1    
              info.stepType  = '[Ran|'; % Random direction
           else
              info.stepType  = 'Ran|'; % Random direction
           end
          info.nstep.nRan = info.nstep.nRan+1;
          end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%