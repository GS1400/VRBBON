
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DS.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,par,step,info]=DS(fun,point,step,par,tune,info)
% Repeatedly improve the function value
%
% for details of input and output structures see VRBBON.m
%

function [point,par,step,info]=DS(fun,point,step,par,tune,info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RLS with R-random direction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if info.prt>=1
   disp(' ')
   disp('***********************************************')
   disp('start of Phases I & II:')
   disp('***********************************************')
   disp(' ')
end
  


for ii = 1:tune.T0
    
    par.totalngood= 0;
    
    
  

     switch  tune.comBound
         case 2
             
           par.t=0; par.T=tune.C; par.dir = 1;  
           % MLS using random approximate coordinate directions
          [point,par,step,info] = MLS(fun,point,step,par,tune,info);
          if info.prt>=2 &  par.ngood>0
            disp(['fnoise=',num2str(point.fm),' improved at nf=',...
                       num2str(info.nf)]) 
          end
          
           if info.prt>=1
              disp(' ')
              disp('***********************************************')
              disp('end of Phase I:')
              disp('***********************************************')
              disp(' ')
           end
          
          
           if info.done, break; end
          
         case 0
            par.t=0; par.T = tune.R; par.dir = 0;  
           % MLS using scaled random  directions
          [point,par,step,info] = MLS(fun,point,step,par,tune,info);
          
           if info.prt>=2 &  par.ngood>0
             disp(['fnoise=',num2str(point.fm),' improved at nf=',...
                       num2str(info.nf)]) 
          end
          
           if info.prt>=1
              disp(' ') 
              disp('***********************************************')
              disp('end of Phase I:')
              disp('***********************************************')
              disp(' ')
           end
          
          
           if info.done, break; end
           
         case 1
           par.t=0; par.T=tune.C; par.dir = 1;  
           % MLS using random  coordunate directions
          [point,par,step,info] = MLS(fun,point,step,par,tune,info);  
          
          if info.prt>=2 &  par.ngood>0
             disp(['fnoise=',num2str(point.fm),' improved at nf=',...
                       num2str(info.nf)]) 
          end
          
           if info.done, 
               if info.prt>=1
                disp(' ') 
                disp('***********************************************')
                disp('end of Phase I:')
                disp('***********************************************')
                disp(' ')
              end
               break; 
           end
          
           par.t=0; par.T=tune.R; par.dir = 0;  
           % MLS using scaled random directions
          [point,par,step,info] = MLS(fun,point,step,par,tune,info);
          
          if info.prt>=1
              disp(' ') 
              disp('***********************************************')
              disp('end of Phase I:')
              disp('***********************************************')
              disp(' ')
           end
          
          if info.done, break; end
          
     end
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MLS with random subspace direction if m>=3 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if point.m>=3
       if info.prt>=1
         disp('start of MLS using subspace direction')
       end
       
       while 1
            par.t=0; par.T=1; 
            par.dir  = 2;   % pick random subspace direction
            % MLS using random subspace direction
            [point,par,step,info] = MLS(fun,point,step,par,tune,info);
            if info.prt>=2 & par.ngood>0
               disp(['fnoise=',num2str(point.fm),' improved at nf=',...
                       num2str(info.nf)]) 
            end
            par.totalngood =  par.totalngood+par.ngood;
            if info.done, break; end  
            ok = (par.ngood ==0);
            if ok, break; end
       end
       if info.prt>=1, disp('end of MLS using subspace direction'); end 

    end
    
    
    if info.done, break; end  
    
  
    if info.prt>=1
       disp(' ')
       disp('***********************************************')
       disp('end of Phases III')
       disp('***********************************************')
       disp(' ')
    end
   

    if ~info.done


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct a robust quadratic model by using %
    % the best points then solve trust region     %
    % subproblem by minq8.m                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if tune.model

    ok = (point.m>=5); % I checked 

    if ok
       info.model = info.model+1;

        if info.prt>=1
            disp(' ');
            disp('***********************************************')
            disp('start of Phase III'),
            disp('***********************************************')
            disp(' ');
        end

        par.ngoodold = par.ngood;

        % construct a reduced quadratic model
        [point,step] = fit(point,step,tune);

        if ~ point.randg
            info.nquadmodel = info.nquadmodel+1;
            if info.prt>=1
              fprintf(['fit is successful;\nthe gradient & Hessain', ...
                    ' are estimated\n'])
            end
            [step]=genRadius(point,step,tune); % generate the radius
            while 1
                par.t = 0; par.T =1; 
                par.dir  = 3;   % pick TR direction
                if info.prt>=1
                    disp('start of MLS using trust region direction')
                end 
                % MLS using random subspace direction
                [point,par,step,info] = ...
                        MLS(fun,point,step,par,tune,info);
                 par.totalngood =  par.totalngood+par.ngood;   
                if info.prt>=1 
                    disp('end of MLS using trust region direction')
                end     

                 if info.prt>=2 & par.ngood>0
                    disp(['fnoise=',num2str(point.fm),...
                          ' improved at nf=',num2str(info.nf)]) 
                 end  

                if info.done, break; end  
                ok = (par.ngood ==0);
                if ok, break; end

                step.d = (0.5+rand)*step.d;

            end
        else
            info.nlinearmodel = info.nlinearmodel+1;
            if info.prt>=1
              disp('fit is successful only for estimating the gradient')
            end

            while 1
                 par.t=0; par.T=1; 
                 par.dir = 4; % pick  random scaled gradient direction
                if info.prt>=1
                    disp('start of MLS using scaled gradient direction')
                end 
                % MLS using scaled gradient direction
                [point,par,step,info] = ...
                        MLS(fun,point,step,par,tune,info);
                 par.totalngood =  par.totalngood+par.ngood;   
                 if info.prt>=1
                    disp('end of MLS using scaled gradient direction')
                 end     

                if info.prt>=2 & par.ngood>0
                    disp(['fnoise=',num2str(point.fm),...
                          ' improved at nf=',num2str(info.nf)]) 
                end

                if info.done, break; end  
                ok = (par.ngood ==0);

                if ok, break; end

            end


        end
        

         if info.prt>=1
           disp(' ')
           disp('***********************************************')
           disp('end of Phase III')
           disp('***********************************************')
           disp(' ')
         end
         
         if info.done, break; end  
    end


    end % the model used

    end % end of Phase III
 
    
    
    % restart the interval
    res = (par.totalngood==0 & point.m>=3);
    % restate the interval since there is no improvement
    if res 
       if info.prt>=1
          disp('the interval has been restarted by resInterval')
       end
       [step]=resInterval(point,step,tune);
    end
  
end % for loop




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
