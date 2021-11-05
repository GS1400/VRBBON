

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% driverVRBBON.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:3
    for j=1:100
        fprintf('=')
    end
    fprintf('\n')
end

fprintf(['VSBBON solves a unconstrained noisy black box',...
    ' optimization of a not necessarily smooth function of \n',...
    ' many continuous arguments. No derivatives are needed.',...
    ' A limited amount of noise is tolerated.\n\n']);
disp('===============================================================')

clear;close all


solverPath = input(['Insert the VSBBON path \n',...
    '>> solverPath='],'s');


if ~exist(solverPath, 'dir')
    disp('the directory does not exist')
    return
end





disp('===============================================================')



 
n=5000;

 

% To solve your own problem, simply replace in the previous line 
% the expression after @(x) by your expression or function call. 
% Parameters in this expression are passed by value, hence are 
% fixed during minimization.

% start and stop info
x0      = rand(n,1)-0.5;  % starting point

% problem definition complete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('===============================================================')


solverlist = {'VRBBON' 'VRBBO' 'LMMAES'  };

problemlist = {'Driven_cavity' };

problemlistfig = { 'Driven cavity' };

levellist = [0.0001 0.0001 0.001];

levelfig = {'low' 'medium' 'large'};

     
color     = {'b' 'r' 'g'};
lines     = {'-' ':' '-.' };
mark      = {'<' '+' '*'};
linewid   = [0.1,0.1,0.1];   

 eval(['addpath ',solverPath,'/fun'])
 eval(['addpath ',solverPath,'/RealLife/testproblems'])
 

 
for pr=1:length(problemlist)
 
for nl=1:length(levellist) 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create noise

    fprintf(['noise.level: 0.0001/0.01/0.1 \n',...
             'noise.type:  1 (absolute) or  2 (relative)\n',...
             'noise.distr: 1 (uniform)  or  2 (Gauss)\n'])

    noise = struct('noisefun',1,'level',levellist(nl),'type',1,'distr',1);

     
     % create real-life objective function
    fun0  = str2num(['@(x)',problemlist{pr},'(x)']);
    fun  = @(x) funf(x,fun0,noise);

     st = struct('secmax',inf,'nfmax',200*n,'finit',fun(x0),...
     'fbest',0.01,'accf',0.01,'prt',-1,'reallife',1);
    tune=[];
     
    for so=1:length(solverlist)
       
        switch solverlist{so}
            case 'VRBBON'
                eval(['addpath ',solverPath,'/solvers/VRBBON'])
                eval(['addpath ',solverPath,'/solvers/VRBBON/minq8'])
                [~,~,info] = VRBBON(fun,x0,st,tune);
            case 'VRBBO'
                eval(['addpath ',solverPath,'/solvers/VRBBO'])
                [~,~,info] = VRBBO(fun,x0,st,tune);
            case 'LMMAES'
                  eval(['addpath ',solverPath,'/solvers/LMMAES'])
                  [~,~,~,info] = LM_MA_ES(fun,x0,st,tune);
        end
                
       
            
        
        [arrayqf,I] = sort(info.arrayqf,'descend');
        arraynf     = info.arraynf;

        loglog(arraynf,arrayqf,[lines{so},color{so},mark{so}],...
                 'LineWidth',linewid(so),'markersize',7);
             
         legendname{so} = solverlist{so};
      
         hold on
         
         
       psid=[problemlist{pr},levelfig{nl},'.',solverlist{pr}] 
       %eval([psid,'.TotalTime=showtime0;']) 
       eval([psid,'.arrayqf=arrayqf;']);     
       eval([psid,'.arraynf=arraynf;']);     
   
         
        
         
         
         
         
         
    end
    
    com=['save ', solverPath,'/data/',[problemlist{pr},levelfig{nl}],'.mat;'];
    eval(com)
    
   
    
    xxtrict =  [2 5 10 20 50 200 1e3 1e4 1e5 1e6 1e7];
    
    yytrict = [0.00001 0.00003 0.00005 0.0001 0.0003 ...
        0.0005 0.001 0.003 0.005 0.01 0.03 0.05 0.1 0.3 0.5 1];
    
    I = find(yytrict>=st.accf); yytrict = yytrict(I);
    
    
    lgd = legend(legendname,'Location','southwest');
    lgd.FontSize = 20;
    lgd.FontWeight = 'bold';
    set(gca, 'XTick',xxtrict)
    set(gca, 'YTick', yytrict)
    
    xlabel('{\bf nf}','fontsize',14)
    ylabel('{\bf qs}','fontsize',14)
    title(['Problem = ',problemlistfig{pr},', n = ',...
            num2str(n)],'fontsize',10)
   com=['print -depsc2 ',solverPath,'/eps/',[problemlist{pr},...
         levelfig{nl}],'.eps'];
   disp(com)
   eval(com)
   com=['print -dpng ',solverPath,'/eps/',[problemlist{pr},...
         levelfig{nl}],'.png'];
   disp(com)
   eval(com)    
   close
   
   
    
end



end


for i=1:3
    for j=1:100
        fprintf('=')
    end
    fprintf('\n')
end
