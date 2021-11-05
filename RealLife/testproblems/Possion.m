%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Possion
% written by Morteza Kimiaei Oct. 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function fk = Possion(x)

n=length(x);
if n<16, error('n must be >=16'); end

M=floor(sqrt((n)));
PAR=1/(M+1)^2;
n=M*M; fk=0;
for k=1:n
      FA=4*x(k);
      J=floor((k-1)/M)+1;
      I=k-(J-1)*M;
      FA=FA+PAR*x(k)^3/(1+PAR*(I)^2+PAR*(J)^2);
      if(I==1) 
        FA=FA-1;
      end
      if(I>1) 
        FA=FA-x(k-1);
      end
      if(I<M) 
         FA=FA-x(k+1);
      end
      if(I==M) 
         FA=FA-2+exp((J)/(M+1));
      end
      if(J==1) 
         FA=FA-1;
      end
      if(J>1) 
        FA=FA-x(k-M);
      end
      if(J<M)
         FA=FA-x(k+M);
      end
      if(J==M) 
         FA=FA-2+exp((I)/(M+1));
      end
      
       fk=fk+FA^2;
end  
fk=0.5*fk;
% ---------------------------------------------------
end % end of function
% --------------------------------------------------