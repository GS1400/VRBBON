%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convection-diffusion
% written by Morteza Kimiaei Oct. 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function fk = Convection_diffusion(x)

n=length(x);
if n<16, error('n must be >=16'); end

M=floor(sqrt((n)));
PAR=1/(M+1);
n=M*M; fk=0;
for k=1:n
      FA=4*x(k);
      A3=0;
      J=floor((k-1)/M)+1;
      I=k-(J-1)*M;
      A1=PAR*(I);
      A2=PAR*(J);
      FA=FA-2000*A1*A2*(1-A1)*(1-A2)*PAR^2;
      if(I>1) 
         FA=FA-x(k-1);
         A3=A3-x(k-1);
      end
      if(I<M) 
         FA=FA-x(k+1);
         A3=A3+x(k+1);
      end
      if(J>1) 
          FA=FA-x(k-M);
          A3=A3-x(k-M);
      end
      if(J<M) 
          FA=FA-x(k+M);
          A3=A3+x(k+M);
      end
      FA=FA+20*PAR*A3*x(k);
      
       fk=fk+FA^2;
end  
fk=0.5*fk;
% ---------------------------------------------------
end % end of function
% --------------------------------------------------