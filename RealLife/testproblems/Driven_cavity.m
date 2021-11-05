%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Driven cavity 
% written by Morteza Kimiaei Oct. 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function fk = Driven_cavity(x)

n=length(x);
if n<16, error('n must be >=16'); end

M=floor(sqrt((n)));
PAR=500/(M+2)^4;
n=M*M; fk=0;
for k=1:n
      H=0.5/(M+2);
      J=floor((k-1)/M)+1;
      I=k-(J-1)*M;
      FA=20*x(k);
      A1=0;
      A2=0;
      A3=0;
      A4=0;
      if(J>2) 
        FA=FA+x(k-M-M);
        A4=A4+x(k-M-M);
      end
      if(J>1) 
         if (I>1)
             FA=FA+2*x(k-M-1);
             A3=A3+x(k-M-1);
             A4=A4+x(k-M-1);
         end
         FA=FA-8*x(k-M);
         A1=A1-x(k-M);
         A4=A4-4*x(k-M);
         if(I<M)
          FA=FA+2*x(k-M+1);
          A3=A3-x(k-M+1);
          A4=A4+x(k-M+1);
         end
      end
      if(J>1) 
         if (I>2)
           FA=FA+x(k-2);
            A3=A3+x(k-2);
         end
         FA=FA-8*x(k-1);
         A2=A2-x(k-1);
         A3=A3-4*x(k-1);
      end
      if(I<M) 
          FA=FA-8*x(k+1);
          A2=A2+x(k+1);
          A3=A3+4*x(k+1);
         if (I<M-1)
            FA=FA+x(k+2);
            A3=A3-x(k+2);
         end
      end
      if(J<M) 
         if (I>1)
           FA=FA+2*x(k+M-1);
           A3=A3+x(k+M-1);
           A4=A4-x(k+M-1);
         end
         FA=FA-8*x(k+M);
         A1=A1+x(k+M);
         A4=A4+4*x(k+M);
         if (I>M)
            FA=FA+2*x(k+M+1);
            A3=A3-x(k+M+1);
            A4=A4-x(k+M+1);
         end   
      end
      if (J<M-1)
        FA=FA+x(k+M+M);
        A4=A4-x(k+M+M);
      end
      if (J==M)
        if (I>1)
          FA=FA-H-H;
          A3=A3-H;
          A4=A4+H;
        end
        FA=FA+8*H;
        A1=A1-H;
        A4=A4-4*H;
        if (I<M)
          FA=FA-2*H;
          A3=A3+H;
          A4=A4+H;
        end
        FA=FA+H;
        A4=A4-H;
      end
      if (J==M-1)
        FA=FA-H;
        A4=A4+H;
      end
      FA=FA+0.25*PAR*(A1*A3-A2*A4);
       fk=fk+FA^2;
end  
fk=0.5*fk;
% ---------------------------------------------------
end % end of function
% --------------------------------------------------