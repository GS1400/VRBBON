
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adjustX.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point]=adjustX(point,tune);
% remove NaN or inf from X
%
% for details of input and output structures see VRBBON.m
%

function [point]=adjustX(point,tune);

for i=1:size(point.X,1)
    Xi        = point.X(i,:);
    Inan      = isnan(Xi);
    if any(Inan), Xi(Inan) = tune.gammaX; point.X(i,:)=Xi; end
    Iinf= (Xi==inf);
    if any(Iinf), Xi(Iinf) = tune.gammaX; point.X(i,:)=Xi; end
    Iinf= (Xi==-inf);
    if any(Iinf), Xi(Iinf) = -tune.gammaX; point.X(i,:)=Xi; end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

