
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% scaleModel.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [point,step]=scaleModel(point,step,tune) 
% do scaling for the model to insure the affine invarriant of fitting
% procedure
%
% for details of input and output structures see VRBBON.m
%

function [point,step]=scaleModel(point,step,tune) 

xm = point.X(:,point.b)'; mmax = tune.mmax; X = point.X';


[point]=sizeSample(point);


if point.msim < point.n
  I   = randperm(point.n,point.msim);
  X   = X(:,I);
  xm  = xm(I);
end

K         = min(size(X,1)-1,mmax);
distance  = sum((X-ones(size(X,1),1)*xm).^2,2);
[~,ind]   = sort(distance);
ind       = ind(2:K+1);
S         = X(ind,:) -ones(K,1)*xm;
R         = triu(qr(S,0));
R         = R(1:point.msim,:);

warning off
L = inv(R)';
warning on
step.sc=sum((S*L').^2,2).^(3/2);

point.ind=ind; point.K=K; step.distance=distance;
point.xxm=xm; point.XX=X;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%