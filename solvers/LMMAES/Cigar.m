%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copyright by Hans-Georg Beyer (HGB)
%% For teaching use only! It is not allowed to use 
%% this program without written permission by HGB 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = Cigar(y)
  f = y(1)*y(1) + 1e6*sum( y(2:end).^2 );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
