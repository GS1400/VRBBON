%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copyright by Hans-Georg Beyer (HGB)
%% For teaching use only! It is not allowed to use 
%% this program without written permission by HGB 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = diffPowers(y)
  nm1 = length(y)-1;
  f = sum( (y.^2).^(1 + 5/nm1*(0:nm1)') );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
