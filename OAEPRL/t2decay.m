%% T2 Decay Function
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.13
% Last Modified 02/21/2011

function f = t2decay(x,b,spec,params) % Sets up the t2decay function (goes with t2fit)
f = ((sum(abs(spec-((x(1)*exp(-b/x(2))))))).^2/params.XPTS).^0.5; % The T2 decay function
