%% T1 Decay Function
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.13
% Last Modified 02/21/2011

function f = t1decay(x,b,spec,params) % Sets up the t1decay function (goes with t1fit)
f = ((sum(abs(spec-(x(1)*(1-2*exp(-b/x(2)))+x(3))))).^2/params.XPTS).^0.5; % The T1 decay function
