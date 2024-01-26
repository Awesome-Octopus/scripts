%% SRT Optimization Function
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.13
% Last Modified 02/21/2011

function f = t1opt(x,b,spec,params); % Sets up the t1opt function (goes with t1fit)
opt = (x(1)*(1-2*exp(-b/x(2)))+x(3)).*(1./b).^0.5; % Signal intensity as a function of SRT
[maxopt index] = max(opt); % Maximum value of the T1 optimization function
f = b(index); % Time point of the maximum value of the T1 optimization function
plot(b,opt)
