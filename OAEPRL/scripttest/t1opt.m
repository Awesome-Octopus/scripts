%% SRT Optimization Function
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.05
% Last Modified 07/15

function f = t1opt(x,b,spec,params); % Sets up the t1opt function (goes with t1fit)
opt = (x(1)*(1-2*exp(-b/x(2)))+x(3)).*(1./b).^0.5; % Signal intensity as a function of SRT
[m n] = size(opt); % Size of the opt array
for i = 1:n; % sets up a loop over the values of opt
    if (b(i) > 510); % If opt is greater than 510 us (the minimum SRT on the ELEXSYS)
        newopt(i) = opt(i); % Creates a new array with only attainable values of the SRT
    end; % ends the if structure
end; % ends the for loop
[maxopt index] = max(newopt); % Maximum value of the T1 optimization function
f = b(index); % Time point of the maximum value of the T1 optimization function
plot(b,newopt)
