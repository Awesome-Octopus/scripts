%% Power Saturation Linear Fit
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.6
% Last Modified 8/18/2020

function f = linearfit(x,sqrtpower,height) % Establishes the powerfit function
for i = 1:4; % for loop structure to create the fit
fit(i) = x(1)*sqrtpower(i)+x(2);
end; % ends the for loop
f = (((sum(abs(height(1:4)-fit))).^2)/4)^0.5; % root mean squared deviation between data and fit
