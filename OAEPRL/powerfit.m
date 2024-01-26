%% Power Saturation Function
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.6
% Last Modified 8/18/2020

function f = powerfit(x,power,sqrtpower,height) % Establishes the powerfit function
[m n] = size(sqrtpower); % gets the size of the sqrtpower array
for i = 1:n; % for loop structure to create the fit
% fit(i) = x(1)*sqrtpower(i)/((1+(power(i)/x(2)))^(0.5*x(3))); % the fit function
fit(i) = x(1)*sqrtpower(i)*(1 + ((2^(1/x(3))-1) * sqrtpower(i)^2)/x(2))^-x(3);
end; % ends the for loop
f = ((sum(abs(height-fit)).^2)/n)^0.5; % root mean squared deviation between data and fit
