clear all

[x y params] = eprload();
x = x';

signal = max(y(420:615))-min(y(420:615))
noise = mean((y(93:247) - mean(y(93:247))).^2)^0.5
SN = signal/noise

plot(x,y)