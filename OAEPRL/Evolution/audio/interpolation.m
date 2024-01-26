clear all

freq = 500
x = linspace(0,10,21);
y = sin(x);

xi = linspace(0,10,101)
yi = interp1(x,y,xi)

subplot(2,1,1)
plot(x,y)
subplot(2,1,2)
plot(xi,yi)