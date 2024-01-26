clear all

freq = 500
x = linspace(0,10,10*16000);
y = sin(x*2*pi*500);

plot(x,y)

sound(y,16000)