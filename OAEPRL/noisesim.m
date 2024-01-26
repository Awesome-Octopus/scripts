clear all

x = linspace(0,10,1001);
y = exp(-((x-2.5).^2)/(2*0.4.^2))+rand(1,1001)-0.5;



loops = 1000;
snx = linspace(1,loops,loops);
sny = zeros(1,loops);
signal = zeros(1,loops);
noise = zeros(1,loops);

for i=1:loops
    y = y+exp(-((x-2.5).^2)/(2*0.4.^2))+rand(1,1001)-0.5;
    signal(i) = y(251);
    noise(i) = mean((y(501:1001) - mean(y(501:1001))).^2)^0.5;
    sny(i) = signal(i)/noise(i);
    subplot(2,2,1)
    plot(x,y)
    axis([0 10 -50 1100])
    title('Accumulated Data Through 1000 Scans')
    subplot(2,2,2)
    plot(snx,signal)
    axis([-10 1000 -50 1100])
    title('Signal as a Function of Scans')
    xlabel('number of scans')
    ylabel('signal intensity')
    subplot(2,2,3)
    plot(snx,noise)
    axis([-10 1000 -1 11])
    title('Noise Level as a Function of Scans')
    xlabel('number of scans')
    ylabel('noise level')
    subplot(2,2,4)
    plot(snx,sny)
    axis([-10 1000 0 110])
    title('S/N as a Function of Scans')
    xlabel('number of scans')
    ylabel('S/N')
pause(0.01)
end