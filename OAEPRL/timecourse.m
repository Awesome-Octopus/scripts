clear all;
[x y params] = eprload('Sample1A.spc')
x = x'
ysmooth = smooth(y,10)
[m n] = size(y)

for i=1:params.SSY
int(1,i) = max(cumtrapz(x,cumtrapz(x,ysmooth(:,i))))
time(1,i) = params.JNS*params.SSX*params.RTC*0.001/60*i
end
int = int'
time = time'
plot(time,int)

data = [x y]
processed = [time int]

save('Sample1A.txt','data','-ascii');
save('Sample1A_proc.txt','processed','-ascii');
