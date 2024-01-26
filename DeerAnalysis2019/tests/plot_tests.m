data=load('test_bckg.dat');
[m,n]=size(data);

figure(1); clf;
plot(data(:,1),data(:,2),'k');
hold on;
plot(data(:,1),data(:,3),'r');
if n>3,
    plot(data(:,1),data(:,4),'m');
end;

data=load('test_fit.dat');
[m,n]=size(data);

figure(2); clf;
plot(data(:,1),data(:,2),'k');
hold on;
plot(data(:,1),data(:,3),'r');

data=load('test_spc.dat');
[m,n]=size(data);

figure(3); clf;
plot(data(:,1),data(:,2),'k');
hold on;
plot(data(:,1),data(:,3),'r');

data=load('test_distr.dat');
[m,n]=size(data);

figure(4); clf;
plot(data(:,1),data(:,2),'k');

data=load('test_Lcurve.dat');
[m,n]=size(data);

figure(5); clf;
plot(data(:,1),data(:,2),'ko');
