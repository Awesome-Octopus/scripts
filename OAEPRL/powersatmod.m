%% Power Saturation Analysis Package
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.6
% Last Modified 8/18/2020

clear all;  % Clears all variables and structures

%% Read in the Bruker EPR File

[x y params,fname] = eprload(); % Uses EasySpin eprload function to read in data
fname = fname(1:end-4); % Truncates the file extension
x = x{1};
[m n] = size(y);
% mmin = round(m*0.4);
% mmax = round(m*0.6);
ysmooth = datasmooth(y,10); % Smoothes the data for analysis

for i=1:params.YPTS; % for structure to loop over all scans
    height(i) = max(ysmooth(:,i))-min(ysmooth(:,i)); % detemines the peak to trough height for each scan
end

powerminatt = 10*log10(200/(params.YMIN));
powermaxatt = 10*log10(200/(params.YWID-params.YMIN));
poweratt = linspace(powerminatt,powermaxatt,n);
power = 200./(10.^(poweratt/10));
sqrtpower = power.^0.5;

if height(1) > height(n)
    height = flip(height);
end


% 
% a = max(height)/2
% b = 10
% c = 1.2
% xdata = linspace(min(sqrtpower),max(sqrtpower),m); % Creates an x data array for the plotted fit
% for i = 1:m % for loop to create an array for the fit
% ydata(i) = a*xdata(i)/((1+(xdata(i)^2/b))^(0.5*c)); % the fit function
% % ydata(i) = x(1)*xdata(i)*(1 + ((2^(1/x(3))-1) * xdata(i)^2)/x(2))^-x(3);
% end; % ends the for loop
% 
% plot(sqrtpower,height,xdata,ydata)

%% Fit the Data to obtain P1/2
options = optimset('MaxFunEvals',10000'); % Sets the Maximum number of function evaluations

[x,fvalA] = fmincon(@powerfit,[max(height)/2,10,1.2],[],[],[],[],[0,0,0.5],[1e6,1e6,1.5],[],options,power,sqrtpower,height); % Uses the Matlab builtin fminsearch to minimize the t1fit function
%[x,fvalA] = fmincon(@powerfit,[max(height),5,0.5],[],[],[],[],[],[],[],options,power,sqrtpower,height);
% [x,fvalA] = fminsearch(@powerfit,[max(height),5,0.5],options,power,sqrtpower,height); % Uses the Matlab builtin fminsearch to minimize the t1fit function

[y,fvalB] = fminsearch(@linearfit,[(height(4)-height(1))/(sqrtpower(4)-sqrtpower(1)),height(1)],options,sqrtpower,height);

xdata = linspace(min(sqrtpower),max(sqrtpower),1001); % Creates an x data array for the plotted fit
[m n] = size(xdata); % gets the size of the x data array
for i = 1:n % for loop to create an array for the fit
ydata(i) = x(1)*xdata(i)*(1 + ((2^(1/x(3))-1) * xdata(i)^2)/x(2))^-x(3);
end; % ends the for loop
plot(xdata,xdata*y(1)+y(2),xdata,xdata*y(1)/2+y(2),sqrtpower,height,'o',xdata,ydata,'b','MarkerEdgeColor','k','MarkerFaceColor','b'); % plots the data and the fit
axis tight; % makes the axis fit the data tightly
xlabel('square root power (mw)','FontSize',14), ylabel('peak to trough height','FontSize',14); % Constructs the axes labels
title('Power Saturation Curve','FontSize',16); % Constructs the title
powerlabel = [sprintf('P_1_/_2 = %.1f ',x(2)) 'mw']; % Assigns a string to the obtained value for T1 on the plot
text(x(2)^0.5+max(xdata)/30,max(ydata)/2,powerlabel,'FontSize',18,'color','k'); % Prints the T1 string on the plot
line([x(2)^0.5;x(2)^0.5],[0;max(height)],'Color','r','LineWidth',2); % Draws a line at the T1 value
data = [sqrtpower' height']; % creates a data array
fname = input('Choose a File Name: ','s'); % Asks for a filename input
fnamepower = sprintf('%s.txt',fname); % Constructs a filename string from the inputed filename
save(fnamepower,'data','-ascii'); % Saves and ASCII file of the data array
%data array
data = [xdata' ydata'];
fnamepowerf=sprintf('%s_fit.txt',fname);
save(fnamepowerf,'data','-ascii');
%second curve
[x1 y1 params1,fname] = eprload(); % Uses EasySpin eprload function to read in data
fname = fname(1:end-4); % Truncates the file extension
x1 = x1{1};
[m n] = size(y1);
% mmin = round(m*0.4);
% mmax = round(m*0.6);
ysmooth1 = datasmooth(y1,10); % Smoothes the data for analysis

for i=1:params1.YPTS; % for structure to loop over all scans
    height1(i) = max(ysmooth1(:,i))-min(ysmooth1(:,i)); % detemines the peak to trough height for each scan
end

powerminatt1 = 10*log10(200/(params1.YMIN));
powermaxatt1 = 10*log10(200/(params1.YWID-params1.YMIN));
poweratt1 = linspace(powerminatt1,powermaxatt1,n);
power1 = 200./(10.^(poweratt1/10));
sqrtpower1 = power1.^0.5;

if height1(1) > height1(n)
    height1 = flip(height1);
end
% 
% a = max(height)/2
% b = 10
% c = 1.2
% xdata = linspace(min(sqrtpower),max(sqrtpower),m); % Creates an x data array for the plotted fit
% for i = 1:m % for loop to create an array for the fit
% ydata(i) = a*xdata(i)/((1+(xdata(i)^2/b))^(0.5*c)); % the fit function
% % ydata(i) = x(1)*xdata(i)*(1 + ((2^(1/x(3))-1) * xdata(i)^2)/x(2))^-x(3);
% end; % ends the for loop
% 
% plot(sqrtpower,height,xdata,ydata)

%% Fit the Data to obtain P1/2
options = optimset('MaxFunEvals',10000'); % Sets the Maximum number of function evaluations

[x1,fvalA] = fmincon(@powerfit,[max(height1)/2,10,1.2],[],[],[],[],[0,0,0.5],[1e6,1e6,1.5],[],options,power1,sqrtpower1,height1); % Uses the Matlab builtin fminsearch to minimize the t1fit function
%[x,fvalA] = fmincon(@powerfit,[max(height),5,0.5],[],[],[],[],[],[],[],options,power,sqrtpower,height);
% [x,fvalA] = fminsearch(@powerfit,[max(height),5,0.5],options,power,sqrtpower,height); % Uses the Matlab builtin fminsearch to minimize the t1fit function

[y1,fvalB] = fminsearch(@linearfit,[(height1(4)-height1(1))/(sqrtpower1(4)-sqrtpower1(1)),height1(1)],options,sqrtpower1,height1);

xdata1 = linspace(min(sqrtpower1),max(sqrtpower1),1001); % Creates an x data array for the plotted fit
[m n] = size(xdata); % gets the size of the x data array
for i = 1:n % for loop to create an array for the fit
ydata1(i) = x1(1)*xdata1(i)*(1 + ((2^(1/x1(3))-1) * xdata1(i)^2)/x1(2))^-x1(3);
end; % ends the for loop
figure;
plot(xdata1,xdata1*y1(1)+y1(2),xdata1,xdata1*y1(1)/2+y1(2),sqrtpower1,height1,'o',xdata1,ydata1,'b','MarkerEdgeColor','k','MarkerFaceColor','b'); % plots the data and the fit
axis tight; % makes the axis fit the data tightly
xlabel('square root power (mw)','FontSize',14), ylabel('peak to trough height','FontSize',14); % Constructs the axes labels
title('Power Saturation Curve','FontSize',16); % Constructs the title
powerlabel1 = [sprintf('P_1_/_2 = %.1f ',x1(2)) 'mw']; % Assigns a string to the obtained value for T1 on the plot
text(x1(2)^0.5+max(xdata1)/30,max(ydata1)/2,powerlabel1,'FontSize',18,'color','k'); % Prints the T1 string on the plot
line([x1(2)^0.5;x1(2)^0.5],[0;max(height1)],'Color','r','LineWidth',2); % Draws a line at the T1 value
data = [sqrtpower1' height1']; % creates a data array
fname1 = input('Choose a File Name: ','s'); % Asks for a filename input
fnamepower1 = sprintf('%s.txt',fname1); % Constructs a filename string from the inputed filename
save(fnamepower1,'data','-ascii'); % Saves and ASCII file of the data array
% data array
data1 = [xdata1' ydata1'];
fnamepowerf1=sprintf('%s_fit.txt',fname1);
save(fnamepowerf1,'data1','-ascii');

% Third Curve
[x2 y2 params2,fname] = eprload(); % Uses EasySpin eprload function to read in data
fname = fname(1:end-4); % Truncates the file extension
x2 = x2{1};
[m n] = size(y2);
% mmin = round(m*0.4);
% mmax = round(m*0.6);
ysmooth2 = datasmooth(y2,10); % Smoothes the data for analysis

for i=1:params2.YPTS; % for structure to loop over all scans
    height2(i) = max(ysmooth2(:,i))-min(ysmooth2(:,i)); % detemines the peak to trough height for each scan
end

powerminatt2 = 10*log10(200/(params2.YMIN));
powermaxatt2 = 10*log10(200/(params2.YWID-params2.YMIN));
poweratt2 = linspace(powerminatt2,powermaxatt2,n);
power2 = 200./(10.^(poweratt2/10));
sqrtpower2 = power2.^0.5;

if height2(1) > height2(n)
    height2 = flip(height2);
end

%%% Scaling the value of NiEDDA Intensity to match same starting point of
%%% Nitrogen data assuming that the data collected from the lowest dB
%%% (highest power) to highest dB (lowest power)
  
height2 = height2* height2(19)/height1 (19);
% 
% a = max(height)/2
% b = 10
% c = 1.2
% xdata = linspace(min(sqrtpower),max(sqrtpower),m); % Creates an x data array for the plotted fit
% for i = 1:m % for loop to create an array for the fit
% ydata(i) = a*xdata(i)/((1+(xdata(i)^2/b))^(0.5*c)); % the fit function
% % ydata(i) = x(1)*xdata(i)*(1 + ((2^(1/x(3))-1) * xdata(i)^2)/x(2))^-x(3);
% end; % ends the for loop
% 
% plot(sqrtpower,height,xdata,ydata)

%% Fit the Data to obtain P1/2
options = optimset('MaxFunEvals',10000'); % Sets the Maximum number of function evaluations

[x2,fvalA] = fmincon(@powerfit,[max(height2)/2,10,1.2],[],[],[],[],[0,0,0.5],[1e6,1e6,1.5],[],options,power2,sqrtpower2,height2); % Uses the Matlab builtin fminsearch to minimize the t1fit function
%[x,fvalA] = fmincon(@powerfit,[max(height),5,0.5],[],[],[],[],[],[],[],options,power,sqrtpower,height);
% [x,fvalA] = fminsearch(@powerfit,[max(height),5,0.5],options,power,sqrtpower,height); % Uses the Matlab builtin fminsearch to minimize the t1fit function

[y2,fvalB] = fminsearch(@linearfit,[(height2(4)-height2(1))/(sqrtpower2(4)-sqrtpower2(1)),height2(1)],options,sqrtpower2,height2);

xdata2 = linspace(min(sqrtpower2),max(sqrtpower2),1001); % Creates an x data array for the plotted fit
[m n] = size(xdata2); % gets the size of the x data array
for i = 1:n % for loop to create an array for the fit
ydata2(i) = x2(1)*xdata2(i)*(1 + ((2^(1/x2(3))-1) * xdata2(i)^2)/x2(2))^-x2(3);
end; % ends the for loop
figure;
plot(xdata2,xdata2*y2(1)+y2(2),xdata2,xdata2*y2(1)/2+y2(2),sqrtpower2,height2,'o',xdata2,ydata2,'b','MarkerEdgeColor','k','MarkerFaceColor','b'); % plots the data and the fit
axis tight; % makes the axis fit the data tightly
xlabel('square root power (mw)','FontSize',14), ylabel('peak to trough height','FontSize',14); % Constructs the axes labels
title('Power Saturation Curve','FontSize',16); % Constructs the title
powerlabel2 = [sprintf('P_1_/_2 = %.1f ',x2(2)) 'mw']; % Assigns a string to the obtained value for T1 on the plot
text(x2(2)^0.5+max(xdata2)/30,max(ydata2)/2,powerlabel2,'FontSize',18,'color','k'); % Prints the T1 string on the plot
line([x2(2)^0.5;x2(2)^0.5],[0;max(height2)],'Color','r','LineWidth',2); % Draws a line at the T1 value
data = [sqrtpower2' height2']; % creates a data array
fname2 = input('Choose a File Name: ','s'); % Asks for a filename input
fnamepower2 = sprintf('%s.txt',fname2); % Constructs a filename string from the inputed filename
save(fnamepower2,'data','-ascii'); % Saves and ASCII file of the data array
% data array
data12 = [xdata2' ydata2'];
fnamepowerf2=sprintf('%s_fit.txt',fname2);
save(fnamepowerf2,'data12','-ascii');

%%%%%%%%%
figure;
plot(sqrtpower,height,'or',xdata,ydata,'r',sqrtpower1,height1,'ob',xdata1,ydata1,'b',sqrtpower2,height2,'ok',xdata2,ydata2,'k');%,'MarkerEdgeColor','k','MarkerFaceColor','b'); % plots the data and the fit
xlabel('P^1^/^2,mW')
ylabel('Central peak height')