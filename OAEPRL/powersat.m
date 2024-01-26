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

powerminatt = 10*log10(200/(params.YMIN))
powermaxatt = 10*log10(200/(params.YWID-params.YMIN))
poweratt = linspace(powerminatt,powermaxatt,n)
power = 200./(10.^(poweratt/10))
sqrtpower = power.^0.5

if height(1) > height(n)
    height = flip(height)
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

[x,fvalA] = fmincon(@powerfit,[max(height)/4,10,1.2],[],[],[],[],[0,0,0.5],[1e6,1e6,1.5],[],options,power,sqrtpower,height); % Uses the Matlab builtin fminsearch to minimize the t1fit function
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
fnamepower = sprintf('%s.txt',fname); % Constructs a filename string from the inputed filename
save(fnamepower,'data','-ascii'); % Saves and ASCII file of the data array
