%% Power Saturation Analysis Package
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.7
% Last Modified 11/11/2019

clear all;  % Clears all variables and structures

%% Read in the Bruker EPR File

filenum = input('Number of files: ','s'); % Prompts the user for the number of files
filenum = str2num(filenum); % converts the input string into a floating point number
fname = input('File Name: ','s'); % Prompts the user for the filename

for i = 1:filenum % loops over the number of files and reads in the data
    [x,y,params,fnamedata] = eprload();
    y = y/filenum;
    arrayy(:,:,i) = y;
end

[p,q] = size(y); % gets the dimensions of the power sat data
ysum = zeros(p,q); % creates a array for summing the data

for i = 1:filenum
    ysum = arrayy(:,:,i) + ysum; % sums the data
end

y = ysum; % writes the averaged data to the y array for processing

% fname = fname(1:end-4); % Truncates the file extension
[m n] = size(y);
mmin = round(m*0.4);
mmax = round(m*0.6);
ysmooth = datasmooth(y,10); % Smoothes the data for analysis

for i=1:params.SSY; % for structure to loop over all scans
    power(i) = 200/(10^((params.MPD+((i-1)*params.MPS))/10)); % calculates the power based on the attenuation
    sqrtpower(i) = power(i).^0.5; % takes the square root of the power
    height(i) = max(ysmooth(mmin:mmax,i))-min(ysmooth(mmin:mmax,i)); % detemines the peak to trough height for each scan
end

if sqrtpower(params.SSY) < sqrtpower(1)
    power = fliplr(power)
    sqrtpower = fliplr(sqrtpower)
    height = fliplr(height)
end

%% Fit the Data to obtain P1/2
options = optimset('MaxFunEvals',10000); % Sets the Maximum number of function evaluations
[x,fvalA] = fmincon(@powerfit,[max(height)/10,30,0.5],[],[],[],[],[0,0,0.5],[1e6,1e6,1.5],[],options,power,sqrtpower,height); % Uses the Matlab builtin fminsearch to minimize the t1fit function
%[x,fvalA] = fmincon(@powerfit,[max(height),5,0.5],[],[],[],[],[],[],[],options,power,sqrtpower,height);
%[x,fvalA] = fminsearch(@powerfit,[max(height),5,0.5],options,power,sqrtpower,height); % Uses the Matlab builtin fminsearch to minimize the t1fit function

[y,fvalB] = fminsearch(@linearfit,[(height(4)-height(1))/(sqrtpower(4)-sqrtpower(1)),height(1)],options,sqrtpower,height);

xdata = linspace(min(sqrtpower),max(sqrtpower),1001); % Creates an x data array for the plotted fit
[m n] = size(xdata); % gets the size of the x data array
for i = 1:n % for loop to create an array for the fit
% ydata(i) = x(1)*xdata(i)/((1+(xdata(i)^2/x(2)))^(0.5*x(3))); % the fit function
ydata(i) = x(1)*xdata(i)*(1 + ((2^(1/x(3))-1) * xdata(i)^2)/x(2))^-x(3);
end; % ends the for loop
plot(xdata,xdata*y(1)+y(2),xdata,xdata*y(1)/2+y(2),sqrtpower,height,'o',xdata,ydata,'b','MarkerEdgeColor','k','MarkerFaceColor','b'); % plots the data and the fit
axis tight; % makes the axis fit the data tightly
xlabel('square root power (mw)','FontSize',14), ylabel('peak to trough height','FontSize',14); % Constructs the axes labels
title('Power Saturation Curve','FontSize',16); % Constructs the title
powerlabel = [sprintf('P_1_/_2 = %.1f ',x(2)) 'mw']; % Assigns a string to the obtained value for T1 on the plot
text(x(2)^0.5+max(xdata)/30,max(ydata)/2,powerlabel,'FontSize',18,'color','k'); % Prints the T1 string on the plot
line([x(2)^0.5;x(2)^0.5],[0;max(height)],'Color','r','LineWidth',2); % Draws a line at the T1 value
data = [sqrtpower' height'] % creates a data array
fnamepower = sprintf('%s.txt',fname) % Constructs a filename string from the inputed filename
save(fnamepower,'data','-ascii') % Saves and ASCII file of the data array
