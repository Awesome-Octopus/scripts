
%% Concentration Determination Package
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 2.0
% Last Modified 11/14/2019

clear all; % clears all of the variables and structures

%% Read in Bruker Binary File and Establish Constants

slope = 6.1667 % Standard from experiment
offset = 1;
[x y params] = eprload(); % Uses EasySpin eprload function to read in data

y = y*1e4;

%% Smooth and Integrate the Data

[m,n] = size(y)
first = mean(y(1:50))
last = mean(y(m-50:m))
fit = linspace(first,last,m)
y = y-fit';

smoothy = datasmooth(y,20);
inty = cumtrapz(x,smoothy);
intfirst = mean(inty(1:50))
intlast = mean(inty(m-50:m))
intfit = linspace(intfirst,intlast,m)
inty = inty-intfit';
dblinty = cumtrapz(x,inty);
maxint = (((max(dblinty)-min(dblinty))/(10^(params.RCAG/10)))/params.AVGS)-offset

%% Determine the Concentation

conc = (maxint*slope);

%% Plot the Data

subplot(2,2,1)
plot(x,y); % plots the spectrum
xlabel('Field (G)','FontSize',14), ylabel('d\chi"/dB','FontSize',14); % Constructs the axes labels
title('CW EPR Spectrum','FontSize',16); % Construacts the title
conclabel = [sprintf('Concentration = %.0f ',conc) '\muM']; % Assigns a string to the obtained value for T1 on the plot
text(min(x)+(max(x)-min(x))/20,min(y)+(max(y)-min(y))/200,conclabel,'FontSize',18,'color','r'); % Prints the T1 string on the plot

subplot(2,2,2)
plot(x,smoothy)

subplot(2,2,3)
plot(x,inty)

subplot(2,2,4)
plot(x,dblinty)


% plot(x,y); % plots the spectrum
% xlabel('Field (G)','FontSize',14), ylabel('d\chi"/dB','FontSize',14); % Constructs the axes labels
% title('CW EPR Spectrum','FontSize',16); % Construacts the title
% conclabel = [sprintf('Concentration = %.0f ',conc) '\muM']; % Assigns a string to the obtained value for T1 on the plot
% text(min(x)+(max(x)-min(x))/20,min(y)+(max(y)-min(y))/200,conclabel,'FontSize',18,'color','r'); % Prints the T1 string on the plot
% 
% 
% % 1.6778 1.4563 1.3062 1.6753 120 128 127 59 62 54