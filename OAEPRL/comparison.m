%% TEMPO Comparison Package
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.13
% Last Modified 02/21/2011

clear all; % clears all of the variables and structures

%% Read in Bruker Binary File and Establish Constants

[x y params] = eprload(); % Uses EasySpin eprload function to read in data
y = real(y); % Eliminates the Imaginary Data
[gain db] = strread(params.VideoGain,'%f %s'); % Reads in the video Gain setting
conc = input('Sample Concentration (uM): ','s'); % Prompts the user for the sample concentration
conc = str2num(conc); % converts the input string into a floating point number

%% Smooth and Integrate the Data

smoothy = datasmooth(y,10); % smoothes the data for better integration
inty = cumtrapz(x,smoothy); % integral of data
maxint = max(inty); % max value of double integral

%%  Adjustment for Video Gain

maxint = maxint/(10^((gain-27)/20)); % adjusts the intensity for the video gain setting

%% Compare with the TEMPO Standard

compare = maxint/(conc*106237-163298)*100; % Compares the integrated data to the TEMPO curve

%%  Plot the Data

plot(x,y); % plots the data
xlabel('Field (G)','FontSize',14), ylabel('intensity','FontSize',14); % constructs the axes labels
title('ESE EPR Spectrum','FontSize',16); % constructs the title
label = [sprintf('TEMPO Comparison = %.1f ',compare) '%']; % creates the figure label
text(min(x)+(max(x)-min(x))*0.05,max(y)*0.95,label,'FontSize',16,'color','k'); % plots the figure label
