
%% Polymer Spin Concentration Determination Package
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.13
% Last Modified 06/09/2017

clear all; % clears all of the variables and structures

%% Read in Bruker Binary File and Establish Constants

standardint = 62.0867 % Standard from experiment
standardconc = 150; % Standard from experiment
[x y params] = eprload('TEMPO in water 003uM gain E4.spc'); % Uses EasySpin eprload function to read in data
% height = input('Sample Height (mm): ','s'); % Prompts the user for the sample height
% height = str2num(height); % converts the input string into a floating point number
height = 40;
%% Smooth anhelpd Integrate the Data

smoothy = filter(2,1,y); % Smoothes the data for better integration
inty = cumtrapz(x,smoothy); % integral of data
dblinty = cumtrapz(x,inty); % double integral of data
maxint = (max(dblinty)/params.JSD)/params.RRG; % max value of double integral

%% Determine the Concentation

heightfactor = (1-exp(-height/7.5)) % Scaling for sample volume
conc = standardconc*maxint/standardint/heightfactor; % concentration based on standard

%% Plot the Data

subplot(2,2,1)
plot(x,y,x,smoothy); % plots the spectrum
xlabel('Field (G)','FontSize',14), ylabel('d\chi"/dB','FontSize',14); % Constructs the axes labels
title('CW EPR Spectrum','FontSize',16); % Constructs the title
conclabel = [sprintf('Concentration = %.0f ',conc) '\muM']; % Assigns a string to the obtained value for T1 on the plot
text(min(x)+(max(x)-min(x))/20,min(y)+(max(y)-min(y))/200,conclabel,'FontSize',18,'color','r'); % Prints the T1 string on the plot
subplot(2,2,2)
plot(x,dblinty)
subplot(2,2,3)
plot(x,inty)