%% Concentration Determination Package - ELEXSYS Revision
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.13
% Last Modified 02/21/2011

clear all; % clears all of the variables and structures

%% Read in Bruker Binary File and Establish Constants

standardint = 147.9562; % Standard from experiment
standardconc = 100; % Standard from experiment
[x y params] = eprload(); % Uses EasySpin eprload function to read in data

%% Smooth and Integrate the Data

smoothy = smooth(y,10); % Smoothes the data for better integration
inty = cumtrapz(x,smoothy); % integral of data
dblinty = cumtrapz(x,inty); % double integral of data
maxint = max(dblinty); % max value of double integral

%% Determine the Concentation

conc = standardconc*maxint/standardint; % concentration based on standard

%% Plot the Data

plot(x,y); % plots the spectrum
xlabel('Field (G)','FontSize',14), ylabel('d\chi"/dB','FontSize',14); % Constructs the axes labels
title('CW EPR Spectrum','FontSize',16); % Constructs the title
conclabel = [sprintf('Concentration = %.0f ',conc) '\muM']; % Assigns a string to the obtained value for T1 on the plot
text(min(x)+(max(x)-min(x))/20,min(y)+(max(y)-min(y))/200,conclabel,'FontSize',18,'color','r'); % Prints the T1 string on the plot
