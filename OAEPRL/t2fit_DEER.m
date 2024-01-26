%% T2-DEER Fitting Package
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.13
% Last Modified 02/21/2011

clear all;  % Clears all variables and structures

%% Read in the Bruker EPR File

[b spec params fname] = eprload(); % Uses EasySpin eprload function to read in data
fname = fname(1:end-4); % Truncates the file extension
[shift ns] = strread(params.FTEzDelay1,'%f %s'); % Reads in the value of tau (d1)
b = (b' + shift)/1000; % Shifts the x axis values by tau (d1)
spec = real(spec'/max(spec)); % Removes the imaginary component of the y data

%% Minimize the Exponential Fit

options = optimset('MaxFunEvals',10000); % Sets the Maximum number of function evaluations
[x,fval] = fminsearch(@t2decay,[max(spec)/2,max(b)/4],options,b,spec,params); % Uses the Matlab builtin fminsearch to minimize the t2fit function

%% Plot the Data

plot(b,spec,b,x(1)*exp(-b/x(2))); % Plots the experimental data and the optained fit
axis tight; % Shrinks the axis to fit the data tightly
title('T_2 Determination','FontSize',20); % Constructs the title
xlabel('time(\mus)','FontSize',18), ylabel('intensity','FontSize',18); % Constructs the axes labels
label = [sprintf('T_2 = %0.2f ',x(2)) '\mus']; % Assigns a string to the obtained value for T2 on the plot
text(max(b)*0.6,max(spec)*0.8,label,'FontSize',18,'Color','r'); % Prints the T2 string on the plot
line([x(2);x(2)],[0;max(spec)],'Color','r','LineWidth',2); % Draws a line at the T2 value
line([2.548;2.548],[0;max(spec)],'Color','b','LineWidth',2); % Draws a line at the 1 us value
line([4.548;4.548],[0;max(spec)],'Color','b','LineWidth',2); % Draws a line at the 2 us value
line([6.548;6.548],[0;max(spec)],'Color','b','LineWidth',2); % Draws a line at the 3 us value
line([8.548;8.548],[0;max(spec)],'Color','b','LineWidth',2); % Draws a line at the 4 us value
text(2.65,max(spec)*0.5,'1 \mus','FontSize',18,'Color','b'); % Prints the 1 us string on the plot
text(4.65,max(spec)*0.4,'2 \mus','FontSize',18,'Color','b'); % Prints the 2 us string on the plot
text(6.65,max(spec)*0.3,'3 \mus','FontSize',18,'Color','b'); % Prints the 3 us string on the plot
text(8.65,max(spec)*0.2,'4 \mus','FontSize',18,'Color','b'); % Prints the 4 us string on the plot

oneus = x(1)*exp(-2.548/x(2))*100; % % max echo intensity at 1 us
twous = x(1)*exp(-4.548/x(2))*100; % % max echo intensity at 2 us
threeus = x(1)*exp(-6.548/x(2))*100; % % max echo intensity at 3 us
fourus = x(1)*exp(-8.548/x(2))*100; % % max echo intensity at 4 us

oneuslabel = ['1 \mus' sprintf(' = %0.1f ',oneus) '%']; % 1 us % label
text(max(b)*0.6,max(spec)*0.7,oneuslabel,'FontSize',18,'Color','b'); % Prints the 1 us string on the plot
twouslabel = ['2 \mus' sprintf(' = %0.1f ',twous) '%'];% 2 us % label
text(max(b)*0.6,max(spec)*0.6,twouslabel,'FontSize',18,'Color','b'); % Prints the 2 us string on the plot
threeuslabel = ['3 \mus' sprintf(' = %0.1f ',threeus) '%'];% 3 us % label
text(max(b)*0.6,max(spec)*0.5,threeuslabel,'FontSize',18,'Color','b'); % Prints the 3 us string on the plot
fouruslabel = ['4 \mus' sprintf(' = %0.1f ',fourus) '%'];% 4 us % label
text(max(b)*0.6,max(spec)*0.4,fouruslabel,'FontSize',18,'Color','b'); % Prints the 4 us string on the plot


%% Save and ASCII File of the Results

fit = x(1)*exp(-b/x(2)); % Consrtucts an array from the best fit to the data
data = [b' spec' fit']; % Cosntructs an array with the x, y and fit data
fnamet2 = sprintf('%s.txt',fname); % Constructs a filename string from the params structure
save(fnamet2,'data','-ascii'); % Saves the data array to an ASCII file
