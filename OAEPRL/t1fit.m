%% T1 Fitting Package
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.13
% Last Modified 02/21/2011

clear all;  % Clears all variables and structures

%% Read in the Bruker EPR File

[b spec params,fname] = eprload(); % Uses EasySpin eprload function to read in data
fname = fname(1:end-4); % Truncates the file extension
[shift ns] = strread(params.FTEzDelay1,'%f %s'); % Reads in the value of tau (d1)
[shift2 ns2] = strread(params.FTEzDelay2,'%f %s'); % Reads in the value of T (d2)
b = (b'+shift*2+shift2)/1000; % Shifts the x axis values by tau (d1) + T (d2)
spec = real(spec'/max(spec)); % Removes the imaginary component of the y data

%% Minimize the Exponential Fit

options = optimset('MaxFunEvals',10000); % Sets the Maximum number of function evaluations
[x,fval] = fminsearch(@t1decay,[min(spec),max(b)/4,0],options,b,spec,params); % Uses the Matlab builtin fminsearch to minimize the t1fit function
opt = t1opt(x,b,spec,params); % Uses the t1opt function to find the optimum T1

%% Plot the Data

plot(b,spec,b,x(1)*(1-2*exp(-b/x(2)))+x(3)); % Plots the experimental data and the optained fit
axis tight; % Shrinks the axis to fit the data tightly
xlabel('time(\mus)','FontSize',18), ylabel('intensity','FontSize',18); % Constructs the axes labels
title('T_1 Determination','FontSize',20); % Constructs the title
T1label = [sprintf('T_1 = %.0f ',x(2)) '\mus']; % Assigns a string to the obtained value for T1 on the plot
text(x(2)+max(b)/30,0,T1label,'FontSize',18,'color','r'); % Prints the T1 string on the plot
optlabel = [sprintf('Optimum SRT = %.0f ',opt) '\mus']; % Assigns a string to the obtained optimum SRT
text(opt+max(b)/30,max(spec)/2,optlabel,'FontSize',18,'color','b'); % Prints the Optimum SRT string on the plot
line([x(2);x(2)],[min(spec);max(spec)],'Color','r','LineWidth',2); % Draws a line at the T1 value
line([opt;opt],[min(spec);max(spec)],'Color','b','LineWidth',2); % Draws a line 

%% Save and ASCII File of the Results

fit = x(1)*(1-2*exp(-b/x(2)))+x(3); % Consrtucts an array from the best fit to the data
data = [b' spec' fit']; % Cosntructs an array with the x, y and fit data
fnamet1 = sprintf('%s.txt',fname); % Constructs a filename string from the params structure
save(fnamet1,'data','-ascii'); % Saves the data array to an ASCII file
