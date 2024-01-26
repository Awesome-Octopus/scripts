
%% Bruker ASCII Conversion Package
% Ohio Advanced EPR Laboratory
% Rob McCarrick
% OAEPRL EPR Processings Package, Version 1.13
% Last Modified 02/21/2011

clear all;  % Clears all variables and structures

%% Read in the Bruker EPR File

[x y params fname] = eprload(); % Uses EasySpin eprload function to read in data
yreal = real(y); % Creates an array with the real y data (green channel)
yimag = imag(y); % Creates an array with the imaginary y data (yellow channel)
fname = fname(1:end-4); % Truncates the file extension

%% if/else structure that determines which type of file 

if (iscell(x) == 0); % if stucture that checks to see if the SSY field is absent (find 1D data)
if (sum(yimag) == 0); % if structure that checks for an imaginary component to the data
    if (isfield(params,'SSY') == 0); % if stucture that checks to see if the SSY field is absent (find 1D data)
        data = [x yreal]; % Cosntructs an array with the x and y data
        fnamecw = sprintf('%s.txt',fname); % Constructs a filename string from the inputed filename
        save(fnamecw,'data','-ascii'); % Saves and ASCII file of the data array
        plot(x,yreal); % Plots the data
        xlabel('Field (G)','FontSize',14), ylabel('d\chi"/dB','FontSize',14); % Constructs the axes labels
        title('CW EPR Spectrum','FontSize',16); % Constructs the title
    elseif (isfield(params,'SSY') == 1); % elseif structure that checks for a 2D data set (for power saturation)
        for i = 1:params.SSY; % for loop that writes an ASCII file for each slice of data
          datax = x; % Creates an array of the x data
          datay = y(:,i); % Creates an array for one slice of the y data
          power = linspace(200/(10^((params.MPD+((i-1)*params.MPS))/10)),200/(10^((params.MPD+((i-1)*params.MPS))/10)),params.SSX); % Creates an array with the power value for the scan
          data = [power' datax' datay]; % Creates an array with the data
          fnamecwPOWER.i = sprintf('%s_%i.txt',fname,i); % Creates a string with the filename and index
          save(fnamecwPOWER.i,'data','-ascii'); % Saves and ASCII file of the data array
        end; % ends the for loop
    end; % ends the if/elseif structure
elseif (isfield(params,'YPTS') == 1); % if stucture that checks to see if the YPTS field is present (find 2D data)
    [delay1 ns] = strread(params.FTEzDelay1,'%f %s'); % Pulls d1 from the params structure
    [delay2 ns] = strread(params.FTEzDelay2,'%f %s'); % Pulls d2 from the params structure
    [deltax ns] = strread(params.FTEzDeltaX,'%f %s'); % Pulls dx from the params structure
    [deltay ns] = strread(params.FTEzDeltaY,'%f %s'); % Pulls dy from the params structure
    if (isvector(findstr(params.TITL,'3P')) || isvector(findstr(params.TITL,'3p'))); % if structure that checks for a 3P vs. tau data set
        for i = 1:params.YPTS; % for loop that writes an ASCII file for each slice of data
        datax = linspace(delay2,delay2+((params.XPTS-1)*deltax),params.XPTS); % Creates an array of the x data from the experimental parameters
        datay = linspace(delay1+(deltay*(i-1)),delay1+(deltay*(i-1)),params.XPTS); % Creates an array for one slice of the y data
        datayreal = yreal(:,i); % Creates an array for one slice of the y real data
        datayimag = yimag(:,i); % Creates an array for one slice of the y imaginary data
        data = [datax' datay' datayreal datayimag]; % Creates an array with the data
        fname3P.i = sprintf('%s_%i.txt',fname,i); % Creates a string with the filename and index
        save(fname3P.i,'data','-ascii'); % Saves and ASCII file of the data array
        end; % ends the for loop
    elseif (isvector(findstr(params.TITL,'HYSCORE')) || isvector(findstr(params.TITL,'hyscore')) || isvector(findstr(params.TITL,'Hyscore'))); % if structure that checks for a HYSCORE data set
        datax = linspace(delay1,delay1+((params.XPTS-1)*deltax),params.XPTS); % Creates an array of the x data from the experimental parameters
        data = [datax' yreal']; % Creates an array with the data
        fname = sprintf('%s.txt',fname); % Creates a string with the filename and index
        save(fname,'data','-ascii'); % Saves and ASCII file of the data array
    end; % ends the if/elseif structure
elseif (isfield(params,'YPTS') == 0);
    data = [x yreal yimag]; % Creates an array with the data
    fnamepulse = sprintf('%s.txt',fname); % Creates a string with the filename
    save(fnamepulse,'data','-ascii'); % Saves and ASCII file of the data array
        if (isvector(findstr(params.TITL,'DEER')) || isvector(findstr(params.TITL,'deer'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('time(ns)','FontSize',14), ylabel('refocused echo intensity','FontSize',14); % Constructs the axes labels
        title('DEER Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'T2')) || isvector(findstr(params.TITL,'t2'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('time(ns)','FontSize',14), ylabel('echo intensity','FontSize',14); % Constructs the axes labels
        title('T2 Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'T1')) || isvector(findstr(params.TITL,'t1'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('time(ns)','FontSize',14), ylabel('inverted echo intensity','FontSize',14); % Constructs the axes labels
        title('T1 Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'FS')) || isvector(findstr(params.TITL,'fs'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('Field (G)','FontSize',14), ylabel('echo intensity','FontSize',14); % Constructs the axes labels
        title('ESE-EPR Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'2P')) || isvector(findstr(params.TITL,'2p'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('time (ns)','FontSize',14), ylabel('echo intensity','FontSize',14); % Constructs the axes labels
        title('2P ESEEM Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'3P')) || isvector(findstr(params.TITL,'3p'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('time (ns)','FontSize',14), ylabel('stimulated echo intensity','FontSize',14); % Constructs the axes labels
        title('3P ESEEM Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'MIMS')) || isvector(findstr(params.TITL,'mims')) || isvector(findstr(params.TITL,'Mims'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('frequency (MHz)','FontSize',14), ylabel('stimulated echo intensity','FontSize',14); % Constructs the axes labels
        title('Mims ENDOR Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'DAVIES')) || isvector(findstr(params.TITL,'davies')) || isvector(findstr(params.TITL,'Davies'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('frequency (MHz)','FontSize',14), ylabel('stimulated echo intensity','FontSize',14); % Constructs the axes labels
        title('Davies ENDOR Spectrum','FontSize',16); % Constructs the title
        else; % if structure grab everything else that doesn't fit
        plot(x,yreal,x,yimag); % Plots the data
        end; % ends the if/elseif structure
end; % ends the if/elseif structure
elseif (iscell(x) == 1); % if stucture that checks to see if the SSY field is absent (find 1D data)
x = x{1}
if (sum(yimag) == 0); % if structure that checks for an imaginary component to the data
    if (isfield(params,'SSY') == 0); % if stucture that checks to see if the SSY field is absent (find 1D data)
        if (isfield(params,'TITL') == 1);
        data = [x' yreal']; % Creates an array with the data
        fname = sprintf('%s.txt',fname); % Creates a string with the filename and index
        save(fname,'data','-ascii'); % Saves and ASCII file of the data array
        elseif (isfield(params,'TITL') == 0);
        data = [x' yreal]; % Cosntructs an array with the x and y data
        fnamecw = sprintf('%s.txt',fname); % Constructs a filename string from the inputed filename
        save(fnamecw,'data','-ascii'); % Saves and ASCII file of the data array
        plot(x,yreal); % Plots the data
        xlabel('Field (G)','FontSize',14), ylabel('d\chi"/dB','FontSize',14); % Constructs the axes labels
        title('CW EPR Spectrum','FontSize',16); % Constructs the title
        end    
    elseif (isfield(params,'SSY') == 1); % elseif structure that checks for a 2D data set (for power saturation)
        fname = input('Choose a File Name: ','s'); % Asks for a filename input
        for i = 1:params.SSY; % for loop that writes an ASCII file for each slice of data
          datax = x; % Creates an array of the x data
          datay = y(:,i); % Creates an array for one slice of the y data
          power = linspace(200/(10^((params.MPD+((i-1)*params.MPS))/10)),200/(10^((params.MPD+((i-1)*params.MPS))/10)),params.SSX); % Creates an array with the power value for the scan
          data = [power' datax' datay]; % Creates an array with the data
          fnamecwPOWER.i = sprintf('%s_%i.txt',fname,i); % Creates a string with the filename and index
          save(fnamecwPOWER.i,'data','-ascii'); % Saves and ASCII file of the data array
        end; % ends the for loop
    end; % ends the if/elseif structure
elseif (isfield(params,'YPTS') == 1); % if stucture that checks to see if the YPTS field is present (find 2D data)
    [delay1 ns] = strread(params.FTEzDelay1,'%f %s'); % Pulls d1 from the params structure
    [delay2 ns] = strread(params.FTEzDelay2,'%f %s'); % Pulls d2 from the params structure
    [deltax ns] = strread(params.FTEzDeltaX,'%f %s'); % Pulls dx from the params structure
    [deltay ns] = strread(params.FTEzDeltaY,'%f %s'); % Pulls dy from the params structure
    if (isvector(findstr(params.TITL,'3P')) || isvector(findstr(params.TITL,'3p'))); % if structure that checks for a 3P vs. tau data set
        for i = 1:params.YPTS; % for loop that writes an ASCII file for each slice of data
        datax = linspace(delay2,delay2+((params.XPTS-1)*deltax),params.XPTS); % Creates an array of the x data from the experimental parameters
        datay = linspace(delay1+(deltay*(i-1)),delay1+(deltay*(i-1)),params.XPTS); % Creates an array for one slice of the y data
        datayreal = yreal(:,i); % Creates an array for one slice of the y real data
        datayimag = yimag(:,i); % Creates an array for one slice of the y imaginary data
        data = [datax' datay' datayreal datayimag]; % Creates an array with the data
        fname3P.i = sprintf('%s_%i.txt',fname,i); % Creates a string with the filename and index
        save(fname3P.i,'data','-ascii'); % Saves and ASCII file of the data array
        end; % ends the for loop
    elseif (isvector(findstr(params.TITL,'HYSCORE')) || isvector(findstr(params.TITL,'hyscore')) || isvector(findstr(params.TITL,'Hyscore'))); % if structure that checks for a HYSCORE data set
        datax = linspace(delay1,delay1+((params.XPTS-1)*deltax),params.XPTS); % Creates an array of the x data from the experimental parameters
        data = [datax' yreal']; % Creates an array with the data
        fname = sprintf('%s.txt',fname); % Creates a string with the filename and index
        save(fname,'data','-ascii'); % Saves and ASCII file of the data array
    end; % ends the if/elseif structure
elseif (isfield(params,'YPTS') == 0);
    data = [x yreal yimag]; % Creates an array with the data
    fnamepulse = sprintf('%s.txt',fname); % Creates a string with the filename
    save(fnamepulse,'data','-ascii'); % Saves and ASCII file of the data array
        if (isvector(findstr(params.TITL,'DEER')) || isvector(findstr(params.TITL,'deer'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('time(ns)','FontSize',14), ylabel('refocused echo intensity','FontSize',14); % Constructs the axes labels
        title('DEER Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'T2')) || isvector(findstr(params.TITL,'t2'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('time(ns)','FontSize',14), ylabel('echo intensity','FontSize',14); % Constructs the axes labels
        title('T2 Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'T1')) || isvector(findstr(params.TITL,'t1'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('time(ns)','FontSize',14), ylabel('inverted echo intensity','FontSize',14); % Constructs the axes labels
        title('T1 Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'FS')) || isvector(findstr(params.TITL,'fs'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('Field (G)','FontSize',14), ylabel('echo intensity','FontSize',14); % Constructs the axes labels
        title('ESE-EPR Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'2P')) || isvector(findstr(params.TITL,'2p'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('time (ns)','FontSize',14), ylabel('echo intensity','FontSize',14); % Constructs the axes labels
        title('2P ESEEM Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'3P')) || isvector(findstr(params.TITL,'3p'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('time (ns)','FontSize',14), ylabel('stimulated echo intensity','FontSize',14); % Constructs the axes labels
        title('3P ESEEM Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'MIMS')) || isvector(findstr(params.TITL,'mims')) || isvector(findstr(params.TITL,'Mims'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('frequency (MHz)','FontSize',14), ylabel('stimulated echo intensity','FontSize',14); % Constructs the axes labels
        title('Mims ENDOR Spectrum','FontSize',16); % Constructs the title
        elseif (isvector(findstr(params.TITL,'DAVIES')) || isvector(findstr(params.TITL,'davies')) || isvector(findstr(params.TITL,'Davies'))); % if structure to test for experiment type
        plot(x,yreal,x,yimag); % Plots the data
        xlabel('frequency (MHz)','FontSize',14), ylabel('stimulated echo intensity','FontSize',14); % Constructs the axes labels
        title('Davies ENDOR Spectrum','FontSize',16); % Constructs the title
        else; % if structure grab everything else that doesn't fit
        plot(x,yreal,x,yimag); % Plots the data
        end; % ends the if/elseif structure
end
end