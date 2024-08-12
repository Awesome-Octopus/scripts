
%prompt for the names of the condition groups
group_titles = strsplit(input(...
    'Enter the names of the groups of samples to be plotted, separated by |: ',...
    's'), '|');
n_groups = numel(group_titles);

% here we will make a struct to hold all our data for a set of grouped
% experimental conditions, containing ensemble data, as well as structs
% containing data for the individual replicates
group = struct('title', [], 'replicates', ...
    struct('filename', [], 'params', [], 'ft', [], ...
        'freq', [], 'subtracted_echo', [], 'scaled_echo', [], ...
        'raw_echo', [], 'time', [], 'bkgd_fit', []), ...
        'average', [], 'nreps', [], 'color', []);

pallet = {'red', 'blue', 'black', 'green', 'magenta', 'cyan', 'yellow'};

% this the agglomeration of all the experimental conditions
data = repmat(group, 1, n_groups);

for i=1:n_groups
    filenames = {};
    data(i).nreps = 0;
    data(i).color = pallet{i};
    data(i).title = group_titles{i};
    % while true
        txt = sprintf("Select one or more spectra for the group %s", group_titles{i});
        [new_filenames,new_path] = uigetfile("*.DTA", txt, "MultiSelect", "on");
    
        % provide instructions the first time
    
        % ITERATIVE SELECTION TEMPORARILY DISABLED
        % if isnumeric(new_filenames) == 1
        %     % user cancelled -> done selecting files
        %     break
        % end
    
        % this appears to be necessary if you only have one entry in your list of 
        % files or multiple files in the same path
        % because in that case the datatype for the filename is just a string,
        % whereas multiples are a cell array containing strings
        if iscell(new_filenames) == 0
            new_filenames = {new_filenames};
        end
        filenames = [filenames, fullfile(new_path,new_filenames)]
        

    % end
    data(i).nreps = numel(filenames)
    for j=1:numel(filenames)
        data(i).replicates(j).filename = filenames{j}
    end

end

for i=1:n_groups
    for j=1:data(i).nreps
        tmpname = data(i).replicates(j).filename;
        data(i).replicates(j) = eseem_analysis(tmpname);
    end
end    

function dataObj = eseem_analysis(file)
    
    [b,spec,params] = eprload(file); % Uses EasySpin eprload function to read in data
    
    if (isfield(params,'PlsSPELGlbTxt') == 0) % Check for absense of PulseSPEL information
        [shift] = 400;
        [shift2] = 100;
        [pw] = 12;
        [inc] = 8;
    
        %[shift ns] = strread(params.FTEzDelay1,'%f %s'); % Reads in the value of tau (d1)
        %[shift2 ns2] = strread(params.FTEzDelay2,'%f %s'); % Reads in the value of T (d2)
        %[pw ns] = strread(params.FTEzMWPiHalf,'%f %s'); % Reads in the value of the pulse width (p0)
        b = (b' + shift + shift2 + pw*2); % Shifts the x axis values by tau (d1) + T (d2)
        spec = real(spec'); % Removes the imaginary component of the y data
    
    elseif (isfield(params,'PlsSPELGlbTxt') == 1); % Checks for presence of PulseSPEL Data
        pulseparams = textscan(params.PlsSPELGlbTxt,'%s %s %d %*[^\n]',50,'HeaderLines',7); % Reads the parameters text
        [d1arrayx d1arrayy] = find(strcmp(pulseparams{1},'d1')); % Finds the index of the d1 parameter
        d1array = pulseparams{3}; % Writes the d1 parameter from the cell to a variable
        shift = double(d1array(d1arrayx)); % Converts the d1 value to a floating 
        [d2arrayx d2arrayy] = find(strcmp(pulseparams{1},'d2')); % Finds the index of the d3 parameter
        d2array = pulseparams{3}; % Writes the d3 parameter from the cell to a variable
        shift2 = double(d2array(d2arrayx));
        [d30arrayx d30arrayy] = find(strcmp(pulseparams{1},'d30')); % Finds the index of the d30 parameter
        d30array = pulseparams{3}; % Writes the d30 parameter from the cell to a variable
        inc = double(d30array(d30arrayx));
        
        [p0arrayx p0arrayy] = find(strcmp(pulseparams{1},'p0')); % Finds the index of the d1 parameter
        p0array = pulseparams{3}; % Writes the d1 parameter from the cell to a variable
        pw = double(p0array(p0arrayx)); % Converts the d1 value to a floating 
        
        b = (b' + shift + shift2 + pw*2); % Shifts the x axis values by tau (d1) + T (d2)
        spec = real(spec'); % Removes the imaginary component of the y data
    
        % Background Subtraction and Scaling to 1
        
        [k, c, yfit] = exponfit(b,spec,2); % Uses the EasySpin function expnfit to fit the background
        scale = max(yfit); % Establishes the maximum of the exponential fit as a scaling factor
        spec_scaled = spec/max(scale); % Scales the experimental data
        yfit_scaled = yfit/max(scale); % Scales the exponential background
        ysub = spec_scaled-yfit_scaled; % Subtracts the background function
        
        % Calculate Fouier Transform
    
        %[inc ns] = strread(params.FTEzDeltaX,'%f %s'); % Reads the time increment from the paramaters structure
        xf = fdaxis(8/1000,512); % Uses the EasySpin fdaxis function to generate an x axis for the FT
        ft = real(ctafft(ysub,params.XPTS-1,params.XPTS)); % Constructs the Cross Term Averaged FFT of the experimental data
        ft = fftshift(ft); % Shifts the FFT using the Matlab function fftshift to obtain the correct spectrum

        
    elseif (isfield(params,'YPTS') == 1); % if stucture that checks to see if the YPTS field is present (find 2D data)

    %% Read in Bruker Binary Data
      
    % [inc ns] = strread(params.FTEzDeltaX,'%f %s');
    [incy ns] = strread(params.FTEzDeltaY,'%f %s');
    [shift ns] = strread(params.FTEzDelay1,'%f %s'); % Reads in the value of tau (d1)
    [shift2 ns2] = strread(params.FTEzDelay2,'%f %s'); % Reads in the value of T (d2)
    [pw ns] = strread(params.FTEzMWPiHalf,'%f %s'); % Reads in the value of the pulse width (p0)
    x = b{1}; % Gets the first cell from the x axis
    x = (x' + shift + shift2 + pw*2); % Shifts the x axis values by tau (d1) + T (d2)
    spec = real(spec); % Removes the imaginary component of the y data

    %% Exponential Fit, Subtraction and Crossterm Averaged FFT

    for i = 1:params.YPTS; % for loop that writes an ASCII file for each slice of data
            time(:,i) = x; % Creates a mutlidimensional time array for plotting
            [k,c,yfit(:,i)] = exponfit(x,spec(:,i),2); % Uses the EasySpin function expnfit to fit the background
            scale = max(yfit(:,i)); % Establishes the maximum of the exponential fit as a scaling factor
            datay_scaled(:,i) = spec(:,i)/max(scale); % Scales the experimental data
            yfit_scaled(:,i) = yfit(:,i)/max(scale); % Scales the exponential background
            ysub(:,i) = datay_scaled(:,i)-yfit_scaled(:,i); % Subtracts the background function
            xf(:,i) = fdaxis(inc/1000,params.XPTS*2); % Uses the EasySpin fdaxis function to generate an x axis for the FT
            ft(:,i) = real(ctafft(ysub(:,i),params.XPTS-1,params.XPTS*2)); % Constructs the Cross Term Averaged FFT of the experimental data
            ft(:,i) = fftshift(ft(:,i)); % Shifts the FFT using the Matlab function fftshift to obtain the correct spectrum
            for j = 1:params.XPTS; % For loop to create a y offset for the data
            stack(j,i) = shift+ (i-1)*incy; % Creates a y dimension of tau
            end; % Cnds the for loop
    end; % Ends the for loop
    xf(1:params.XPTS,:) = []; % gets rid of the negative FFT values
    ft(1:params.XPTS,:) = []; % gets rid of the negative FFT Values

    %% return the data structure
    dataObj = struct('filename', file, 'params', params, 'ft', ft, ...
        'freq', xf, 'subtracted_echo', ysub, 'scaled_echo', spec_scaled, ...
        'raw_echo', spec, 'time', b, 'bkgd_fit', yfit_scaled);
end

fig = figure;
t = tiledlayout(fig, 3, 1);
scaled_plot_ax = nexttile(t, 1);
title(scaled_plot_ax, 'Time Domain');
scaled_plot_ax.XLabel.String = 'Time (ns)'
scaled_plot_ax.YLabel.String = 'Intensity (a.u.)'
hold (scaled_plot_ax, 'on');

subtracted_plot_ax = nexttile(t, 2);
title(subtracted_plot_ax, 'Background Subtracted Time Domain');
subtracted_plot_ax.XLabel.String = 'Time (ns)'
subtracted_plot_ax.YLabel.String = 'Intensity (a.u.)'
hold (subtracted_plot_ax, 'on');

ft_plot_ax = nexttile(t, 3);
title(ft_plot_ax, 'Fourier Transform');
ft_plot_ax.XLabel.String = 'Frequency (MHz)'
ft_plot_ax.YLabel.String = 'Intensity (a.u.)'
hold (ft_plot_ax, 'on');

for i=1:n_groups % for each sample group
    for j=1:data(i).nreps % for each rep
        curr = data(i).replicates(j);
        plot(scaled_plot_ax, curr.time, curr.scaled_echo, curr.time, ...
            curr.bkgd_fit, 'Marker', 'none', 'LineStyle', '-', ...
            'Color', data(i).color);
        plot(subtracted_plot_ax, curr.time, curr.subtracted_echo, ...
            'Marker', 'none', 'LineStyle', '-', 'Color', data(i).color);
        plot(ft_plot_ax, curr.freq, curr.ft, 'Marker', 'none', ...
            'LineStyle','-', 'Color', data(i).color);
    end
end