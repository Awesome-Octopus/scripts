% a script to plot lovely CW spectra stacked on top of each other for easy
% visualization, and making publication figures.

%% Limitations
% - this is only intended to nitroxide radical spectra and wasn't tested on
% anything else
% - the spectra are normalized when plotting based on the maximum of their
% integrated spectra, not total integrated area
% - if your center peak isnt the tallest peak, which should be a vey rare
% occurance, this program will not allign that spectrum correctly with
% others where the center peak is the tallest.



% all relevent plotting parameters and the GUI elements that select them
% are saved in fig.UserData so that it can be passed to callback functions



clc;

clear all;
close all;

%% Load the Data
filenames = {};

i = 1;

while true

    [new_filenames,new_path] = uigetfile("*.DTA", "Select one or more spectra", "MultiSelect", "on");

    % provide instructions the first time


    if isnumeric(new_filenames) == 1
        % user cancelled -> done selecting files
        break
    end

    % this appears to be necessary if you only have one entry in your list of 
    % files or multiple files in the same path
    % because in that case the datatype for the filename is just a string,
    % whereas multiples are a cell array containing strings
    if iscell(new_filenames) == 0
        new_filenames = {new_filenames};
    end
    filenames = [filenames fullfile(new_path,new_filenames)];
    i = i + 1;
end

if isempty(filenames)
    return
end    

for i=1:numel(filenames);
   [spectra(i).field, spectra(i).ampl, spectra(i).params] ...
       = eprload(filenames{i});
   
end
clear new_filenames filenames i new_path;

%% make GUI elements
fig = uifigure( "Position", [30, 100, 1200, 700],...
    "Visible", "on", "Name", "CW Plotting");

ax = uiaxes(fig,'Position',[350, 40, 820, 650]);
ax.XLabel.String = 'Magnetic Field (G)';

fig.UserData.sel_ndx_1 = 1;
fig.UserData.sel_ndx_2 = 1;

% we need to move the data to be a child of the figure so that it can be
% passed between functions smoothly
fig.UserData.spectra = spectra;


%% color pallette
fig.UserData.colorSchemes(1) = {[0 0.45 0.74;0.85 0.33 0.1;0.93 0.69 0.13;...
    0.49 0.18 0.56;0.47 0.67 0.19;0.3 0.75 0.93;0.64 0.08 0.18]}; % <--- default
fig.UserData.colorSchemes(2) = {[0 0 0]}; % <-- all black
fig.UserData.colorSchemes(3) = {[0.8 0.6 0.79;0.62 0.75 0.81;...
    0.62 0.88 0.62;1 0.69 0.27; 1 0.4 0.39]}; % <---- pastel
fig.UserData.colorSchemes(4) = {[0.21 0.12 0.4;0.91 0.24 0.67;...
    0.58 0.79 0.06;1 0.66 0;0.68 0.59 0.86; 0.16 0.82 0.79]}; % <-- retro
fig.UserData.colorSchemes(5) = {[0.65 0.82 0.8;0.95 0.83 0.53;...
    0.79 0.52 0.45;0.47 0.3 0.38]}; % <--- earth tones
fig.UserData.colorSchemes(6) = {[0 0 0; 1 0 0; 0 1 0; 0 0 1]}; % <-- bold
fig.UserData.colorSchemes(7) = {[0 0 0; 1 0 0]}; % <-- alternating black and red

fig.UserData.colorNames = {'Default', 'Black', 'Pastel', 'Retro', 'Earth Tones', 'Bold', 'Alternating'};


% save plot options and dialog box
fig.UserData.file_filter = {'*.fig';'*.png';'*.jpg';'*.tif';'*.ps';'*.pdf';'*.eps';'*.ai'};
fig.UserData.saveBtn = uibutton(fig, 'push', 'Position', [100, 610, 120, 60], ...
    'Text', 'Save Plot', 'ButtonPushedFcn', @(src, event)save_spectra(ax, fig));

%% ---------------------!!!!!!!!!!!!!!!!!!!!!!---------------------------
% these ui elements will not work no matter what for linux. UIedit field
% and uidropdown
% will all not let you pick anything (but only at this one particular label)
% An improvised work arround is to use 2 uibuttons to increase/decrease the selection index


fig.UserData.selected = uilabel (fig, ... %'Limits', [1 numel(raw_spectra)], ...
      'Position', [35 540 200 55], 'text', spectra(1).params.TITL, 'WordWrap','on');

fig.UserData.sel1decr = uibutton (fig, 'push', 'Position', [250, 550, 20, 20], ...
    'Text', '<', 'ButtonPushedFcn', @(src, event) select_spectra(fig, src));
fig.UserData.sel1incr = uibutton (fig, 'push', 'Position', [275, 550, 20, 20], ...
    'Text', '>', 'ButtonPushedFcn', @(src, event) select_spectra(fig, src));

%uilabel (fig, 'Position', [35 375 250 22], 'Text', 'Spectrum 2: ');
fig.UserData.selected2 = uilabel (fig, ... %'Limits', [1 numel(raw_spectra)], ...
     'Position', [35 485 200 55], 'text', spectra(1).params.TITL, 'WordWrap','on');

fig.UserData.sel2decr = uibutton (fig, 'push', 'Position', [250, 495, 20, 20], ...
    'Text', '<', 'ButtonPushedFcn', @(src, event) select_spectra(fig, src));
fig.UserData.sel2incr = uibutton (fig, 'push', 'Position', [275, 495, 20, 20], ...
    'Text', '>', 'ButtonPushedFcn', @(src, event) select_spectra(fig, src));

fig.UserData.swapBtn = uibutton(fig, 'push', 'Position', [75, 395, 170, 50], ...
    'Text', 'Swap Spectra Positions', 'ButtonPushedFcn', ...
    @(src, event)swap_spectra(ax, fig), 'Enable','on');

uilabel (fig,'Position', [35 355 250 22], 'Text', 'Font Size:');
fig.UserData.fontSize = uispinner(fig, "Limits", [1 100], "Position", [225 355 65 26], ...
        'Value', 14, 'ValueChangedFcn', @(src, event)update(ax, fig));

uilabel (fig,'Position', [35 300 250 50], 'Text', 'Number of points for smoothing:');
fig.UserData.smoothingPts = uispinner(fig, "Limits", [0 15], "Position", [225 312 65 26], ...
        'Value', 0, 'ValueChangedFcn', @(src, event)update(ax, fig));

uilabel (fig,'Position', [35 270 250 22], 'Text', 'Horizontal staggering (G):');
fig.UserData.horizStagger = uispinner(fig, "Position",[225, 270, 65 26],...
        'Value', 0, 'Step', 0.1, 'Editable', 'on', ...
        'ValueChangedFcn', @(src, event)update(ax, fig));

uilabel (fig,'Position', [35 220 250 22], 'Text', 'Vertical staggering:');
fig.UserData.vertStagger = uispinner(fig,"Position",[225, 220, 65 26],...
        'Value', 0.2, 'Step', 0.05, 'Editable', 'on', ...
        'ValueChangedFcn', @(src, event)update(ax, fig));

uilabel (fig,'Position', [35 170 250 22], 'Text', 'Line thickness:');
fig.UserData.LineThickness = uispinner(fig,"Position",[130, 170, 65 26],...
        'Value', 2, 'Step', 0.1, 'Editable', 'on', ...
        'ValueChangedFcn', @(src, event)update(ax, fig));

uilabel (fig, 'Position', [35 120 150 22], 'Text', 'Color Scheme:');
if ispc
    fig.UserData.pallete = uidropdown(fig, 'position', [130 120 120 22], ...
        'Value', fig.UserData.colorSchemes{1}, 'Items', fig.UserData.colorNames, ...
        'ItemsData', fig.UserData.colorSchemes, ...
        'ValueChangedFcn', @(src, event)update(ax, fig));
else
    fig.UserData.pallete = uispinner(fig, 'position', [130 120 65 22], ...
        'Value', 1, 'Step', 1, 'Limits', [1 numel(fig.UserData.colorSchemes)], ...
    'ValueChangedFcn', @(src, event)update(ax, fig));
end

uilabel(fig, 'Position', [237 192 100 22], 'Text', 'show axes:');
fig.UserData.axes = uiswitch(fig, 'slider', 'position', [240 170 200 22], ...
    'ItemsData', [false true], 'Items', {'off', 'on'}, ...
    'ValueChangedFcn', @(src, event)update(ax, fig));

uilabel(fig, 'Position', [237 140 100 22], 'Text', 'Normalize:');
fig.UserData.norm = uiswitch(fig, 'slider', 'position', [240 118 200 22], ...
    'ItemsData', [false true], 'Items', {'off', 'on'}, 'Value', true, ...
    'ValueChangedFcn', @(src, event)update(ax, fig));

if ispc % <-- dropdowns only work on PC not linux .
    uilabel(fig, 'Position', [70 80 100 22], 'Text', 'Annotation style:');
    fig.UserData.annotation = uidropdown(fig,"Value", 'labels', 'position', ...
    [60 60 100 22], 'Items', {'labels', 'legend', 'none'}, ...
    'ValueChangedFcn', @(src, event)update(ax, fig));

else
    uilabel(fig, 'Position', [70 80 100 22], 'Text', 'Annotation style:');
    fig.UserData.annotation = uiswitch(fig, 'slider', 'position', [80 60 200 22], ...
    'ItemsData', {'labels', 'legend'}, 'Items', {'labels', 'legend'}, ...
    'ValueChangedFcn', @(src, event)update(ax, fig));
end 

uilabel(fig, 'Position', [220 80 100 22], 'Text', 'label offset');
fig.UserData.labelOffset = uispinner(fig, 'Position', [210 60 85 22], 'Value', 0.1, ...
    'Step', 0.01, 'ValueChangedFcn', @(src, event)update(ax, fig));


%% plot the spectra once intially as is
for i =1:numel(spectra)
        plot(ax, spectra(i).field, ...
            spectra(i).ampl, 'LineWidth', 1);
        hold (ax, "on");

        %this makes the text labels for spectra as invisible placeholders for
        %now
        ax.UserData(i) = text(ax,0, 0,'','Visible','off');
end

hold(ax, "off");


% this creates a button down function that gets a callback whenever a click
% is made in the plotting area
%ax.ButtonDownFcn = @(src, event) ParseClick (src, event, ax, fig);


% call update automatically for initial plotting with spectra in an
% unmodified order, future calls to swap will reorder the array of raw
% spectra, but then all processing is done over by the update function
update(ax, fig);

%% opens save dialog
function save_spectra (ax, fig)
    [fname, path, ndx] = uiputfile(fig.UserData.file_filter, ...
        'Save the plotted spectra');
    fullpath = strcat(path,fname);
    if ndx == 1
         saveas(fig,fullpath);
    elseif ndx > 1
         % if we want to save a picture only save the axis and not the GUI elements
         % so convert axes to image
         exportgraphics(ax,fullpath);
    else    
        % otherwise if dialog cancelled
    end

end

%% Changes the plotting parameters after reordering is done. 
function update (ax, fig)        

    processed_spectra = process_spectra(fig);
    staggered_spectra = stagger (processed_spectra, ...
        fig.UserData.horizStagger.Value, fig.UserData.vertStagger.Value); 

    % change the colors to the currently selected color theme
    if ispc
        ax.ColorOrder = fig.UserData.pallete.Value;
    else
        ax.ColorOrder = fig.UserData.colorSchemes{fig.UserData.pallete.Value};
    end

    if strcmp(fig.UserData.annotation.Value, 'labels') == 1
    % with data labels, ax.Children is arranged like:
    % 1 - labeltext 1
    % 2 - line for spectrum 1
    % 3 - labeltext 2
    % 4 - line for spectrum 2
        for i=1:numel(staggered_spectra)
        % change the values stored in x and y data in the ax object
        ax.Children(i*2).XData = staggered_spectra(i).field;
        ax.Children(i*2).YData = staggered_spectra(i).ampl;
        ax.Children(i*2).LineWidth = fig.UserData.LineThickness.Value;
        

       % now adjust position of label and give it text
        
        ax.Children(2*i-1).Position = [staggered_spectra(i).field(end), ...
            (mean(staggered_spectra(i).ampl(end-20:end)) + ...
            fig.UserData.labelOffset.Value), 0];
        ax.Children(2*i-1).HorizontalAlignment = 'right';
        ax.Children(2*i-1).String = staggered_spectra(i).params.TITL;
        ax.Children(2*i-1).FontSize = fig.UserData.fontSize.Value;
        ax.Children(2*i-1).Interpreter = "none";
        ax.Children(2*i-1).Visible = "on";

        end
        legend(ax,'off');
    elseif strcmp(fig.UserData.annotation.Value, 'legend') == 1
    % if the user wants a legend instead
        l = legend(ax);
        l.Interpreter = "none";
        for i=1:numel(staggered_spectra)
        % change the values stored in x and y data in the ax object
            ax.Children(i*2).XData = staggered_spectra(i).field;
            ax.Children(i*2).YData = staggered_spectra(i).ampl;
            ax.Children(i*2).LineWidth = fig.UserData.LineThickness.Value;
            ax.Children(2*i-1).Visible = "off";
            l.String{end-i+1} = staggered_spectra(i).params.TITL;
            l.FontSize = fig.UserData.fontSize.Value;
        end
    else
        for i=1:numel(staggered_spectra)
            ax.Children(2*i-1).Visible = "off";
        end    
        legend(ax,'off');    
    end

    if fig.UserData.axes.Value == false
        set(ax, 'Visible', 'off');
    else
        set(ax, 'Visible', 'on');
    end
    
end

%% process spectrum
function processed_spectra = process_spectra (fig)

    input_spectra = fig.UserData.spectra;
    smoothingPts = fig.UserData.smoothingPts.Value;

    processed_spectra = input_spectra; % << initialize by copying
    for s = 1:numel(input_spectra)

        
        %perform n point data smoothing
        processed_spectra(s).ampl = datasmooth(input_spectra(s).ampl, ...
            smoothingPts);
        % perform baseline correction
        processed_spectra(s).ampl = processed_spectra(s).ampl - ...
            mean(processed_spectra(s).ampl(1:20));
        
        % in order to accurately allign spectra we should allign not to the
        % max point of the derivative spectrum but to the maximum of the
        % integrated absorbance spectrum. this will correspond to where the
        % baseline corrected derivative spectrum crosses the 0 point on the
        % center line. note: if for some reason your center peak isn't the
        % tallest this will grab the wrong peak. 


        absorption_spectra(s).ampl = cumtrapz(processed_spectra(s).ampl);
        
        %get the index of the max value of the integrated spectrum, and use it to get
        %the field where the value is max.
        [~, ndx] = max(absorption_spectra(s).ampl);

        % find the total range of each derivative spectrum
        peakHeight(s) = max(processed_spectra(s).ampl) - min(processed_spectra(s).ampl);
        centerfield(s) = processed_spectra(s).field(ndx);
    end
    

    scalingfactor = max(peakHeight);
    allignmentpoint = mean(centerfield);

    % if the toggle is off don't normalize each spectrum to max height
    if fig.UserData.norm.Value == true
        for s =1:numel(input_spectra)
            processed_spectra(s).field = input_spectra(s).field - ...
                (centerfield(s) - allignmentpoint);
            processed_spectra(s).ampl = processed_spectra(s).ampl ...
                * scalingfactor / peakHeight(s);
        end
    else
        for s =1:numel(input_spectra)
            processed_spectra(s).field = input_spectra(s).field - ...
                (centerfield(s) - allignmentpoint);
        end
    end

end

function staggered_spectra  = stagger (input_spectra, horizStagger, vertStagger)
    for s = 1:numel(input_spectra)
        %adjust center so that that maxima line up with staggering added for
        %horizontal offset 

        staggered_spectra(s).field = input_spectra(s).field - horizStagger * (s-1);

        %scale relative to the largest spectrum height and apply vertical
        %offset
        staggered_spectra(s).ampl = input_spectra(s).ampl - vertStagger * (s-1);
        staggered_spectra(s).params = input_spectra(s).params;
    end
end

function swap_spectra (ax, fig)

    % swap elements for the list of spectra
    tmp = fig.UserData.spectra;

    fig.UserData.spectra(fig.UserData.sel_ndx_2) = tmp(fig.UserData.sel_ndx_1);
    fig.UserData.spectra(fig.UserData.sel_ndx_1) = tmp(fig.UserData.sel_ndx_2);
    
    % swap what they are pointing to to reflect their swapped positions
    a = fig.UserData.sel_ndx_1;
    b = fig.UserData.sel_ndx_2;
    fig.UserData.sel_ndx_1 = b;
    fig.UserData.sel_ndx_2 = a;

    update(ax, fig);

end

function select_spectra(fig, src)


    % find what is currently in the selected spectrum displayed for #1 and
    % #2 this is the necessary workarround snce we can't use dropdown menus
    % on linux


    i = 1;
    while strcmp(fig.UserData.spectra(i).params.TITL, fig.UserData.selected.Text) ~= 1
         i = i + 1;
    end
    fig.UserData.sel_ndx_1 = i;
    
    i = 1;
    while strcmp(fig.UserData.spectra(i).params.TITL, fig.UserData.selected2.Text) ~= 1
         i = i + 1;
    end
    fig.UserData.sel_ndx_2 = i;
    
    % then find which button was pressed and whether a spectrum that is
    % before or after it can be selected
    if src == fig.UserData.sel1decr
        if strcmp(fig.UserData.selected.Text, fig.UserData.spectra(1).params.TITL) == 1
        % if you attempt to select an option before the first one set it to
        % 1 instead
            fig.UserData.sel_ndx_1 = 1;
            
        else
            % set the value to that found at index - 1
            fig.UserData.selected.Text = fig.UserData.spectra(fig.UserData.sel_ndx_1 -1).params.TITL;
            % decrement the index and keep record of the index for use with swap spectra function
            fig.UserData.sel_ndx_1 = fig.UserData.sel_ndx_1 - 1;
            
        end
    elseif src == fig.UserData.sel1incr
        if strcmp(fig.UserData.selected.Text, fig.UserData.spectra(end).params.TITL) == 1
        % if you attempt to select an option to be after the end of the
        % list set it to the end instead
            fig.UserData.sel_ndx_1 = numel(fig.UserData.spectra);
            
        else
            % set the value to that found at index + 1
            fig.UserData.selected.Text = fig.UserData.spectra(fig.UserData.sel_ndx_1 +1).params.TITL;
            % increment the index and keep record of the index for use with swap spectra function
            fig.UserData.sel_ndx_1 = fig.UserData.sel_ndx_1 + 1;
            
        end
    elseif src == fig.UserData.sel2decr
        if strcmp(fig.UserData.selected2.Text, fig.UserData.spectra(1).params.TITL) == 1
        % if you attempt to select an option before the first one set it to
        % 1 instead
            fig.UserData.sel_ndx_2 = 1;
           
        else
            % set the value to that found at index - 1
            fig.UserData.selected2.Text = fig.UserData.spectra(fig.UserData.sel_ndx_2-1).params.TITL;
            % decrement the index and keep record of the index for use with swap spectra function
            fig.UserData.sel_ndx_2 = fig.UserData.sel_ndx_2 - 1;
            
        end
    elseif src == fig.UserData.sel2incr
        if strcmp(fig.UserData.selected2.Text, fig.UserData.spectra(end).params.TITL) == 1
        % if you attempt to select an option to be after the end of the
        % list set it to the end instead
            fig.UserData.sel_ndx_2 = numel(fig.UserData.spectra);
            
        else
            % set the value to that found at index + 1
            fig.UserData.selected2.Text = fig.UserData.spectra(fig.UserData.sel_ndx_2+1).params.TITL;
            % increment the index and keep record of the index for use with swap spectra function
            fig.UserData.sel_ndx_2 = fig.UserData.sel_ndx_2 + 1;
            
        end
    end
end

%% wrapper function for update
% whenever a cliuck is made in the plotting area get the cursor position,
% then check if it corresponds to a graphics object, if it does, pass that
% handle to update