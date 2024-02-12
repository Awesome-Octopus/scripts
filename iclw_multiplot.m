clear all;
if license('test','Signal_Processing_Toolbox') == 1
   warndlg('You need to download and install the Signal Processing Toolbox for Matlab for this script to be able to run.')  
end

fprintf("Do you want to specify a text file containing a list of residues being selected?\nThis will be used to plot ICLW against residue number if so (yes / no):\n")
userIn = input('', 's');
validResponses = {'yes', 'no', 'NO', 'YES', 'No', 'Yes', 'y', 'n', 'Y', 'N'};
while ~ismember(userIn, validResponses)
    disp("Enter 'yes' or 'no'")
    userIn = input('', 's');
end
positive_response = {'yes', 'YES', 'Yes', 'y', 'Y'};
do_residues = false;
if ismember(userIn, positive_response)
    do_residues = true;
end
if do_residues
    [listfile, listpath] = uigetfile('*.txt', 'Select Residue List File');
    residue_list = importdata(strcat(listpath, listfile));
    [rows, cols] = size(residue_list);
    if cols ~= 1
        warndlg(fprintf("The residue list file appears to be incorrectly formatted.\nIt should be a plain text file containing a column of integers, one per line, each indicating a residue number."));
    end
    limit = rows;
    
else
    limit = -1;
end


if limit ~= -1
      filenames = cell(1,limit);
      for i=1:limit
        title = sprintf("Select spectrum for residue %d", residue_list(i));
        [new_file,new_path] = uigetfile("*.DTA", title, "MultiSelect", "off");
        filenames{i} = fullfile(new_path,new_file);
      end
else
    filenames = {};
    while true
    
            
            [new_filenames,new_path] = uigetfile("*.DTA", "Select one or more spectra, press cancel to exit", "MultiSelect", "on");
        
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
            
    
    end
end

for i=1:numel(filenames)
    [x, y, params] = eprload(filenames{i});
    spectra(i).x = x;
    spectra(i).y = y;

    spectra(i).params = params;
    spectra(i).interpol = griddedInterpolant(spectra(i).x, spectra(i).y);
    spectra(i).xq = linspace(min(spectra(i).x),max(spectra(i).x),10000);
    spectra(i).yq = spectra(i).interpol(spectra(i).xq);
    %adjust baseline based on first few points (this assumes that your
    %spectra start with a flat baseline)
    spectra(i).yq = spectra(i).yq - mean(spectra(i).yq(1:100));
    [spectra(i).peaks, spectra(i).field]  ...
        = findpeaks(spectra(i).yq, spectra(i).xq, 'MinPeakDistance', ...
        10, 'MinPeakProminence', max(spectra(i).yq)/10, 'Annotate',...
        'extents','WidthReference','halfheight');

    %flip the spectrum upside down and find the max
    spectra(i).inverseyq = -1*spectra(i).yq;

    [spectra(i).inversepeaks, spectra(i).inversefield] ...
        = findpeaks(spectra(i).inverseyq, spectra(i).xq, 'MinPeakDistance', ...
        10, 'MinPeakProminence', max(spectra(i).inverseyq)/10, 'Annotate',...
        'extents','WidthReference','halfheight');

    % if the spectrum has 2 or more detectable peaks, the correct peak is the 2nd
    % one
    if numel(spectra(i).peaks) >1 && numel(spectra(i).inversepeaks >1)
        spectra(i).centerpeakx = spectra(i).field(2);
        spectra(i).centertroughx = spectra(i).inversefield(2);
        spectra(i).iclw = 1/abs((spectra(i).centertroughx - spectra(i).centerpeakx));

    % otherwise it is the only one
    elseif numel(spectra(i).peaks) < 2 && numel(spectra(i).inversepeaks >1)
        spectra(i).centerpeakx = spectra(i).field(1);
        spectra(i).centertroughx = spectra(i).inversefield(1);
        spectra(i).iclw = 1/abs((spectra(i).centertroughx - spectra(i).centerpeakx));

    else
         warndlg('A spectrum appears to not have an assignable peak, so its width has been set to 1. Go through the code and play with the parameters for the findpeaks function probably.');
         spectra(i).iclw = 1;
    end

    if do_residues
        spectra(i).resNum = residue_list(i);
        residueNumbers(i) = spectra(i).resNum;
    end
     iclws(i) = spectra(i).iclw;
end

if do_residues
    scatter(residueNumbers, iclws, 'Filled');
    ax = gca;
    ax.YLabel.String = 'Gauss^-^1';
    ax.XLabel.String = 'Residue Number';
    ax.Title.String = 'Inverse Central Linewidths';
    set(ax,'Color','None');
else
    
    %prune off the path and extension from the filename
    disp(filenames)
    for i = 1:numel(filenames)
        cutoff = 0;
        new = '';
        for s = 1:numel(filenames{i})
            if strcmp(filenames{i}(s), "\") || strcmp(filenames{i}(s), "/")
                cutoff = s;
                if cutoff ~= numel(filenames{i})
                    new = filenames{i}(cutoff+1:end);
                end
            end
        end
        if ~isempty(new)
            filenames{i} = new;
            new = '';
        end
        for s = 1:numel(filenames{i})
            if strcmp(filenames{i}(s), ".")
                cutoff = s;
                if cutoff ~= 1
                    new = filenames{i}(1:cutoff-1);
                end
            end
        end
        if ~isempty(new)
            filenames{i} = new;
        end

    end    
    filenames = categorical(filenames);
    % Create a horizontal bar plot with customized colors
    bar(filenames, iclws, 'FaceColor', [0.5, 0.7, 0.9]);
    ax = gca;
    % Add title and axis labels
    ax.Title.String = 'Inverse Central Line Widths';
    ax.YLabel.String = 'Gauss^-^1';

end

