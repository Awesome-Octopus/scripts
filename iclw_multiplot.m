if license('test','Signal_Processing_Toolbox') == 1
   warndlg('You need to download and install the Signal Processing Toolbox for Matlab for this script to be able to run.')  
end

warndlg('To use the iclw_multiplot script, you first need to write a plain .txt file with the residue numbers for your desired spectra. One number per line. Select the text file when the script starts, then load each spectrum individually from the selection dialog the in the order listed in your text file.')
[file, path] = uigetfile({'*.txt'}, 'Select the text file with your list of residue numbers')
resNum = importdata(strcat(path, file));


for i=1:numel(resNum)
    [x, y, params] = eprload();
    spectra(i).x = x;
    spectra(i).y = y;
    spectra(i).resNum = resNum(i);
    spectra(i).params = params;
    spectra(i).interpol = griddedInterpolant(spectra(i).x, spectra(i).y);
    spectra(i).xq = linspace(min(spectra(i).x),max(spectra(i).x),10000);
    spectra(i).yq = spectra(i).interpol(spectra(i).xq);
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
        spectra(i).iclw = 1/(spectra(i).centertroughx - spectra(i).centerpeakx);
    
    % otherwise it is the only one
    elseif numel(spectra(i).peaks) < 2 && numel(spectra(i).inversepeaks >1)
        spectra(i).centerpeakx = spectra(i).field(1);
        spectra(i).centertroughx = spectra(i).inversefield(1);
        spectra(i).iclw = 1/(spectra(i).centertroughx - spectra(i).centerpeakx);
                  
    else
         warndlg('A spectrum appears to not have an assignable peak, so its width has been set to 1. Go through the code and play with the parameters for the findpeaks function probably.');
         spectra(i).iclw = 1;
    end

    
     
     residueNumbers(i) = spectra(i).resNum;
     iclws(i) = spectra(i).iclw;
end

scatter(residueNumbers, iclws, 'Filled');
ax = gca;
ax.YLabel.String = 'Gauss^-^1';
ax.XLabel.String = 'Residue Number';
ax.Title.String = 'Inverse Central Linewidths';
set(ax,'Color','None');